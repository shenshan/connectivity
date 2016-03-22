import datajoint as dj
import numpy as np
import matplotlib.pyplot as plt
from cortical_map.utils import rasterize, plot_cells, spin_shuffle, compute_overlap_density
from . import connectivity as conn
from itertools import product

schema = dj.schema('shan_reconstruction', locals())


def bugfix_reshape(a):
    return a.ravel(order='C').reshape(a.shape, order='F')


@schema
class CellRegion(dj.Lookup):
    definition = """
    # table for cell region labels

    cell_region_id      : tinyint   # cell region identifier
    ---

    cell_region_name    : char(8) # name of the cell region
    """

    contents = [
        (1, 'axon'), (2, 'cellbody'), (3, 'dendrite')
    ]


@schema
class DensityParameters(dj.Lookup):
    definition = """
    density_param_id        : tinyint            # identifier of the parameter set
    ---
    bin_width               : float              # bin width for the 3d histogram in mu
    raster_resolution       : float              # resolution of the rasterization in mu
    cell_offset             : float              # assumed offset between cells into the slice
    multiplication_factor   : int                # multiply the number of points by this factor when spin shuffling
    cut_compensation        : tinyint unsigned   # try to account for the fact that parts of the neuron might be missing
    """

    contents = [
        dict(density_param_id=0, bin_width=15, raster_resolution=.5, cell_offset=0, multiplication_factor=2,
             cut_compensation=1)
    ]

@schema
class CellTypePair(dj.Lookup):
    definition = """
    # table with cell type name pairs fetched from CellType

    cell_type_from      : varchar(256)
    cell_type_to        : varchar(256)
    ---
    """

    @property
    def contents(self):
        ct = CellType().fetch['cell_type_name']
        yield from product(ct,ct)

@schema
class OverlapGroup(dj.Computed):
    definition = """
    # Overlap densities grouping table
    -> CellTypePair
    -> DensityParameters
    ---
    """

    class OrthoOverlapDensity(dj.Part):
        definition = """
        # overlap densities for a particular cell type pair and distance
        -> OverlapGroup
        -> conn.Distance
        ---
        d           : longblob # distance to the connecting line between two neurons
        p           : longblob # overlap density along that line
        """


    def _make_tuples(self, key):
        # get identifiers for different regions
        AXON = (CellRegion() & dict(cell_region_name='axon')).fetch1['cell_region_id']
        DENDRITE = (CellRegion() & dict(cell_region_name='dendrite')).fetch1['cell_region_id']

        # get parameters for density estimation
        raster_resolution, bin_width, cell_offset, multiplication_factor, cut_compensation = \
            (DensityParameters() & key).fetch1['raster_resolution','bin_width',
                                               'cell_offset','multiplication_factor', 'cut_compensation']

        # load all cells from the particular pair of cell types
        X, region, cells = {}, {}, {}
        for role in ('from','to'):
            try:
                layer, morph = (CellType() & dict(cell_type_name=key['cell_type_' + role])).fetch1['layer', 'cell_type_morph']
            except:
                #----------------------------------
                # TODO: Remove this later
                from IPython import embed
                embed()
                exit()
                #----------------------------------

            cells[role] = (conn.ConnectMembership() * conn.Cell()
                          & dict(role=role, cell_layer=layer, cell_type_morph=morph)
                          ).project(**{('cell_id_' +role):'cell_id'})
            X[role], region[role] = OverlapGroup.load_cells(Tree() * CellType()
                                                            & dict(layer=layer, cell_type_morph=morph),
                    resolution=raster_resolution, correct_cut=cut_compensation)


        # when there are no distances
        if not conn.Distance() * cells['from'] * cells['to']:
            self.insert1(key)
            return
        
        # get deltas
        keys, delta_x, delta_y = (conn.Distance() * cells['from'] * cells['to']).fetch[dj.key, 'delta_x', 'delta_y']
        offsets = np.c_[delta_x, -cell_offset * np.ones_like(delta_x), delta_y]

        
        # save a control figure
        fig, ax = plot_cells(X['from'], X['to'], offsets[0], dict(color='silver', label=key['cell_type_from']),
                             dict(color='skyblue', label=key['cell_type_to']))
        fig.savefig('{cell_type_from}_to_{cell_type_to}.png'.format(**key))

        plt.close(fig)
        # shuffle by rotation around z-axis and increase data volume
        X['from'] = spin_shuffle(X['from'], copy=multiplication_factor)
        X['to'] = spin_shuffle(X['to'], copy=multiplication_factor)
        for role in region:
            region[role] = np.tile(region[role], multiplication_factor)
            
        # insert grouping key
        self.insert1(key)
        part = OverlapGroup.OrthoOverlapDensity()

        # insert an element in part (sub) table for each distance between cells of that type
        for k, delta in zip(keys, offsets):

            H, E = compute_overlap_density(X['from'][region['from'] == AXON], X['to'][region['to'] == DENDRITE], bin_width, delta)

            # compute bin centers along y-axis and only get positive side
            y = E[1]
            y = 0.5 * (y[1:] + y[:-1])
            idx = y > 0

            # marginalize density
            p = H[:,idx,:].sum(axis=(0,2))
            key.update({f:k[f] for f in k.dtype.fields})
            key.pop('cell_id_from')
            key.pop('cell_id_to')
            key['p'] = p
            key['d'] = y[idx]
            part.insert1(key)


    @staticmethod
    def load_cells(trees, resolution=.5, correct_cut=False):
        """
        Loads cells from the database
        :param cell_type: string containing layer and cell type morphology
        :param resolution: resolution for rasterization (default = .5 mu)
        :return: rasterized 3d points, region labels
        """
        cells = []
        regions = []

        CELLBODY = (CellRegion() & dict(cell_region_name='cellbody')).fetch1['cell_region_id']

        for k, node_coords, connected_pairs, region in \
                zip(*trees.fetch[dj.key, 'node_coords', 'connected_pairs', 'node_region']):
            node_coords = bugfix_reshape(node_coords)
            node_coords[:, [1, 2]] = node_coords[:,
                                     [2, 1]]  # TODO: check with Shan and Xiaolong whether cells got rotated again

            connected_pairs = bugfix_reshape(connected_pairs.astype(int)) - 1  # matlab to python index shift
            region = region.squeeze()
            if not np.any(region == CELLBODY):
                print(k, 'has no cellbody')
                continue
            mu = np.mean(node_coords[region == CELLBODY, :], axis=0)
            node_coords -= mu


            X, region = rasterize(node_coords, connected_pairs, region, resolution)
            if correct_cut:
                ymax = X[:,1].max()
                idx = X[:,1] <= -ymax
                tmp  = X[idx, :]
                tmp[:,1] *= -1
                X = np.vstack((X, tmp))
                region = np.hstack((region, region[idx]))

            cells.append(X)
            regions.append(region)
        return np.vstack(cells), np.hstack(regions)
