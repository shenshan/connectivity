import datajoint as dj
import numpy as np
import matplotlib.pyplot as plt
from cortical_map.utils import rasterize, plot_cells, spin_shuffle, compute_overlap_density
from . import connectivity as conn

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
        (1, 'dendrite'), (2, 'cellbody'), (3, 'axon')
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
class OrthoOverlap(dj.Computed):
    definition = """
    -> conn.Distance
    -> DensityParameters
    ---
    d           : longblob # distance to the connecting line between two neurons
    p           : longblob # overlap density along that line
    """

    @property
    def populated_from(self):
        return CellType().project(ctn_from='cell_type_name') * CellType().project(
            ctn_to='cell_type_name') * DensityParameters()

    def _make_tuples(self, key):
        AXON = (CellRegion() & dict(cell_region_name='axon')).fetch1['cell_region_id']
        DENDRITE = (CellRegion() & dict(cell_region_name='dendrite')).fetch1['cell_region_id']

        raster_resolution, bin_width, cell_offset, multiplication_factor, cut_compensation = \
            (DensityParameters() & key).fetch1['raster_resolution','bin_width',
                                               'cell_offset','multiplication_factor', 'cut_compensation']
        layer_from, morph_from = (CellType() & dict(cell_type_name=key['ctn_from'])).fetch1['layer', 'cell_type_morph']
        cells_from = (conn.ConnectMembership() * conn.Cell()
                      & dict(role='from', cell_layer=layer_from, cell_type_morph=morph_from)
                      ).project(cell_id_from='cell_id')
        X_from, region_from = OrthoOverlap.load_cells(
            Tree() * CellType() & dict(layer=layer_from, cell_type_morph=morph_from),
            resolution=raster_resolution, correct_cut=cut_compensation)

        layer_to, morph_to = (CellType() & dict(cell_type_name=key['ctn_to'])).fetch1['layer', 'cell_type_morph']
        cells_to = (conn.ConnectMembership() * conn.Cell()
                    & dict(role='to', cell_layer=layer_to, cell_type_morph=morph_to)
                        ).project( cell_id_to='cell_id')
        X_to, region_to = OrthoOverlap.load_cells(Tree() * CellType() & dict(layer=layer_to, cell_type_morph=morph_to),
                                                  resolution=raster_resolution, correct_cut=cut_compensation)


        keys, offset_x, offset_y = (conn.Distance() * cells_from * cells_to).fetch[dj.key, 'delta_x', 'delta_y']
        offsets = np.c_[offset_x, -cell_offset * np.ones_like(offset_x), offset_y]


        fig, ax = plot_cells(X_from, X_to, offsets[0], dict(color='silver', label=key['ctn_from']),
                             dict(color='skyblue', label=key['ctn_to']))
        fig.savefig('{ctn_from}_to_{ctn_to}.png'.format(**key))


        X_from = spin_shuffle(X_from, copy=multiplication_factor)
        X_to = spin_shuffle(X_to, copy=multiplication_factor)

        for k, delta in zip(keys, offsets):
            H, E = compute_overlap_density(X_from[region_from == AXON], X_to[region_to == DENDRITE], bin_width, delta)
            y = E[1]
            y = 0.5 * (y[1:] + y[:-1])
            idx = y > 0
            p = H[:,idx,:].sum(axis=(0,2))
            ins_key = (conn.Distance() * DensityParameters() & k).project().fetch1()
            ins_key['p'] = p
            ins_key['d'] = y[idx]
            self.insert1(ins_key)


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
                print("correcting cut")
                ymax = X[:,1].max()
                idx = X[:,1] <= -ymax
                tmp  = X[idx, :]
                tmp[:,1] *= -1
                X = np.vstack((X, tmp))
                region = np.hstack((region, region[idx]))

            cells.append(X)
            regions.append(region)
        return np.vstack(cells), np.hstack(regions)
