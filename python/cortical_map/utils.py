import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable
from matplotlib.pyplot import fill, cm, Rectangle
import os

def rasterize(node_coords, connected_pairs, region, resolution=1):
    """
    Rasterizes a neuron graph

    :param node_coords: coordinates of the single nodes
    :param connected_pairs: index pairs of connected nodes
    :param region: label vector for regions (1=DENDRITE, 2=CELLBODY, 3=AXON)
    :param resolution: raster step width along the edges
    :return: rasterized points in 3d, new region labels
    """
    X = []
    reg = []
    for fro, to in connected_pairs:
        node = node_coords[fro]
        conn = node_coords[to]
        v = (conn - node)
        l = np.sqrt((v ** 2).sum())
        if l > 0:
            v /= l
            for delta in np.arange(resolution, l, resolution):
                X.append(node + delta * v)
                reg.append(region[fro])
    return np.vstack((X, node_coords)), np.hstack((reg, region))


def spin_shuffle(X, copy=1):
    """
    Shuffles the rows of X by random rotations in the xy plane.
    :param X: data
    :param copy: multiplication factor of the data; there will be len(X)*mult datapoints in the return value
    :return: multiplied and rotation shuffled data
    """
    if copy > 1:
        X = np.vstack([np.array(X) for _ in range(copy)])
    return np.vstack([spin(x, th) for x, th in zip(X, np.random.rand(len(X)) * 2 * np.pi)])


def spin(X, theta):
    """
    Rotates the rows of X by theta.

    :param X: data matrix
    :param theta: angle in radians
    :return: rotated data
    """
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return np.dot(X, R.T)


def plot_cells(X, Y, delta, param1, param2, threed=False):

    Y = Y+delta
    sns.set_style('whitegrid')
    fig = plt.figure(figsize=(12, 12), dpi=400)

    if threed:
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot(X[:, 0], X[:, 1], X[:, 2], '.', ms=.5, **param1)
        ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], '.', ms=.5, **param2)
        ax.set_zlabel(r'relative cortical depth [$\mu$m]')
    else:
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(X[:, 1], X[:, 2], '.', ms=.5, **param1)
        ax.plot(Y[:, 1], Y[:, 2], '.', ms=.5, **param2)

    ax.set_xlabel(r'parallel to cut [$\mu$m]')
    ax.set_ylabel(r'perpendicular to cut [$\mu$m]')

    lgnd = ax.legend()
    for h in lgnd.legendHandles:
        h._legmarker.set_markersize(15)

    return fig, ax

def compute_overlap_density(X, Y, bin_width, delta):
    """
    Computes the overlap density between X and Y shifted by dist along the x axis.

    :param X: data points for neuron 1
    :param Y: data points for neuron 2
    :param bin_width: bin width. Bins will be bin_width x bin_width x bin_width
    :param delta: amount by which Y will be shifted relative to X
    :return: overlap density, bin borders of the histograms
    """
    Y = Y + delta
    ma = np.vstack((X, Y)).max(axis=0)
    mi = np.vstack((X, Y)).min(axis=0)
    bins = tuple((np.hstack((np.arange(0, low, -bin_width)[::-1], np.arange(bin_width, high + bin_width, bin_width)))
                  for low, high in zip(mi, ma)))

    H1, _ = np.histogramdd(X, bins)
    H2, E = np.histogramdd(Y, bins)
    H1 /= bin_width ** 3
    H2 /= bin_width ** 3
    return H1 * H2, E

def extended_hinton(ax, V, C, vmax=None, cmin=None, cmax=None, cmap=None, matrix_style=False, alpha=1,
                    enforce_box=False):
    if cmap is None:
        cmap = cm.jet

    if vmax is None:  vmax = np.amax(np.abs(V))
    if cmax is None:  cmax = np.amax(C)
    if cmin is None:  cmin = np.amin(C)

    cnorm = Normalize(vmin=cmin, vmax=cmax, clip=True)
    cmapable = ScalarMappable(norm=cnorm, cmap=cmap)

    if matrix_style:
        V, C = V.T, C.T

    ax.patch.set_facecolor([0, 0, 0, 0])

    for (x, y), w in np.ndenumerate(V):
        s = C[x, y]
        color = cmap(cnorm(s))  # cmap(s / cmax)
        size = np.abs(w / vmax)
        rect = Rectangle([x - size / 2, y - size / 2], size, size,
                         facecolor=color, edgecolor=color, alpha=alpha)
        ret = ax.add_patch(rect)

    if enforce_box:
        ax.axis('tight')
        try:
            ax.set_aspect('equal', 'box')
        except:
            pass
    #ax.autoscale_view()
    #ax.invert_yaxis()
    return cnorm

def plot_connections(P, labels, vmax=None):
    colors = ["black", "azure","apple","golden yellow",   "neon pink"]
    cmap = sns.blend_palette(sns.xkcd_palette(colors), as_cmap=True)

    sns.set_style('whitegrid')
    sns.set_context('paper',rc={"lines.linewidth": 2, 'font.size':12, 'font.family':'Helvetica'})
    fig = plt.figure(figsize=(4.6,4.6), dpi=300)
    gs = plt.GridSpec(5,20)

    ax = fig.add_subplot(gs[:,:17])
    n = len(labels)

    cnorm = extended_hinton(ax, P, P, matrix_style=True, cmap=cmap,
                            vmax=vmax if vmax is not None else P.max(),
                            cmin=0, cmax=P.max(),enforce_box=True)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)
#     ax.set_title(r'full connection probability matrix', fontweight='bold')
    ax.set_xlim((-1,n))
    ax.set_ylim((-1,n))
    ax.set_ylabel('postsynaptic')
    ax.set_xlabel('presynaptic')
    ax_color = fig.add_subplot(gs[1:4,18:])
    cbar = ColorbarBase(ax_color, cmap=cmap, norm=cnorm)
    fig.tight_layout()
    ax.invert_yaxis()

    return fig, {'matrix':ax, 'color':ax_color}

def layer(name):
    if 'L1' in name:
        return 1
    elif 'L23' in name:
        if 'Pyr' in name:
            return 3
        else:
            return 2
    else:
        if 'Pyr' in name:
            return 5
        else:
            return 4


def load_data(transpose=False, collapse=None, remove=None):
    path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    with open(path + '/files/matrix09112015.csv', 'r') as fid:
        labels = [e.strip() for e in fid.readline().split(',')[1:]]
        K, N = [], []
        for l in fid.readlines():
            K.append([list(map(float, e.strip().split('/')))[0] for e in l.split(',')[1:]])
            N.append([list(map(float, e.strip().split('/')))[1] for e in l.split(',')[1:]])
    layers = [layer(name) for name in labels]
    K = np.asarray(K)
    N = np.asarray(N)

    if collapse is not None:
        for k, v in collapse.items():
            idx = list(sorted([labels.index(e) for e in v]))
            i = idx[0]
            labels[i] = k
            K[i, :] = K[idx, :].sum(axis=0)
            K[:, i] = K[:, idx].sum(axis=1)
            N[i, :] = N[idx, :].sum(axis=0)
            N[:, i] = N[:, idx].sum(axis=1)
            K = np.delete(K, idx[1:], axis=0)
            K = np.delete(K, idx[1:], axis=1)
            N = np.delete(N, idx[1:], axis=0)
            N = np.delete(N, idx[1:], axis=1)
            labels = [labels[i] for i in range(len(labels)) if i not in idx[1:]]
            layers = [layers[i] for i in range(len(layers)) if i not in idx[1:]]
    if remove is not None:
        idx = list(sorted([labels.index(e) for e in remove]))
        K = np.delete(K, idx, axis=0)
        K = np.delete(K, idx, axis=1)
        N = np.delete(N, idx, axis=0)
        N = np.delete(N, idx, axis=1)
        labels = [labels[i] for i in range(len(labels)) if i not in idx]
        layers = [layers[i] for i in range(len(layers)) if i not in idx]

    if transpose:
        return labels, layers, K.T, N.T
    else:
        return labels, layers, K, N

