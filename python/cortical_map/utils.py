import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


def plot_cells(X, Y, delta, param1, param2):
    fig = plt.figure(figsize=(12, 12), dpi=400)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    Y = Y+delta
    ax.plot(X[:, 0], X[:, 1], X[:, 2], '.', ms=.5, **param1)
    ax.plot(Y[:, 0], Y[:, 1], Y[:, 2], '.', ms=.5, **param2)
    ax.set_xlabel(r'parallel to cut [$\mu$m]')
    ax.set_ylabel(r'perpendicular to cut [$\mu$m]')
    ax.set_zlabel(r'relative cortical depth [$\mu$m]')

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
