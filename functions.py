import igraph as ig
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cairocffi
from scipy import stats
import scanpy as sc
import pandas as pd


def read_and_filter_grn(grn_path, var=0, ci=0, net=None):
    grn = pd.read_table(grn_path, header=0)
    grn['id'] = grn[['regulator', 'target']].apply(lambda x: '{}_{}'.format(x[0], x[1]), axis=1)
    if var > 0:
        grn = grn[grn['var.exp.median'] >= var]
    if ci > 0:
        grn = grn[grn['combined_confidences.values'] >= ci]
    if net is not None:
        grn['net'] = net

    return grn


def shared_edges(grns):
    edges = [grns[x].id for x in grns]
    ret = set(edges[0])
    for x in range(len(edges)-1, 0, -1):
        ret = ret & set(edges[x]) & set(edges[x-1])
    return ret


def shared_edges_variance(grns, var):
    out = []
    for v_ in var:
        sub = {}
        for x in grns:
            grn = grns[x]
            sub[x] = grn[grn['var.exp.median'] >= v_]

        out.append(len(shared_edges(sub)))

    return out


def fraction_edges_variance(grns, var):
    pass
