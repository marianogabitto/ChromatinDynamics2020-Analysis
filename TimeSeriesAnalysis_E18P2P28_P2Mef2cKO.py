# ###################################################################################################
# ###################################################################################################
# Script to read Gene Regulatory Network information and Generate Panels in the figures.
# Joint Work: Giuseppe Saldi and Mariano Gabitto
# ###################################################################################################
# ###################################################################################################
import functions as fnc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# READING NETS
nets = {'P28_PV': 'Data/P28/P28PV.tsv', 'P28_SST': 'Data/P28/P28SST.tsv', 'E18_PV': 'Data/E18/E18PV.tsv',
        'E18_SST': 'Data/E18/E18SST.tsv', 'P2_PV': 'Data/nov11P2-PVP-mtl/PV/network.tsv',
        'P2_SST': 'Data/nov11P2-PVP-mtl/SST/network.tsv'}
grns = {x: fnc.read_and_filter_grn(nets[x], var=.0, net=x) for x in nets}

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL B: NUMBER OF EXPRESSED GENES
values = [len(set(grns['E18_SST'].target)), len(set(grns['E18_PV'].target)),
          len(set(grns['E18_SST'].target) & set(grns['P2_SST'].target)),
          len(set(grns['E18_PV'].target) & set(grns['P2_PV'].target)),
          len(set(grns['P2_SST'].target)), len(set(grns['P2_PV'].target)),
          len(set(grns['P2_SST'].target) & set(grns['P28_SST'].target)),
          len(set(grns['P2_PV'].target) & set(grns['P28_PV'].target)),
          len(set(grns['P28_SST'].target)), len(set(grns['P28_PV'].target))]

pp = PdfPages("panelA_number_expressed_genes.pdf")
plt.xticks(rotation=45)
plt.bar(['E18_SST {}'.format(values[0]), 'E18_PV {}'.format(values[1]), 'P2-E18_SST {}'.format(values[2]),
         'P2-E18_PV {}'.format(values[3]), 'P2_SST {}'.format(values[4]), 'P2_PV {}'.format(values[5]),
         'P2-P28_SST {}'.format(values[6]), 'P2-P28_PV {}'.format(values[7]), 'P28_SST {}'.format(values[8]),
         'P28_PV {}'.format(values[9])],
        values)
plt.ylabel("Number of Expressed Genes")
plt.tight_layout()
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL C: TOTAL NUMBER OF EDGES
tot = {}
for gr in grns.keys():
    list_tot = []
    edges_vem = grns[gr]['var.exp.median'].to_numpy()
    for i_ in np.arange(0.01, 0.41, 0.01):
        list_tot.append(np.sum(edges_vem > i_))
    tot[gr] = list_tot

pp = PdfPages("panelB_total_number_edges.pdf")
markers = ['o', 'o', '^', '^', 's', 's']
for i_, gr in enumerate(tot.keys()):
    plt.scatter(np.arange(0.01, 0.4, 0.01), tot[gr], alpha=0.6, label=gr, marker=markers[i_])
plt.legend(loc='upper right', shadow=False)
plt.xlabel('variance explained')
plt.ylabel('# edges')
plt.tight_layout()
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL D: SHARD EDGES PV AND SST AT EACH TIME POINT
variance_range = np.arange(0.01, 0.41, 0.01)
e18 = fnc.shared_edges_variance({'E18_SST': grns['E18_SST'], 'E18_PV': grns['E18_PV']}, variance_range)
p2 = fnc.shared_edges_variance({'P2_SST': grns['P2_SST'], 'P2_PV': grns['P2_PV']}, variance_range)
p28 = fnc.shared_edges_variance({'P28_SST': grns['P28_SST'], 'P28_PV': grns['P28_PV']}, variance_range)

pp = PdfPages("panelC_shared_edges_pv_sst.pdf")
plt.scatter(variance_range, p28, alpha=0.6, marker='o')
plt.scatter(variance_range, p2, alpha=0.6, marker='s')
plt.scatter(variance_range, e18, alpha=0.6, marker='^')
plt.legend(loc='upper right', shadow=False)
plt.xlabel('variance explained')
plt.ylabel('# shared edges')
plt.tight_layout()
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL E: FRACTION OF UNIQUE EDGES PV AND SST AT EACH TIME POINT
pp = PdfPages("panelD_fraction_unique_edges_pv_sst.pdf")
plt.scatter(variance_range, (np.array(tot['E18_SST']) - np.array(e18)) / np.array(tot['E18_SST']),
            alpha=0.6, marker='^', color='b', label="E18 SST")
plt.scatter(variance_range, (np.array(tot['E18_PV']) - np.array(e18)) / np.array(tot['E18_PV']),
            alpha=0.6, marker='^', color='r', label="E18 PV")
plt.scatter(variance_range, (np.array(tot['P2_SST']) - np.array(p2)) / np.array(tot['P2_SST']),
            alpha=0.6, marker='s', color='b', label="P2 SST")
plt.scatter(variance_range, (np.array(tot['P2_PV']) - np.array(p2)) / np.array(tot['P2_PV']),
            alpha=0.6, marker='s', color='r', label="P2 PV")
plt.scatter(variance_range, (np.array(tot['P28_SST']) - np.array(p2)) / np.array(tot['P28_SST']),
            alpha=0.6, marker='o', color='b', label="P28 SST")
plt.scatter(variance_range, (np.array(tot['P28_PV']) - np.array(p2)) / np.array(tot['P28_PV']),
            alpha=0.6, marker='o', color='r', label="P28 PV")
plt.legend(loc='upper right', shadow=False)
plt.xlabel('variance explained')
plt.ylabel('Fraction of unique edges')
plt.axis([0, 0.4, 0, 1.0])
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL F: NUMBER OF EXPRESSED GENES

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL G: NUMBER OF EXPRESSED GENES
