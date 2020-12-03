# ###################################################################################################
# ###################################################################################################
# Script to read Gene Regulatory Network information and Generate Panels in the
# figures.
# Joint Work: Giuseppe Saldi and Mariano Gabitto
# ###################################################################################################
# ###################################################################################################
from functions import *

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# READING inputs

nets = {'P28_PV':  'Data/P28PV.tsv',
        'P28_SST': 'Data/P28SST.tsv',
        'E18_PV':  'Data/E18PV.tsv',
        'E18_SST': 'Data/E18SST.tsv',
        'P2_PV':   'Data/16NOVP2PVnet.tsv',
        'P2_SST':  'Data/16NOVP2SSTnet.tsv'}
# RNA
paths = {'P28_exp': 'Data/P28_dwk.h5ad',
         'E18_exp': 'Data/E18_Lhx6.h5ad',
         'P2_exp':  'Data/nov15_P2_noftt.h5ad',
         'tf_list': 'Data/NOV6_aug_mus_tf.txt'}

grns = {x: read_and_filter_grn(nets[x], var=.0, net=x) for x in nets}
adatas = load_expression(paths)
fraction = .01

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
plt.bar(['E18_SST {}'.format(values[0]), 'E18_PV {}'.format(values[1]),
         'P2-E18_SST {}'.format(values[2]),'P2-E18_PV {}'.format(values[3]),
         'P2_SST {}'.format(values[4]), 'P2_PV {}'.format(values[5]),
         'P2-P28_SST {}'.format(values[6]), 'P2-P28_PV {}'.format(values[7]),
         'P28_SST {}'.format(values[8]),'P28_PV {}'.format(values[9])],
        values)
plt.ylabel("Number of Expressed Genes")
plt.tight_layout()
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL C: TOTAL NUMBER OF EDGES
v_range = np.arange(0.01, 0.41, 0.01)
tot = {}
for gr in grns:
    list_tot = []
    edges_vem = grns[gr]['var.exp.median'].to_numpy()
    for i_ in v_range:
        list_tot.append(np.sum(edges_vem > i_))
    tot[gr] = list_tot

pp = PdfPages("panelB_total_number_edges.pdf")
markers = ['o', 'o', '^', '^', 's', 's']
colors = ['#5A738B', '#B8767E', '#5A738B', '#B8767E', '#5A738B', '#B8767E']
for i_, gr in enumerate(tot.keys()):
    plt.scatter(v_range, tot[gr], alpha=0.6, label=gr, marker=markers[i_],
                c=colors[i_])
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
e18 = shared_edges_variance({'E18_SST': grns['E18_SST'],
                            'E18_PV': grns['E18_PV']}, variance_range)
p2 = shared_edges_variance({'P2_SST':  grns['P2_SST'],
                            'P2_PV': grns['P2_PV']}, variance_range)
p28 = shared_edges_variance({'P28_SST': grns['P28_SST'],
                             'P28_PV': grns['P28_PV']}, variance_range)

pp = PdfPages("panelC_shared_edges_pv_sst.pdf")
plt.scatter(variance_range, p28, alpha=0.6, marker='o', c='#999')
plt.scatter(variance_range, p2,  alpha=0.6, marker='s', c='#999')
plt.scatter(variance_range, e18, alpha=0.6, marker='^', c='#999')
plt.legend(loc='upper right', shadow=False, labels=['P28', 'P2', 'E18'])
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
plt.scatter(variance_range, (np.array(tot['E18_SST']) - np.array(e18)) /
            np.array(tot['E18_SST']),alpha=0.6, marker='^', color='#B8767E',
            label="E18 SST")
plt.scatter(variance_range, (np.array(tot['E18_PV']) - np.array(e18)) /
            np.array(tot['E18_PV']), alpha=0.6, marker='^', color='#5A738B',
            label="E18 PV")
plt.scatter(variance_range, (np.array(tot['P2_SST']) - np.array(p2)) /
            np.array(tot['P2_SST']), alpha=0.6, marker='s', color='#B8767E',
            label="P2 SST")
plt.scatter(variance_range, (np.array(tot['P2_PV']) - np.array(p2)) /
            np.array(tot['P2_PV']), alpha=0.6, marker='s', color='#5A738B',
            label="P2 PV")
plt.scatter(variance_range, (np.array(tot['P28_SST']) - np.array(p2)) /
            np.array(tot['P28_SST']), alpha=0.6, marker='o', color='#B8767E',
            label="P28 SST")
plt.scatter(variance_range, (np.array(tot['P28_PV']) - np.array(p2)) /
            np.array(tot['P28_PV']), alpha=0.6, marker='o', color='#5A738B',
            label="P28 PV")

plt.legend(loc='best', shadow=False)
plt.xlabel('variance explained')
plt.ylabel('Fraction of unique edges')
plt.axis([0, 0.4, 0, 1.0])
pp.savefig()
pp.close()

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL F: TARGETS AMONGST NETWORK
tbv = {}
variance_range = np.arange(0.01, 0.41, 0.01)

for net in grns:
    tbv[net] = [targets_by_var(grns[net], var=_v) for _v in variance_range]

pv_shared_tbv = [len(shared_targets_by_var({k: v for k, v in grns.items() if
                 'PV' in k}, _v))for _v in variance_range]
sst_shared_tbv = [len(shared_targets_by_var({k: v for k, v in grns.items() if
                  'SST' in k}, _v)) for _v in variance_range]
all_shared_tbv = [len(shared_targets_by_var(grns, _v)) for _v in
                  variance_range]

ese = {}
ese['PV'] = pd.DataFrame(zip(tbv['P28_PV'], tbv['P2_PV'], tbv['E18_PV'],
                         pv_shared_tbv, all_shared_tbv),
                         columns=['P28', 'P2', 'E18', 'pv_shared',
                         'all_shared'], index=variance_range)
ese['SST'] = pd.DataFrame(zip(tbv['P28_SST'], tbv['P2_SST'], tbv['E18_SST'],
                          pv_shared_tbv, all_shared_tbv),
                          columns=['P28', 'P2', 'E18', 'sst_shared',
                          'all_shared'],index=variance_range)

fig, axs = plt.subplots(ncols=2)

ax, ax1 = axs.flatten()
pla = ax.scatter(ese['PV'].index.values, ese['PV']['P28'], alpha=0.6,
                 marker='o', c='#5A738B')
plb = ax.scatter(ese['PV'].index.values, ese['PV']['P2'], alpha=0.6,
                 marker='s', c='#5A738B')
plc = ax.scatter(ese['PV'].index.values, ese['PV']['E18'], alpha=0.6,
                 marker='^', c='#5A738B')
plg = ax.scatter(ese['PV'].index.values, ese['PV']['pv_shared'], alpha=0.6,
                 marker='*', c='#5a8b77')
plh = ax.scatter(ese['PV'].index.values, ese['PV']['all_shared'], alpha=0.6,
                 marker='*', c='#658b5a')
ax.legend((pla, plb, plc, plg, plh), ('P28', 'P2', 'E18', 'PV shared',
                                      'all shared'),loc='upper right',
          shadow=False)
ax.set_title('PV')

pld = ax1.scatter(ese['SST'].index.values, ese['SST']['P28'], alpha=0.6,
                  marker='o', c='#B8767E')
ple = ax1.scatter(ese['SST'].index.values, ese['SST']['P2'], alpha=0.6,
                  marker='s', c='#B8767E')
plf = ax1.scatter(ese['SST'].index.values, ese['SST']['E18'], alpha=0.6,
                  marker='^', c='#B8767E')
pli = ax1.scatter(ese['SST'].index.values, ese['SST']['sst_shared'], alpha=0.6,
                  marker='*', c='#5a8b77')
plj = ax1.scatter(ese['SST'].index.values, ese['SST']['all_shared'], alpha=0.6,
                  marker='*', c='#658b5a')

ax1.legend((pld, ple, plf, pli, plj), ('P28', 'P2', 'E18', 'SST shared',
                                       'all shared'), loc='upper right',
           shadow=False)
ax1.set_title('SST')
plt.tight_layout()
# plt.ylabel('# edges')
# plt.xlabel('variance explained')
# plt.show()
plt.savefig('panelF_Targets_among_networks.pdf')

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL G: UNIQUE CELL TYPE EDGES FOR NON SPECIFIC TFS

dedf = {}
seltf = {}

var = .05
f_grns = {}
for grn in grns:
    f_grns[grn] = grns[grn][grns[grn]['var.exp.median'] >= var]

seltf['P28'] = set(set(f_grns['P28_PV'].regulator) |
                   set(f_grns['P28_SST'].regulator))
seltf['P2'] = set(set(f_grns['P2_PV'].regulator) |
                  set(f_grns['P2_SST'].regulator))
seltf['E18'] = set(set(f_grns['E18_PV'].regulator) |
                   set(f_grns['E18_SST'].regulator))


dedf = get_de_info(adatas, seltf)

seltf['P28'] = get_non_specific_tf(dedf['P28'])
seltf['P2']  = get_non_specific_tf(dedf['P2'])
seltf['E18'] = get_non_specific_tf(dedf['E18'])

rpl = compute_edge_fractions(f_grns)

fig, axs = plt.subplots(ncols=3)
ax, ax1, ax2 = axs.flatten()
t = 'E18'
ax.hist(rpl[t][rpl[t].index.isin(seltf['E18'])].T*100, alpha=0.6, density=True,
        label=['E18 PV', 'E18 SST'], bins=10, color=['#5A738B', '#B8767E'])
ax.set_title('{}'.format(t))
t = 'P2'
ax1.hist(rpl[t][rpl[t].index.isin(seltf['P2'])].T*100, alpha=0.6, density=True,
         label=['P2 PV', 'P2 SST'], bins=10, color=['#5A738B', '#B8767E'])
ax1.set_title('{}'.format(t))
t = 'P28'
ax2.hist(rpl[t][rpl[t].index.isin(seltf['P28'])].T*100, alpha=0.6, density=True,
         label=['P28 PV', 'P28 SST'], color=['#5A738B', '#B8767E'])
ax2.set_title('{}'.format(t))
plt.tight_layout()
plt.savefig('panelG_unique_celltype_edges_for_non_specific_TFS.pdf')


t = ['E18', 'P2', 'P28']
tmp = {x: np.mean(rpl[x][rpl[x].index.isin(seltf[x])]) for x in t}
tmp = (pd.DataFrame(tmp)*100).T
fig, ax = plt.subplots(1, 1)
title = 'Average proportion of unique edges for cell type non specific TFs'
ax.set_title(title)
plt.box(False)
plt.tick_params(left=False, labelleft=False, bottom=False,labelbottom=False)
table(ax,tmp,loc='center', colWidths=[0.2, 0.2, 0.2])
plt.subplots_adjust(wspace=0.2, hspace=0.2)
plt.savefig('panelG_stats.pdf')

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# networks files for rendering
var = .05
f_grns = {}
cols = ['target', 'regulator', 'var.exp.median']
for grn in grns:
    f_grns[grn] = grns[grn][grns[grn]['var.exp.median'] >= var]
    tmp = f_grns[grn][cols]
    file_n = 'ipynb/{}.tsv'.format(grn)
    tmp[cols].to_csv(file_n, sep='\t', index=None)


e18s = f_grns['E18_PV'][f_grns['E18_PV'].id.isin(f_grns['E18_SST'].id)][cols]
e18s.to_csv('ipynb/{}.tsv'.format('E18_shared'), sep='\t', index=None)

p2s = f_grns['P2_PV'][f_grns['P2_PV'].id.isin(f_grns['P2_SST'].id)][cols]
p2s.to_csv('ipynb/{}.tsv'.format('P2_shared'), sep='\t', index=None)

p28s = f_grns['P28_PV'][f_grns['P28_PV'].id.isin(f_grns['P28_SST'].id)][cols]
p28s.to_csv('ipynb/{}.tsv'.format('P28_shared'), sep='\t', index=None)
