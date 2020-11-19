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
from scipy.sparse import isspmatrix

# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# READING inputs

nets = {'P28_PV' :'Data/P28PV.tsv',
        'P28_SST':'Data/P28SST.tsv',
        'E18_PV' :'Data/E18PV.tsv',
        'E18_SST':'Data/E18SST.tsv',
        'P2_PV'  :'Data/16NOVP2PVnet.tsv',
        'P2_SST' :'Data/16NOVP2SSTnet.tsv'}

grns = {x:read_and_filter_grn(nets[x],var=.0,net=x) for x in nets}

# RNA
paths = {'P28_exp':'Data/P28_dwk.h5ad',
         'E18_exp':'Data/E18_Lhx6.h5ad',
         'P2_exp' :'Data/nov15_P2_noftt.h5ad',
         'tf_list':'Data/NOV6_aug_mus_tf.txt'}


adatas = load_expression(paths)
fraction = .01

# TF = pd.read_table(paths['tf_list'], header=None)
# TF_per_TP = {k:{} for k in adatas}
# for k in TF_per_TP:
#     for cardinal in ['PV','SST']:
#         TF_per_TP[k][cardinal] = extract_tf(adatas[k],TF, obs_key = 'cardinal_class',obs_value = cardinal, my_layer='normed')

# ### Sort TF on decreasing TMM expression
# tmms = {cl:{} for cl in TF_per_TP}
# bps = {}
# for ct in tmms:
#     for cl in TF_per_TP[ct]:
#         tmms[ct][cl] = pd.DataFrame(stats.trim_mean(abs(TF_per_TP[ct][cl]),  0.025, axis=0),index=TF_per_TP[ct][cl].columns.values)
#     bps[ct] = pd.concat([tmms[ct][x] for x in tmms[ct]], axis=1) 
#     bps[ct].columns = list(tmms[ct].keys())

# sel_tfs = {}
# for ct in bps:
#     sel_tfs[ct] = bps[ct][abs(bps[ct]).sum(axis=1)>0]

# del tmms
# del bps
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
# PANEL F: TARGETS AMONGST NETWORK


# ############################################################################################################
# ############################################################################################################
# ############################################################################################################
# PANEL G: UNIQUE CELL TYPE EDGES FOR NON SPECIFIC TFS

dedf ={}
seltf = {}

var = .05
f_grns = {}
for grn in grns:
    f_grns[grn] = grns[grn][grns[grn]['var.exp.median'] >= var]

seltf['P28'] = set(set(f_grns['P28_PV'].regulator) | set(f_grns['P28_SST'].regulator))
seltf['P2']  = set(set(f_grns['P2_PV' ].regulator) | set(f_grns['P2_SST' ].regulator))
seltf['E18'] = set(set(f_grns['E18_PV'].regulator) | set(f_grns['E18_SST'].regulator))

for tp in adatas:
    dedf[tp] = sc.get.rank_genes_groups_df(adatas[tp], group="PV",key='DE')
    gene_ids = adatas[tp].var.index.values
    if(isspmatrix(adatas[tp].layers['normed'])):
        obs = adatas[tp].layers['normed'].todense()
    else: obs = adatas[tp].layers['normed'].copy()
    obs = pd.DataFrame(obs,columns=gene_ids,index=adatas[tp].obs['cardinal_class'])
    obs_bool = obs.astype(bool)
    fraction_obs = obs_bool.groupby(level=0).sum()/obs_bool.groupby(level=0).count()
    mean_obs = obs.groupby(level=0).mean()
    mean_obs.loc['SST'][mean_obs.loc['SST'] == 0] = 1
    mean_obs = pd.DataFrame(np.log2(mean_obs.loc['PV']/mean_obs.loc['SST']))
    mean_obs.columns=['l2']
    dedf[tp] = dedf[tp].merge(fraction_obs.T, left_on='names', right_index=True)
    dedf[tp] = dedf[tp].merge(mean_obs, left_on='names', right_index=True)
    dedf[tp] = dedf[tp][dedf[tp].names.isin(seltf[tp])]





seltf['P28'] = dedf['P28'][((abs(dedf['P28'].l2)<=.25) | 
                            ((abs(dedf['P28'].l2)>.25) & 
                             (dedf['P28'].pvals>.05))) &
                             ((dedf['P28'].PV >=.1) | 
                             (dedf['P28'].SST>=.1))].names.values.tolist()

seltf['P2']  = dedf['P2'][((abs(dedf['P2'].l2)<=.25) |
                          ((abs(dedf['P2'].l2)>.25)  &
                            (dedf['P2'].pvals>.05))) &
                            ((dedf['P2'].PV >=.1) |
                             (dedf['P2'].SST>=.1))].names.values.tolist()

seltf['E18'] = dedf['E18'][((abs(dedf['E18'].l2)<=.25) |
                           ((abs(dedf['E18'].l2)>.25)  &
                            (dedf['E18'].pvals>.05)))  &
                            ((dedf['E18'].PV >=.1) |
                            (dedf['E18'].SST>=.1))].names.values.tolist()

out_reg_s = {}
tp = set([x.split('_')[0] for x in f_grns])
for t in tp:
    report = {}
    network_a, network_b = [x for x in f_grns.keys() if t == x.split('_')[0]]
    networks = [network_a, network_b]
    cell_type_a, cell_type_b = [x.split('_')[1] for x in networks]
    f_net_a  = f_grns[network_a]
    f_net_b  = f_grns[network_b]
    for regulator in set(set(f_net_a.regulator) | set(f_net_b.regulator)):
        n=0
        edge_id = []
        if regulator in f_net_a.regulator.values:
            n = f_net_a.regulator.value_counts()[regulator]
            edge_id = list(f_net_a[f_net_a.regulator==regulator].id.values)
        report[regulator] = {cell_type_a:{'n':n,'edge_id':edge_id}}
        n=0
        edge_id = []
        if regulator in f_net_b.regulator.values:
            n = f_net_b.regulator.value_counts()[regulator]
            edge_id = list(f_net_b[f_net_b.regulator==regulator].id.values)

        report[regulator][cell_type_b] = {'n': n , 'edge_id':edge_id}
        
        s_edges = set(report[regulator][cell_type_b]['edge_id']) & set(report[regulator][cell_type_a]['edge_id'])
        report[regulator]['shared_id'] = len(s_edges)
        divider = report[regulator][cell_type_a]['n']
        if divider == 0: divider = 1
        report[regulator][cell_type_a]['prop_s'] = report[regulator]['shared_id'] / divider
        report[regulator][cell_type_a]['prop_u'] = None
        if report[regulator][cell_type_a]['n'] > 0:
            report[regulator][cell_type_a]['prop_u'] = abs(1 - report[regulator][cell_type_a]['prop_s'])
        divider = report[regulator][cell_type_b]['n'] 
        if divider == 0: divider = 1
        report[regulator][cell_type_b]['prop_s'] = report[regulator]['shared_id'] / divider
        report[regulator][cell_type_b]['prop_u'] = None
        if report[regulator][cell_type_b]['n'] > 0:
            report[regulator][cell_type_b]['prop_u'] = abs(1 - report[regulator][cell_type_b]['prop_s'])
    out_reg_s[t] = report

rpl = {}
for t in tp:
    rpl[t] = pd.DataFrame(zip([out_reg_s[t][x]['PV']['prop_u'] for x in out_reg_s[t].keys()],[out_reg_s[t][x]['SST']['prop_u'] for x in out_reg_s[t]]))
    rpl[t].columns = ['PV','SST']
    rpl[t].index = out_reg_s[t].keys()

fig, axs = plt.subplots(ncols=3)
ax, ax1,ax2 = axs.flatten()
t = 'E18'
ax.hist(rpl[t][rpl[t].index.isin(seltf['E18'])]*100, alpha=0.6,density=True,label=['E18 PV','E18 SST'],bins=10, color = ['#5A738B','#B8767E'])
ax.set_title('{}'.format(t))
t = 'P2'
ax1.hist(rpl[t][rpl[t].index.isin(seltf['P2'])]*100, alpha=0.6,density=True,label=['P2 PV','P2 SST'],bins=10, color = ['#5A738B','#B8767E'])
ax1.set_title('{}'.format(t))
t = 'P28'
ax2.hist(rpl[t][rpl[t].index.isin(seltf['P28'])]*100, alpha=0.6,density=True,label=['P28 PV','P28 SST'], color = ['#5A738B','#B8767E'])
ax2.set_title('{}'.format(t))
plt.tight_layout()
plt.savefig('panelG_unique_celltype_edges_for_non_specific_TFS.pdf')
