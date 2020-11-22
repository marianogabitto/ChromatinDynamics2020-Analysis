import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cairocffi
import scanpy as sc
from scipy.sparse import isspmatrix
from matplotlib.backends.backend_pdf import PdfPages


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


def extract_tf(data, TF, obs_key = None, obs_value = None,my_layer='Ms'):
    sub = -1
    if(obs_key is not None and obs_key in data.obs.keys()):
        data = data[data.obs[obs_key].isin([obs_value])]
        if(isspmatrix(data.layers[my_layer])):
            data.layers[my_layer] = data.layers[my_layer].todense()
        sub = pd.DataFrame(data.layers[my_layer][:,data.var.index.isin(TF[0])])
        sub.columns = data.var.index[data.var.index.isin(TF[0])]
        sub.index=data.obs.index
    return(sub)


def load_expression(paths):
    adatas = {k.split('_')[0]:sc.read_h5ad(paths[k]) for k in paths if '_exp' in k}
    adatas['P2'].obs['cardinal_class'] = adatas['P2'].obs['CellType']
    for tp in adatas:
        adatas[tp] = adatas[tp][:,adatas[tp].var.pct_dropout_by_counts != 100]
        adatas[tp] = adatas[tp][adatas[tp].obs.cardinal_class.isin(['PV', 'SST'])]
        sc.tl.rank_genes_groups(adatas[tp], groupby='cardinal_class',rankby_abs=True,
                                pts = True,use_raw=False, layers='normed',
                                method='wilcoxon',n_genes = adatas[tp].shape[1], 
                                key_added = 'DE')

    return adatas


def targets_by_var(grn, var=.05, tfs=None):
    tmpg = grn
    if tfs is not None:
        tmpg = grn[grn.regulator.isin(tfs.index.values)] 
    tmpg = tmpg[tmpg['var.exp.median']>=var]
    return(len(set(tmpg.target)))

def shared_targets(grns):
    edges = [grns[x].target for x in grns]
    ret = set(edges[0])
    for x in range(len(edges)-1,0,-1):
        ret = ret & set(edges[x]) & set(edges[x-1])
    return(ret)

def shared_targets_by_var(grns, by_var=.05):
    sub = {}
    for x in grns:
        grn = grns[x]
        sub[x] = grn[grn['var.exp.median'] >= by_var]
    return(shared_targets(sub))


def get_de_info(adatas, seltf):
    dedf = {}
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
    return dedf

def get_non_specific_tf(dedf, lfc=.25, pval=.05,fpv=.1,fsst=.1):
    return dedf[((abs(dedf.l2)<=lfc) | 
                ((abs(dedf.l2)>lfc) & 
                 (dedf.pvals>pval))) &
                 ((dedf.PV >=fpv) | 
                 (dedf.SST>=fsst))].names.values.tolist()

def compute_edge_fractions(f_grns):
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
        rpl[t] = pd.DataFrame(zip([out_reg_s[t][x]['PV']['prop_u'] for x in out_reg_s[t].keys()],
                                  [out_reg_s[t][x]['SST']['prop_u'] for x in out_reg_s[t]]))
        rpl[t].columns = ['PV','SST']
        rpl[t].index = out_reg_s[t].keys()
    
    return rpl