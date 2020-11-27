# ###################################################################################################
# ###################################################################################################
# Script to read Gene Regulatory Network information and create simulated cells.
# Joint Work: Giuseppe Saldi and Mariano Gabitto
# ###################################################################################################
# ###################################################################################################
from functions import *
import glob


# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data
print("Building Cells from Average Expression Data")
TFA = {'PV': pd.read_table('Data/nov11P2-PVP-mtl/nov11-p2_PV_TFA.tsv', index_col=0, header=0),
       'SST': pd.read_table('Data/nov11P2-PVP-mtl/nov11-p2_SST_TFA.tsv', index_col=0, header=0)}

betas = {}
path = 'Data/nov11P2-PVP-mtl/nov11_betas_task_{t}_boot_{b}.tsv'
betas['PV'] = read_betas(0, 'PV', path, n_boot=10)
betas['SST'] = read_betas(1, 'SST', path, n_boot=10)

TFA_s = TFA['SST'][TFA['SST'].columns[TFA['SST'].columns.isin(betas['SST']['boot_1'].columns)]]
TFA_p = TFA['PV'][TFA['PV'].columns[TFA['PV'].columns.isin(betas['PV']['boot_1'].columns)]]

sst_avg = B_dot_A(betas['SST'], TFA_s, 10)
pv_avg = B_dot_A(betas['PV'], TFA_p, 10)

sst_cells = list(sst_avg.columns.values)
sst_genes = list(sst_avg.index)
paths = {'mtx': 'Out/sst_avg.mtx', 'genes': 'Out/sst_avg_genes.tsv', 'cells': 'Out/sst_avg_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_avg)

pv_cells = list(pv_avg.columns.values)
pv_genes = list(pv_avg.index)
paths = {'mtx': 'Out/pv_avg.mtx', 'genes': 'Out/pv_avg_genes.tsv', 'cells': 'Out/pv_avg_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_avg)
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data after Zeroing MEF2C Activity
print("Building Cells from Average Expression Data after Zeroing MEF2C Activity")

TFA_sk = TFA_s.copy()
TFA_pk = TFA_p.copy()
TFA_sk.loc[:, 'Mef2c'] = 0
TFA_pk.loc[:, 'Mef2c'] = 0

sst_k_avg = B_dot_A(betas['SST'], TFA_sk, 10)
pv_k_avg = B_dot_A(betas['PV'], TFA_pk, 10)

paths = {'mtx': 'Out/sst_avg_ko.mtx', 'genes': 'Out/sst_avg_ko_genes.tsv', 'cells': 'Out/sst_avg_ko_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_k_avg)

paths = {'mtx': 'Out/pv_avg_ko.mtx', 'genes': 'Out/pv_avg_ko_genes.tsv', 'cells': 'Out/pv_avg_ko_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_k_avg)
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data TRANSFAC
print("Building Cells from Average Expression Data TRANSFAC")
TFA_TRANSFAC = {'PV': pd.read_table('Data/nov12_P2-tfac-mtl/nov12-p2tfac_PV_TFA.tsv', index_col=0, header=0),
                'SST': pd.read_table('Data/nov12_P2-tfac-mtl/nov12-p2tfac_SST_TFA.tsv', index_col=0, header=0)}

betas_transfac = {}
path = 'Data/nov12_P2-tfac-mtl/nov12_itfac_betas_task_{t}_boot_{b}.tsv'
betas_transfac['PV'] = read_betas(0, 'PV', path, n_boot=10)
betas_transfac['SST'] = read_betas(1, 'SST', path, n_boot=10)

TFA_TRANSFAC_s = TFA_TRANSFAC['SST'][TFA_TRANSFAC['SST'].columns[TFA_TRANSFAC['SST'].columns.isin(
    betas_transfac['SST']['boot_1'].columns)]]
TFA_TRANSFAC_p = TFA_TRANSFAC['PV'][TFA_TRANSFAC['PV'].columns[TFA_TRANSFAC['PV'].columns.isin(
    betas_transfac['PV']['boot_1'].columns)]]

sst_transfac_avg = B_dot_A(betas_transfac['SST'], TFA_TRANSFAC_s, 10)
pv_transfac_avg = B_dot_A(betas_transfac['PV'], TFA_TRANSFAC_p, 10)

sst_cells = list(sst_transfac_avg.columns.values)
sst_genes = list(sst_transfac_avg.index)
paths = {'mtx': 'Out/sst_transfac_avg.mtx', 'genes': 'Out/sst_transfac_avg_genes.tsv', 'cells':
         'Out/sst_transfac_avg_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_transfac_avg)

pv_cells = list(pv_transfac_avg.columns.values)
pv_genes = list(pv_transfac_avg.index)
paths = {'mtx': 'Out/pv_transfac_avg.mtx', 'genes': 'Out/pv_transfac_avg_genes.tsv',
         'cells': 'Out/pv_transfac_avg_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_transfac_avg)
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data after Zeroing MEF2C Activity
print("Building Cells from Average Expression Data after Zeroing MEF2C Activity TRANSFAC")
TFA_TRANSFAC_sk = TFA_TRANSFAC_s.copy()
TFA_TRANSFAC_pk = TFA_TRANSFAC_p.copy()
TFA_TRANSFAC_sk.loc[:, 'Mef2c'] = 0
TFA_TRANSFAC_pk.loc[:, 'Mef2c'] = 0

kbetas_transfac = betas_transfac.copy()
sst_transfac_k_avg = B_dot_A(kbetas_transfac['SST'], TFA_TRANSFAC_sk, 10)
pv_transfac_k_avg = B_dot_A(kbetas_transfac['PV'], TFA_TRANSFAC_pk, 10)

paths = {'mtx': 'Out/sst_transfac_avg_ko.mtx', 'genes': 'Out/sst_transfac_avg_ko_genes.tsv', 'cells':
         'Out/sst_transfac_avg_ko_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_transfac_k_avg)

paths = {'mtx': 'Out/pv_transfac_avg_ko.mtx', 'genes': 'Out/pv_transfac_avg_ko_genes.tsv', 'cells':
         'Out/pv_transfac_avg_ko_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_transfac_k_avg)
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data
# NO FTT NORMALIZATION, 1LOGP
print("Building Cells from Average Expression Data. No FTT Norm.")
TFA_1lp = {'PV': pd.read_table('Data/nov16P2-PVP-noftt-mtl/nov15-p2_PV_TFA.tsv', index_col=0, header=0),
           'SST': pd.read_table('Data/nov16P2-PVP-noftt-mtl/nov15-p2_SST_TFA.tsv', index_col=0, header=0)}

betas_1lp = {}
path = 'Data/nov16P2-PVP-noftt-mtl/nov15-noftt_betas_task_{t}_boot_{b}.tsv'
betas_1lp['PV'] = read_betas(0, 'PV', path, n_boot=10)
betas_1lp['SST'] = read_betas(1, 'SST', path, n_boot=10)

TFA_1lp_s = TFA_1lp['SST'][TFA_1lp['SST'].columns[TFA_1lp['SST'].columns.isin(betas_1lp['SST']['boot_1'].columns)]]
TFA_1lp_p = TFA_1lp['PV'][TFA_1lp['PV'].columns[TFA_1lp['PV'].columns.isin(betas_1lp['PV']['boot_1'].columns)]]

sst_1lp_avg = B_dot_A(betas_1lp['SST'], TFA_1lp_s, 10)
pv_1lp_avg = B_dot_A(betas_1lp['PV'], TFA_1lp_p, 10)

sst_cells = list(sst_1lp_avg.columns.values)
sst_genes = list(sst_1lp_avg.index)
paths = {'mtx': 'Out/sst_1lp_avg.mtx', 'genes': 'Out/sst_1lp_avg_genes.tsv', 'cells': 'Out/sst_1lp_avg_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_1lp_avg)

pv_cells = list(pv_1lp_avg.columns.values)
pv_genes = list(pv_1lp_avg.index)
paths = {'mtx': 'Out/pv_1lp_avg.mtx', 'genes': 'Out/pv_1lp_avg_genes.tsv', 'cells': 'Out/pv_1lp_avg_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_1lp_avg)
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data after Zeroing MEF2C Activity
# NO FTT, 1logp
print("Building Cells from Average Expression Data after Zeroing MEF2C Activity. NoFTT Norm")

TFA_1lp_sk = TFA_1lp_s.copy()
TFA_1lp_pk = TFA_1lp_p.copy()
TFA_1lp_sk.loc[:, 'Mef2c'] = 0
TFA_1lp_pk.loc[:, 'Mef2c'] = 0

sst_1lp_k_avg = B_dot_A(betas_1lp['SST'], TFA_1lp_sk, 10)
pv_1lp_k_avg = B_dot_A(betas_1lp['PV'], TFA_1lp_pk, 10)

sst_cells = list(sst_1lp_k_avg.columns.values)
sst_genes = list(sst_1lp_k_avg.index)
paths = {'mtx': 'Out/sst_1lp_avg.mtx', 'genes': 'Out/sst_1lp_avg_genes.tsv', 'cells': 'Out/sst_1lp_avg_cells.tsv'}
write_reconstructed_exp(sst_cells, sst_genes, paths, sst_1lp_k_avg)

pv_cells = list(pv_1lp_k_avg.columns.values)
pv_genes = list(pv_1lp_k_avg.index)
paths = {'mtx': 'Out/pv_1lp_avg.mtx', 'genes': 'Out/pv_1lp_avg_genes.tsv', 'cells': 'Out/pv_1lp_avg_cells.tsv'}
write_reconstructed_exp(pv_cells, pv_genes, paths, pv_1lp_k_avg)
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Turn mtx into andata
mtx = {}
for x in glob.glob("Out/*.mtx"):
    f = x.split('/')[-1].split('.')[0]
    mtx[f] = pd.DataFrame(scio.mmread(x).todense())
    mtx[f].columns = pd.read_table(x.replace('.mtx', '_cells.tsv')).columns
    mtx[f].index = pd.read_table(x.replace('.mtx', '_genes.tsv')).columns

adatas = {x: sc.AnnData(mtx[x].T) for x in mtx}
for x in adatas:
    adatas[x].obs['cardinal'] = x.split('_')[0].upper()

rec = {}
batches = ['P2_rWT_PV', 'P2_rWT_SST']
rec['P2_rWT'] = adatas['pv_avg'].concatenate(adatas['sst_avg'], batch_categories=batches)
rec['P2_rKO'] = adatas['pv_avg_ko'].concatenate(adatas['sst_avg_ko'], batch_categories=['P2_rKO_PV', 'P2_rKO_SST'])
rec['P2_rTfacKO'] = adatas['pv_transfac_avg_ko'].concatenate(adatas['sst_transfac_avg_ko'],
                                                             batch_categories=['P2_rTfacKO_PV', 'P2_rTfacKO_SST'])
rec['P2_rTfacWT'] = adatas['pv_transfac_avg'].concatenate(adatas['sst_transfac_avg'],
                                                          batch_categories=['P2_rTfacWT_PV', 'P2_rTfacWT_SST'])
rec['P2_rT1lpWT'] = adatas['pv_1lp_avg'].concatenate(adatas['sst_1lp_avg'],
                                                     batch_categories=['P2_r1lpWT_PV', 'P2_rT1lpWT_SST'])
