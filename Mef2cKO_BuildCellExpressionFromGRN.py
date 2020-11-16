import pandas as pd
import scipy.io as scio
import csv
import scipy.sparse as sp

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data
print("Building Cells from Average Expression Data")
TFA = {'PV': pd.read_table('Data/nov11P2-PVP-mtl/2020-11-13_01-22-52/nov11-p2_PV_TFA.tsv', index_col=0, header=0),
       'SST': pd.read_table('Data/nov11P2-PVP-mtl/2020-11-13_01-22-56/nov11-p2_SST_TFA.tsv', index_col=0, header=0)}

t = 0
betas = {'PV': {'boot_{}'.format(b): pd.read_table(
    'Data/nov11P2-PVP-mtl/nov11_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b), index_col=0,
    header=0) for b in range(1, 11)}}
t = 1
betas['SST'] = {'boot_{}'.format(b): pd.read_table(
    'Data/nov11P2-PVP-mtl/nov11_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b), index_col=0,
    header=0) for b in range(1, 11)}

TFA_s = TFA['SST'][TFA['SST'].columns[TFA['SST'].columns.isin(betas['SST']['boot_1'].columns)]]
sst_avg = (betas['SST']['boot_1'].dot(TFA_s.T) +
           betas['SST']['boot_2'].dot(TFA_s.T) +
           betas['SST']['boot_3'].dot(TFA_s.T) +
           betas['SST']['boot_4'].dot(TFA_s.T) +
           betas['SST']['boot_5'].dot(TFA_s.T) +
           betas['SST']['boot_6'].dot(TFA_s.T) +
           betas['SST']['boot_7'].dot(TFA_s.T) +
           betas['SST']['boot_8'].dot(TFA_s.T) +
           betas['SST']['boot_9'].dot(TFA_s.T) +
           betas['SST']['boot_10'].dot(TFA_s.T)) / 10
sst_cells_names = list(sst_avg.columns.values.tolist())
sst_cells_genes = list(sst_avg.index)

TFA_p = TFA['PV'][TFA['PV'].columns[TFA['PV'].columns.isin(betas['PV']['boot_1'].columns)]]
pv_avg = (betas['PV']['boot_1'].dot(TFA_p.T) +
          betas['PV']['boot_2'].dot(TFA_p.T) +
          betas['PV']['boot_3'].dot(TFA_p.T) +
          betas['PV']['boot_4'].dot(TFA_p.T) +
          betas['PV']['boot_5'].dot(TFA_p.T) +
          betas['PV']['boot_6'].dot(TFA_p.T) +
          betas['PV']['boot_7'].dot(TFA_p.T) +
          betas['PV']['boot_8'].dot(TFA_p.T) +
          betas['PV']['boot_9'].dot(TFA_p.T) +
          betas['PV']['boot_10'].dot(TFA_p.T)) / 10
pv_cells_names = list(pv_avg.columns.values.tolist())
pv_cells_genes = list(pv_avg.index)

scio.mmwrite("Out/pv_avg.mtx", sp.csr_matrix(pv_avg.to_numpy()))
scio.mmwrite("Out/sst_avg.mtx", sp.csr_matrix(sst_avg.to_numpy()))
with open('Out/sst_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_cells_genes)
with open('Out/sst_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_cells_names)
with open('Out/pv_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_cells_genes)
with open('Out/pv_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_cells_names)
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

kbetas = betas.copy()
for ct in kbetas:
    for boot in kbetas[ct]:
        kbetas[ct][boot].loc['Mef2c'] = 0
sst_k_avg = (kbetas['SST']['boot_1'].dot(TFA_sk.T) +
             kbetas['SST']['boot_2'].dot(TFA_sk.T) +
             kbetas['SST']['boot_3'].dot(TFA_sk.T) +
             kbetas['SST']['boot_4'].dot(TFA_sk.T) +
             kbetas['SST']['boot_5'].dot(TFA_sk.T) +
             kbetas['SST']['boot_6'].dot(TFA_sk.T) +
             kbetas['SST']['boot_7'].dot(TFA_sk.T) +
             kbetas['SST']['boot_8'].dot(TFA_sk.T) +
             kbetas['SST']['boot_9'].dot(TFA_sk.T) +
             kbetas['SST']['boot_10'].dot(TFA_sk.T)) / 10
sst_cells_names = list(sst_k_avg.columns.values.tolist())
sst_cells_genes = list(sst_k_avg.index)

pv_k_avg = (kbetas['PV']['boot_1'].dot(TFA_pk.T) +
            kbetas['PV']['boot_2'].dot(TFA_pk.T) +
            kbetas['PV']['boot_3'].dot(TFA_pk.T) +
            kbetas['PV']['boot_4'].dot(TFA_pk.T) +
            kbetas['PV']['boot_5'].dot(TFA_pk.T) +
            kbetas['PV']['boot_6'].dot(TFA_pk.T) +
            kbetas['PV']['boot_7'].dot(TFA_pk.T) +
            kbetas['PV']['boot_8'].dot(TFA_pk.T) +
            kbetas['PV']['boot_9'].dot(TFA_pk.T) +
            kbetas['PV']['boot_10'].dot(TFA_pk.T)) / 10
pv_cells_names = list(pv_k_avg.columns.values.tolist())
pv_cells_genes = list(pv_k_avg.index)

scio.mmwrite("Out/pv_avg_ko.mtx", sp.csr_matrix(pv_avg.to_numpy()))
scio.mmwrite("Out/sst_avg_ko.mtx", sp.csr_matrix(sst_avg.to_numpy()))
with open('Out/sst_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_cells_genes)
with open('Out/sst_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_cells_names)
with open('Out/pv_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_cells_genes)
with open('Out/pv_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_cells_names)
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data TRANSFAC
print("Building Cells from Average Expression Data TRANSFAC")
TFA_TRANSFAC = {'PV': pd.read_table('Data/nov12_P2-tfac-mtl/2020-11-13_01-51-36/nov12-p2tfac_PV_TFA.tsv', index_col=0,
                                    header=0),
                'SST': pd.read_table('Data/nov12_P2-tfac-mtl/2020-11-13_01-51-42/nov12-p2tfac_SST_TFA.tsv',
                                     index_col=0, header=0)}

t = 0
betas_transfac = {'PV': {'boot_{}'.format(b): pd.read_table(
    'Data/nov12_P2-tfac-mtl/nov12_itfac_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b),
    index_col=0, header=0) for b in range(1, 11)}}
t = 1
betas_transfac['SST'] = {'boot_{}'.format(b): pd.read_table(
    'Data/nov12_P2-tfac-mtl/nov12_itfac_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b),
    index_col=0, header=0) for b in range(1, 11)}

TFA_TRANSFAC_s = TFA_TRANSFAC['SST'][
    TFA_TRANSFAC['SST'].columns[TFA_TRANSFAC['SST'].columns.isin(betas_transfac['SST']['boot_1'].columns)]]
sst_transfac_avg = (betas_transfac['SST']['boot_1'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_2'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_3'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_4'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_5'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_6'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_7'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_8'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_9'].dot(TFA_TRANSFAC_s.T) +
                    betas_transfac['SST']['boot_10'].dot(TFA_TRANSFAC_s.T)) / 10
TFA_TRANSFAC_p = TFA_TRANSFAC['PV'][
    TFA_TRANSFAC['PV'].columns[TFA_TRANSFAC['PV'].columns.isin(betas_transfac['PV']['boot_1'].columns)]]
pv_transfac_avg = (betas_transfac['PV']['boot_1'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_2'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_3'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_4'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_5'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_6'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_7'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_8'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_9'].dot(TFA_TRANSFAC_p.T) +
                   betas_transfac['PV']['boot_10'].dot(TFA_TRANSFAC_p.T)) / 10

scio.mmwrite("Out/pv_transfac_avg.mtx", sp.csr_matrix(pv_transfac_avg.to_numpy()))
scio.mmwrite("Out/sst_transfac_avg.mtx", sp.csr_matrix(sst_transfac_avg.to_numpy()))
with open('Out/sst_transfac_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(sst_transfac_avg.index))
with open('Out/sst_transfac_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(sst_transfac_avg.columns.values.tolist()))
with open('Out/pv_transfac_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(pv_transfac_avg.index))
with open('Out/pv_transfac_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(pv_transfac_avg.columns.values.tolist()))
# ###################################################################################################
# ###################################################################################################

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data after Zeroing MEF2C Activity TRANSFAC
print("Building Cells from Average Expression Data after Zeroing MEF2C Activity TRANSFAC")
TFA_TRANSFAC_sk = TFA_TRANSFAC_s.copy()
TFA_TRANSFAC_pk = TFA_TRANSFAC_p.copy()
TFA_TRANSFAC_sk.loc[:, 'Mef2c'] = 0
TFA_TRANSFAC_pk.loc[:, 'Mef2c'] = 0

kbetas_transfac = betas_transfac.copy()
for ct in kbetas_transfac:
    for boot in kbetas_transfac[ct]:
        kbetas_transfac[ct][boot].loc['Mef2c'] = 0
sst_transfac_k_avg = (kbetas_transfac['SST']['boot_1'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_2'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_3'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_4'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_5'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_6'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_7'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_8'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_9'].dot(TFA_TRANSFAC_sk.T) +
                      kbetas_transfac['SST']['boot_10'].dot(TFA_TRANSFAC_sk.T)) / 10
pv_transfac_k_avg = (kbetas_transfac['PV']['boot_1'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_2'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_3'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_4'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_5'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_6'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_7'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_8'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_9'].dot(TFA_TRANSFAC_pk.T) +
                     kbetas_transfac['PV']['boot_10'].dot(TFA_TRANSFAC_pk.T)) / 10

scio.mmwrite("Out/pv_transfac_avg_ko.mtx", sp.csr_matrix(pv_transfac_k_avg.to_numpy()))
scio.mmwrite("Out/sst_transfac_avg_ko.mtx", sp.csr_matrix(sst_transfac_k_avg.to_numpy()))
with open('Out/sst_transfac_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(sst_transfac_avg.index))
with open('Out/sst_transfac_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(sst_transfac_avg.columns.values.tolist()))
with open('Out/pv_transfac_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(pv_transfac_avg.index))
with open('Out/pv_transfac_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(list(pv_transfac_avg.columns.values.tolist()))
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
