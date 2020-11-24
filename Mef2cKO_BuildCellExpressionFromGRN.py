# ###################################################################################################
# ###################################################################################################
# Script to read Gene Regulatory Network information and create simulated cells.
# Joint Work: Giuseppe Saldi and Mariano Gabitto
# ###################################################################################################
# ###################################################################################################
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
TFA_p = TFA['PV'][TFA['PV'].columns[TFA['PV'].columns.isin(betas['PV']['boot_1'].columns)]]

acc = betas['SST']['boot_1']
for i in range(2, 11):
    acc += betas['SST']['boot_{}'.format(i)]
acc = acc.dot(TFA_s.T) / 10
sst_avg = acc 

acc = betas['PV']['boot_1']
for i in range(2, 11):
    acc += betas['PV']['boot_{}'.format(i)]
acc = acc.dot(TFA_p.T) / 10
pv_avg = acc 

sst_cells_names = list(sst_avg.columns.values.tolist())
sst_cells_genes = list(sst_avg.index)
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

acc = betas['SST']['boot_1']
for i in range(2, 11):
    acc += betas['SST']['boot_{}'.format(i)]
acc = acc.dot(TFA_sk.T) / 10
sst_k_avg = acc 

acc = betas['PV']['boot_1']
for i in range(2, 11):
    acc += betas['PV']['boot_{}'.format(i)]
acc = acc.dot(TFA_pk.T) / 10
pv_k_avg = acc 

pv_cells_names = list(pv_k_avg.columns.values.tolist())
pv_cells_genes = list(pv_k_avg.index)

scio.mmwrite("Out/pv_avg_ko.mtx", sp.csr_matrix(pv_k_avg.to_numpy()))
scio.mmwrite("Out/sst_avg_ko.mtx", sp.csr_matrix(sst_k_avg.to_numpy()))
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
# Building Cells from Average Expression Data after Zeroing MEF2C Activity
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

# ###################################################################################################
# ###################################################################################################
# Building Cells from Average Expression Data
# NO FTT NORMALIZATION, 1LOGP
print("Building Cells from Average Expression Data. No FTT Norm.")
TFA_1lp = {'PV': pd.read_table('Data/nov16P2-PVP-noftt-mtl/2020-11-16_06-44-41/nov15-p2_PV_TFA.tsv', index_col=0,
                               header=0),
           'SST': pd.read_table('Data/nov16P2-PVP-noftt-mtl/2020-11-16_06-44-44/nov15-p2_SST_TFA.tsv', index_col=0,
                                header=0)}

t = 0
betas_1lp = {'PV': {'boot_{}'.format(b): pd.read_table(
    'Data/nov16P2-PVP-noftt-mtl/nov15-noftt_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b), index_col=0,
    header=0) for b in range(1, 11)}}
t = 1
betas_1lp['SST'] = {'boot_{}'.format(b): pd.read_table(
    'Data/nov16P2-PVP-noftt-mtl/nov15-noftt_betas_task_{t}_boot_{b}.tsv'.format(t=t, b=b), index_col=0,
    header=0) for b in range(1, 11)}

TFA_1lp_s = TFA_1lp['SST'][TFA_1lp['SST'].columns[TFA_1lp['SST'].columns.isin(betas_1lp['SST']['boot_1'].columns)]]
TFA_1lp_p = TFA_1lp['PV'][TFA_1lp['PV'].columns[TFA_1lp['PV'].columns.isin(betas_1lp['PV']['boot_1'].columns)]]

acc = betas_1lp['SST']['boot_1']
for i in range(2, 11):
    acc += betas_1lp['SST']['boot_{}'.format(i)]
acc = acc.dot(TFA_1lp_s.T) / 10
sst_1lp_avg = acc

acc = betas_1lp['PV']['boot_1']
for i in range(2, 11):
    acc += betas_1lp['PV']['boot_{}'.format(i)]
acc = acc.dot(TFA_1lp_p.T) / 10
pv_1lp_avg = acc

sst_1lp_cells_names = list(sst_1lp_avg.columns.values.tolist())
sst_1lp_cells_genes = list(sst_1lp_avg.index)
pv_1lp_cells_names = list(pv_1lp_avg.columns.values.tolist())
pv_1lp_cells_genes = list(pv_1lp_avg.index)

scio.mmwrite("Out/pv_1lp_avg.mtx", sp.csr_matrix(pv_1lp_avg.to_numpy()))
scio.mmwrite("Out/sst_1lp_avg.mtx", sp.csr_matrix(sst_1lp_avg.to_numpy()))
with open('Out/sst_1lp_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_1lp_cells_genes)
with open('Out/sst_1lp_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_1lp_cells_names)
with open('Out/pv_1lp_avg_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_1lp_cells_genes)
with open('Out/pv_1lp_avg_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_1lp_cells_names)
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

acc = betas_1lp['SST']['boot_1']
for i in range(2, 11):
    acc += betas_1lp['SST']['boot_{}'.format(i)]
acc = acc.dot(TFA_1lp_sk.T) / 10
sst_1lp_k_avg = acc

acc = betas_1lp['PV']['boot_1']
for i in range(2, 11):
    acc += betas_1lp['PV']['boot_{}'.format(i)]
acc = acc.dot(TFA_1lp_pk.T) / 10
pv_1lp_k_avg = acc

pv_1lp_cells_names = list(pv_1lp_k_avg.columns.values.tolist())
pv_1lp_cells_genes = list(pv_1lp_k_avg.index)

scio.mmwrite("Out/pv_1lp_avg_ko.mtx", sp.csr_matrix(pv_1lp_k_avg.to_numpy()))
scio.mmwrite("Out/sst_1lp_avg_ko.mtx", sp.csr_matrix(sst_1lp_k_avg.to_numpy()))
with open('Out/sst_1lp_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_1lp_cells_genes)
with open('Out/sst_1lp_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(sst_1lp_cells_names)
with open('Out/pv_1lp_avg_ko_genes.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_1lp_cells_genes)
with open('Out/pv_1lp_avg_ko_cells.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter='\t')
    spamwriter.writerow(pv_1lp_cells_names)
# ###################################################################################################
# ###################################################################################################
# ###################################################################################################
