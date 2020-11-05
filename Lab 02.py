from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import scipy.signal as sig
import pandas as pd
import numpy as np

B = (10, 20, 25, 30)
X = np.zeros(4)
Y = np.zeros(4)
R = 8.1345

# DSC Example
df_a = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Example - 5K min.csv')
df_b = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Example - 10K min.csv')
df_c = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Example - 15K min.csv')
df_d = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Example - 20K min.csv')

fig_1,ax_1 = plt.subplots(4, sharex=True, sharey=True, figsize=(5, 6), dpi=300)
ax_1[0].plot(df_a['Temperature'], df_a['Heat Flow'], color='C0', linewidth=1, label='5 K/min')
ax_1[1].plot(df_b['Temperature'], df_b['Heat Flow'], color='C1', linewidth=1, label='10 K/min')
ax_1[2].plot(df_c['Temperature'], df_c['Heat Flow'], color='C2', linewidth=1, label='15 K/min')
ax_1[3].plot(df_d['Temperature'], df_d['Heat Flow'], color='C3', linewidth=1, label='20 K/min')
    
fig_1.suptitle('DSC Scan Example Data')
ax_1[3].set(xlabel=r'Temperature ($\degree$C)')
ax_1[1].set(ylabel='Heat Flow (mW/mg) - Endo Down')
for i in range(0,4):
    ax_1[i].legend(loc=3)

fig_1.savefig('DSC Example.png')


# DSC Sample D
df_e = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Sample D - 10K min.csv')
df_f = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Sample D - 20K min.csv')
df_g = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Sample D - 25K min.csv')
df_h = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\DSC ' \
                    '- Sample D - 30K min.csv')

fig_2,ax_2 = plt.subplots(4, sharex=True, sharey=True, figsize=(5, 6), dpi=300)
ax_2[0].plot(df_e['Temperature'], df_e['Heat Flow'], color='C4', linewidth=1, \
           label='10 K/min')
ax_2[1].plot(df_f['Temperature'], df_f['Heat Flow'], color='C5', linewidth=1, \
           label='20 K/min')
ax_2[2].plot(df_g['Temperature'], df_g['Heat Flow'], color='C6', linewidth=1, \
           label='25 K/min')
ax_2[3].plot(df_h['Temperature'], df_h['Heat Flow'], color='C7', linewidth=1, \
           label='30 K/min')
    
fig_2.suptitle('DSC Scan Sample D')
ax_2[3].set(xlabel=r'Temperature ($\degree$C)')
ax_2[1].set(ylabel='Heat Flow (mW/mg) - Endo Down')
for i in range(0,4):
    ax_2[i].legend(loc=3)

fig_2.savefig('DSC Sample.png')

peaks_1 = sig.find_peaks(df_e['Heat Flow'], prominence=0.2, width=10, distance=10)
peaks_2 = sig.find_peaks(df_f['Heat Flow'], prominence=0.2, width=10, distance=10)
peaks_3 = sig.find_peaks(df_g['Heat Flow'], prominence=0.2, width=10, distance=10)
peaks_4 = sig.find_peaks(df_h['Heat Flow'], prominence=0.2, width=10, distance=10)

e = df_e['Temperature'][peaks_1[0]].to_numpy()
f = df_f['Temperature'][peaks_2[0]].to_numpy()
g = df_g['Temperature'][peaks_3[0]].to_numpy()
h = df_h['Temperature'][peaks_4[0]].to_numpy()

T_p = e[0], f[0], g[0], h[0]
T_p = tuple(i + 273.15 for i in T_p)

print(T_p)


for i in range(0,4):
    X[i] = 1 / T_p[i]
for i in range(0,4):
    Y[i] = np.log(B[i] * X[i])
    
lr = LinearRegression().fit(X.reshape(-1, 1), Y.reshape(-1, 1))
E_c = lr.coef_[0,0] * -R

print(E_c,' J/mol')

x_vals = np.linspace(0.00095, 0.00105, 2)
y_vals = lr.coef_[0,0] * x_vals + lr.intercept_[0]

fig_3 = plt.figure(dpi=300)
plt.scatter(X, Y, color='k')
plt.plot(x_vals, y_vals, color='C7', label=r'$E_c = 115$ kJ/mol')
    
plt.title('Crystallization Energy')
plt.xlabel(r'1 / $T_p$ (K$^{-1}$)')
plt.ylabel(r'$\ln{(B/T_P)}$ (A.U.)')
plt.legend()

fig_3.savefig('Crystallization.png')

# FTIR
df_i = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\FTIR ' \
                    '- 2 cm - 16 scans.csv')
df_j = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\FTIR ' \
                    '- 2 cm - 32 scans.csv')
df_k = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\FTIR ' \
                    '- 8 cm - 16 scans.csv')
df_l = pd.read_csv (r'C:\Users\andre\Documents\GitHub\MSE-353-Lab-02\FTIR ' \
                    '- 8 cm - 32 scans.csv')

fig_4,ax_4 = plt.subplots(4, sharex=True, sharey=True, figsize=(5, 6), dpi=300)
ax_4[0].plot(df_i['Wavenumber'], df_i['Transmittance'], color='C0', linewidth=1, \
           label=r'2 cm$^{-1}$ | 16 scans')
ax_4[1].plot(df_j['Wavenumber'], df_j['Transmittance'], color='C1', linewidth=1, \
           label=r'2 cm$^{-1}$ | 32 scans')
ax_4[2].plot(df_k['Wavenumber'], df_k['Transmittance'], color='C2', linewidth=1, \
           label=r'8 cm$^{-1}$ | 16 scans')
ax_4[3].plot(df_l['Wavenumber'], df_l['Transmittance'], color='C3', linewidth=1, \
           label=r'8 cm$^{-1}$ | 32 scans')
    
fig_4.suptitle('FTIR')
ax_4[3].set(xlabel=r'Wavenumber (cm$^{-1}$)')
ax_4[1].set(ylabel='Transmittance (%)')
for i in range(0,4):
    ax_4[i].legend()

fig_4.savefig('FTIR.png')