#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
a = 13
folder = "Full_nve_data/"
subfolder1 = "Points/"
subfolder2 = "Speed/"
subfolder3 = "System/"
#%%
file0 = "Parameters_"+str(a)+".txt"
params = np.loadtxt(folder+subfolder3+file0)

N = int(params[0])
time = params[4]
tau = params[6]
#%%
file1 = "Speed_data_"+str(a)+".txt"
Data1 = np.loadtxt(folder+subfolder2+file1)

V = Data1
Vx = Data1[:, 0]
Vy = Data1[:, 1]
Vz = Data1[:, 2]
print(np.size(Vx))
#%%
file2 = "System_data_"+str(a)+".txt"
Data2 = np.loadtxt(folder+subfolder3+file2)

t = Data2[:, 0]
E = Data2[:, 1]
Q = Data2[:, 2]
T = Data2[:, 3]
r2 = Data2[:, 4]
#%%
n_bins = 20
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
ax1.hist(Vx, bins = n_bins, density = False, log = False)
ax1.set_title('Vx histogram', fontsize = 12)
ax2.hist(Vy, bins = n_bins, density = False, log = False)
ax2.set_title('Vy histogram', fontsize = 12)
ax3.hist(Vz, bins = n_bins, density = False, log = False)
ax3.set_title('Vz histogram', fontsize = 12)
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), bins = n_bins, density = False)
ax4.set_title('V histogram', fontsize = 12)
plt.tight_layout()
plt.show()
#%%
# =============================================================================
# Vxyz = np.hstack((Vx, Vy, Vz))
# Vmax = np.sqrt(np.max(Vxyz))
# Vmin = 0.0
# dV = 1.0
# Vh = np.arange(Vmin, Vmax, dV)**2
# =============================================================================
#%%
# =============================================================================
# fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
# ax1.hist(Vx**2, density = True, log = True)
# ax1.set_title('Vx histogram', fontsize = 12)
# ax2.hist(Vy**2, density = True, log = True)
# ax2.set_title('Vy histogram', fontsize = 12)
# ax3.hist(Vz**2, density = True, log = True)
# ax3.set_title('Vz histogram', fontsize = 12)
# #ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = True)
# h, bins, p = ax4.hist(np.hstack((Vx, Vy, Vz))**2, bins = 50, density = True, log = True)
# ax4.set_title('V histogram', fontsize = 12)
# plt.tight_layout()
# plt.show()
# =============================================================================
#%%
n_bins2 = 40
h, bins = np.histogram(np.abs(np.hstack((Vx, Vy, Vz))), bins = n_bins2, density=False)

if_b = 1
x = bins[np.nonzero(h)]**2
y = np.log(h[np.nonzero(h)])
A = np.vstack([x, np.full(len(x), if_b)]).T
m, c = np.linalg.lstsq(A, y, rcond=None)[0]

fig, ax1 = plt.subplots()
ax1.scatter(x, y, label='Эксп. точки')
ax1.plot(x, m*x+c, label='Аппркосимация')
ax1.legend(loc='upper right', fontsize=12)
ax1.set_title('Линеаризорованная плотность распределения скоростей')
plt.show()

print('Аппроксимационная температура: ', np.round(-1/2/m, 3))

#%%
fig, ax1 = plt.subplots()
ax1.plot(t[1:], E[1:])
#ax1.set_ylim([0.999*np.min(E[1:]), 1.001*np.max(E[1:])])
ax1.set_xlabel('t')
plt.title('E from t', fontsize = 12)
plt.show()
#%%
fig, ax2 = plt.subplots()
ax2.plot(t, Q)
#ax2.set_ylim([0.999*np.min(Q), 1.001*np.max(Q)])
ax2.set_xlabel('t')
plt.title('Impulse from t', fontsize = 12)
plt.show()
#%%
fig, ax3 = plt.subplots()
ax3.plot(t[1:], T[1:])
plt.title('T from t', fontsize = 12)
ax3.set_xlabel('t')
plt.show()

print('Средняя температура: ', np.round(np.mean(T), 3))

#%%

#%%
# =============================================================================
# cut2 = 0
# v = np.reshape(V, [int(time/tau), N, 3])
# Z = np.sum(v[:, :, :] * v[0, :, :]/N, axis=(1, 2))
# D2 = np.round(np.sum(Z[cut2:])*tau/3, 4)
# print('Коэф. диффузии D2 = ' + str(D2))
# 
# #%%
# fig, ax = plt.subplots()
# ax.plot(t[cut2:], Z[cut2:])
# ax.set_title('АКФС1', fontsize=12)
# plt.show()
# =============================================================================

