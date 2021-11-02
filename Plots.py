#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
#%%
folder = "Data/"
file0 = "Params.txt"
params = np.loadtxt(file0)
print(params)
N = int(params[0])
time = params[4]
tau = params[6]

#%%
file1 = "Speed_data.txt"
Data1 = np.loadtxt(folder+file1)


#%%
V = Data1
Vx = Data1[:, 0]
Vy = Data1[:, 1]
Vz = Data1[:, 2]
print(np.size(Vx))

#%%
file2 = "System_data.txt"
Data2 = np.loadtxt(folder+file2)

#%%
t = Data2[:, 0]
E = Data2[:, 1]
Q = Data2[:, 2]
T = Data2[:, 3]
r2 = Data2[:, 4]

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
ax1.hist(Vx, density = False, log = False)
ax1.set_title('Vx histogram', fontsize = 12)
ax2.hist(Vy, density = False, log = False)
ax2.set_title('Vy histogram', fontsize = 12)
ax3.hist(Vz, density = False, log = False)
ax3.set_title('Vz histogram', fontsize = 12)
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = False)
ax4.set_title('V histogram', fontsize = 12)
plt.tight_layout()
plt.show()

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
ax1.hist(Vx**2, density = True, log = True)
ax1.set_title('Vx histogram', fontsize = 12)
ax2.hist(Vy**2, density = True, log = True)
ax2.set_title('Vy histogram', fontsize = 12)
ax3.hist(Vz**2, density = True, log = True)
ax3.set_title('Vz histogram', fontsize = 12)
#ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = True)
h, bins, p = ax4.hist(np.hstack((Vx, Vy, Vz))**2, density = True, log = True)
ax4.set_title('V histogram', fontsize = 12)
plt.tight_layout()
plt.show()
#%%
fig, (ax1, ax2) = plt.subplots(1, 2)
h, bins, p = ax1.hist(np.hstack((Vx, Vy, Vz)), density = False, log = False)
ax2.scatter(bins[:-1]**2, np.log(h))
plt.show()

#%%
fig, ax1 = plt.subplots()
ax1.plot(t[1:], E[1:])
ax1.set_ylim([0.999*np.min(E[1:]), 1.001*np.max(E[1:])])
plt.title('E from t', fontsize = 12)
plt.show()
#%%
fig, ax2 = plt.subplots()
ax2.plot(t, Q)
#ax2.set_ylim([0.999*np.min(Q), 1.001*np.max(Q)])
plt.title('Impulse from t', fontsize = 12)
plt.show()
#%%
fig, ax3 = plt.subplots()
ax3.plot(t[1:], T[1:])
plt.title('T from t', fontsize = 12)
plt.show()

#%%
cut1 = 2
x1 = t[cut1:]
y1 = r2[cut1:]
if_b = 1
A = np.vstack([x1, np.full(len(x1), if_b)]).T
m, c = np.linalg.lstsq(A, y1, rcond=None)[0]
D1 = np.round(m/6, 4)
print('Коэф. диффузии D1 = ' + str(D1))
#%%
fig, ax4 = plt.subplots()
ax4.plot(x1, y1, 'b', label='model data')
ax4.plot(x1, m*x1+c, 'r', label='LSF approx')
ax4.set_xlabel('t')
ax4.set_ylabel('r^2')
ax4.legend(loc='lower right', fontsize=12)
plt.title('r^2 from t', fontsize = 12)
plt.show()
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

#%%
vx = np.reshape(Vx, (int(time/tau), N))
vy = np.reshape(Vy, (int(time/tau), N))
vz = np.reshape(Vz, (int(time/tau), N))
autocorr0 = np.zeros(2*int(time/tau)-1)
for i in range(N):
    sigx = vx[:,i]
    sigy = vy[:,i]
    sigz = vz[:,i]
    autocorr0 += fftconvolve(sigx, sigx[::-1], mode='full')
    autocorr0 += fftconvolve(sigy, sigy[::-1], mode='full')
    autocorr0 += fftconvolve(sigy, sigy[::-1], mode='full')
autocorr0 /= N*int(time/tau)
autocorr = autocorr0[len(sigx)-1:]
D3 = np.round(np.sum(autocorr)*tau/3, 4)
print('Коэф. диффузии D3 = ' + str(D3))
#%%
fig, ax = plt.subplots()
ax.plot(np.arange(len(sigx)), autocorr)
ax.plot(np.arange(len(sigx)), Z[cut2:])
ax.set_title('АКФС2', fontsize=12)
plt.show()
