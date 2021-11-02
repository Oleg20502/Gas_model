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
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = True)
ax4.set_title('V histogram', fontsize = 12)
plt.tight_layout()
plt.show()

#%%
fig, ax1 = plt.subplots()
ax1.plot(t[1:], E[1:])
ax1.set_ylim([0.999*np.min(E), 1.001*np.max(E)])
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
if_b = 0
A = np.vstack([t, np.full(len(t), if_b)]).T
m, c = np.linalg.lstsq(A, r2, rcond=None)[0]
D1 = np.round(m/6, 4)
print('Коэф. диффузии D1 = ' + str(D1))
#%%
fig, ax4 = plt.subplots()
ax4.plot(t, r2, 'b')
ax4.plot(t, m*t+c, 'r')
plt.title('r^2 from t', fontsize = 12)
plt.show()
#%%
v = np.reshape(V, [int(time/tau), N, 3])
Z = np.sum(v[:, :, :] * v[0, :, :]/N, axis=(1, 2))
D2 = np.round(np.sum(Z)*tau/3, 4)
print('Коэф. диффузии D2 = ' + str(D2))

#%%
fig, ax = plt.subplots()
ax.plot(t, Z)
ax.set_title('АКФС', fontsize=12)
plt.show()

#%%
# =============================================================================
# sig = Vx[::N]
# autocorr0 = fftconvolve(sig, sig[::-1], mode='full')
# autocorr = autocorr0[len(sig)-1:]
# D2 = np.round(np.sum(autocorr)*tau/3/4, 4)
# print('Коэф. диффузии D2 = ' + str(D2))
# #%%
# fig, ax = plt.subplots()
# ax.plot(np.arange(len(sig)), autocorr)
# #ax.plot(np.arange(len(sig)), autocorr)
# plt.show()
# =============================================================================
