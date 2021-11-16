#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
#%%
a = 12
folder = "Full_nve_data/"
subfolder1 = "Points/"
subfolder2 = "Speed/"
subfolder3 = "System/"
#%%
file0 = "Parameters_"+str(a)+".txt"
params = np.loadtxt(folder+subfolder3+file0)
print(params)
N = int(params[0])
x_size = int(params[1])
y_size = int(params[2])
z_size = int(params[3])
time = params[4]
tau = params[6]
sigma = params[8]

#%%
file1 = "Speed_data_"+str(a)+".txt"
Data1 = np.loadtxt(folder+subfolder2+file1)

V = Data1
Vx = Data1[:, 0]
Vy = Data1[:, 1]
Vz = Data1[:, 2]
print(np.size(Vx))

#%%
v_mean = np.mean(np.sqrt(Vx**2 + Vy**2 + Vz**2))
print('Средняя скорость: ', np.round(v_mean, 2))

#%%
file2 = "System_data_"+str(a)+".txt"
Data2 = np.loadtxt(folder+subfolder3+file2)

t = Data2[:, 0]
E = Data2[:, 1]
Q = Data2[:, 2]
T = Data2[:, 3]
r2 = Data2[:, 4]
#%%
a = 0.4
S = np.pi*(2*sigma*a)**2
n = N/x_size/y_size/z_size
mfp1 = np.round(1/n/S, 5)
print('Длина свободного пробега 1: ', mfp1)

#%%
cut1 = 2
x1 = t[cut1:]
y1 = r2[cut1:]
if_b = 1
A = np.vstack([x1, np.full(len(x1), if_b)]).T
m, c = np.linalg.lstsq(A, y1, rcond=None)[0]
D1 = np.round(m/6, 4)
print('Коэф. диффузии D1 = ' + str(D1))

mfp2 = np.round(3*D1/v_mean, 5)
print('Длина свободного пробега 2: ', mfp2)

a_exp = np.round(np.sqrt(v_mean/3/np.pi/n/D1)/2/sigma, 4)
print('Эксп. параметр радиуса: ', a_exp)
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
ax.plot(t, autocorr)
ax.set_xlabel('t')
ax.set_title('АКФС', fontsize=12)
plt.show()
