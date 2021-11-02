#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
folder = "Data/"
file0 = "Params.txt"
params = np.loadtxt(file0)
print(params)
N = int(params[0])

#%%
file1 = "Speed_data.txt"
Data1 = np.loadtxt(folder+file1)


#%%
Vx = Data1[:, 0]
Vy = Data1[:, 1]
Vz = Data1[:, 2]
print(np.size(Vx))

#%%
file2 = "System_data.txt"
Data2 = np.loadtxt(folder+file2)

#%%
t = Data2[1:, 0]
E = Data2[1:, 1]
Q = Data2[1:, 2]
T = Data2[1:, 3]
r2 = Data2[1:, 4]

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
plt.tight_layout(1.14)
ax1.hist(Vx, density = False, log = False)
plt.title(str(np.size(Vx)) + 'points')
#plt.show()
#fig, ax = plt.subplot()
ax2.hist(Vy, density = False, log = False)
#plt.show()
#fig, ax = plt.subplot()
ax3.hist(Vz, density = False, log = False)
#plt.show()
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = False)
plt.title('v histogram', fontsize = 12)
plt.show()

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
plt.tight_layout(1.14)
ax1.hist(Vx**2, density = True, log = True)
plt.title(str(np.size(Vx)) + 'points')
#plt.show()
#fig, ax = plt.subplot()
ax2.hist(Vy**2, density = True, log = True)
#plt.show()
#fig, ax = plt.subplot()
ax3.hist(Vz**2, density = True, log = True)
#plt.show()
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = True)
plt.title('v histogram', fontsize = 12)
plt.show()

#%%
fig, ax1 = plt.subplots()
ax1.plot(t, E)
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
ax3.plot(t, T)
plt.title('T from t', fontsize = 12)
plt.show()

#%%
if_b = 0
A = np.vstack([t, np.full(len(t), if_b)]).T
m, c = np.linalg.lstsq(A, r2, rcond=None)[0]
D = np.round(m/6, 4)
print('Коэф. диффузии D = ' + str(D))

fig, ax4 = plt.subplots()
ax4.plot(t, r2, 'b')
ax4.plot(t, m*t+c, 'r')
plt.title('r^2 from t', fontsize = 12)
plt.show()
#%%



