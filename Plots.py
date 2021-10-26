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
Vx = Data1[:N, 0]
Vy = Data1[:N, 1]
Vz = Data1[:N, 2]
print(np.size(Vx))

#%%
file2 = "System_data.txt"
Data2 = np.loadtxt(folder+file2)

#%%
t = Data2[1:, 0]
E = Data2[1:, 1]
Q = Data2[1:, 2]
T = Data2[1:, 3]

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
plt.tight_layout(1.14)
h1, b, p = ax1.hist(Vx**2, density = True, log = True)
plt.title(str(np.size(Vx)) + 'points')
#plt.show()
#fig, ax = plt.subplot()
ax2.hist(Vx**2, density = True, log = True)
#plt.show()
#fig, ax = plt.subplot()
ax3.hist(Vx**2, density = True, log = True)
#plt.show()
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2), density = True)
plt.title('v histogram', fontsize = 12)
plt.show()

#%%
fig, ax1 = plt.subplots()
ax1.plot(t, E)
#plt.ylim(0.9999*np.min(E), 1.00001*np.max(E))
plt.title('E from t', fontsize = 12)
plt.show()
#%%
fig, ax2 = plt.subplots()
ax2.plot(t, Q)
plt.title('Impulse from t', fontsize = 12)
plt.show()
#%%
fig, ax3 = plt.subplots()
ax3.plot(t, T)
plt.title('T from t', fontsize = 12)
plt.show()

#%%




