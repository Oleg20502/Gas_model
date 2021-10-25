#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
file = "Speed_data.txt"
Data = np.loadtxt(file)

Vx = Data[:, 1]
Vy = Data[:, 2]
Vz = Data[:, 3]

#%%
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
ax1.hist(Vx)
#plt.show()
#fig, ax = plt.subplot()
ax2.hist(Vy)
#plt.show()
#fig, ax = plt.subplot()
ax3.hist(Vz)
#plt.show()
#%%
ax4.hist(np.sqrt(Vx**2 + Vy**2 + Vz**2))
plt.show()