#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
def vec_mod(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def get_distance(Pos, i, j):
    dx = Pos[i][0] - Pos[j][0]
    dy = Pos[i][1] - Pos[j][1]
    dz = Pos[i][2] - Pos[j][2]
    return vec_mod(dx, dy, dz)

def get_pgu_distance(Pos, i, j):
    global x_size, y_size, z_size
    dx = np.abs(Pos[i][0] - Pos[j][0])
    dy = np.abs(Pos[i][1] - Pos[j][1])
    dz = np.abs(Pos[i][2] - Pos[j][2])
    return vec_mod(np.min([dx, x_size-dx]), np.min([dy, y_size-dy]), np.min([dz, z_size-dz]))

#%%
a = 8
folder = "Full_nve_data/"
subfolder1 = "Points/"
subfolder2 = "Speed/"
subfolder3 = "System/"
file0 = "Parameters_"+str(a)+".txt"
params = np.loadtxt(folder+subfolder3+file0)
N = int(params[0])
x_size = int(params[1])
y_size = int(params[2])
z_size = int(params[3])
time = params[4]
tau = params[6]

miss = 'Lattice='+str(x_size)+' 0 0 0 '+str(y_size)+' 0 0 0 '+str(z_size)

#%%
Pos = np.zeros((int(time/tau), N, 3))
n = N+2
file1 = "Points_data_"+str(a)+".txt"
with open(folder+subfolder1+file1) as f:
    i = 0
    j = 0
    for line in f.readlines():
        if(i%n != 0 and i%n != 1):
            Pos[j//N, j%N, :] = np.array(list(map(float, line.split(' '))))[1:]
            j += 1
        i += 1
        
# =============================================================================
# for i in range(1, 9):
#     Pos[:, i*N:(i+1)*N, :] = Pos[:, :N, :]
# 
# Pos[:, N:2*N, 0] += x_size
# Pos[:, 2*N:3*N, 1] += y_size
# Pos[:, 3*N:4*N, 2] += z_size
# 
# Pos[:, 4*N:5*N, 0] -= x_size
# Pos[:, 5*N:6*N, 1] -= y_size
# Pos[:, 6*N:7*N, 2] -= z_size
# =============================================================================
#%%
# =============================================================================
# t_number = 10
# step = int(time/tau)//(t_number+1)
# distances = np.zeros(t_number * N * (N-1))
# I = 0
# for k in range(1, t_number+1):
#     for i in range(N):
#         for j in range(N):
#             if i != j:
#                 distances[I] = get_pgu_distance(Pos[k*step], i, j)
#                 I += 1
# distances /= N*t_number
# #%%
# fig, ax = plt.subplots()
# ax.hist(distances[:])
# plt.show()
# 
# #%%
# fig, ax = plt.subplots()
# ax.plot(np.arange(len(distances)), distances[:])
# plt.show()
# =============================================================================
#%%
Rmax = min(x_size, y_size, z_size)/2
Rmin = 0.3
dR = 0.1
t_number2 = 5
step2 = int(time/tau)//(t_number2+1)
distances2 = np.zeros(t_number2 * N * (N-1))
R = np.arange(Rmin, Rmax+dR, dR)
h2 = np.zeros(np.size(R)-1)
I = 0
for k in range(1, t_number2+1):
    for i in range(N):
        for j in range(N):
            if i != j:
                r = get_pgu_distance(Pos[k*step2], i, j)
                if Rmin < r < Rmax:
                    distances2[I] = r
                    h2[int((r - Rmin)//dR)] += 1
                    I += 1

#%%
#Rh = h2/np.sum(h2)/dR
# =============================================================================
# Rh = h2
# fig, ax = plt.subplots()
# ax.scatter(R[:-1], Rh, color = 'r')
# plt.show()
# =============================================================================
#%%
fig, ax = plt.subplots()
D = distances2[:I]
h3, b3, p = ax.hist(D, bins = R, density=True)
ax.set_title('Плотность распределения растояний между молекулами')
plt.show()

#%%
fig, ax = plt.subplots()
D = distances2[:I]
h4, b4, p = ax.hist(D, bins = 25, density=True, cumulative=True)
ax.set_title('Распределение растояний между молекулами')
plt.show()

#%%
n = N/x_size/y_size/z_size
g = h3/(4*np.pi*R[:-1]**2)/n * 50
fig, ax = plt.subplots()
ax.plot(R[:-1], g[:])
ax.set_title('Парная корреляционная функция')
plt.show()
