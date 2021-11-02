#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
def vec_mod(x, y, z):
    return np.sqrt(x** + y**2 + z**2)

def get_distance(Pos):
    global x_size, y_size, z_size
    dx = np.abs(Pos[i][0]-Pos[j][0])
    dy = np.abs(Pos[i][1]-Pos[j][1])
    dz = np.abs(Pos[i][2]-Pos[j][2])
    return vec_mod(np.min(dx, x_size-dx), np.min(dy, y_size-dy), np.min(dz, z_size-dz))

#%%
folder = "Data/"
file0 = "Params.txt"
params = np.loadtxt(file0)
N = int(params[0])
x_size = int(params[1])
y_size = int(params[2])
z_size = int(params[3])
time = params[4]
tau = params[6]

miss = 'Lattice='+str(x_size)+' 0 0 0 '+str(y_size)+' 0 0 0 '+str(z_size)

#%%
Pos = np.zeros((int(time/tau)*N, 3))
n = N+2
file1 = "Points_data.txt"
with open(folder+file1) as f:
    i = 0
    j = 0
    for line in f.readlines():
        if(i%n != 0 and i%n != 1):
            #print(np.array(list(map(float, line.split(' '))))[1:])
            Pos[j, :] = np.array(list(map(float, line.split(' '))))[1:]
            j += 1
        i += 1
#%%
Rmax = np.sqrt(x_size**2)/2
Rmin = 0.0
dR = 1.0
distances = np.zeros(int(time/tau)*N*(N-1)/2)
R  = np.arrange(Rmin, Rmax+dR, dR)
h = np.zeros((int(time/tau)*N, np.size(R)))
for k in range(int(time/tau)):
    for i in range(N):
        for j in range(i, N):
            distances[0] = get_distance(Pos[k])