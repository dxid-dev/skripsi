import numpy as np
import matplotlib.pyplot as plt

# Inisiasi dan pembacaan data VACF
N = 868
Nc = 1000 
dt = 0.001
omega = np.arange(1, 380.5, 0.5)
nu = omega / 2 / np.pi
Nf = Nc*1
num_frame = Nf
fileName = "pdos.lammpstrj"

def find_pdos(v_all, Nc, dt, omega):
    Nf = v_all.shape[0] 
    M  = Nf - Nc 
    vacf = np.zeros(Nc)
    for nc in range(Nc):
        ratio = (nc+1)/Nc * 100    
        print("Calculate PDOS Progress %s%%" %ratio)
        for m in range(M+1):
            delta = np.sum(v_all[m + 0]*v_all[m + nc])
            vacf[nc] = vacf[nc] + delta
    vacf = vacf / vacf[0] 
    vacf_output = vacf 
    vacf = vacf*(np.cos(np.pi*np.arange(Nc)/Nc)+1)*0.5
    vacf = vacf*np.append(np.ones(1), 2*np.ones(Nc-1))/np.pi
    pdos = np.zeros(len(omega))
    for n in range(len(omega)):
        pdos[n] = dt * sum(vacf * np.cos(omega[n] * np.arange(Nc) * dt))
    return(vacf_output, pdos)

# Kalkulasi PDOS
v_all = np.zeros((num_frame, N, 3))
fin = open(fileName, "r")
for i in range(num_frame):
    ratio = (i+1)/num_frame * 100 
    print("Read Data Progress %s%%" %ratio)
    initial = i * (9 + N)
    for j in range(9):
        fin.readline()
    for k in range(N):
        line = fin.readline().split()[2:]
        line = [float(l) for l in line]
        v_all[i, k] = line

vacf, pdos = find_pdos(v_all, Nc, dt, omega)
t = np.arange(Nc)*dt