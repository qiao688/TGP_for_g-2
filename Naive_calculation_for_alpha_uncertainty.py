import numpy as np
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
root_s = dataset[:,0]
R_ratio = dataset[:,1]
root_var = dataset[:,2]

w_alpha = 0
Calpha = []
for i in range(859):
    if i == 0:
      w_alpha = 0.5 * (root_s[1] - root_s[0])
      Calpha.append((2*alpha*(Mz**2)/(3*np.pi))*(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2))))
    elif i == 858:
      w_amuon = 0.5 * (root_s[858] - root_s[857])
      Calpha.append((2*alpha*(Mz**2)/(3*np.pi))*(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2))))
    if not i==0 and not i==858 :
      w_alpha = 0.5 * (root_s[i+1] - root_s[i-1])
      Calpha.append((2*alpha*(Mz**2)/(3*np.pi))*(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2))))
Calpha_naive = np.array(Calpha).flatten()

alpha_var = 0
for i in range(859):
    alpha_var += Calpha_naive[i]**2*(root_var[i]**2)
    
alpha_uncertainty = alpha_var**0.5
print(alpha_uncertainty)
