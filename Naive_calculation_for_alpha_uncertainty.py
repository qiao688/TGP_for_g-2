import numpy as np
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
root_s = dataset[:,0]
R_ratio = dataset[:,1]
root_var = dataset[:,2]

Kalpha = 2 * (root_s[1:] - root_s[:-1])/ (root_s[1:]*(Mz**2-root_s[1:]**2))
print(Kalpha.shape)

alpha_had_var = 0
for i in range(858):
    alpha_had_var += Kalpha[i]**2*(root_var[i]**2)
    
alpha_uncertainty = alpha*(Mz**2)/(3*np.pi)*(alpha_had_var**0.5)
print(alpha_uncertainty)
