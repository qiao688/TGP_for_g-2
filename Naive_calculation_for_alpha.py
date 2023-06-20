import numpy as np
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
root_s = dataset[:,0]
R_ratio = dataset[:,1]

def alpha_had(R_ratio, root_s):
    K = 2*R_ratio/(root_s*(Mz**2-root_s**2))
    return alpha*(Mz**2)/(3*np.pi)*np.trapz(K, x=root_s)
    
print(alpha_had(R_ratio, root_s))
