import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
Mpion = 139.57039e-3 #Mass pion = 139.57039 ± 0.00018 MeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
root_s = dataset[:,0]
R_ratio = dataset[:,1]
root_var = dataset[:,2]

def beita(root_s):
    return np.sqrt(np.abs(1-4*Mmuon**2 / root_s**2))

def Kernel(root_s):
    beita1 = beita(root_s)
    x = (1-beita1) / (1+beita1)
    return (0.5*x**2 * (2-x**2) + (1+x**2)*(1+x)**2/(x**2) * (np.log1p(x)-x+0.5*x**2)+(1+x)/(1-x)*(x**2)*np.log(x))
    
K = Kernel(root_s)
w_amuon = 0
Camuon = []
for i in range(859):
    if i == 0:
      w_amuon = 0.5 * (root_s[1] - root_s[0])
      Camuon.append(2*alpha**2/(3*np.pi**2)*(K[i] / root_s[i] * w_amuon))
    elif i == 858:
      w_amuon = 0.5 * (root_s[858] - root_s[857])
      Camuon.append(2*alpha**2/(3*np.pi**2)*(K[i] / root_s[i] * w_amuon))
    if not i==0 and not i==858 :
      w_amuon = 0.5 * (root_s[i+1] - root_s[i-1])
      Camuon.append(2*alpha**2/(3*np.pi**2)*(K[i] / root_s[i] * w_amuon))
Camuon_naive = np.array(Camuon).flatten()

amuon_var = 0
for i in range(859):
    amuon_var += Camuon_naive[i]**2 * root_var[i]**2
 
amuon_uncertainty = amuon_var**0.5
print(amuon_uncertainty)
