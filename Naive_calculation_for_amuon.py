import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
Mpion = 139.57039e-3 #Mass pion = 139.57039 ± 0.00018 MeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
root_s = dataset[:,0]
R_ratio = dataset[:,1]

def beita(root_s):
    return np.sqrt(np.abs(1-4*Mmuon**2 / root_s**2))

def Kernel(root_s):
    beita1 = beita(root_s)
    x = (1-beita1) / (1+beita1)
    return (0.5*x**2 * (2-x**2) + (1+x**2)*(1+x)**2/(x**2) * (np.log1p(x)-x+0.5*x**2)+(1+x)/(1-x)*(x**2)*np.log(x))

def amuon_hvp(R_ratio, root_s):
    Kernel1 = Kernel(root_s)
    K1 = 2*Kernel1*R_ratio/root_s
    return alpha**2/(3*np.pi**2)*np.trapz(K1, x=root_s)

print(amuon_hvp(R_ratio, root_s))
