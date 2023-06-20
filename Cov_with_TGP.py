from Compute_with_TGP import amuon_mean
from Compute_with_TGP import alpha_mean
from Compute_with_TGP import amuon_uncertainty
from Compute_with_TGP import alpha_uncertainty

import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
Mpion = 139.57039e-3 #Mass pion = 139.57039 ± 0.00018 MeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
root_s = np.load("x_predict_TGP.npy")

def Kernel(root_s):
    a = 4*Mmuon**2/root_s**2
    beita = np.sqrt(1-a)
    x = (1-beita)/(1+beita)
    return 0.5*x**2*(2-x**2) + (1+x**2)*(1+x)**2/x**2 * (np.log1p(x)-x+0.5*x**2) + (1+x)/(1-x)*x**2*np.log(x)

R_mean = np.load("mean_TGP.npy")
R_cov = np.load("cov_TGP.npy")
E_X = amuon_mean
E_Y = alpha_mean
Kamuon = Kernel(root_s)[1:] *2 * (root_s[1:] - root_s[:-1])/ root_s[1:]
Kalpha = 2 * (root_s[1:] - root_s[:-1])/ (root_s[1:]*(Mz**2-root_s[1:]**2))

def Cov(root_s):
    E = 0
    for i in range(2574):
        for j in range (2574):
            E += Kamuon[i]*Kalpha[j]*(R_cov[i,j])
    return (alpha**2/(3*np.pi**2))*(alpha*(Mz**2)/(3*np.pi))*E

    
def Coe(root_s):
    Cov1 = Cov(root_s)
    sigma_X = amuon_uncertainty
    sigma_Y = alpha_uncertainty
    return Cov1/(sigma_X*sigma_Y)
    
print(Cov(root_s))
print(Coe(root_s))
