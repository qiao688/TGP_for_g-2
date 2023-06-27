from TGP_compute import amuon_uncertainty
from TGP_compute import alpha_uncertainty

import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
Mpion = 139.57039e-3 #Mass pion = 139.57039 ± 0.00018 MeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
root_s = np.load("x_predict_TGP.npy")
R_mean = np.load("mean_TGP.npy")
R_cov = np.load("cov_TGP.npy")
C_HVP = np.load("C_HVP.npy")
C_had = np.load("C_had.npy")

# Calculate the covariance between amuon_HVP and alpha_had

Cov = 0
for i in range(2575):
    for j in range(2575):
        Cov += C_HVP[i] * C_had[j] * R_cov[i,j]
print(Cov)

# Calculate the correlation coefficient

def Coe(root_s):
    Cov1 = Cov
    sigma_X = amuon_uncertainty
    sigma_Y = alpha_uncertainty
    return Cov1/(sigma_X*sigma_Y)
print(Coe(root_s))
