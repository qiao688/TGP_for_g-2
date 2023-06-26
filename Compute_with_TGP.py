import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
Mpion = 139.57039e-3 #Mass pion = 139.57039 ± 0.00018 MeV
alpha = 1/137.035999084 #fine-structure constant: 7.297 352 5693(11)×10−3 = 1/137.035 999 084(21) 0.15(uncertainty)
root_s = np.load("x_predict_TGP.npy")
R_mean = np.load("mean_TGP.npy")
R_cov = np.load("cov_TGP.npy")

# Calculate C_HVP

def Kernel(root_s):
    a = 4*Mmuon**2/root_s**2
    beita = np.sqrt(1-a)
    x = (1-beita)/(1+beita)
    return 0.5*x**2*(2-x**2) + (1+x**2)*(1+x)**2/x**2 * (np.log1p(x)-x+0.5*x**2) + (1+x)/(1-x)*x**2*np.log(x)

K = Kernel(root_s)
w_amuon = 0
Camuon_HVP = []
for i in range(2575):
    if i == 0:
      w_amuon = 0.5 * (root_s[1] - root_s[0])
      Camuon_HVP.append(K[i] / root_s[i] * w_amuon)
    elif i == 2574:
      w_amuon = 0.5 * (root_s[2574] - root_s[2573])
      Camuon_HVP.append(K[i] / root_s[i] * w_amuon)
    if not i==0 and not i==2574 :
      w_amuon = 0.5 * (root_s[i+1] - root_s[i-1])
      Camuon_HVP.append(K[i] / root_s[i] * w_amuon)
C_HVP = np.array(Camuon_HVP).flatten()
np.save("C_HVP",C_HVP)

# Calculate amuon_HVP 

amuon_mean = 0
for i in range(2575):
    amuon_mean += C_HVP[i] * R_mean[i]
amuon_HVP = 2*alpha**2/(3*np.pi**2)*amuon_mean
print(amuon_HVP) 

# Calculate amuon_uncertainty

amuon_var = 0
for i in range(2575):
    for j in range(2575):
        amuon_var += C_HVP[i] * C_HVP[j] * R_cov[i,j]
amuon_uncertainty = 2*alpha**2/(3*np.pi**2)*amuon_var**0.5
print(amuon_uncertainty)

# Calculate C_had

w_alpha = 0
Calpha_had = []
for i in range(2575):
    if i == 0:
      w_alpha = 0.5 * (root_s[1] - root_s[0])
      Calpha_had.append(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2)))
    elif i == 2574:
      w_alpha = 0.5 * (root_s[2574] - root_s[2573])
      Calpha_had.append(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2)))
    if not i==0 and not i==2574 :
      w_alpha = 0.5 * (root_s[i+1] - root_s[i-1])
      Calpha_had.append(w_alpha / (root_s[i] * (Mz**2 - root_s[i]**2)))
C_had = np.array(Calpha_had).flatten()
np.save("C_had",C_had)

# Calculate alpha_had

alpha_mean = 0
for i in range(2575):
    alpha_mean += C_had[i] * R_mean[i]
alpha_had = 2*alpha*(Mz**2)/(3*np.pi)*alpha_mean
print(alpha_had) 

# Calculate alpha_uncertainty

alpha_var = 0
for i in range(2575):
    for j in range(2575):
        alpha_var += C_had[i] * C_had[j] * R_cov[i,j]
alpha_uncertainty = 2*alpha*(Mz**2)/(3*np.pi)*alpha_var**0.5
print(alpha_uncertainty)
