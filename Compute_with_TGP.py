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
    
def Kernel_hat(root_s):
    K = Kernel(root_s)
    return 3*root_s**2/Mmuon**2*K
         
    
R_mean = np.load("mean_TGP.npy")
R_cov = np.load("cov_TGP.npy")
K_mean = Kernel(root_s)
Kamuon_var = Kernel(root_s)[1:] *2 * (root_s[1:] - root_s[:-1])/ root_s[1:]
Kalpha_var = 2 * (root_s[1:] - root_s[:-1])/ (root_s[1:]*(Mz**2-root_s[1:]**2))

amuon_hvp_mean = 0
for i in range(2574):
    amuon_hvp_mean += (root_s[i+1]-root_s[i])*(2*K_mean[i]*R_mean[i]/(root_s[i]))
    
alpha_had_mean = 0
for i in range(2574):
    alpha_had_mean += (root_s[i+1]-root_s[i])*(2*R_mean[i]/(root_s[i]*(Mz**2-root_s[i]**2)))
    
amuon_hvp_var = 0
for j in range(2574):
    amuon_hvp_var += Kamuon_var[j]**2*R_cov[j,j]
    for m in range(j+1,2574):
        amuon_hvp_var += 2*Kamuon_var[j]*Kamuon_var[m]*R_cov[j,m]
    
alpha_had_var = 0
for j in range(2574):
    alpha_had_var += Kalpha_var[j]**2*R_cov[j,j]
    for m in range(j+1,2574):
        alpha_had_var += 2*Kalpha_var[j]*Kalpha_var[m]*R_cov[j,m]      

amuon_mean = alpha**2/(3*np.pi**2)*amuon_hvp_mean    
alpha_mean = alpha*(Mz**2)/(3*np.pi)*alpha_had_mean
amuon_uncertainty = alpha**2/(3*np.pi**2)*amuon_hvp_var**0.5
alpha_uncertainty = alpha*(Mz**2)/(3*np.pi)*alpha_had_var**0.5
print(amuon_mean)
print(alpha_mean)
print(amuon_uncertainty)
print(alpha_uncertainty)
