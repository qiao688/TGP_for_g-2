import numpy as np
Mz = 91.1876 #Mass Z = 91.1876 ± 0.0021 GeV
root_s = np.load("x_predict_TGP.npy")
#print(root_s[1300])

Cov = np.load("cov_TGP.npy")
C_hvp = (root_s[1:] - root_s[:-1])/ (root_s[1:]*(Mz**2-root_s[1:]**2))

Vmax = 0
n = 0
for i in range(2574):
    Var = 0
    for j in range(2574):
        Var += C_hvp[i] * C_hvp[j] * Cov[i,j]
    #print(Var)
    if Var > Vmax:
        Vmax = Var
        n = i
    else:
        continue
print(Vmax,n)
