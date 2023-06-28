import numpy as np
root_s = np.load("x_predict_TGP.npy")
Cov = np.load("cov_TGP.npy")
C_hvp = np.load("C_HVP.npy")

Vmax = 0
n = 0
for i in range(2575):
    Var = 0
    for j in range(2575):
        Var += C_hvp[i] * C_hvp[j] * Cov[i,j]
    #print(Var)
    if Var > Vmax:
        Vmax = Var
        n = i
    else:
        continue
print(Vmax,n)
#print(root_s[704])
