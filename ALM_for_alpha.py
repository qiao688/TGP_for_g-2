import numpy as np
root_s = np.load("x_predict_TGP.npy")
Cov = np.load("cov_TGP.npy")
C_had = np.load("C_had.npy")

Vmax = 0
n = 0
for i in range(2575):
    Var = 0
    for j in range(2575):
        Var += C_had[i] * C_had[j] * Cov[i,j]
    #print(Var)
    if Var > Vmax:
        Vmax = Var
        n = i
    else:
        continue
print(Vmax,n)
#print(root_s[704])
