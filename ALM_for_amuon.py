import numpy as np
Mmuon = 105.6583755e-3 #Mass muon = 105.6583755 ± 0.0000023 MeV
root_s = np.load("x_predict_TGP.npy")

def Kernel(root_s):
    a = 4*Mmuon**2/root_s**2
    beita = np.sqrt(1-a)
    x = (1-beita)/(1+beita)
    return 0.5*x**2*(2-x**2) + (1+x**2)*(1+x)**2/x**2 * (np.log1p(x)-x+0.5*x**2) + (1+x)/(1-x)*x**2*np.log(x)
#print(Kernel(root_s))

Cov = np.load("cov_TGP.npy")
#print(Cov)
C_hvp = Kernel(root_s)[1:] * (root_s[1:] - root_s[:-1])/ root_s[1:]

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
