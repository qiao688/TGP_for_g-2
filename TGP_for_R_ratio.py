import numpy as np
import kingpin as kp
import matplotlib.pyplot as plt

dataset = np.loadtxt("CM_energy(0.3-1.937GeV).txt")
x = dataset[:,0]
y = dataset[:,1]
noise = dataset[:,2]

p = []
for i in range(858):
    p.append(np.linspace(start = x[i], stop = x[i+1], num = 3, endpoint = False))
p1 = np.array(p).flatten()
p2 = np.append(p1,1.937)

# Bespoke choices

model = kp.Celerite2(x, y, noise = noise, x_predict=p2)

# Priors and proposals for mean, sigma, length and nugget

mean = kp.Uniform(0., 150.)
sigma = kp.Uniform(0., 500.)
length = kp.Uniform(0., 5.)
params = kp.Independent(mean, sigma, length)

# Run RJ-MCMC

rj = kp.TGP(model, params)
rj.walk(thin=4, n_threads=8, n_iter=100000, n_burn=5000, n_iter_params=1)

print(rj.acceptance)
print(rj.arviz_summary())
print(rj.mean)
print(rj.cov)
print(rj.mean.shape)
print(rj.cov.shape)
mean = rj.mean
cov = rj.cov
np.save("x_predict_TGP",p2)
np.save("mean_TGP", mean)
np.save("cov_TGP", cov)
np.save("edge_counts", rj.edge_counts)
rj.plot()
plt.xlabel(r"$\sqrt{s} [GeV]$")
plt.ylabel("R")
plt.savefig("plot_TGP.pdf")
