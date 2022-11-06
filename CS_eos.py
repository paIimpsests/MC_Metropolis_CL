import matplotlib.pyplot as plt
import numpy as np


def CS_p(rho):
    p = rho * (1 + (np.pi / 6 * rho) + (np.pi / 6 * rho)**2 - (np.pi / 6 * rho)**3) / (1 - np.pi / 6 * rho)**3
    return p

rho = np.linspace(0,1.2,100)

i = np.empty([])
rho2 = np.empty([])
i, rho2 = np.loadtxt("density.txt", delimiter="\t", unpack=True)

fig, ax = plt.subplots()
ax.plot(rho, CS_p(rho))
ax.scatter(rho2[-1], 5, marker="d")
plt.xlabel(r"$\rho \sigma^3$")
plt.ylabel(r"$P \sigma^3 \varepsilon^{-1}$")
plt.title(r"$\textrm{Carnahan-Starling equation of state}$")
plt.savefig("./figures/CS_eos.pdf")
plt.show()


fig, ax = plt.subplots()
ax.plot(i, rho2, label=r"$P \sigma^3 \varepsilon^{-1} = 5$")
plt.xlabel(r"$Cycles$")
plt.ylabel(r"$\rho \sigma^3$")
plt.title(r"$\textrm{MC Metropolis implementation}$")
plt.legend()
plt.show()

