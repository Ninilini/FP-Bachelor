# coding=utf-8
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
import sympy
from uncertainties import correlated_values, correlation_matrix
import scipy.integrate as int
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as sdevs
import scipy.constants as con
from scipy.constants import physical_constants as pcon

T, alpha = np.loadtxt('data/data_alpha.txt', unpack = True)

def f(x, m, b):
    return m*1/x+b

params, cov = curve_fit(f, T, alpha)

m = params[0]
b = params[1]

dm = np.sqrt(cov[0][0])
db = np.sqrt(cov[1][1])

print(m)
print(dm)
print(b)
print(db)

x = np.linspace(70, 300, 1000)

plt.plot(1/T, alpha, 'rx', label='Messwerte')
plt.plot(1/x, f(x, m, b), 'k-', label='Regressionsgerade')
plt.grid()
plt.xlabel(r"$T^{-1}$ in $\mathrm{K}^{-1}$")
plt.ylabel(r"$\alpha\cdot 10^{-6}$ in $\mathrm{K}^{-1}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plots/plot_alpha.pdf")
plt.close()

mcu = 342
Mcu = 63.546
k = 140
V = 7.092*10**(-6)
dt, U, I, T1, T2 = np.loadtxt('data/data.txt', unpack = True)

cp = (Mcu/mcu) * (U * I * dt)/(T2-T1)

cv = cp - (9 * f(T2, m, b)**2 * k * V * T2)
