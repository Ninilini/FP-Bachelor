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

def linear(x, m, b):
    return m*x + b

t = np.array([-20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
n = np.array([225, 260, 244, 271, 333, 308, 360, 318, 382, 345, 360, 304, 317, 357, 307, 333, 326, 291, 247,
172, 101])

hoehe = np.mean(n[4:16])
links = n[0:5]
rechts = n[16:]
dt_rechts = t[16:]

# linearer Fit links und rechts
params_links, cov_links = curve_fit(linear, t[0:5], links)
errors_links = np.sqrt(np.diag(cov_links))
m = ufloat(params_links[0], errors_links[0])
b = ufloat(params_links[1], errors_links[1])
print("Steigung links: ", m)
print("y-Achsenabschnitt rechts: ", b)
print("Halbe Höhe: ", hoehe/2)

params_rechts, cov_rechts = curve_fit(linear, dt_rechts, rechts)
errors_rechts = np.sqrt(np.diag(cov_rechts))
m = ufloat(params_rechts[0], errors_rechts[0])
b = ufloat(params_rechts[1], errors_rechts[1])
print("Steigung rechts: ", m)
print("y-Achsenabschnitt rechts: ", b)

# Berechnen des Schnittpunktes
x_links = np.linspace(-25, -17)
x_rechts = np.linspace(6, 20)

links_w = linear(x_links, *params_links)
rechts_w = linear(x_rechts, *params_rechts)

plt.figure(1)
plt.ylabel(r"$N(t) \, / \, (\mathrm{s})^{-1}$")
plt.xlabel(r"$\mathrm{d}t \, / \, \mathrm{ns}$")
plt.errorbar(t, n, yerr=np.sqrt(n), fmt='kx', label="Messwerte")
plt.plot(x_links, linear(x_links, *params_links), 'r',
         label="Regression links")
plt.plot(x_rechts, linear(x_rechts, *params_rechts), 'r',
         label="Regression rechts")
plt.axhline(y=hoehe, xmin=0.15, xmax=0.75, label="Plateau")
plt.axhline(y=hoehe/2, xmin=0.11, xmax=0.87, color="green",
            label="Halbwertsbreite")
plt.ylim(0, 410)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plateau.pdf")
plt.clf()
