# coding=utf-8
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as err
import scipy.integrate as int
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

m = ufloat(m, dm)
b = ufloat(b, db)

print('m = ', m)
print('b = ', b)


x = np.linspace(70, 300, 1000)

plt.plot(1/T, alpha, 'rx', label='Messwerte')
plt.plot(1/x, noms(f(x, m, b)), 'k-', label='Regressionsgerade')
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
dt2, U2, I2, R12, R22, T12, T22 = np.loadtxt('data/data.txt', unpack = True)
dt = unumpy.uarray(dt2, 1)
U = unumpy.uarray(U2, 0.01)
I = unumpy.uarray(I2, 0.1)
R1 = unumpy.uarray(R12, 0.1)
R2 = unumpy.uarray(R22, 0.1)

def t(R):
    return 0.00134*R**2 + 2.296*R + 30.13

T1 = t(R1)
T2 = t(R2)

print('''
alpha:
--------------------------
''',f(T2, m, b))

print('''
T1:
--------------------------
''',T1)

print('''
T2:
--------------------------
''',T2)

print('''
dT:
--------------------------
''',T2-T1)

cp = (Mcu/mcu) * (U * I/1000 * dt)/(T2-T1)

print('''
cp:
--------------------------
''',cp)

cv = cp - (9 * (f(T2, m, b)*10**(-5))**2 * k * V * T2)

print('''
cv:
--------------------------
''',cv)

cpf = err(cp)
cvf = err(cv)
cp1 = noms(cp)
cv1 = noms(cv)

# plt.errorbar(noms(T2),cp1,yerr=cpf,fmt='bx',label='Werte für $C_{\mathrm{P}}$')
plt.errorbar(noms(T2),cv1,yerr=cvf,fmt='rx',label='Werte für $C_{\mathrm{V}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.grid()
plt.xlabel(r"$T$ in $\mathrm{K}$")
plt.ylabel(r"$C_{\mathrm{P,V}}$ in $\mathrm{JK}^{-1}\mathrm{mol}^{-1}$")
plt.tight_layout()
plt.savefig("plots/plot_Cv.pdf")
plt.close()
