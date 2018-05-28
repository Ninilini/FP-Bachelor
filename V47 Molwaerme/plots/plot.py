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

# x = np.linspace(70, 300, 1000)
#
# plt.plot(1/T, alpha, 'rx', label='Messwerte')
# plt.plot(1/x, f(x, m, b), 'k-', label='Regressionsgerade')
# plt.grid()
# plt.xlabel(r"$T^{-1}$ in $\mathrm{K}^{-1}$")
# plt.ylabel(r"$\alpha\cdot 10^{-6}$ in $\mathrm{K}^{-1}$")
# plt.legend(loc='best')
# plt.tight_layout()
# plt.savefig("plots/plot_alpha.pdf")
# plt.close()

m = ufloat(m, dm)
b = ufloat(b, db)

print('m = ', m)
print('b = ', b)

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

D1 = T2[1]*3.9
D2 = T2[2]*3.1
D3 = T2[3]*3
D4 = T2[4]*2.2
D5 = T2[5]*2.2
D6 = T2[6]*2.7
D7 = T2[7]*0
D8 = T2[8]*3.6
D9 = T2[9]*3.9
D10 = T2[10]*4.0

DM = 1/10 * (D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10)

print(D1)
print(D2)
print(D3)
print(D4)
print(D5)
print(D6)
print(D7)
print(D8)
print(D9)
print(D10)

print(DM)
w = con.k/con.hbar * DM
print(w)
# # plt.errorbar(noms(T2),cp1,yerr=cpf,fmt='bx',label='Werte für $C_{\mathrm{P}}$')
# plt.errorbar(noms(T2),cv1,yerr=cvf,fmt='rx',label='Werte für $C_{\mathrm{V}}$')
# plt.legend(loc='best')
# plt.tight_layout()
# plt.grid()
# plt.xlabel(r"$T$ in $\mathrm{K}$")
# plt.ylabel(r"$C_{\mathrm{V}}$ in $\mathrm{JK}^{-1}\mathrm{mol}^{-1}$")
# plt.tight_layout()
# plt.savefig("plots/plot_Cv.pdf")
# plt.close()
