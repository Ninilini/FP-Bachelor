# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat

z, B = np.genfromtxt('data/B.txt', unpack=True)
l, th1, th2, th3 = np.genfromtxt('data/GaAs.txt', unpack=True)
d1 = 0.00511
d2 = 0.00136
d3 = 0.001296

# Normierung auf Einheitslänge:
th1 = th1/d1
th2 = th2/d2
th3 = th3/d3

####################################################
def f1(x, a, b, c):
	return a*(x**2) + b*x + c
params1, cov1 = curve_fit(f1, z, B)
z_fit = np.linspace(-20, 20)
plt.plot(z, B, 'b.', label='Messwerte')
plt.plot(z_fit, f1(z_fit, *params1), 'k-', label='Regression')
plt.xlabel(r'$z \:/\: \si{mm}$')
plt.ylabel(r'$B \:/\: \si{mT}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plots/plot_B.pdf')
plt.close()
print('Magnetfeld-Fit:')
print('Parameter: ', params1, '\nFehler: ', np.sqrt(np.diag(cov1)))

#####################################################

plt.plot(l**2, th1, 'b.', label='hochreines GaAs')
plt.plot(l**2, th2, 'r.', label=r'n-dotiertes GaAs mit $N = \SI{1.2e18}{cm^{-3}}$')
plt.plot(l**2, th3, 'g.', label=r'n-dotiertes GaAs mit $N = \SI{2.8e18}{cm^{-3}}$')
plt.xlabel(r'$\lambda^2 \:/\: \si{\micro\meter^{2}}$')
plt.ylabel(r'$\frac{\theta}{d} \:/\: \si{\degree\meter^{-1}}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plots/plot_GaAs1.pdf')
plt.close()

########################################################

y1 = (th2-th1)*2*np.pi/360
y2 = (th3-th1)*2*np.pi/360

def f2(x, a, b):
	return a*x + b
params2, cov2 = curve_fit(f2, l**2, y1)
params3, cov3 = curve_fit(f2, l**2, y2)
l2 = np.linspace(1, 6.5)
plt.plot(l**2, y1, 'r.', label=r'$N = \SI{1.2e18}{cm^{-3}}$ (Differenz der Messwerte)')
plt.plot(l**2, y2, 'g.', label=r'$N = \SI{2.8e18}{cm^{-3}}$ (Differenz der Messwerte)')
plt.plot(l2, f2(l2, *params2), 'r-', label=r'$N = \SI{1.2e18}{cm^{-3}}$ (lineare Regression)')
plt.plot(l2, f2(l2, *params3), 'g-', label=r'$N = \SI{2.8e18}{cm^{-3}}$ (lineare Regression)')
plt.xlabel(r'$\lambda^2 \:/\: \si{\micro\meter^{2}}$')
plt.ylabel(r'$\frac{\theta}{d} \:/\: \si{rad\,\meter^{-1}}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plots/plot_GaAs2.pdf')
plt.close()

print('Fit für N=1.2e18:')
print('Parameter: ', params2, '\nFehler: ', np.sqrt(np.diag(cov2)))
print('Fit für N=2.8e18:')
print('Parameter: ', params3, '\nFehler: ', np.sqrt(np.diag(cov3)))

########################################################

a = (1.602176462**3)/(8*(np.pi**2)*8.854187817*(2.99792458**3))
N1 = 1.2
N2 = 2.8
n = 3.354
m_e = 0.910938188
B_max = ufloat(408, 10)
A1 = ufloat(8.75936855, 9.34901912)
A2 = ufloat(9.35898797, 6.57753851)
print('N = 1.2: m*=', unp.sqrt(a*N1*B_max/(n*A1))*m_e, ' m_e')
print('N = 2.8: m*=', unp.sqrt(a*N2*B_max/(n*A2))*m_e, ' m_e')

# Magnetfeld-Fit:
# Parameter:  [ -1.31382591   2.99072852 405.30603237]
# Fehler:  [ 0.07652639  0.90317205 10.12785969]

# Fit für N=1.2e18:
# Parameter:  [8.75936855 9.66626756]
# Fehler:  [ 9.34901912 36.69455966]
# Fit für N=2.8e18:
# Parameter:  [ 9.35898797 13.2345978 ]
# Fehler:  [ 6.57753851 25.81659948]

# N = 1.2: m*= 0.055+/-0.029  m_e
# N = 2.8: m*= 0.081+/-0.029  m_e
