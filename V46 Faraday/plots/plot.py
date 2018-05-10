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
l = l*(10**(-6))

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
plt.plot(z_fit, f1(z_fit, *params1), 'k-', label='Fit')
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
plt.xlabel(r'$\lambda^2 \:/\: \si{\meter^{2}}$')
plt.ylabel(r'$\frac{\theta}{d} \:/\: \si{\degree\meter^{-1}}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plots/plot_GaAs1.pdf')
plt.close()

########################################################

y1 = (th2-th1)*2*np.pi/360
y2 = (th3-th1)*2*np.pi/360

u1 = np.array([y1[0], y1[1], y1[2], y1[4], y1[5]])
v1 = np.array([y1[3], y1[6], y1[7]])
lu1 = np.array([l[0], l[1], l[2], l[4], l[5]])
lv1 = np.array([l[3], l[6], l[7]])
u2 = np.array([y2[1], y2[2], y2[3], y2[4]])
v2 = np.array([y2[0], y2[5], y2[6], y2[7]])
lu2 = np.array([l[1], l[2], l[3], l[4]])
lv2 = np.array([l[0], l[5], l[6], l[7]])
def f2(x, a, b):
	return a*x + b
params2, cov2 = curve_fit(f2, lu1**2, u1)
params3, cov3 = curve_fit(f2, lu2**2, u2)
l2 = np.linspace(10**(-12), 6.5*(10**(-12)))
plt.plot(lu1**2, u1, 'r.', label=r'$N = \SI{1.2e18}{cm^{-3}}$ (Differenz der Messwerte)')
plt.plot(lu2**2, u2, 'g.', label=r'$N = \SI{2.8e18}{cm^{-3}}$ (Differenz der Messwerte)')
plt.plot(lv1**2, v1, 'rx')
plt.plot(lv2**2, v2, 'gx')
plt.plot(l2, f2(l2, *params2), 'r-', label=r'$N = \SI{1.2e18}{cm^{-3}}$ (linearer Fit)')
plt.plot(l2, f2(l2, *params3), 'g-', label=r'$N = \SI{2.8e18}{cm^{-3}}$ (linearer Fit)')
plt.xlabel(r'$\lambda^2 \:/\: \si{\meter^{2}}$')
plt.ylabel(r'$\frac{\theta}{d} \:/\: \si{rad\,\meter^{-1}}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plots/plot_GaAs2.pdf')
plt.close()

print('Fit für N=1.2e18:')
print('Parameter: ', params2, '\nFehler: ', np.sqrt(np.diag(cov2)))
print('Fit für N=2.8e18:')
print('Parameter: ', params3, '\nFehler: ', np.sqrt(np.diag(cov3)))

# x = np.linspace(0, 10, 1000)
# y = x ** np.sin(x)
#
# plt.subplot(1, 2, 1)
# plt.plot(x, y, label='Kurve')
# plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
# plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
# plt.legend(loc='best')
#
# plt.subplot(1, 2, 2)
# plt.plot(x, y, label='Kurve')
# plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
# plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
# plt.legend(loc='best')
