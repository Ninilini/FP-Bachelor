# coding: utf-8
                       
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython import embed
counts = np.genfromtxt('data/Alu_leer/alu_l_2.Spe', skip_header=12,
                       skip_footer=17)
                       
leer = counts
error_leer = np.sqrt(counts)
x = np.linspace(0,1922056,511)
plt.bar(x, leer, yerr=error_leer)
plt.xlim(0, 250)
plt.title('Messung des leeren Würfels')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('Alu_leer.pdf')
plt.close()
