import numpy as np
import matplotlib.pyplot as plt
import sys

input_string = sys.argv[1]

name_string = input_string
if name_string.endswith('.asc'):
    name_string = name_string[:-4]

x, density = np.loadtxt(input_string, usecols=(0, 1), unpack=True)

fig, ax1 = plt.subplots(1, 1)

ax1.set_title(name_string)

ax1.plot(x, density, color='blue', label='Density')

ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'Density')
ax1.set_autoscaley_on(True)
ax1.set_autoscalex_on(True)

ax1.legend()

xmin = np.amin(x)
xmax = np.amax(x)

ax1.set_xlim([xmin, xmax])

plt.savefig(name_string + '.png', dpi=150.)
