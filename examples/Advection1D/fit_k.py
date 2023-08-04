import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Given lists of variables
errors = [ 0.676131, 0.66194, 0.64546, 0.635339, 0.629921, 0.62718]
dx_values = [0.03125, 0.015625, 0.0078125, 0.00390625, 0.00195312, 0.000976562]

# Define the power law function
def power_law(dx, A, k):
    return A * (dx ** k)

# Perform curve fitting
params, _ = curve_fit(power_law, dx_values, errors)

# Extract the fitted parameters
A = params[0]
k = params[1]

# Print the values of A and k
print("A:", A)
print("k:", k)

# Plotting the results
plt.scatter(dx_values, errors, color='blue', label='Data')
plt.plot(dx_values, power_law(dx_values, A, k), color='red', label='Fit')
plt.xlabel('dx_values')
plt.ylabel('Errors')
plt.legend()

# Save the plot as a PNG file
plt.savefig('fit_plot.png')

# Show the plot
plt.show()

