'''
Curve fitting with a customized tanh function
'''

import sys
import numpy as np
from scipy.optimize import curve_fit

#############################################################################
# Functions
#############################################################################
def hyperTan(x, a, b, c, d):
    return a * (b - np.tanh((x + d) / c))

# Sample data
def tanh_fit(x_data, y_data, initial_guess):
    # Perform the curve fitting
    fit_params, _ = curve_fit(hyperTan, x_data, y_data, p0=initial_guess)
    
    # Extract the optimized parameters
    a_opt, b_opt, c_opt, d_opt = fit_params
    
    # Generate the fitted y values using the optimized parameters
    y_fit = hyperTan(x_data, a_opt, b_opt, c_opt, d_opt)
    
    # R_squared analysis
    residuals = y_data - y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return a_opt, b_opt, c_opt, d_opt, y_fit, r_squared


#############################################################################
# User-input parameters
#############################################################################

if len(sys.argv) != 9:
    print("Usage: python3 this_code data_file a_ini b_ini c_ini d_ini skipped_bins_front skipped_bins_end output_file")
    exit()

# Read user input
profile = sys.argv[1]
initial_guesses        = [float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
skipped_bins_front     = int(sys.argv[6]) 
skipped_bins_end       = int(sys.argv[7])
output_file = sys.argv[8]


#############################################################################
# Read avg coordination number profile
#############################################################################
center_bin_coords = []
quantity = []

with open(profile, 'r') as file:
    for line in file:
        if not line.strip() or line.startswith("#"):
            continue
        data = line.strip().split()
        center_bin_coords.append(float(data[0]))
        quantity.append(float(data[1]))

center_bin_coords = np.array(center_bin_coords)
quantity = np.array(quantity)

# Hyperbolic tangent fitting 
x_for_fit = center_bin_coords[skipped_bins_front : -skipped_bins_end]
y_for_fit = quantity[skipped_bins_front : -skipped_bins_end]
a_opt, b_opt, c_opt, d_opt, y_fit, r_square = tanh_fit(x_for_fit, y_for_fit, initial_guesses)

# Write fitted tanh function
with open(output_file, "w") as f:
    f.write("# center_bin_coords tanh_fit \n")
    for x, y in zip(x_for_fit, y_fit):
        f.write("{:.4f} {:.4f}\n".format(x, y))
    f.write("# a_opt:             {:.4f}\n".format(a_opt))
    f.write("# b_opt:             {:.4f}\n".format(b_opt))
    f.write("# c_opt:             {:.4f}\n".format(c_opt))
    f.write("# d_opt:             {:.4f}\n".format(d_opt))
    f.write("# R_square:          {:.4f}\n\n".format(r_square))
    f.write("# Interfacial width: {:.4f}\n".format(abs(c_opt)))
    