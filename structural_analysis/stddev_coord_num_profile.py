'''

Calculate the standard deviation of the coordination number profile and fit it with tanh function

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

if len(sys.argv) != 11:
    print("Usage: python3 this_code coord_num_profile cutoff binWidth nbins a_ini b_ini c_ini d_ini skipped_bins_front skipped_bins_end")
    exit()

# Read user input
coord_num_profile_file = sys.argv[1]
cutoff_distance        = float(sys.argv[2]) 
binWidth               = float(sys.argv[3])
num_bins               = int(sys.argv[4])
initial_guesses        = [float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])]
skipped_bins_front     = int(sys.argv[9]) 
skipped_bins_end       = int(sys.argv[10])


#############################################################################
# Read coordination number profile
#############################################################################
coord_num_data     = {}                    # See sample at the end of the section
tstep   = None

with open(coord_num_profile_file, 'r') as file:
    for line in file:
        line = line.strip()  # Remove extra spaces at the beginning and end of each line

        if not line:         # skip blank lines
            continue

        if line.startswith("# timestep:"):
            tstep                 = int(line.split(":")[1].strip())
            coord_num_data[tstep] = {
                    "CENTER_X_COORDS": [],
                    "COORD_NUMS": []
            }
        
            center_xs  = []
            coord_nums = []
            
            next(file) # This line is "# center_x_coord  coord_num"
            
            for bin in range(num_bins):
                data = next(file).strip().split()
                
                center_x  = float(data[0])
                coord_num = float(data[1])

                center_xs.append(center_x)
                coord_nums.append(coord_num)

            coord_num_data[tstep]["CENTER_X_COORDS"]  = center_xs
            coord_num_data[tstep]["COORD_NUMS"]       = coord_nums

######## Sample data for demonstration ########
# data = {
#     150000: {
#         "CENTER_X_COORDS": [1, 2, 3, 4, ...],
#         "COORD_NUMS": [3.9, 3.9, ...]
#     },
#     
#     150050: {
#         "CENTER_X_COORDS": [1, 2, 3, 4, ...],
#         "COORD_NUMS": [3.9, 3.9, ...]
#     }
# }

center_x_coords   = []
coord_num_profile = []

for tstep, info in coord_num_data.items():
    center_x_coords.append(info["CENTER_X_COORDS"])
    coord_num_profile.append(info["COORD_NUMS"])

center_x_coords   = np.array(center_x_coords) 
coord_num_profile = np.array(coord_num_profile)
    
# calculate average center_x_coord
avg_center_x_coords = np.mean(center_x_coords, axis=0)

# calculate standard deviation of corodination number in each bin
stddev_coord_num_profile = np.std(coord_num_profile, axis=0)

# Hyperbolic tangent fitting 
x_for_fit = avg_center_x_coords[skipped_bins_front : -skipped_bins_end]
y_for_fit = stddev_coord_num_profile[skipped_bins_front : -skipped_bins_end]
a_opt, b_opt, c_opt, d_opt, y_fit, r_square = tanh_fit(x_for_fit, y_for_fit, initial_guesses)

with open("stddev_coord_num_profile.dat", "w") as f:
    f.write("# center_x_coords stddev_coord_num \n")
    for bin, stddev_coord_num in zip(avg_center_x_coords, stddev_coord_num_profile):
        f.write("{:.4f} {:.4f}\n".format(bin, stddev_coord_num))

with open("stddev_coord_num_profile_tanh_fit.dat", "w") as f:
    f.write("# center_x_coords stddev_coord_num_tanh_fit \n")
    for xfit, yfit in zip(x_for_fit, y_fit):
        f.write("{:.4f} {:.4f}\n".format(xfit, yfit))

with open("stddev_coord_num_profile_tanh_fit_paras.dat", "w") as f:
    f.write("# a_opt b_opt c_opt d_opt r_square \n")
    f.write("{:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(a_opt, b_opt, c_opt, d_opt, r_square))