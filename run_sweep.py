import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import xml.etree.ElementTree as ET

elastin_injury_rate_vals =  np.array([0.01, 0.05]) #np.array([0.01, 0.015, 0.025, 0.05])
elastin_injury_vals = np.array([0.40, 0.80]) #np.array([0.0, 0.12, 0.36, 0.60, 0.72])
Kp_injury_vals = np.array([0.0, 0.82])
Kd_injury_vals = np.array([0.0, 0.82])

def run_parameter_sweep_folders():
    args_list = [
        (e_val, kp_val, kd_val, ke_val)
        for e_val in elastin_injury_vals
        for kp_val in Kp_injury_vals
        for kd_val in Kd_injury_vals
        for ke_val in elastin_injury_rate_vals
    ]

    for arg in args_list:
        print(arg)
        e_val, kp_val, kd_val, ke_val = arg
        os.chdir(f"FEFSG_e{e_val:.2f}_Kp{kp_val:.3f}_Kd{kd_val:.3f}_ke{ke_val:.3f}")
        os.system("cp ../FEFSG/expanse_driver_script .")
        os.system('sbatch expanse_driver_script')
        #os.system('python create_vessel.py TAA_aneurysm.feb')
        os.chdir("../")
    return

# Your plot_heatmap function remains unchanged

if __name__ == "__main__":
    run_parameter_sweep_folders()
    #plot_heatmap(elastin_injury_vals, elastin_injury_rate_vals)  # Adjust x-label if needed
