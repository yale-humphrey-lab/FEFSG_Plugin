import numpy as np
import os
import shutil
import xml.etree.ElementTree as ET


# --------------------------
# Parameter Ranges
# --------------------------
# These arrays define the parameter values to sweep through.
elastin_injury_rate_vals = np.array([0.01, 0.05])  # Example: Injury rate values
elastin_injury_vals = np.array([0.40, 0.80])       # Example: Elastin injury values
Kp_injury_vals = np.array([0.0, 0.82])             # Production mechanosensitivity injury values
Kd_injury_vals = np.array([0.0, 0.82])             # Degradation mechanosensitivity injury values


def update_feb(elastin_val, kp_val, kd_val, rate_val, filename):
    """
    Update parameter values in a FEBio XML file.

    Args:
        elastin_val (float): New elastin injury value.
        kp_val (float): New production mechanosensitivity injury value.
        kd_val (float): New degradation mechanosensitivity injury value.
        rate_val (float): New injury rate value for the load controller.
        filename (str): Path to the FEBio file to modify.

    Returns:
        None. The file is updated in place.
    """
    # Parse XML
    tree = ET.parse(filename)
    root = tree.getroot()

    # Update material parameters
    for material in root.iter('material'):
        # Elastin injury value
        elastin = material.find('elastin_injury_val')
        if elastin is not None and elastin.text:
            elastin.text = elastin.text.replace("VVV", str(elastin_val))

        # Production mechanosensitivity
        K_p = material.find('production_mechanosensitivity_injury_val')
        if K_p is not None and K_p.text:
            K_p.text = K_p.text.replace("VVV", str(kp_val))

        # Degradation mechanosensitivity
        K_d = material.find('degradation_mechanosensitivity_injury_val')
        if K_d is not None and K_d.text:
            K_d.text = K_d.text.replace("VVV", str(kd_val))

    # Update load controller math expression (id="2")
    for controller in root.findall('.//load_controller'):
        if controller.get('id') == "2":
            math_elem = controller.find('math')
            if math_elem is not None and math_elem.text:
                math_elem.text = math_elem.text.replace("VVV", str(rate_val))

    # Save updated XML
    tree.write(filename, encoding='utf-8', xml_declaration=True)


def create_parameter_sweep_folders():
    """
    Create simulation folders for all parameter combinations.

    Each combination of elastin, Kp, Kd, and rate values will get its own folder
    containing a copy of the baseline "FEFSG" folder, with parameters updated
    in the FEBio input file.

    Returns:
        None
    """
    # Generate all parameter combinations
    param_combinations = [
        (e_val, kp_val, kd_val, rate_val)
        for e_val in elastin_injury_vals
        for kp_val in Kp_injury_vals
        for kd_val in Kd_injury_vals
        for rate_val in elastin_injury_rate_vals
    ]

    # Loop through each combination
    for e_val, kp_val, kd_val, rate_val in param_combinations:
        folder_name = (
            f"FEFSG_e{e_val:.2f}_Kp{kp_val:.3f}_Kd{kd_val:.3f}_ke{rate_val:.3f}"
        )

        # Only create if folder doesn't exist
        if not os.path.isdir(folder_name):
            print(f"Creating: {folder_name}")
            shutil.copytree("FEFSG", folder_name)

            # Update the FEBio XML file inside the new folder
            feb_file = os.path.join(folder_name, "TAA_aneurysm.feb")
            update_feb(e_val, kp_val, kd_val, rate_val, feb_file)


if __name__ == "__main__":
    create_parameter_sweep_folders()