import re
import numpy as np


def load_hexstrings_simple(file_path_input):
    """
    Loads only valid hex strings (0-9, A-F) from a text file.
    Lines starting with '*' or empty lines are ignored.
    Returns a list of hex strings.
    """
    hexstrings = []
    
    with open(file_path_input, 'r') as file:
        for line in file:
            line = line.strip()
            # Skip empty or commented lines
            if not line or line.startswith('*'):
                continue
            
            # Check if the line contains only valid hex characters (0-9, A-F, a-f)
            if all(c in '0123456789ABCDEFabcdef' for c in line):
                hexstrings.append(line.upper())  # Normalize to uppercase
    
    return hexstrings 

def load_header_hex(file_path_input, file_path_output):
    """
    Reads a hexadecimal file header and writes everything until '*END*' to output.
    Creates the output file if it does not already exist.
    """
    with open(file_path_input, 'r') as input_file, open(file_path_output, 'w') as output_file:
        for line in input_file:
            if line.startswith('*END*'):
                break
            if not line[0].isdigit() and not line[0].isalpha():
                output_file.write(line)

def load_header_XMLCON(file_path_input, file_path_output):
    """
    Extracts sensor blocks from XML and writes them as commented lines into the output file.
    Appends to the existing file.
    """
    with open(file_path_input, 'r') as xml_file, open(file_path_output, 'a') as cnv_file:
        content = xml_file.read()
        sensors = re.findall(r'<Sensor index="\d+" SensorID="\d+" >.*?</Sensor>', content, re.DOTALL)
        for sensor in sensors:
            for line in sensor.splitlines():
                cnv_file.write(f"# {line}\n")



def write_arrays_to_file(*arrays, file_path_output, decimals=6, total_width=12):
    """
    Writes multiple NumPy arrays side-by-side into a formatted text file.

    Parameters:
    -----------
    *arrays : variable number of numpy arrays
        All arrays must have equal length (otherwise trimmed to shortest).
    filename : str
        Name of the output file (no directory).
    output_dir : str
        Directory where the file will be written. Creates directory if necessary.
    decimals : int
        Number of decimal digits for floats.
    total_width : int
        Width of each printed value.
    """
    if not arrays:
        return
    
    # Ensure directory exists
    lengths = [len(arr) for arr in arrays]
    if len(set(lengths)) > 1:
        # Finde die minimale Länge und schneide alle Arrays darauf zu
        min_length = min(lengths)
        arrays = [arr[:min_length] for arr in arrays]
        print(f"Warning: Arrays have different lengths. Using minimum length: {min_length}")
    
    n_rows = len(arrays[0])
    n_cols = len(arrays)
    
    with open(file_path_output, 'a') as file:
        for i in range(n_rows):
            line = ""
            for j in range(n_cols):
                value = arrays[j][i]
                
                if isinstance(value, (int, np.integer)):
                    formatted_value = f"{value:>{total_width}d}"
                else:
                    formatted_value = f"{value:>{total_width}.{decimals}f}"
                
                line += formatted_value
                
                if j < n_cols - 1:
                    line += " "
            
            file.write(line + "\n")



def append_ctd_info(file_path_output, *arrays):
    """
    Appends CTD variable names and span (min/max) information to an existing file.
    The order of variables is determined by the order of the input arrays.
    """
    ordered_names = [
        ("prDM: Pressure, Digiquartz [db]", "pressure"),
        ("t090C: Temperature [ITS-90, deg C]", "temperatures1"),
        ("t190C: Temperature, 2 [ITS-90, deg C]", "temperatures2"),
        ("c0mS/cm: Conductivity [mS/cm]", "conductivitys1"),
        ("c1mS/cm: Conductivity, 2 [mS/cm]", "conductivitys2"),
        ("sbox0Mm/Kg: Oxygen, SBE 43 [umol/kg]", "oxygens1"),
        ("sbox1Mm/Kg: Oxygen, SBE 43, 2 [umol/kg]", "oxygens2"),
        ("sbox2Mm/Kg: Oxygen, SBE 43, 3 [umol/kg]", "oxygens3"),
        ("sbox3Mm/Kg: Oxygen, SBE 43, 4 [umol/kg]", "oxygens4"),
        ("flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]", "flurometer1"),
        ("turbWETntu0: Turbidity, WET Labs ECO [NTU]", "turbidity"),
        ("par: PAR/Irradiance, Biospherical/Licor", "par"),
        ("spar: SPAR, Biospherical/Licor", "spar"),
        ("timeS: Time, Elapsed [seconds]", "time_elapsed"),
        ("dz/dtM: Descent Rate [m/s]", "descent_rate_smooth2"),
        ("sal00: Salinity, Practical [PSU]", "salinitys1"),
        ("sal11: Salinity, Practical, 2 [PSU]", "salinitys2"),
        ("flag: 0.000e+00", "flag"),
    ]

    used = ordered_names[:len(arrays)]

    with open(file_path_output, "a") as f:

        # Write variable names
        for i, (label, _) in enumerate(used):
            f.write(f"# name {i} = {label}\n")

        # Writing spans
        for i, arr in enumerate(arrays):
            arr = np.asarray(arr, dtype=float)
            f.write(f"# span {i} = {np.nanmin(arr):12.4g}, {np.nanmax(arr):12.4g}\n")

