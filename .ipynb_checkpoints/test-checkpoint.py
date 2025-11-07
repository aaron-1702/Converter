import numpy as np
from typing import Literal, List
import gsw
from scipy import stats
from math import floor

# Konstanten
DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44.660
ITS90_TO_IPTS68 = 1.00024
UMNO3_TO_MGNL = 0.014007
R = 8.3144621
F = 96485.365

# Kalibrierungskoeffizienten (hier Ihre bestehenden Coefficient-Klassen und Werte einfügen)
# ...

class TemperatureCoefficients:
    """
    :param a0: coefficient
    :param a1: coefficient
    :param a2: coefficient
    :param a3: coefficient
    """

    a0: float
    a1: float
    a2: float
    a3: float

class PressureCoefficients:
    """
    :param pa0: coefficient
    :param pa1: coefficient
    :param pa2: coefficient
    :param ptca0: coefficient
    :param ptca1: coefficient
    :param ptca2: coefficient
    :param ptcb0: coefficient
    :param ptcb1: coefficient
    :param ptcb2: coefficient
    :param ptempa0: coefficient
    :param ptempa1: coefficient
    :param ptempa2: coefficient
    """

    pa0: float
    pa1: float
    pa2: float
    ptca0: float
    ptca1: float
    ptca2: float
    ptcb0: float
    ptcb1: float
    ptcb2: float
    ptempa0: float
    ptempa1: float
    ptempa2: float

class ConductivityCoefficients:
    """
    :param g: coefficient
    :param h: coefficient
    :param i: coefficient
    :param j: coefficient
    :param cpcor: compressibility coefficient
    :param ctcor: coefficient
    :param wbotc: bridge oscillator temperature coefficient see the
        37 Manual: https://www.seabird.com/asset-get.download.jsa?id=54627862348
    """

    g: float
    h: float
    i: float
    j: float
    cpcor: float
    ctcor: float
    wbotc: float

class Oxygen43Coefficients:
    """
    :param soc: linear scaling coefficient
    :param v_offset: voltage at zero oxygen signal
    :param tau_20: sensor time constant tau(T,P) at 20 C, 1 atmosphere, 0 PSU;
        slope term in calculation of tau(T,P)
    :param a: coefficient
    :param b: coefficient
    :param c: coefficient
    :param e: coefficient
    :param d0: tau(T,P) coefficient
    :param d1: tau(T,P) coefficient
    :param d2: tau(T,P) coefficient
    :param h1: hysteresis correction coefficient
    :param h2: hysteresis correction coefficient
    :param h3: hysteresis correction coefficient
    """

    soc: float
    v_offset: float
    tau_20: float
    a: float
    b: float
    c: float
    e: float
    d0: float  # TODO: the sensor lists a d0 coefficient, but it doesn't seem to be used?
    d1: float
    d2: float
    h1: float
    h2: float
    h3: float

temp1_coefs = TemperatureCoefficients()
temp1_coefs.a0 = 4.33459406e-003
temp1_coefs.a1 = 6.28735680e-004
temp1_coefs.a2 = 2.05782221e-005
temp1_coefs.a3 = 1.68624541e-006

temp2_coefs = TemperatureCoefficients()
temp2_coefs.a0 = 4.34471154e-003
temp2_coefs.a1 = 6.38146193e-004
temp2_coefs.a2 = 2.04525401e-005
temp2_coefs.a3 = 1.60292340e-006


pres_coefs = PressureCoefficients()
pres_coefs.pa0 = -4.212464e+004
pres_coefs.pa1 = -1.972709e-001
pres_coefs.pa2 = 1.365930e-002

pres_coefs.ptca0 = 3.604100e-002
pres_coefs.ptca1 = 0.000000e+000
pres_coefs.ptca2 = 3.046128e+001

pres_coefs.ptcb0 = -3.791930e-004
pres_coefs.ptcb1 = 4.202260e-006
pres_coefs.ptcb2 = 3.386190e-00

pres_coefs.ptempa0 = 0.000000e+000
pres_coefs.ptempa1 = 1.279410e-002
pres_coefs.ptempa2 = -9.215740e+000


cond1_coefs = ConductivityCoefficients()
cond1_coefs.g = -9.94900904e+000
cond1_coefs.h = 1.42164741e+000
cond1_coefs.i = 1.26008747e-004
cond1_coefs.j = 6.16820694e-005
cond1_coefs.cpcor = -9.57000000e-008
cond1_coefs.ctcor = 3.2500e-006
cond1_coefs.wbotc = 0.00000000e+000

cond2_coefs = ConductivityCoefficients()
cond2_coefs.g = -9.77824825e+000
cond2_coefs.h = 1.31877668e+000
cond2_coefs.i = 2.26556605e-004
cond2_coefs.j = 6.59127902e-005
cond2_coefs.cpcor = -9.57000000e-008
cond2_coefs.ctcor = 3.2500e-006
cond2_coefs.wbotc = 0.00000000e+000


oxy1_coefs = Oxygen43Coefficients()
oxy1_coefs.soc = 4.7612e-001
oxy1_coefs.v_offset = -0.4891
oxy1_coefs.tau_20 = 1.7400
oxy1_coefs.a = -3.2939e-003
oxy1_coefs.b = 1.4617e-004
oxy1_coefs.c = -2.1408e-006
oxy1_coefs.e = 3.6000e-002
oxy1_coefs.d0 = 2.5826e+000
oxy1_coefs.d1 = 1.92634e-004
oxy1_coefs.d2 = -4.64803e-002
oxy1_coefs.h1 = -3.3000e-002
oxy1_coefs.h2 = 5.0000e+003
oxy1_coefs.h3 = 1.4500e+003

oxy2_coefs = Oxygen43Coefficients()
oxy2_coefs.soc = 5.4771e-001
oxy2_coefs.v_offset = -0.5020
oxy2_coefs.tau_20 = 1.2900
oxy2_coefs.a = -3.9126e-003
oxy2_coefs.b = 1.5415e-004
oxy2_coefs.c = -2.2575e-006
oxy2_coefs.e = 3.6000e-002
oxy2_coefs.d0 = 2.5826e+000
oxy2_coefs.d1 = 1.92634e-004
oxy2_coefs.d2 = -4.64803e-002
oxy2_coefs.h1 = -3.3000e-002
oxy2_coefs.h2 = 5.0000e+003
oxy2_coefs.h3 = 1.4500e+003




def hex_to_voltage(hex_values: np.ndarray):
    """Konvertiert Hex-Werte zu Spannungswerten (array-fähig)"""
    # Hex zu Dezimal
    N = np.array([int(hex_val, 16) for hex_val in hex_values])
    # Konvertiere zu Volt
    V = 5 * (1 - (N / 4095))
    return V

def extract_bytes_from_hexstrings(hex_strings: List[str], byte_positions: List[int]):
    """Extrahiert Bytes aus mehreren Hexstrings an bestimmten Positionen"""
    results = []
    for hex_string in hex_strings:
        bytes_list = [hex_string[i:i+2] for i in range(0, len(hex_string), 2)]
        extracted_bytes = [bytes_list[pos] for pos in byte_positions if pos < len(bytes_list)]
        results.append(extracted_bytes)
    return np.array(results)

def convert_temperature_array(
    temperature_counts_in: np.ndarray,
    coefs: TemperatureCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
    use_mv_r: bool = False,
):
    """Array-fähige Temperaturkonvertierung"""
    
    if use_mv_r:
        mv = (temperature_counts_in - 524288) / 1.6e007
        r = (mv * 2.900e009 + 1.024e008) / (2.048e004 - mv * 2.0e005)
        temperature_counts = r
    else:
        temperature_counts = temperature_counts_in

    log_t = np.log(temperature_counts)
    temperature = (
        1 / (coefs.a0 + coefs.a1 * log_t + coefs.a2 * log_t**2 + coefs.a3 * log_t**3)
    ) - KELVIN_OFFSET_0C

    if standard == "IPTS68":
        temperature *= ITS90_TO_IPTS68
    if units == "F":
        temperature = temperature * 9 / 5 + 32

    return temperature

def convert_pressure_array(
    pressure_counts: np.ndarray,
    compensation_voltages: np.ndarray,
    coefs: PressureCoefficients,
    units: Literal["dbar", "psia"] = "psia",
):
    """Array-fähige Druckkonvertierung"""
    sea_level_pressure = 14.7

    t = (
        coefs.ptempa0
        + coefs.ptempa1 * compensation_voltages
        + coefs.ptempa2 * compensation_voltages**2
    )
    x = pressure_counts - coefs.ptca0 - coefs.ptca1 * t - coefs.ptca2 * t**2
    n = x * coefs.ptcb0 / (coefs.ptcb0 + coefs.ptcb1 * t + coefs.ptcb2 * t**2)
    pressure = coefs.pa0 + coefs.pa1 * n + coefs.pa2 * n**2 - sea_level_pressure

    if units == "dbar":
        pressure *= PSI_TO_DBAR

    return pressure

def convert_conductivity_array(
    conductivity_counts: np.ndarray,
    temperatures: np.ndarray,
    pressures: np.ndarray,
    coefs: ConductivityCoefficients,
):
    """Array-fähige Leitfähigkeitskonvertierung"""
    f = conductivity_counts * np.sqrt(1 + coefs.wbotc * temperatures) / 1000
    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperatures + coefs.cpcor * pressures
    return numerator / denominator


def convert_sbe43_oxygen(
    voltage: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: Oxygen43Coefficients,
    apply_tau_correction: bool = False,
    apply_hysteresis_correction: bool = False,
    window_size: float = 1,
    sample_interval: float = 1,
):
    """Returns the data after converting it to ml/l.

    voltage is expected to be in volts, temperature in deg c, pressure
    in dbar, and salinity in practical salinity (PSU). All equation
    information comes from the June 2013 revision of the SBE43 manual

    :param voltage: SBE43 voltage
    :param temperature: temperature value converted to deg C
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs: calibration coefficients for the SBE43 sensor
    :param apply_tau_correction: whether or not to run tau correction
    :param apply_hysteresis_correction: whether or not to run hysteresis
        correction
    :param window_size: size of the window to use for tau correction, if
        applicable, in seconds
    :param sample_interval: sample rate of the data to be used for tau
        correction, if applicable. In seconds.

    :return: converted Oxygen values, in ml/l
    """
    # start with all 0 for the dvdt
    dvdt_values = np.zeros(len(voltage))
    if apply_tau_correction:
        # Calculates how many scans to have on either side of our median
        # point, accounting for going out of index bounds
        scans_per_side = floor(window_size / 2 / sample_interval)
        for i in range(scans_per_side, len(voltage) - scans_per_side):
            ox_subset = voltage[i - scans_per_side : i + scans_per_side + 1]

            time_subset = np.arange(
                0, len(ox_subset) * sample_interval, sample_interval, dtype=float
            )

            result = stats.linregress(time_subset, ox_subset)

            dvdt_values[i] = result.slope

    correct_ox_voltages = voltage.copy()
    if apply_hysteresis_correction:
        # Hysteresis starts at 1 because 0 can't be corrected
        for i in range(1, len(correct_ox_voltages)):
            # All Equation info from APPLICATION NOTE NO. 64-3
            d = 1 + coefs.h1 * (np.exp(pressure[i] / coefs.h2) - 1)
            c = np.exp(-1 * sample_interval / coefs.h3)
            ox_volts = correct_ox_voltages[i] + coefs.v_offset

            prev_ox_volts_new = correct_ox_voltages[i - 1] + coefs.v_offset
            ox_volts_new = ((ox_volts + prev_ox_volts_new * c * d) - (prev_ox_volts_new * c)) / d
            ox_volts_final = ox_volts_new - coefs.v_offset
            correct_ox_voltages[i] = ox_volts_final

    oxygen = _convert_sbe43_oxygen(
        correct_ox_voltages,
        temperature,
        pressure,
        salinity,
        coefs,
        dvdt_values,
    )
    return oxygen

def _convert_sbe43_oxygen(
    voltage: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: Oxygen43Coefficients,
    dvdt_value: np.ndarray,
):
    """Returns the data after converting it to ml/l.

    voltage is expected to be in volts, temperature in deg c, pressure
    in dbar, and salinity in practical salinity (PSU). All equation
    information comes from the June 2013 revision of the SBE43 manual.
    Expects that hysteresis correction is already performed on the
    incoming voltage, if desired.

    :param voltage: SBE43 voltage
    :param temperature: temperature value converted to deg C
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs: calibration coefficients for the SBE43 sensor
    :param dvdt_value: derivative value of voltage with respect to time
        at this point. Expected to be 0 if not using Tau correction

    :return: converted Oxygen value, in ml/l
    """

    # Oxygen Solubility equation constants, From SBE43 Manual Appendix A
    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.0501
    a3 = 4.94457
    a4 = -0.256847
    a5 = 3.88767
    b0 = -0.00624523
    b1 = -0.00737614
    b2 = -0.010341
    b3 = -0.00817083
    c0 = -0.000000488682

    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    a_term = a0 + a1 * ts + a2 * ts**2 + a3 * ts**3 + a4 * ts**4 + a5 * ts**5
    b_term = salinity * (b0 + b1 * ts + b2 * ts**2 + b3 * ts**3)
    c_term = c0 * salinity**2
    solubility = np.exp(a_term + b_term + c_term)

    # Tau correction
    tau = coefs.tau_20 * np.exp(coefs.d1 * pressure + coefs.d2 * (temperature - 20)) * dvdt_value

    soc_term = coefs.soc * (voltage + coefs.v_offset + tau)
    temp_term = 1.0 + coefs.a * temperature + coefs.b * temperature**2 + coefs.c * temperature**3
    oxygen = (
        soc_term
        * solubility
        * temp_term
        * np.exp((coefs.e * pressure) / (temperature + KELVIN_OFFSET_0C))
    )
    return oxygen

def convert_oxygen_to_umol_per_kg(ox_values: np.ndarray, potential_density: np.ndarray):
    """Converts given oxygen values to milligrams/kg.

    Note: Sigma-Theta is expected to be calculated via gsw_sigma0,
    meaning is it technically potential density anomaly. Calculating
    using gsw_rho(SA, CT, p_ref = 0) results in actual potential
    density, but this function already does the converison, so values
    will need to have 1000 subtracted from them before being passed into
    this function. The function is done this way to stay matching to
    Application Note 64, but the results of either method are identical.

    :param ox_values: oxygen values, already converted to ml/L
    :param potential_density: potential density (sigma-theta) values.
        Expected to be the same length as ox_values

    :return: oxygen values converted to milligrams/Liter
    """

    oxygen_umolkg = (ox_values * OXYGEN_MLPERL_TO_UMOLPERKG) / (potential_density + 1000)
    return oxygen_umolkg

def potential_density_from_t_s_p(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
    reference_pressure=0.0,
):
    """Derive potential density from measured temperature, salinity, and
    pressure.

    :param temperature: Measure temperature, in degrees C
    :param salinity: Measured salinity, in practical salinity units
    :param pressure: Measured pressure, in decibars
    :param lon: Longitude
    :param lat: Latitude
    :param reference_pressure: Reference pressure in decibars. Defaults
        to 0.0.

    :return: Potential density in kg/m^3
    """

    absolute_salinity = gsw.SA_from_SP(salinity, pressure, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature, pressure)
    potential_density = (
        gsw.rho(absolute_salinity, conservative_temperature, reference_pressure) - 1000
    )
    return potential_density

def depth_from_pressure(
    pressure_in: np.ndarray,
    latitude: float,
    depth_units: Literal["m", "ft"] = "m",
    pressure_units: Literal["dbar", "psi"] = "dbar",
):
    """Derive depth from pressure and latitude.

    :param pressure: Numpy array of floats representing pressure, in
        dbar or psi
    :param latitude: Latitude (-90.0 to 90.0)
    :param depth_units: 'm' for meters, 'ft' for feet. Defaults to 'm'.
    :param pressure_units: 'dbar' for decibars, 'psi' for PSI. Defaults
        to 'dbar'.

    :return: A numpy array representing depth in meters or feet
    """
    pressure = pressure_in.copy()
    if pressure_units == "psi":
        pressure /= DBAR_TO_PSI

    depth = -gsw.z_from_p(pressure, latitude)

    if depth_units == "ft":
        depth *= 3.28084

    return depth


def process_hex_strings_array(hex_strings: List[str]):
    """Verarbeitet alle Hexstrings auf einmal mit Array-Operationen"""
    
    # Extrahiere alle benötigten Bytes
    all_bytes = extract_bytes_from_hexstrings(hex_strings, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
    
    # Konvertiere Bytes zu numerischen Werten
    byte_arrays = {}
    for i in range(all_bytes.shape[1]):
        byte_arrays[f'byte_{i}'] = np.array([int(byte, 16) for byte in all_bytes[:, i]])
    
    # Temperatur 1 (Bytes 0, 1, 2)
    temp1_counts = 1000 / (byte_arrays['byte_0'] * 256 + byte_arrays['byte_1'] + byte_arrays['byte_2'] / 256)
    temperatures1 = convert_temperature_array(temp1_counts, temp1_coefs)
    
    # Temperatur 2 (Bytes 9, 10, 11)
    temp2_counts = 1000 / (byte_arrays['byte_9'] * 256 + byte_arrays['byte_10'] + byte_arrays['byte_11'] / 256)
    temperatures2 = convert_temperature_array(temp2_counts, temp2_coefs)
    
    '''    # Druck (Bytes 6, 7, 8)
    pressure_counts = byte_arrays['byte_6'] * 256 + byte_arrays['byte_7'] + byte_arrays['byte_8'] / 256
    compensation_voltages = np.ones_like(pressure_counts)  # Annahme: immer 1
    pressures = convert_pressure_array(pressure_counts, compensation_voltages, pres_coefs) '''
    
    pressures = np.array([1.149,1.197,1.197,1.149,1.253])

    # Leitfähigkeit 1 (Bytes 3, 4, 5)
    cond1_counts = byte_arrays['byte_3'] * 256 + byte_arrays['byte_4'] + byte_arrays['byte_5'] / 256
    conductivities1 = convert_conductivity_array(cond1_counts, temperatures1, pressures, cond1_coefs)
    
    # Leitfähigkeit 2 (Bytes 12, 13, 14)
    cond2_counts = byte_arrays['byte_12'] * 256 + byte_arrays['byte_13'] + byte_arrays['byte_14'] / 256
    conductivities2 = convert_conductivity_array(cond2_counts, temperatures2, pressures, cond2_coefs)
    
    # Sauerstoff-Spannungen (Bytes 30-33 und 33-36)
    oxy_voltages1 = hex_to_voltage([hex_str[30:33] for hex_str in hex_strings])
    oxy_voltages2 = hex_to_voltage([hex_str[33:36] for hex_str in hex_strings])
    
    # Salinität
    salinities1 = gsw.SP_from_C(conductivities1, temperatures1, pressures)
    salinities2 = gsw.SP_from_C(conductivities2, temperatures2, pressures)
    
    # Sauerstoff (vereinfacht - für echte Array-Implementierung müssten die Sauerstofffunktionen auch array-fähig gemacht werden)
    oxygen1_ml = np.array([convert_sbe43_oxygen(np.array([v]), t, p, s, oxy1_coefs)[0] 
                          for v, t, p, s in zip(oxy_voltages1, temperatures1, pressures, salinities1)])
    oxygen2_ml = np.array([convert_sbe43_oxygen(np.array([v]), t, p, s, oxy2_coefs)[0] 
                          for v, t, p, s in zip(oxy_voltages2, temperatures2, pressures, salinities2)])
    
    # Tiefe
    deg = 54
    minutes = 9.31
    latitude = deg + minutes / 60
    depths = depth_from_pressure(pressures, latitude)
    
    results = {
        'temperature1': temperatures1,
        'temperature2': temperatures2,
        'pressure': pressures,
        'conductivity1': conductivities1,
        'conductivity2': conductivities2,
        'oxygen1_ml': oxygen1_ml,
        'oxygen2_ml': oxygen2_ml,
        'salinity1': salinities1,
        'salinity2': salinities2,
        'depth': depths
    }
    
    return results

# Hauptprogramm
if __name__ == "__main__":
    hex_strings = [
        "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465",
        "0CD0811076768045960CCA6F1107B27D48A3986C40F0DEF2FFFFBD0000BB29512A089DB7004D638F5C91C465",
        "0CD08710767B8045960CCA6D1107B17D48A3986C40F10EF0FFFFBE0000BB29512A089DB7004D63905D91C465",
        "0CD08710767A80458F0CCA6C1107AE7D58A3986C3FF16EE7FFFFBE0000BB29512A089DB7004D63915D91C465",
        "0CD08A10767680459E0CCA6F1107B27D58A3987C40F18EE3FFFFBE0000BB29512A089DB7004D63925D91C465"
    ]
    
    results = process_hex_strings_array(hex_strings)
    
    print("Ergebnisse als Arrays:")
    for key, value in results.items():
        print(f"{key}: {value}")