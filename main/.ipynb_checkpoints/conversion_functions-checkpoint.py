import numpy as np
import math
import gsw
from datetime import datetime, timezone
from typing import Literal
from scipy import stats
from math import floor

from constants import *
from calibration_coefficients import *


# self-made 
def hex_to_voltage(hex_value: str):
    clean_hex = hex_value.strip().upper()
    
    # Hex to Decimal
    N = int(clean_hex, 16)

    # Convert to volt
    V = 5 * (1 - (N / 4095))

    return V

def Time_elapsed(hexstrings, freq=24):
    """
    Computes elapsed time for a list of hex strings.
    Assumption: Each string = 1 measurement, starting at 0 s.
    
    :param hexstrings: List of hex strings (each string is one scan).
    :param freq: Sampling frequency in Hz (default 24 Hz).
    :return: List of elapsed times (in seconds).
    """
    dt = 1 / freq  # time interval per measurement
    elapsed = [i * dt for i in range(len(hexstrings))]
    return elapsed

def convert_par(
    raw_par: np.ndarray,
    coefs: PARCoefficients,
):
    """Converts a raw voltage value for PAR to µmol photons/m2*s.

    All equation information comes from application note 96

    :param raw_par: raw output voltage from PAR sensor
    :param coefs: calibration coefficients for the PAR sensor

    :return: converted PAR in µmol photons/m2*s
    """
    
    par = coefs.multiplier * ((10**9 * 10**(raw_par/coefs.im))/coefs.a0) + coefs.a1

    return par

def convert_spar(
    raw_spar: np.ndarray,
    coefs: SPAR_Coefficients,
):
    """Converts a raw voltage value for SPAR to µmol photons/m2*s.

    All equation information comes from application note 96

    :param raw_par: raw output voltage from SPAR sensor
    :param coefs: calibration coefficients for the SPAR sensor

    :return: converted SPAR in µmol photons/m2*s
    """
    
    spar = raw_spar * coefs.conversionfactor * coefs.ratiomultiplier

    return spar

def convert_alimeter(
    voltage: np.array,
    coefs: Altimeter_Coefficients,
    units: Literal["m", "ft"] = "m"
):
    """Converts a raw voltage value for Alimeter to m.
    
    :param voltage: raw output voltage from Alimeter sensor
    :param coefs: calibration coefficients for the Alimeter sensor
    """
    alt_m = (voltage * 300)/ coefs.scalefactor + coefs.offset
    
    if units == "ft":
        alt_m *= 3.28084
    
    return alt_m

def pressure_from_frequency(
    pressure_counts: np.ndarray, 
    temp_counts: np.ndarray, 
    coefs: PressureCoefficients,
    units: Literal["dbar", "psia"] = "psia"
):
    """
    Berechnet den Druck [psia] aus Digiquartz-Frequenz und Temperaturcount (arrays erlaubt).

    :param freq_hz: gemessene Druckfrequenz [Hz]
    :param temp_freq_hz: gemessene Temperaturcount
    :param coefs: PressureCoefficients
    :return: Druck in psia
    """

    U = coefs.AD590M * temp_counts + coefs.AD590B
    C = coefs.C1 + coefs.C2 * U + coefs.C3 * U **2
    D = coefs.D1 + coefs.D2 * U
    T0 = (
        coefs.T1 
        + coefs.T2 * U 
        + coefs.T3 * U **2 
        + coefs.T4 * U **3 
        + coefs.T5 * U **4
    )
    # Pressure Periode in µs
    T = 10**6 / pressure_counts
    
    pressure = C * (1 - (T0**2/T**2)) * (1 - D * (1 -(T0**2/T**2)))
    pressure = pressure * coefs.Slope + (coefs.Offset * DBAR_TO_PSI) - sea_level_pressure
    
    if units == "dbar":
        pressure *= PSI_TO_DBAR
    
    return pressure

# Seabird-functions
def convert_temperature(
    temperature_counts_in: np.ndarray,
    coefs: TemperatureCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
    use_mv_r: bool = False,
):
    """Returns the value after converting it to degrees C, ITS-90.

    Data is expected to be raw data from an instrument in A/D counts

    :param temperature_counts_in: temperature value to convert in A/D
        counts
    :param coefs: calibration coefficients for the temperature sensor
    :param standard: whether to use ITS90 or to use IPTS-68 calibration
        standard
    :param units: whether to use celsius or to convert to fahrenheit
    :param use_mv_r: true to perform extra conversion steps required by
        some instruments (check the cal sheet to see if this is required)

    :return: temperature val converted to ITS-90 degrees C
    """

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
        temperature = temperature * 9 / 5 + 32  # Convert C to F

    return temperature



def convert_pressure(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: PressureCoefficients,
    units: Literal["dbar", "psia"] = "psia",
):
    """Converts pressure counts to PSIA (pounds per square inch, abolute).

    pressure_count and compensation_voltage are expected to be raw data
    from an instrument in A/D counts

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature compensation
        voltage, in counts or volts depending on the instrument
    :param coefs: calibration coefficients for the pressure sensor
    :param units: whether or not to use psia or dbar as the returned
        unit type

    :return: pressure val in PSIA
    """
    sea_level_pressure = 14.7

    t = (
        coefs.ptempa0
        + coefs.ptempa1 * compensation_voltage
        + coefs.ptempa2 * compensation_voltage**2
    )
    x = pressure_count - coefs.ptca0 - coefs.ptca1 * t - coefs.ptca2 * t**2
    n = x * coefs.ptcb0 / (coefs.ptcb0 + coefs.ptcb1 * t + coefs.ptcb2 * t**2)
    pressure = coefs.pa0 + coefs.pa1 * n + coefs.pa2 * n**2 - sea_level_pressure

    if units == "dbar":
        pressure *= PSI_TO_DBAR

    return pressure


def convert_conductivity(
    conductivity_count: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    coefs: ConductivityCoefficients,
):
    """Converts raw conductivity counts to S/m.

    Data is expected to be raw data from instrument in A/D counts

    :param conductivity_count: conductivity value to convert, in A/D
        counts
    :param temperature: reference temperature, in degrees C
    :param pressure: reference pressure, in dbar
    :param coefs: calibration coefficient for the conductivity sensor

    :return: conductivity val converted to S/m
    """
    f = conductivity_count * np.sqrt(1 + coefs.wbotc * temperature) / 1000
    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperature + coefs.cpcor * pressure
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

def convert_eco(
    raw: np.ndarray,
    coefs: ECOCoefficients,
):
    """Converts a raw value for any ECO measurand.

    :param raw: raw counts for digital, raw volts for analog
    :param coefs (ChlorophyllACoefficients): calibration coefficients

    :return: converted ECO measurement in calibration units
    """
    converted = coefs.slope * (raw - coefs.offset)

    return converted

