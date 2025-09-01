import numpy as np
from typing import Literal

DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
# taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/ # pylint: disable=line-too-long
ITS90_TO_IPTS68 = 1.00024
# micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
# [J K^{-1} mol^{-1}] Gas constant, NIST Reference on Constants retrieved 10-05-2015
R = 8.3144621
# [Coulombs mol^{-1}] Faraday constant, NIST Reference on Constants retrieved 10-05-2015
F = 96485.365

hex_string = "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465"

byte_0 = hex_string[0:2]
byte_1 = hex_string[2:4]
byte_2 = hex_string[4:6]
byte_3 = hex_string[6:8]
byte_4 = hex_string[8:10]
byte_5 = hex_string[10:12]
byte_6 = hex_string[12:14]
byte_7 = hex_string[14:16]
byte_8 = hex_string[16:18]
byte_9 = hex_string[18:20]
byte_10 = hex_string[20:22]
byte_11 = hex_string[22:24]
byte_12 = hex_string[24:26]
byte_13 = hex_string[26:28]
byte_14 = hex_string[28:30]


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











temperature_counts_in1 = np.array(1000/(int(byte_0,16)*256 + int(byte_1,16) + int(byte_2,16)/256))
temperature1 = convert_temperature(temperature_counts_in1,coefs= temp1_coefs)
print("Temp_1:",temperature1)

temperature_counts_in2 = np.array(1000/(int(byte_9,16)*256 + int(byte_10,16) + int(byte_11,16)/256))
temperature2 = convert_temperature(temperature_counts_in2,coefs= temp2_coefs)
print("Temp_2:",temperature2)

pressure_counts = np.array(1000/(int(byte_6,16)*256 + int(byte_7,16) + int(byte_8,16)/256))
pressure = convert_pressure(pressure_counts,np.array(0),coefs= pres_coefs)
print("Pressure:", pressure)

pressure = 1.149
conductivity_count = np.array(1000/(int(byte_3,16)*256 + int(byte_4,16) + int(byte_5,16)/256))
conductifity1 = convert_conductivity(conductivity_count,temperature1,pressure,coefs= cond1_coefs)
print("Conductivity_1:", conductifity1)

conductivity_count = np.array(1000/(int(byte_12,16)*256 + int(byte_13,16) + int(byte_14,16)/256))
conductifity2 = convert_conductivity(conductivity_count,temperature2,pressure,coefs= cond2_coefs)
print("Conductivity_2:", conductifity2)