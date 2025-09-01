import numpy as np
from typing import Literal

# ------------------ Konstanten ------------------
KELVIN_OFFSET_0C = 273.15
ITS90_TO_IPTS68 = 1.00024
PSI_TO_DBAR = 0.6894759

# ------------------ Hilfsfunktionen ------------------
def get_24bit_value(b0, b1, b2):
    return (int(b0, 16) << 16) + (int(b1, 16) << 8) + int(b2, 16)

def counts_to_frequency(raw_counts):
    return 2560000 / raw_counts  # Sea-Bird Formel

# ------------------ Temperatur ------------------
def convert_temperature(temperature_counts_in: np.ndarray, coefs, standard="ITS90", units="C"):
    log_t = np.log(temperature_counts_in)
    temperature = (1 / (coefs.a0 + coefs.a1 * log_t + coefs.a2 * log_t**2 + coefs.a3 * log_t**3)) - KELVIN_OFFSET_0C
    if standard == "IPTS68":
        temperature *= ITS90_TO_IPTS68
    if units == "F":
        temperature = temperature * 9 / 5 + 32
    return temperature

# ------------------ Druck ------------------
def convert_pressure(frequency_hz, U, coefs, units="dbar"):
    # Berechne Temperaturkompensation
    t = coefs["AD590M"] * U + coefs["AD590B"]

    # Druckberechnung nach Sea-Bird
    n = frequency_hz - coefs["D1"] - coefs["D2"] * t - coefs["T1"] * t**2
    pressure = coefs["C1"] + coefs["C2"] * n + coefs["C3"] * n**2

    if units == "dbar":
        pressure *= PSI_TO_DBAR

    return pressure

# ------------------ Leitfähigkeit ------------------
def convert_conductivity(frequency_khz, temperature, pressure, coefs):
    f = frequency_khz
    numerator = coefs["G"] + coefs["H"] * f**2 + coefs["I"] * f**3 + coefs["J"] * f**4
    denominator = 1 + coefs["CTcor"] * temperature + coefs["CPcor"] * pressure
    return numerator / denominator

# ------------------ XML-Werte (manuell aus deinem XML) ------------------
class TempCoefs:
    def __init__(self, a0, a1, a2, a3):
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3

temp1_coefs = TempCoefs(4.33459406e-003, 6.28735680e-004, 2.05782221e-005, 1.68624541e-006)
temp2_coefs = TempCoefs(4.34471154e-003, 6.38146193e-004, 2.04525401e-005, 1.60292340e-006)

pressure_coefs = {
    "C1": -4.212464e+004,
    "C2": -1.972709e-001,
    "C3": 1.365930e-002,
    "D1": 3.604100e-002,
    "D2": 0.0,
    "T1": 3.046128e+001,
    "T2": -3.791930e-004,
    "T3": 4.202260e-006,
    "T4": 3.386190e-009,
    "AD590M": 1.279410e-002,
    "AD590B": -9.215740e+000
}

cond1_coefs = {
    "G": -9.94900904e+000,
    "H": 1.42164741e+000,
    "I": 1.26008747e-004,
    "J": 6.16820694e-005,
    "CTcor": 3.25e-006,
    "CPcor": -9.57e-008
}

cond2_coefs = {
    "G": -9.77824825e+000,
    "H": 1.31877668e+000,
    "I": 2.26556605e-004,
    "J": 6.59127902e-005,
    "CTcor": 3.25e-006,
    "CPcor": -9.57e-008
}

# ------------------ Hex-Daten ------------------
hex_string = "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465"

# Bytes extrahieren
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

# ------------------ Berechnungen ------------------
# Temperatur 1
raw_temp1 = get_24bit_value(byte_0, byte_1, byte_2)
temperature1 = convert_temperature(np.array(raw_temp1), temp1_coefs)
print("Temp_1 [°C]:", temperature1)

# Temperatur 2
raw_temp2 = get_24bit_value(byte_9, byte_10, byte_11)
temperature2 = convert_temperature(np.array(raw_temp2), temp2_coefs)
print("Temp_2 [°C]:", temperature2)

# Druck
raw_press = get_24bit_value(byte_6, byte_7, byte_8)
f_press = counts_to_frequency(raw_press)
U = 0.0  # Keine separate Spannung aus Bytes verfügbar → ggf. schätzen
pressure = convert_pressure(f_press, U, pressure_coefs)
print("Pressure [dbar]:", pressure)

# Leitfähigkeit 1
raw_cond1 = get_24bit_value(byte_3, byte_4, byte_5)
f_cond1 = counts_to_frequency(raw_cond1) / 1000.0  # in kHz
conductivity1 = convert_conductivity(f_cond1, temperature1, pressure, cond1_coefs)
print("Conductivity_1 [S/m]:", conductivity1)

# Leitfähigkeit 2
raw_cond2 = get_24bit_value(byte_12, byte_13, byte_14)
f_cond2 = counts_to_frequency(raw_cond2) / 1000.0
conductivity2 = convert_conductivity(f_cond2, temperature2, pressure, cond2_coefs)
print("Conductivity_2 [S/m]:", conductivity2)
