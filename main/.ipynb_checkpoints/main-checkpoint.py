# main.py
import numpy as np
from conversion_functions import *
from calibration_coefficients import *

# FEHLENDE VARIABLEN HINZUFÜGEN:
hexstrings = [
    "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465",
    "0CD0811076768045960CCA6F1107B27D48A3986C40F0DEF2FFFFBD0000BB29512A089DB7004D638F5C91C465",
    "0CD08710767B8045960CCA6D1107B17D48A3986C40F10EF0FFFFBE0000BB29512A089DB7004D63905D91C465",
    "0CD08710767A80458F0CCA6C1107AE7D58A3986C3FF16EE7FFFFBE0000BB29512A089DB7004D63915D91C465",
    "0CD08A10767680459E0CCA6F1107B27D58A3987C40F18EE3FFFFBE0000BB29512A089DB7004D63925D91C465",
    "0CD08410767A8045960CCA6C1107B17D58A3988C41F18EE3FFFFBE0000BB29512A089DB7004D63935D91C465",
    "0CD0811076798045960CCA6A1107B17D58A498AC41F18EE1FFFFBE0000BB29512A089DB7004D63945D91C465",
    "0CD08410767C80458F0CCA6A1107B07D48A398BC40F15ED4FFFFBE0000BB29512A089DB7004D63955D91C465",
    "0CD08710767A8045960CCA6F1107B17D48A498CC41F13ECEFFFFBE0000BB29512A089DB7004D63965D91C465",
    "0CD08110767980459D0CCA6F1107B07D48A398DC41F13ECEFFFFBD0000BB29512A089DB7004D63975D91C465"
]

bytes_list_array = [[hex_string[i:i+2] for i in range(0, len(hex_string), 2)] for hex_string in hexstrings]

pressure_test_array = np.array([1.149,1.197,1.197,1.149,1.253,1.197,1.197,1.149,1.197,1.246])


# name 1 = t090C: Temperature [ITS-90, deg C]
# name 2 = t190C: Temperature, 2 [ITS-90, deg C]

F0 = 1000
temperature_counts1 = []
temperature_counts2 = []

for byte_list in bytes_list_array:
    count1 = int(byte_list[0], 16) * 256 + int(byte_list[1], 16) + int(byte_list[2], 16) / 256
    count2 = int(byte_list[9], 16) * 256 + int(byte_list[10], 16) + int(byte_list[11], 16) / 256
    temperature_counts1.append(count1)
    temperature_counts2.append(count2)

temperature_counts1 = np.array(temperature_counts1)
temperature_counts2 = np.array(temperature_counts2)
print(temperature_counts1)

temperature_counts1 = F0 / temperature_counts1
temperature_counts2 = F0 / temperature_counts2

print(temperature_counts1)

temperatures1 = convert_temperature(temperature_counts1, coefs=temp1_coefs)
temperatures2 = convert_temperature(temperature_counts2, coefs=temp2_coefs)

print("Temperatures1:", temperatures1)
print("Temperatures2:", temperatures2)


# name 0 = prDM: Pressure, Digiquartz [db]   
pressure_counts = []

for byte_list in bytes_list_array:
    count = int(byte_list[6], 16) * 256 + int(byte_list[7], 16) + int(byte_list[8], 16) / 256
    pressure_counts.append(count)

pressure_counts = np.array(pressure_counts)

temperature_counts = []

for hex_string in hexstrings:
    temperature_value = int(hex_string[74:77],16)
    temperature_counts.append(temperature_value)

temperature_counts = np.array(temperature_counts)

pressure = pressure_from_frequency(pressure_counts, temperature_counts, coefs=pres_coefs)


print("Druck in psi", pressure, "gewünscht: 1.667")
pressure *= PSI_TO_DBAR # Druck von psia zu db laut seabird

print("Druck in dbar:", pressure, "gewünscht: 1.149")

#  CTD status
# 8-bit number from CTD
CTD_status_liste = []
bit_liste = []

for hex_string in hexstrings:
    CTD_status = format(int(hex_string[77:78], 16), '04b')
    bit = format(int(hex_string[78:80], 16), '08b')
    CTD_status_liste.append(CTD_status)
    bit_liste.append(bit)

CTD_status_liste = np.array(CTD_status_liste)
bit_liste = np.array(bit_liste)

print("CTD status:", CTD_status_liste)
print("8-Bit Zahl:", bit_liste)


# name 3 = c0mS/cm: Conductivity [mS/cm]
# name 4 = c1mS/cm: Conductivity, 2 [mS/cm]

conductivity_counts1 = []
conductivity_counts2 = []

for byte_list in bytes_list_array:
    count1 = int(byte_list[3], 16) * 256 + int(byte_list[4], 16) + int(byte_list[5], 16) / 256
    count2 = int(byte_list[12], 16) * 256 + int(byte_list[13], 16) + int(byte_list[14], 16) / 256
    conductivity_counts1.append(count1)
    conductivity_counts2.append(count2)

conductivity_counts1 = np.array(conductivity_counts1)
conductivity_counts2 = np.array(conductivity_counts2)

conductifitys1 = convert_conductivity(conductivity_counts1, temperatures1, pressure_test_array, coefs=cond1_coefs)
conductifitys2 = convert_conductivity(conductivity_counts2, temperatures2, pressure_test_array, coefs=cond2_coefs)

print("Conductifitys1:", conductifitys1)
print("Conductifitys2:", conductifitys2)

# name 15 = sal00: Salinity, Practical [PSU]
# name 16 = sal11: Salinity, Practical, 2 [PSU]

salinitys1 = gsw.SP_from_C(conductifitys1, temperatures1, pressure_test_array)
salinitys2 = gsw.SP_from_C(conductifitys2, temperatures2, pressure_test_array)

print("Salinity1:", salinitys1,"PSU")
print("Salinity2:", salinitys2,"PSU")

# name 5 = sbeox0ML/L: Oxygen, SBE 43 [ml/l]
# name 6 = sbeox1ML/L: Oxygen, SBE 43, 2 [ml/l]
# name 7 = sbox0Mm/Kg: Oxygen, SBE 43 [umol/kg]
# name 8 = sbox1Mm/Kg: Oxygen, SBE 43, 2 [umol/kg]

voltages0 = []
voltages1 = []

for hex_string in hexstrings:
    v0 = hex_to_voltage(hex_string[30:33])
    v1 = hex_to_voltage(hex_string[33:36])
    voltages0.append(v0)
    voltages1.append(v1)

voltages0 = np.array(voltages0)
voltages1 = np.array(voltages1)

oxygens1 = convert_sbe43_oxygen(voltages0,temperatures1,pressure_test_array,salinitys1,coefs= oxy1_coefs)
oxygens2 = convert_sbe43_oxygen(voltages1,temperatures2,pressure_test_array,salinitys2,coefs= oxy2_coefs)
print("Oxygens1:", oxygens1,"ml/l")
print("Oxygens2:", oxygens2,"ml/l")


potential_densitys1 = potential_density_from_t_s_p(temperatures1,salinitys1,pressure_test_array)
potential_densitys2 = potential_density_from_t_s_p(temperatures2,salinitys2,pressure_test_array)

oxygens3 = convert_oxygen_to_umol_per_kg(oxygens1, potential_densitys1) *1000
oxygens4 = convert_oxygen_to_umol_per_kg(oxygens2, potential_densitys2) *1000
print("Oxygens1:", oxygens3,"umol/kg")
print("Oxygens2:", oxygens4,"umol/kg")

# name 13 = timeS: Time, Elapsed [seconds]

time_elapsed = Time_elapsed(hexstrings)
print("Time_Elapsed:",time_elapsed,"s")

# name 11 = par: PAR/Irradiance, Biospherical/Licor
voltages2 = []

for hex_string in hexstrings:
    v2 = hex_to_voltage(hex_string[36:39])
    voltages2.append(v2)
    
voltages2 = np.array(voltages2)

par = convert_par(voltages2, coefs = par_coefs)
print("Par:", par, "µmol photons m-2s-1")

# Altimeter
voltages3 = []

for hex_string in hexstrings:
    v3 = hex_to_voltage(hex_string[39:42])
    voltages3.append(v3)
    
voltages3 = np.array(voltages3)

alt_m = 15.0 * voltages3 + 0.0
print("Voltages:", voltages3)
print("Alimeter:",alt_m)

# name 14 = dz/dtM: Descent Rate [m/s] 
pressure_test_array = np.array([1.149,1.197,1.197,1.149,1.253,1.197,1.197,1.149,1.197,1.246])
latitude = 54.50402
depths = depth_from_pressure(pressure, latitude)
print("Deaths: ",depths)



dz = np.diff(depths)
dtM = np.diff(time_elapsed)

descent_rate = dz/dtM


print("Sinkraten (m/s):", descent_rate)

navi_hex = []
latitudes = []
longitudes = []
position = []
for hex_string in hexstrings:
    navi = hex_string[60:74]
    navi_hex.append(navi)
    
navi_hex = np.array(navi_hex)

for hex_str in navi_hex:
        # Bytes extrahieren
        bytes_list = [int(hex_str[i:i+2], 16) for i in range(0, len(hex_str)-2, 2)]
        # Bytes aufsplitten
        b1, b2, b3, b4, b5, b6 = bytes_list[:6]
        last_byte_hex = hex_str[-2:]
        b7 = format(int(last_byte_hex, 16), '08b') 
        # Latitude & Longitude berechnen
        lat = ((b1 * 65536 + b2 * 256 + b3) / 50000) * (-1)**int(b7[6])
        lon = ((b4 * 65536 + b5 * 256 + b6) / 50000) * (-1)**int(b7[7])

        if b7[1] == 1:
            pos = "new"
        else: 
            pos = "old"
        position.append(pos)
        latitudes.append(lat)
        longitudes.append(lon)
print("Position:", position)  
print("Latitudes:", latitudes)
print("Longitudes:",longitudes)

# name 9 = flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
voltages4 = []

for hex_string in hexstrings:
    v4 = hex_to_voltage(hex_string[42:45])
    voltages4.append(v4)
    
voltages4 = np.array(voltages4)

flurometer = convert_eco(voltages4, coefs = fluroeco_coefs)
print("Fluorescence:" , flurometer, "mg/m^3")

# name 10 = turbWETntu0: Turbidity, WET Labs ECO [NTU]
voltages5 = []

for hex_string in hexstrings:
    v5 = hex_to_voltage(hex_string[45:48])
    voltages5.append(v5)
    
voltages5 = np.array(voltages5)

turbidity = convert_eco(voltages5, coefs = turbidity_coefs)
print("Turbidity:", turbidity, "NTU")

# name 12 = spar: SPAR, Biospherical/Licor

voltages9 = []

for hex_string in hexstrings:
    N = int(hex_string[57:60],16)
    v9 = N / 819
    voltages9.append(v9)
    
voltages9 = np.array(voltages9)

spar = convert_spar(voltages9, coefs = spar_coefs)
print("SPar:", spar, "µmol photons m-2s-1")

# ???
fluroeco_coefs2 = ECOCoefficients(slope=1.000, offset=0.000)


voltages7 = []

for hex_string in hexstrings:
    v7 = hex_to_voltage(hex_string[51:54])
    voltages7.append(v7)
    
voltages7 = np.array(voltages7)

flurometer = convert_eco(voltages7, coefs = fluroeco_coefs2)
print("Fluorescence:" , flurometer, "mg/m^3")
print("Why???????")

# Altimeter
voltages3 = []

for hex_string in hexstrings:
    v3 = hex_to_voltage(hex_string[39:42])
    voltages3.append(v3)
    
voltages3 = np.array(voltages3)

alt_m = 15.0 * voltages3 + 0.0

print(alt_m)

# name 17 = flag:  0.000e+00

flag = np.zeros_like(hexstrings, dtype=float)
print(flag)

# Date and Time
swapped = [
    "".join([s[-8:][i:i+2] for i in range(0, 8, 2)][::-1])
    for s in hexstrings
]
timestamps = [int(x, 16) for x in swapped]    
datetimes = [datetime.fromtimestamp(ts, tz=timezone.utc) for ts in timestamps]
datetimes_formatted = np.array([d.strftime("%Y-%m-%d %H:%M:%S") for d in datetimes])
print("Date and Times", datetimes_formatted)

print("Pressures:", pressure,"db\n")

print("Temperatures1:", temperatures1,"°C\n")
print("Temperatures2:", temperatures2,"°C\n")

print("Conductifitys1:", conductifitys1,"mS/cm\n")
print("Conductifitys2:", conductifitys2,"mS/cm\n")

print("Oxygens1:", oxygens1,"ml/l\n")
print("Oxygens2:", oxygens2,"ml/l\n")
print("Oxygens3:", oxygens3,"umol/kg\n")
print("Oxygens4:", oxygens4,"umol/kg\n")

print("Fluorescence:" , flurometer, "mg/m^3\n")

print("Turbidity:", turbidity, "NTU\n")

print("Par:", par, "µmol photons m-2s-1\n")
print("SPar:", spar, "µmol photons m-2s-1\n")

print("Time_Elapsed:", time_elapsed,"s\n")

print("Salinity1:", salinitys1,"PSU\n")
print("Salinity2:", salinitys2,"PSU\n")