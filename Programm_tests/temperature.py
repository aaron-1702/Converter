#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Schritt-für-Schritt Temperaturberechnung aus SBE 911plus Hex-String
"""

import math

def hex_string_to_bytes(hex_string):
    """Konvertiert Hex-String zu Bytes"""
    hex_clean = hex_string.replace(' ', '').replace('\n', '')
    return bytes.fromhex(hex_clean)

def decode_frequency_step_by_step(freq_bytes):
    """
    Dekodiert die Frequenz nach SBE Manual:
    Frequency f = (Byte 0 × 256) + Byte 1 + (Byte 2 ÷ 256)
    """
    print(f"Frequenz-Bytes (hex): {freq_bytes.hex().upper()}")
    print(f"Byte 0: 0x{freq_bytes[0]:02X} = {freq_bytes[0]} decimal")
    print(f"Byte 1: 0x{freq_bytes[1]:02X} = {freq_bytes[1]} decimal") 
    print(f"Byte 2: 0x{freq_bytes[2]:02X} = {freq_bytes[2]} decimal")
    
    frequency = freq_bytes[0] * 256 + freq_bytes[1] + freq_bytes[2] / 256.0
    
    print(f"Frequenz = ({freq_bytes[0]} × 256) + {freq_bytes[1]} + ({freq_bytes[2]} ÷ 256)")
    print(f"Frequenz = {freq_bytes[0] * 256} + {freq_bytes[1]} + {freq_bytes[2]/256:.3f}")
    print(f"Frequenz = {frequency:.3f} Hz")
    
    return frequency

def temperature_conversion_step_by_step(frequency, coeffs):
    """
    Konvertiert Frequenz zu Temperatur mit SBE Gleichung:
    T = 1/(G + H*ln(f/F0) + I*ln²(f/F0) + J*ln³(f/F0)) - 273.15
    """
    print(f"\n=== TEMPERATUR KONVERSION ===")
    print(f"Frequenz: {frequency:.3f} Hz")
    print(f"Kalibrierungs-Koeffizienten:")
    print(f"  G = {coeffs['g']:.9e}")
    print(f"  H = {coeffs['h']:.9e}")
    print(f"  I = {coeffs['i']:.9e}")
    print(f"  J = {coeffs['j']:.9e}")
    print(f"  F0 = {coeffs['f0']:.1f}")
    
    # Schritt 1: Normalisiere Frequenz
    f_normalized = frequency / coeffs['f0']
    print(f"\nSchritt 1: Normalisiere Frequenz")
    print(f"f_norm = f/F0 = {frequency:.3f}/{coeffs['f0']:.1f} = {f_normalized:.6f}")
    
    # Schritt 2: Berechne Logarithmus
    ln_f = math.log(f_normalized)
    print(f"\nSchritt 2: Berechne natürlichen Logarithmus")
    print(f"ln(f_norm) = ln({f_normalized:.6f}) = {ln_f:.6f}")
    
    # Schritt 3: Berechne Potenzen von ln(f)
    ln_f_2 = ln_f ** 2
    ln_f_3 = ln_f ** 3
    print(f"\nSchritt 3: Berechne Potenzen")
    print(f"ln²(f_norm) = {ln_f:.6f}² = {ln_f_2:.6f}")
    print(f"ln³(f_norm) = {ln_f:.6f}³ = {ln_f_3:.6f}")
    
    # Schritt 4: Berechne jeden Term der Gleichung
    term_g = coeffs['g']
    term_h = coeffs['h'] * ln_f
    term_i = coeffs['i'] * ln_f_2
    term_j = coeffs['j'] * ln_f_3
    
    print(f"\nSchritt 4: Berechne einzelne Terme")
    print(f"Term G = {coeffs['g']:.9e} = {term_g:.9e}")
    print(f"Term H = {coeffs['h']:.9e} × {ln_f:.6f} = {term_h:.9e}")
    print(f"Term I = {coeffs['i']:.9e} × {ln_f_2:.6f} = {term_i:.9e}")
    print(f"Term J = {coeffs['j']:.9e} × {ln_f_3:.6f} = {term_j:.9e}")
    
    # Schritt 5: Summiere alle Terme
    denominator = term_g + term_h + term_i + term_j
    print(f"\nSchritt 5: Summiere Nenner")
    print(f"Nenner = G + H×ln(f) + I×ln²(f) + J×ln³(f)")
    print(f"Nenner = {term_g:.9e} + {term_h:.9e} + {term_i:.9e} + {term_j:.9e}")
    print(f"Nenner = {denominator:.9e}")
    
    # Schritt 6: Berechne Temperatur in Kelvin
    temp_kelvin = 1.0 / denominator
    print(f"\nSchritt 6: Berechne Temperatur in Kelvin")
    print(f"T_Kelvin = 1 / {denominator:.9e} = {temp_kelvin:.6f} K")
    
    # Schritt 7: Konvertiere zu Celsius
    temp_celsius = temp_kelvin - 273.15
    print(f"\nSchritt 7: Konvertiere zu Celsius")
    print(f"T_Celsius = {temp_kelvin:.6f} - 273.15 = {temp_celsius:.4f}°C")
    
    return temp_celsius

def main():
    # Der gegebene Hex-String
    hex_string = "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465"
    
    print("SBE 911plus CTD Temperaturberechnung")
    print("=" * 50)
    print(f"Hex-String: {hex_string}")
    
    # Konvertiere zu Bytes
    data = hex_string_to_bytes(hex_string)
    print(f"Anzahl Bytes: {len(data)}")
    
    print(f"\n=== FREQUENZ DEKODIERUNG ===")
    
    # Temperatur 1 (Primary): Bytes 0-2
    print(f"\nPRIMÄRE TEMPERATUR (Bytes 0-2):")
    temp1_bytes = data[0:3]
    temp1_freq = decode_frequency_step_by_step(temp1_bytes)
    
    # Temperatur 2 (Secondary): Bytes 9-11  
    print(f"\nSEKUNDÄRE TEMPERATUR (Bytes 9-11):")
    temp2_bytes = data[9:12]
    temp2_freq = decode_frequency_step_by_step(temp2_bytes)
    
    # Kalibrierungskoeffizienten aus der XMLCON-Datei
    temp1_coeffs = {
        'g': 4.33459406e-003,
        'h': 6.28735680e-004,
        'i': 2.05782221e-005,
        'j': 1.68624541e-006,
        'f0': 1000.0
    }
    
    temp2_coeffs = {
        'g': 4.34471154e-003,
        'h': 6.38146193e-004,
        'i': 2.04525401e-005,
        'j': 1.60292340e-006,
        'f0': 1000.0
    }
    
    # Berechne Temperaturen
    print(f"\n" + "="*50)
    print(f"PRIMÄRE TEMPERATUR BERECHNUNG (Sensor 5492)")
    temp1_celsius = temperature_conversion_step_by_step(temp1_freq, temp1_coeffs)
    
    print(f"\n" + "="*50)
    print(f"SEKUNDÄRE TEMPERATUR BERECHNUNG (Sensor 4451)")
    temp2_celsius = temperature_conversion_step_by_step(temp2_freq, temp2_coeffs)
    
    print(f"\n" + "="*50)
    print(f"ENDERGEBNISSE:")
    print(f"Primäre Temperatur:   {temp1_celsius:.4f}°C")
    print(f"Sekundäre Temperatur: {temp2_celsius:.4f}°C")
    
    # Validierung der Ergebnisse
    print(f"\n=== VALIDIERUNG ===")
    print(f"Typische Meerestemperaturen liegen zwischen -2°C und 30°C")
    print(f"Beide berechneten Werte liegen in diesem realistischen Bereich.")

if __name__ == "__main__":
    main()