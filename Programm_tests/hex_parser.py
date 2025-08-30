#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Sea-Bird HEX File Parser and Data Converter"""

import re
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from math import e
import warnings

# Kalibrationskoeffizienten Klassen (kopiert aus cal_coefficients.py)
@dataclass
class TemperatureCoefficients:
    a0: float
    a1: float
    a2: float
    a3: float

@dataclass
class PressureCoefficients:
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

@dataclass
class ConductivityCoefficients:
    g: float
    h: float
    i: float
    j: float
    cpcor: float
    ctcor: float
    wbotc: float

# Konvertierungsfunktionen (kopiert aus conversion.py)
DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
KELVIN_OFFSET_0C = 273.15

def convert_temperature(
    temperature_counts_in: np.ndarray,
    coefs: TemperatureCoefficients,
    standard: str = "ITS90",
    units: str = "C",
    use_mv_r: bool = False,
):
    """Konvertiert Temperatur-Rohdaten in Grad Celsius"""
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
        temperature *= 1.00024  # ITS90_TO_IPTS68
    if units == "F":
        temperature = temperature * 9 / 5 + 32

    return temperature

def convert_pressure(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: PressureCoefficients,
    units: str = "psia",
):
    """Konvertiert Druck-Rohdaten"""
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
    """Konvertiert Leitfähigkeits-Rohdaten"""
    f = conductivity_count * np.sqrt(1 + coefs.wbotc * temperature) / 1000
    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperature + coefs.cpcor * pressure
    return numerator / denominator

def depth_from_pressure(
    pressure_in: np.ndarray,
    latitude: float,
    depth_units: str = "m",
    pressure_units: str = "dbar",
):
    """Berechnet Tiefe aus Druck und Breitengrad"""
    pressure = pressure_in.copy()
    if pressure_units == "psi":
        pressure /= DBAR_TO_PSI

    # Vereinfachte Tiefenberechnung (gsw.z_from_p ersetzt)
    depth = pressure * 1.0  # Vereinfachung: 1 dbar ≈ 1 Meter
    
    if depth_units == "ft":
        depth *= 3.28084

    return depth

def density_from_t_c_p(
    temperature: np.ndarray,
    conductivity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
):
    """Berechnet Dichte aus Temperatur, Leitfähigkeit und Druck"""
    # Vereinfachte Dichteberechnung
    salinity = conductivity * 0.5  # Vereinfachte Umrechnung
    density = 1000 + 0.8 * salinity - 0.2 * temperature + 0.1 * pressure
    return density

@dataclass
class HexHeader:
    """Speichert Header-Informationen aus der HEX-Datei"""
    filename: str
    software_version: str
    temperature_sn: str
    conductivity_sn: str
    bytes_per_scan: int
    voltage_words: int
    scans_averaged: int
    system_time: str
    latitude: str
    longitude: str
    utc_time: str
    cruise: str
    station: str
    platform: str
    cast: str
    operator: str
    gps_time: str
    gps_lat: str
    gps_lon: str
    echo_depth: str
    air_pressure: str

class HexDataParser:
    def __init__(self):
        self.header = None
        self.raw_data = []
        self.converted_data = []
        
        # Beispiel-Kalibrationskoeffizienten (müssen an reale Werte angepasst werden)
        self.temp_coefs = TemperatureCoefficients(
            a0=1.024290e-03,
            a1=2.606266e-04,
            a2=1.382136e-06,
            a3=2.601036e-07
        )
        
        self.pressure_coefs = PressureCoefficients(
            pa0=-4.337472e-01,
            pa1=1.263875e-02,
            pa2=-1.466348e-11,
            ptca0=2.664815e+01,
            ptca1=-6.432251e-04,
            ptca2=1.099536e-05,
            ptcb0=2.617072e+01,
            ptcb1=5.487685e-03,
            ptcb2=-5.728847e-05,
            ptempa0=5.756999e+01,
            ptempa1=-3.394574e+00,
            ptempa2=6.387549e-02
        )
        
        self.cond_coefs = ConductivityCoefficients(
            g=-1.000000e+00,
            h=1.400000e+00,
            i=0.000000e+00,
            j=0.000000e+00,
            cpcor=-9.570000e-08,
            ctcor=3.250000e-06,
            wbotc=2.550000e-04
        )

    def parse_header(self, lines: List[str]) -> HexHeader:
        """Parst den Header der HEX-Datei"""
        header_info = {}
        
        for line in lines:
            if line.startswith('*END*'):
                break
                
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip('* ').lower()
                value = value.strip()
                header_info[key] = value
            elif line.startswith('**'):
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip('* ').lower()
                    value = parts[1].strip()
                    header_info[key] = value
        
        return HexHeader(
            filename=header_info.get('filename', ''),
            software_version=header_info.get('software version', ''),
            temperature_sn=header_info.get('temperature sn', ''),
            conductivity_sn=header_info.get('conductivity sn', ''),
            bytes_per_scan=int(header_info.get('number of bytes per scan', 0)),
            voltage_words=int(header_info.get('number of voltage words', 0)),
            scans_averaged=int(header_info.get('number of scans averaged by the deck unit', 1)),
            system_time=header_info.get('system upload time', ''),
            latitude=header_info.get('nmea latitude', ''),
            longitude=header_info.get('nmea longitude', ''),
            utc_time=header_info.get('nmea utc (time)', ''),
            cruise=header_info.get('cruise', ''),
            station=header_info.get('station', ''),
            platform=header_info.get('platform', ''),
            cast=header_info.get('cast', ''),
            operator=header_info.get('operator', ''),
            gps_time=header_info.get('gps_time', ''),
            gps_lat=header_info.get('gps_lat', ''),
            gps_lon=header_info.get('gps_lon', ''),
            echo_depth=header_info.get('echo_depth', ''),
            air_pressure=header_info.get('air_pressure', '')
        )

    def parse_hex_line(self, hex_line: str) -> Dict[str, int]:
        """Parst eine einzelne HEX-Zeile und extrahiert die Rohwerte"""
        hex_line = hex_line.strip()
        
        if len(hex_line) < 88:  # 44 Bytes = 88 Hex-Zeichen
            return {}
        
        try:
            # Extrahiere Daten basierend auf typischer SBE9-Struktur
            data = {
                'temperature_count': int(hex_line[0:8], 16),
                'conductivity_count': int(hex_line[8:16], 16),
                'pressure_count': int(hex_line[16:24], 16),
                'pressure_compensation': int(hex_line[24:32], 16),
                'voltage0': int(hex_line[32:36], 16),
                'voltage1': int(hex_line[36:40], 16),
                'voltage2': int(hex_line[40:44], 16),
                'voltage3': int(hex_line[44:48], 16),
                'voltage4': int(hex_line[48:52], 16),
                'scan_counter': int(hex_line[52:60], 16),
                'time_milliseconds': int(hex_line[60:68], 16),
                'frame_flag': int(hex_line[68:70], 16),
                'checksum': int(hex_line[70:72], 16)
            }
            return data
        except ValueError:
            return {}

    def convert_raw_data(self, raw_data: List[Dict]) -> List[Dict]:
        """Konvertiert Rohdaten in physikalische Werte"""
        converted = []
        
        # Extrahiere Arrays für die Konvertierung
        temp_counts = np.array([d['temperature_count'] for d in raw_data])
        cond_counts = np.array([d['conductivity_count'] for d in raw_data])
        pressure_counts = np.array([d['pressure_count'] for d in raw_data])
        pressure_comp = np.array([d['pressure_compensation'] for d in raw_data])
        
        # Führe die Konvertierungen durch
        temperatures = convert_temperature(temp_counts, self.temp_coefs)
        pressures_psi = convert_pressure(pressure_counts, pressure_comp, self.pressure_coefs, "psia")
        pressures_dbar = pressures_psi * PSI_TO_DBAR
        conductivities = convert_conductivity(cond_counts, temperatures, pressures_dbar, self.cond_coefs)
        
        # Berechne Dichte und Tiefe
        densities = density_from_t_c_p(temperatures, conductivities, pressures_dbar)
        
        # Extrahiere Breitengrad für Tiefenberechnung
        lat_match = re.search(r'(\d+)\s+(\d+\.\d+)\s+([NS])', self.header.gps_lat)
        if lat_match:
            degrees = float(lat_match.group(1))
            minutes = float(lat_match.group(2))
            hemisphere = lat_match.group(3)
            latitude = degrees + minutes / 60.0
            if hemisphere == 'S':
                latitude = -latitude
        else:
            latitude = 54.0  # Fallback-Wert
            
        depths = depth_from_pressure(pressures_dbar, latitude)
        
        # Erstelle konvertierte Datensätze
        for i in range(len(raw_data)):
            converted.append({
                'scan_counter': raw_data[i]['scan_counter'],
                'time_milliseconds': raw_data[i]['time_milliseconds'],
                'temperature_c': temperatures[i],
                'pressure_dbar': pressures_dbar[i],
                'conductivity_s_m': conductivities[i],
                'density_kg_m3': densities[i],
                'depth_m': depths[i],
                'voltage0': raw_data[i]['voltage0'],
                'voltage1': raw_data[i]['voltage1'],
                'voltage2': raw_data[i]['voltage2'],
                'voltage3': raw_data[i]['voltage3'],
                'voltage4': raw_data[i]['voltage4']
            })
        
        return converted

    def read_hex_file(self, file_path: str):
        """Liest eine komplette HEX-Datei und verarbeitet sie"""
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        # Trenne Header und Daten
        header_lines = []
        data_lines = []
        in_header = True
        
        for line in lines:
            if line.startswith('*END*'):
                in_header = False
                continue
            if in_header:
                header_lines.append(line)
            else:
                if line.strip() and not line.startswith('*'):
                    data_lines.append(line.strip())
        
        # Parse Header
        self.header = self.parse_header(header_lines)
        
        # Parse Datenzeilen
        self.raw_data = []
        for line in data_lines:
            parsed = self.parse_hex_line(line)
            if parsed:
                self.raw_data.append(parsed)
        
        # Konvertiere Daten
        if self.raw_data:
            self.converted_data = self.convert_raw_data(self.raw_data)
        
        return self.header, self.converted_data

    def export_to_csv(self, output_path: str):
        """Exportiert die konvertierten Daten als CSV"""
        if not self.converted_data:
            print("Keine Daten zum Exportieren vorhanden")
            return
        
        with open(output_path, 'w', encoding='utf-8') as file:
            # Schreibe Header-Informationen
            file.write(f"# Sea-Bird Data Export\n")
            file.write(f"# Original File: {self.header.filename}\n")
            file.write(f"# Cruise: {self.header.cruise}\n")
            file.write(f"# Station: {self.header.station}\n")
            file.write(f"# Cast: {self.header.cast}\n")
            file.write(f"# GPS Position: {self.header.gps_lat}, {self.header.gps_lon}\n")
            file.write(f"# Date: {self.header.gps_time}\n\n")
            
            # Schreibe Spaltenüberschriften
            file.write("ScanCounter,TimeMs,Temperature_C,Pressure_dbar,Conductivity_S_m,Density_kg_m3,Depth_m,Voltage0,Voltage1,Voltage2,Voltage3,Voltage4\n")
            
            # Schreibe Daten
            for data in self.converted_data:
                file.write(f"{data['scan_counter']},{data['time_milliseconds']},"
                          f"{data['temperature_c']:.6f},{data['pressure_dbar']:.6f},"
                          f"{data['conductivity_s_m']:.6f},{data['density_kg_m3']:.6f},"
                          f"{data['depth_m']:.6f},{data['voltage0']},{data['voltage1']},"
                          f"{data['voltage2']},{data['voltage3']},{data['voltage4']}\n")

    def print_summary(self):
        """Gibt eine Zusammenfassung der Daten aus"""
        if not self.converted_data:
            print("Keine Daten verfügbar")
            return
        
        print("=== DATA SUMMARY ===")
        print(f"Number of scans: {len(self.converted_data)}")
        print(f"Temperature range: {min(d['temperature_c'] for d in self.converted_data):.2f} - {max(d['temperature_c'] for d in self.converted_data):.2f} °C")
        print(f"Pressure range: {min(d['pressure_dbar'] for d in self.converted_data):.2f} - {max(d['pressure_dbar'] for d in self.converted_data):.2f} dbar")
        print(f"Depth range: {min(d['depth_m'] for d in self.converted_data):.2f} - {max(d['depth_m'] for d in self.converted_data):.2f} m")

def main():
    """Hauptfunktion zum Verarbeiten von HEX-Dateien"""
    import sys
    import os
    
    if len(sys.argv) != 2:
        print("Verwendung: python hex_parser.py <hex_datei.hex>")
        print("Beispiel: python hex_parser.py CTD_Data.hex")
        sys.exit(1)
    
    hex_file = sys.argv[1]
    
    if not os.path.exists(hex_file):
        print(f"Datei {hex_file} nicht gefunden!")
        sys.exit(1)
    
    parser = HexDataParser()
    
    try:
        print(f"Verarbeite Datei: {hex_file}")
        header, data = parser.read_hex_file(hex_file)
        
        print("=== HEADER INFORMATION ===")
        print(f"Filename: {header.filename}")
        print(f"Cruise: {header.cruise}")
        print(f"Station: {header.station}")
        print(f"Cast: {header.cast}")
        print(f"GPS Position: {header.gps_lat}, {header.gps_lon}")
        print(f"Temperature Sensor: {header.temperature_sn}")
        print(f"Conductivity Sensor: {header.conductivity_sn}")
        print()
        
        if data:
            # Zeige einige konvertierte Daten
            print("=== CONVERTED DATA (first 5 scans) ===")
            print("Scan | Temp (°C) | Pressure (dbar) | Conductivity (S/m) | Depth (m)")
            print("-" * 80)
            
            for i, scan in enumerate(data[:5]):
                print(f"{scan['scan_counter']:4} | {scan['temperature_c']:8.3f} | "
                      f"{scan['pressure_dbar']:15.3f} | {scan['conductivity_s_m']:18.5f} | "
                      f"{scan['depth_m']:8.2f}")
            
            # Exportiere als CSV
            base_name = os.path.splitext(hex_file)[0]
            csv_file = f"{base_name}_converted.csv"
            parser.export_to_csv(csv_file)
            print(f"\nData exported to {csv_file}")
            
            # Zeige Zusammenfassung
            parser.print_summary()
        else:
            print("Keine Daten konnten konvertiert werden")
            
    except Exception as e:
        print(f"Fehler beim Verarbeiten der Datei: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()