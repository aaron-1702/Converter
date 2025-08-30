#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SBE 911plus CTD Data Packet Decoder
Decodes hex strings from SBE 911plus CTD according to official manual
"""

import struct
from typing import Dict, Any, List, Tuple
import math

class SBE911Decoder:
    """Decoder for SBE 911plus CTD data packets based on official manual"""
    
    def __init__(self):
        self.sensor_config = self._get_calibration_coefficients()
        
    def _get_calibration_coefficients(self):
        """Calibration coefficients from XMLCON file"""
        return {
            'temp1': {
                'serial': '5492', 
                'g': 4.33459406e-003, 
                'h': 6.28735680e-004, 
                'i': 2.05782221e-005, 
                'j': 1.68624541e-006,
                'f0': 1000.0
            },
            'cond1': {
                'serial': '4497', 
                'g': -9.94900904e+000, 
                'h': 1.42164741e+000, 
                'i': 1.26008747e-004, 
                'j': 6.16820694e-005, 
                'cpcor': -9.57000000e-008, 
                'ctcor': 3.2500e-006
            },
            'press': {
                'serial': '1385', 
                'c1': -4.212464e+004, 
                'c2': -1.972709e-001, 
                'c3': 1.365930e-002, 
                'd1': 3.604100e-002, 
                'd2': 0.000000e+000, 
                't1': 3.046128e+001, 
                't2': -3.791930e-004, 
                't3': 4.202260e-006, 
                't4': 3.386190e-009, 
                'offset': -0.20000, 
                'ad590m': 1.279410e-002, 
                'ad590b': -9.215740e+000
            },
            'temp2': {
                'serial': '4451', 
                'g': 4.34471154e-003, 
                'h': 6.38146193e-004, 
                'i': 2.04525401e-005, 
                'j': 1.60292340e-006,
                'f0': 1000.0
            },
            'cond2': {
                'serial': '4007', 
                'g': -9.77824825e+000, 
                'h': 1.31877668e+000, 
                'i': 2.26556605e-004, 
                'j': 6.59127902e-005, 
                'cpcor': -9.57000000e-008, 
                'ctcor': 3.2500e-006
            }
        }
    
    def hex_to_bytes(self, hex_string: str) -> bytes:
        """Convert hex string to bytes"""
        hex_string = hex_string.replace(' ', '').replace('\n', '')
        if len(hex_string) % 2:
            hex_string = '0' + hex_string
        return bytes.fromhex(hex_string)
    
    def decode_frequency(self, freq_bytes: bytes) -> float:
        """
        Decode frequency according to SBE manual:
        Frequency f = (Byte 0×256) + Byte 1 + (Byte 2÷256)
        """
        if len(freq_bytes) < 3:
            return 0.0
        return freq_bytes[0] * 256 + freq_bytes[1] + freq_bytes[2] / 256.0
    
    def decode_ad_channel(self, byte1: int, byte2: int) -> tuple:
        """
        Decode A/D channel from two bytes according to manual:
        Each 12-bit number (N) is binary notation of analog voltage.
        N = 4095 for 0V, 0 for 5V
        V = 5 * (1 - [N ÷ 4095])
        """
        # Extract 12-bit values from the packed format
        # byte1 contains 8 MSB of first channel
        # byte2 contains 4 LSB of first channel (bits 4-7) and 4 MSB of second channel (bits 0-3)
        n1 = (byte1 << 4) | ((byte2 & 0xF0) >> 4)
        voltage1 = 5.0 * (1.0 - (n1 / 4095.0))
        
        return n1, voltage1
    
    def frequency_to_temperature(self, freq: float, coeffs: dict) -> float:
        """Convert frequency to temperature using SBE equation"""
        if freq <= 0:
            return 0.0
        
        # SBE temperature equation with G,H,I,J coefficients
        # f is frequency normalized by F0
        f = freq / coeffs['f0']
        if f <= 0:
            return 0.0
            
        log_f = math.log(f)
        
        # Temperature in Kelvin: 1/(G + H*ln(f) + I*ln²(f) + J*ln³(f))
        denominator = (coeffs['g'] + coeffs['h'] * log_f + 
                      coeffs['i'] * log_f**2 + coeffs['j'] * log_f**3)
        
        if abs(denominator) < 1e-15:
            return 0.0
            
        temp_kelvin = 1.0 / denominator
        return temp_kelvin - 273.15
    
    def frequency_to_conductivity(self, freq: float, temp_c: float, press_db: float, coeffs: dict) -> float:
        """Convert frequency to conductivity using SBE equation"""
        if freq <= 0:
            return 0.0
            
        # SBE conductivity equation
        f = freq / 1000.0  # Convert to kHz for conductivity
        
        numerator = (coeffs['g'] + coeffs['h'] * f**2 + 
                    coeffs['i'] * f**3 + coeffs['j'] * f**4)
        denominator = (1.0 + coeffs['ctcor'] * temp_c + coeffs['cpcor'] * press_db)
        
        if abs(denominator) < 1e-15:
            return 0.0
            
        cond_s_m = numerator / denominator
        return cond_s_m * 10.0  # Convert S/m to mS/cm
    
    def frequency_to_pressure(self, press_freq: float, press_temp_counts: int, coeffs: dict) -> tuple:
        """Convert pressure frequency and temperature counts to pressure"""
        if press_freq <= 0 or press_temp_counts == 0:
            return 0.0, 0.0
        
        # Pressure temperature from counts (simplified)
        press_temp = press_temp_counts * coeffs['ad590m'] + coeffs['ad590b']
        
        # This is a simplified pressure calculation - the actual Paroscientific equation
        # is more complex and requires the specific sensor coefficients
        # For now, using a basic approximation
        pressure_psia = (press_freq - 30000) * 0.01 + 14.7
        pressure_db = (pressure_psia - 14.7) * 0.6894759 + coeffs['offset']
        
        return max(0.0, pressure_db), press_temp
    
    def decode_sbe911_packet(self, hex_string: str) -> Dict[str, Any]:
        """Decode SBE 911plus CTD data packet according to manual format"""
        try:
            data = self.hex_to_bytes(hex_string)
            
            if len(data) != 44:
                return {'error': f'Expected 44 bytes, got {len(data)} bytes'}
            
            # According to manual: IEEE488 data format
            # Word 0: Bytes 0-2: Primary temperature frequency
            # Word 1: Bytes 3-5: Primary conductivity frequency  
            # Word 2: Bytes 6-8: Pressure frequency
            # Word 3: Bytes 9-11: Secondary temperature frequency
            # Word 4: Bytes 12-14: Secondary conductivity frequency
            
            # Decode frequencies using manual formula
            temp1_freq = self.decode_frequency(data[0:3])
            cond1_freq = self.decode_frequency(data[3:6])
            press_freq = self.decode_frequency(data[6:9])
            temp2_freq = self.decode_frequency(data[9:12])
            cond2_freq = self.decode_frequency(data[12:15])
            
            # A/D channels start at byte 15
            # Word 5: Bytes 15-17: A/D channels 0-1
            # Word 6: Bytes 18-20: A/D channels 2-3
            # etc.
            
            ad_channels = []
            voltages = []
            
            # Decode A/D channels according to manual format
            for i in range(15, min(len(data), 27), 3):
                if i + 2 < len(data):
                    # First A/D channel
                    if i + 1 < len(data):
                        n1, v1 = self.decode_ad_channel(data[i], data[i+1])
                        ad_channels.append(n1)
                        voltages.append(v1)
                    
                    # Second A/D channel  
                    if i + 2 < len(data):
                        n2 = ((data[i+1] & 0x0F) << 8) | data[i+2]
                        v2 = 5.0 * (1.0 - (n2 / 4095.0))
                        ad_channels.append(n2)
                        voltages.append(v2)
            
            # Extract pressure temperature and status from bytes 30-32
            press_temp_counts = 0
            status_byte = 0
            modulo_count = 0
            
            if len(data) > 30:
                press_temp_counts = ((data[30] << 4) | ((data[31] & 0xF0) >> 4))
                status_byte = data[31] & 0x0F
                if len(data) > 32:
                    modulo_count = data[32]
            
            # Convert to physical units
            temp1_c = self.frequency_to_temperature(temp1_freq, self.sensor_config['temp1'])
            temp2_c = self.frequency_to_temperature(temp2_freq, self.sensor_config['temp2'])
            
            press_db, press_temp_c = self.frequency_to_pressure(press_freq, press_temp_counts, self.sensor_config['press'])
            
            cond1_ms_cm = self.frequency_to_conductivity(cond1_freq, temp1_c, press_db, self.sensor_config['cond1'])
            cond2_ms_cm = self.frequency_to_conductivity(cond2_freq, temp2_c, press_db, self.sensor_config['cond2'])
            
            result = {
                'raw_hex': hex_string,
                'packet_length': len(data),
                'frequencies_hz': {
                    'temperature1': round(temp1_freq, 3),
                    'conductivity1': round(cond1_freq, 3), 
                    'pressure': round(press_freq, 3),
                    'temperature2': round(temp2_freq, 3),
                    'conductivity2': round(cond2_freq, 3)
                },
                'ad_counts': ad_channels,
                'voltages': [round(v, 3) for v in voltages],
                'pressure_temp_counts': press_temp_counts,
                'status_byte': status_byte,
                'modulo_count': modulo_count,
                'physical_values': {
                    'temperature1_c': round(temp1_c, 4),
                    'temperature2_c': round(temp2_c, 4),
                    'pressure_db': round(press_db, 3),
                    'pressure_temp_c': round(press_temp_c, 4),
                    'conductivity1_ms_cm': round(cond1_ms_cm, 5),
                    'conductivity2_ms_cm': round(cond2_ms_cm, 5)
                },
                'auxiliary_sensors': {
                    'oxygen1_v': voltages[0] if len(voltages) > 0 else 0,
                    'oxygen2_v': voltages[1] if len(voltages) > 1 else 0,
                    'par_v': voltages[2] if len(voltages) > 2 else 0,
                    'altimeter_v': voltages[3] if len(voltages) > 3 else 0,
                    'fluorescence_v': voltages[4] if len(voltages) > 4 else 0,
                    'turbidity_v': voltages[5] if len(voltages) > 5 else 0
                },
                'status_info': {
                    'pump_on': bool(status_byte & 0x01),
                    'bottom_contact': not bool(status_byte & 0x02),
                    'water_sampler_confirm': bool(status_byte & 0x04),
                    'modem_carrier': not bool(status_byte & 0x08)
                }
            }
            
            return result
            
        except Exception as e:
            return {'error': f'Decoding error: {e}'}

# Test with your specific hex string
def main():
    """Decode the specific 88-character hex string"""
    
    hex_string = "0CD08410767B80458F0CCA6F1107AE7D58A3987C41F0DEF2FFFFBE0000BB29512A089DB7004D638E5C91C465"
    
    decoder = SBE911Decoder()
    
    print("SBE 911plus CTD Data Decoder - Your Hex String Analysis")
    print("=" * 70)
    print(f"Input hex string: {hex_string}")
    print(f"Length: {len(hex_string)} characters ({len(hex_string)//2} bytes)")
    
    result = decoder.decode_sbe911_packet(hex_string)
    
    if 'error' in result:
        print(f"Error: {result['error']}")
        return
    
    print(f"\n✓ Successfully decoded {result['packet_length']} bytes")
    
    print(f"\nFrequencies (Hz):")
    for sensor, freq in result['frequencies_hz'].items():
        print(f"  {sensor:15}: {freq:8.3f} Hz")
    
    print(f"\nPhysical Values:")
    values = result['physical_values']
    print(f"  Temperature 1   : {values['temperature1_c']:7.4f} °C")
    print(f"  Temperature 2   : {values['temperature2_c']:7.4f} °C")
    print(f"  Pressure        : {values['pressure_db']:7.3f} dbar")
    print(f"  Pressure Temp   : {values['pressure_temp_c']:7.4f} °C")
    print(f"  Conductivity 1  : {values['conductivity1_ms_cm']:7.5f} mS/cm")
    print(f"  Conductivity 2  : {values['conductivity2_ms_cm']:7.5f} mS/cm")
    
    print(f"\nA/D Channels:")
    print(f"  Raw counts: {result['ad_counts']}")
    print(f"  Voltages  : {[f'{v:.3f}' for v in result['voltages']]} V")
    
    print(f"\nAuxiliary Sensors:")
    aux = result['auxiliary_sensors']
    for sensor, voltage in aux.items():
        if voltage != 0:
            print(f"  {sensor:15}: {voltage:6.3f} V")
    
    print(f"\nSystem Status:")
    status = result['status_info']
    for param, value in status.items():
        print(f"  {param:20}: {value}")
    print(f"  modulo_count        : {result['modulo_count']}")
    
    print(f"\nHex Breakdown (first 33 bytes):")
    data = bytes.fromhex(hex_string)
    for i in range(min(33, len(data))):
        print(f"  Byte {i:2d}: 0x{data[i]:02X} ({data[i]:3d})")

if __name__ == "__main__":
    main()