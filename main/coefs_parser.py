import xml.etree.ElementTree as ET
from cal_coefficients import *

def parse_xmlcon_to_coefficients(file_path):
    """
    Extrahiert Kalibrierungskoeffizienten aus XMLCON und gibt sie in der Seabird-kompatiblen Struktur zurück.
    """

    with open(file_path, 'r') as f:
        xml_content = f.read()

    
    # XML parsen
    if isinstance(xml_content, str) and xml_content.strip().startswith('<?xml'):
        root = ET.fromstring(xml_content)
    else:
        tree = ET.parse(xml_content)
        root = tree.get_root()
    
    coefficients = {}
    sensor_counters = {
        'Temperature': 1,
        'Conductivity': 1,
        'Oxygen': 1,
    }
    
    # Durch alle Sensoren iterieren
    sensor_array = root.find('.//SensorArray')
    if sensor_array is None:
        raise ValueError("SensorArray nicht in XML gefunden")
    
    for sensor in sensor_array.findall('Sensor'):
        sensor_type = None
        sensor_element = None
        sensor_index = int(sensor.get('index', -1))
        
        # Sensor-Typ identifizieren
        for child in sensor:
            tag = child.tag
            if 'TemperatureSensor' in tag:
                sensor_type = 'Temperature'
                sensor_element = child
            elif 'ConductivitySensor' in tag:
                sensor_type = 'Conductivity'
                sensor_element = child
            elif 'PressureSensor' in tag:
                sensor_type = 'Pressure'
                sensor_element = child
            elif 'OxygenSensor' in tag:
                sensor_type = 'Oxygen'
                sensor_element = child
            elif 'PAR_BiosphericalLicorChelseaSensor' in tag:
                sensor_type = 'PAR'
                sensor_element = child
            elif 'AltimeterSensor' in tag:
                sensor_type = 'Altimeter'
                sensor_element = child
            elif 'FluoroWetlabECO_AFL_FL_Sensor' in tag:
                sensor_type = 'ECO_Fluoro'
                sensor_element = child
            elif 'TurbidityMeter' in tag:
                sensor_type = 'ECO_Turbidity'
                sensor_element = child
            elif 'FluoroWetlabCDOM_Sensor' in tag:
                sensor_type = 'ECO_CDOM'
                sensor_element = child
            elif 'SPAR_Sensor' in tag:
                sensor_type = 'SPAR'
                sensor_element = child
            elif 'NotInUse' in tag:
                continue  # Nicht verwendete Sensoren überspringen
        
        if sensor_element is None:
            continue
        
        try:
            # Temperatursensoren
            if sensor_type == 'Temperature':
                counter = sensor_counters['Temperature']
                key = f"temp{counter}_coefs"
                
                coefs = TemperatureCoefficients()
                coefs.a0 = float(sensor_element.find('G').text)
                coefs.a1 = float(sensor_element.find('H').text)
                coefs.a2 = float(sensor_element.find('I').text)
                coefs.a3 = float(sensor_element.find('J').text)
                coefs.index = sensor_index
                
                coefficients[key] = coefs
                sensor_counters['Temperature'] += 1
            
            # Leitfähigkeitssensoren
            elif sensor_type == 'Conductivity':
                counter = sensor_counters['Conductivity']
                key = f"cond{counter}_coefs"
                
                # Verwende equation="1" Koeffizienten
                coefs_element = sensor_element.find('Coefficients[@equation="1"]')
                if coefs_element is not None:
                    coefs = ConductivityCoefficients()
                    coefs.g = float(coefs_element.find('G').text)
                    coefs.h = float(coefs_element.find('H').text)
                    coefs.i = float(coefs_element.find('I').text)
                    coefs.j = float(coefs_element.find('J').text)
                    coefs.cpcor = float(coefs_element.find('CPcor').text)
                    coefs.ctcor = float(coefs_element.find('CTcor').text)
                    wbotc_elem = coefs_element.find('WBOTC')
                    coefs.wbotc = float(wbotc_elem.text) if wbotc_elem is not None else 0.0
                    coefs.index = sensor_index
                    
                    coefficients[key] = coefs
                    sensor_counters['Conductivity'] += 1
            
            # Drucksensor
            elif sensor_type == 'Pressure':
                coefs = PressureCoefficients()
                coefs.C1 = float(sensor_element.find('C1').text)
                coefs.C2 = float(sensor_element.find('C2').text)
                coefs.C3 = float(sensor_element.find('C3').text)
                coefs.D1 = float(sensor_element.find('D1').text)
                coefs.D2 = float(sensor_element.find('D2').text)
                coefs.T1 = float(sensor_element.find('T1').text)
                coefs.T2 = float(sensor_element.find('T2').text)
                coefs.T3 = float(sensor_element.find('T3').text)
                coefs.T4 = float(sensor_element.find('T4').text)
                coefs.T5 = float(sensor_element.find('T5').text)
                coefs.Slope = float(sensor_element.find('Slope').text)
                coefs.Offset = float(sensor_element.find('Offset').text)
                coefs.AD590M = float(sensor_element.find('AD590M').text)
                coefs.AD590B = float(sensor_element.find('AD590B').text)
                coefs.index = sensor_index
                
                coefficients['pres_coefs'] = coefs
            
            # Sauerstoffsensoren
            elif sensor_type == 'Oxygen':
                counter = sensor_counters['Oxygen']
                key = f"oxy{counter}_coefs"
                
                # Verwende equation="1" Koeffizienten (Sea-Bird equation 2007+)
                coefs_element = sensor_element.find('CalibrationCoefficients[@equation="1"]')
                if coefs_element is not None:
                    coefs = Oxygen43Coefficients()
                    coefs.soc = float(coefs_element.find('Soc').text)
                    coefs.v_offset = float(coefs_element.find('offset').text)
                    coefs.tau_20 = float(coefs_element.find('Tau20').text)
                    coefs.a = float(coefs_element.find('A').text)
                    coefs.b = float(coefs_element.find('B').text)
                    coefs.c = float(coefs_element.find('C').text)
                    coefs.e = float(coefs_element.find('E').text)
                    coefs.d0 = float(coefs_element.find('D0').text)
                    coefs.d1 = float(coefs_element.find('D1').text)
                    coefs.d2 = float(coefs_element.find('D2').text)
                    coefs.h1 = float(coefs_element.find('H1').text)
                    coefs.h2 = float(coefs_element.find('H2').text)
                    coefs.h3 = float(coefs_element.find('H3').text)
                    coefs.index = sensor_index
                    
                    coefficients[key] = coefs
                    sensor_counters['Oxygen'] += 1
            
            # PAR Sensor
            elif sensor_type == 'PAR':
                coefs = PARCoefficients()
                coefs.im = float(sensor_element.find('M').text)
                coefs.a0 = float(sensor_element.find('CalibrationConstant').text)
                coefs.a1 = float(sensor_element.find('Offset').text)
                coefs.multiplier = float(sensor_element.find('Multiplier').text)
                coefs.index = sensor_index
                
                coefficients['par_coefs'] = coefs
            
            # Altimeter
            elif sensor_type == 'Altimeter':
                coefs = Altimeter_Coefficients()
                coefs.scalefactor = float(sensor_element.find('ScaleFactor').text)
                coefs.offset = float(sensor_element.find('Offset').text)
                coefs.index = sensor_index
                
                coefficients['ali_coefs'] = coefs
            
            # Fluorometer ECO
            elif sensor_type == 'ECO_Fluoro':
                coefs = ECOCoefficients()
                coefs.slope = float(sensor_element.find('ScaleFactor').text)
                coefs.offset = float(sensor_element.find('Vblank').text)
                coefs.index = sensor_index
                
                coefficients['fluroeco_coefs'] = coefs
            
            # Turbidity ECO
            elif sensor_type == 'ECO_Turbidity':
                coefs = ECOCoefficients()
                coefs.slope = float(sensor_element.find('ScaleFactor').text)
                coefs.offset = float(sensor_element.find('DarkVoltage').text)
                coefs.index = sensor_index
                
                coefficients['turbidity_coefs'] = coefs
            
            # CDOM ECO
            elif sensor_type == 'ECO_CDOM':                
                coefs = ECOCoefficients()
                coefs.slope = float(sensor_element.find('ScaleFactor').text)
                coefs.offset = float(sensor_element.find('Vblank').text)
                coefs.index = sensor_index
                
                coefficients['flurocdom_coefs'] = coefs
            
            # SPAR Sensor
            elif sensor_type == 'SPAR':
                coefs = SPAR_Coefficients()
                coefs.conversionfactor = float(sensor_element.find('ConversionFactor').text)
                coefs.ratiomultiplier = float(sensor_element.find('RatioMultiplier').text)
                coefs.index = sensor_index
                
                coefficients['spar_coefs'] = coefs
                
        except Exception as e:
            print(f"Fehler beim Verarbeiten von Sensor {sensor_type} an Index {sensor.get('index')}: {e}")
            continue
    
    return coefficients

# Debug-Funktion um alle verfügbaren Keys anzuzeigen
def print_available_coefficients(coefficients):
    print("Verfügbare Koeffizienten:")
    for key, coefs in coefficients.items():
        # Falls das Objekt einen Index hat, zeige ihn an
        index = getattr(coefs, "index", None)
        if index is not None:
            print(f"  - {key}  -->  Index: {index}")
        else:
            print(f"  - {key}")