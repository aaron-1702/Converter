import re
import numpy as np


def load_hexstrings_simple(file_path):
    """
    Einfache Methode zum Laden von Hexstrings - nimmt alle Zeilen die nur Hex-Zeichen enthalten
    """
    hexstrings = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Überspringe leere Zeilen und Kommentarzeilen
            if not line or line.startswith('*'):
                continue
            
            # Prüfe ob die Zeile nur Hex-Zeichen enthält (0-9, A-F, a-f)
            if all(c in '0123456789ABCDEFabcdef' for c in line):
                hexstrings.append(line.upper())  # Konvertiere zu Großbuchstaben für Konsistenz
    
    return hexstrings 

def load_header_hex(file_path, filename='basic_emb_test.cnv'):
    """
    Lädt den Header aus der Hex-Datei und schreibt ihn in basic_emb_test.cnv
    Erstellt die Datei falls sie nicht existiert
    """
    with open(file_path, 'r') as input_file, open(filename, 'w') as output_file:
        for line in input_file:
            if line.startswith('*END*'):
                break
            if not line[0].isdigit() and not line[0].isalpha():
                output_file.write(line)

def load_header_XMLCON(file_path, filename='basic_emb_test.cnv'):
    with open(file_path, 'r') as xml_file, open(filename, 'a') as cnv_file:
        content = xml_file.read()
        sensors = re.findall(r'<Sensor index="\d+" SensorID="\d+" >.*?</Sensor>', content, re.DOTALL)
        for sensor in sensors:
            for line in sensor.splitlines():
                cnv_file.write(f"# {line}\n")



def write_arrays_to_file(*arrays, filename='basic_emb_test.cnv', decimals=6, total_width=12):
    """
    Schreibt mehrere NumPy-Arrays spaltenweise in eine Datei.
    
    Parameters:
    -----------
    *arrays : variable number of numpy arrays
        NumPy-Arrays als separate Parameter (müssen alle die gleiche Länge haben)
    filename : str
        Name der Ausgabedatei
    decimals : int
        Anzahl der Nachkommastellen
    total_width : int
        Gesamtbreite pro Wert in Zeichen
    """
    if not arrays:
        return
    
    # Überprüfe ob alle Arrays die gleiche Länge haben
    lengths = [len(arr) for arr in arrays]
    if len(set(lengths)) > 1:
        # Finde die minimale Länge und schneide alle Arrays darauf zu
        min_length = min(lengths)
        arrays = [arr[:min_length] for arr in arrays]
        print(f"Warnung: Arrays haben unterschiedliche Längen. Verwende minimale Länge: {min_length}")
    
    n_rows = len(arrays[0])
    n_cols = len(arrays)
    
    with open(filename, 'a') as file:
        # Schreibe Daten zeilenweise
        for i in range(n_rows):
            line = ""
            for j in range(n_cols):
                value = arrays[j][i]
                
                # Entscheide ob Integer oder Float Formatierung
                if isinstance(value, (int, np.integer)):
                    # Integer: rechtsbündig
                    formatted_value = f"{value:>{total_width}d}"
                else:
                    # Float: formatieren mit fester Breite
                    formatted_value = f"{value:>{total_width}.{decimals}f}"
                
                line += formatted_value
                
                # Füge Leerzeichen zwischen Spalten hinzu (außer letzte Spalte)
                if j < n_cols - 1:
                    line += " "
            
            file.write(line + "\n")


