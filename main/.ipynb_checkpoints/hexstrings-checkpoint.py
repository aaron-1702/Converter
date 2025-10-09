import re

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

def load_header_hex(file_path):
    """
    Lädt den Header aus der Hex-Datei und schreibt ihn in basic_emb_test.cnv
    Erstellt die Datei falls sie nicht existiert
    """
    with open(file_path, 'r') as input_file, open('basic_emb_test.cnv', 'w') as output_file:
        for line in input_file:
            if line.startswith('*END*'):
                break
            if not line[0].isdigit() and not line[0].isalpha():
                output_file.write(line)

def load_header_XMLCON(file_path):
    with open(file_path, 'r') as xml_file, open('basic_emb_test.cnv', 'a') as cnv_file:
        content = xml_file.read()
        sensors = re.findall(r'<Sensor index="\d+" SensorID="\d+" >.*?</Sensor>', content, re.DOTALL)
        for sensor in sensors:
            for line in sensor.splitlines():
                cnv_file.write(f"# {line}\n")


