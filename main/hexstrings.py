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