import xml.etree.ElementTree as ET
from calibration_coefficients import (
    TemperatureCoefficients,
    PressureCoefficients,
    ConductivityCoefficients,
    Oxygen43Coefficients,
    PARCoefficients,
    ECOCoefficients,
    Altimeter_Coefficients,
    SPAR_Coefficients,
)


def read_xmlcon(xml_path: str):
    """
    Liest eine XMLCON-Datei (SeaBird CTD) und erzeugt automatisch
    die Kalibrierungskoeffizienten aus calibration_coefficients.py.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Hilfsfunktion: Text aus Tag lesen und zu float konvertieren
    def get_float(element, tag):
        val = element.findtext(tag)
        return float(val) if val is not None else None

    # Temperatur
    temp1_coefs = TemperatureCoefficients()
    temp_tag = root.find(".//TemperatureSensor")
    if temp_tag is not None:
        temp1_coefs.a0 = get_float(temp_tag, "A0")
        temp1_coefs.a1 = get_float(temp_tag, "A1")
        temp1_coefs.a2 = get_float(temp_tag, "A2")
        temp1_coefs.a3 = get_float(temp_tag, "A3")

    # Druck
    pres_coefs = PressureCoefficients()
    pres_tag = root.find(".//PressureSensor")
    if pres_tag is not None:
        pres_coefs.C1 = get_float(pres_tag, "C1")
        pres_coefs.C2 = get_float(pres_tag, "C2")
        pres_coefs.C3 = get_float(pres_tag, "C3")
        pres_coefs.D1 = get_float(pres_tag, "D1")
        pres_coefs.D2 = get_float(pres_tag, "D2")
        pres_coefs.T1 = get_float(pres_tag, "T1")
        pres_coefs.T2 = get_float(pres_tag, "T2")
        pres_coefs.T3 = get_float(pres_tag, "T3")
        pres_coefs.T4 = get_float(pres_tag, "T4")
        pres_coefs.T5 = get_float(pres_tag, "T5")
        pres_coefs.Slope = get_float(pres_tag, "Slope")
        pres_coefs.Offset = get_float(pres_tag, "Offset")
        pres_coefs.AD590M = get_float(pres_tag, "AD590M")
        pres_coefs.AD590B = get_float(pres_tag, "AD590B")

    # Leitfähigkeit
    cond1_coefs = ConductivityCoefficients()
    cond_tag = root.find(".//ConductivitySensor")
    if cond_tag is not None:
        cond1_coefs.g = get_float(cond_tag, "G")
        cond1_coefs.h = get_float(cond_tag, "H")
        cond1_coefs.i = get_float(cond_tag, "I")
        cond1_coefs.j = get_float(cond_tag, "J")
        cond1_coefs.cpcor = get_float(cond_tag, "CPCOR")
        cond1_coefs.ctcor = get_float(cond_tag, "CTCOR")
        cond1_coefs.wbotc = get_float(cond_tag, "WBOTC")

    # Sauerstoff
    oxy_coefs = Oxygen43Coefficients()
    oxy_tag = root.find(".//OxygenSensor")
    if oxy_tag is not None:
        oxy_coefs.soc = get_float(oxy_tag, "Soc")
        oxy_coefs.v_offset = get_float(oxy_tag, "Voffset")
        oxy_coefs.tau_20 = get_float(oxy_tag, "Tau20")
        oxy_coefs.a = get_float(oxy_tag, "A")
        oxy_coefs.b = get_float(oxy_tag, "B")
        oxy_coefs.c = get_float(oxy_tag, "C")
        oxy_coefs.e = get_float(oxy_tag, "E")
        oxy_coefs.d0 = get_float(oxy_tag, "D0")
        oxy_coefs.d1 = get_float(oxy_tag, "D1")
        oxy_coefs.d2 = get_float(oxy_tag, "D2")
        oxy_coefs.h1 = get_float(oxy_tag, "H1")
        oxy_coefs.h2 = get_float(oxy_tag, "H2")
        oxy_coefs.h3 = get_float(oxy_tag, "H3")

    # PAR
    par_coefs = PARCoefficients()
    par_tag = root.find(".//ParSensor")
    if par_tag is not None:
        par_coefs.im = get_float(par_tag, "Im")
        par_coefs.a0 = get_float(par_tag, "A0")
        par_coefs.a1 = get_float(par_tag, "A1")
        par_coefs.multiplier = get_float(par_tag, "Multiplier")

    # ECO (Fluoreszenz, Turbidität etc.)
    eco_coefs = ECOCoefficients()
    eco_tag = root.find(".//EcoSensor")
    if eco_tag is not None:
        eco_coefs.slope = get_float(eco_tag, "Slope")
        eco_coefs.offset = get_float(eco_tag, "Offset")

    # Altimeter
    alt_coefs = Altimeter_Coefficients()
    alt_tag = root.find(".//Altimeter")
    if alt_tag is not None:
        alt_coefs.scalefactor = get_float(alt_tag, "ScaleFactor")
        alt_coefs.offset = get_float(alt_tag, "Offset")

    # SPAR
    spar_coefs = SPAR_Coefficients()
    spar_tag = root.find(".//SparSensor")
    if spar_tag is not None:
        spar_coefs.conversionfactor = get_float(spar_tag, "ConversionFactor")
        spar_coefs.ratiomultiplier = get_float(spar_tag, "RatioMultiplier")

    return {
        "temperature": temp1_coefs,
        "pressure": pres_coefs,
        "conductivity": cond1_coefs,
        "oxygen": oxy_coefs,
        "par": par_coefs,
        "eco": eco_coefs,
        "altimeter": alt_coefs,
        "spar": spar_coefs,
    }
