# dataclasses.py
from dataclasses import dataclass
from typing import Literal, List

@dataclass
class TemperatureCoefficients:
    a0: float
    a1: float
    a2: float
    a3: float

@dataclass
class PressureCoefficients:
    C1: float
    C2: float
    C3: float
    D1: float
    D2: float
    T1: float
    T2: float
    T3: float
    T4: float
    T5: float
    Slope: float
    Offset: float
    AD590M: float
    AD590B: float

@dataclass
class ConductivityCoefficients:
    g: float
    h: float
    i: float
    j: float
    cpcor: float
    ctcor: float
    wbotc: float

@dataclass
class Oxygen43Coefficients:
    soc: float
    v_offset: float
    tau_20: float
    a: float
    b: float
    c: float
    e: float
    d0: float
    d1: float
    d2: float
    h1: float
    h2: float
    h3: float

@dataclass
class PARCoefficients:
    im: float
    a0: float
    a1: float
    multiplier: float

@dataclass
class ECOCoefficients:
    slope: float
    offset: float

@dataclass
class Altimeter_Coefficients:
    scalefactor: float
    offset: float

@dataclass
class SPAR_Coefficients:
    conversionfactor: float
    ratiomultiplier: float