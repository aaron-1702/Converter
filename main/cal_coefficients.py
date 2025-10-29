# Seabird (https://github.com/Sea-BirdScientific/seabirdscientific/tree/main)
class TemperatureCoefficients:
    """
    :param a0: coefficient
    :param a1: coefficient
    :param a2: coefficient
    :param a3: coefficient
    """

    a0: float
    a1: float
    a2: float
    a3: float

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

class ConductivityCoefficients:
    """
    :param g: coefficient
    :param h: coefficient
    :param i: coefficient
    :param j: coefficient
    :param cpcor: compressibility coefficient
    :param ctcor: coefficient
    :param wbotc: bridge oscillator temperature coefficient see the
        37 Manual: https://www.seabird.com/asset-get.download.jsa?id=54627862348
    """

    g: float
    h: float
    i: float
    j: float
    cpcor: float
    ctcor: float
    wbotc: float

class Oxygen43Coefficients:
    """
    :param soc: linear scaling coefficient
    :param v_offset: voltage at zero oxygen signal
    :param tau_20: sensor time constant tau(T,P) at 20 C, 1 atmosphere, 0 PSU;
        slope term in calculation of tau(T,P)
    :param a: coefficient
    :param b: coefficient
    :param c: coefficient
    :param e: coefficient
    :param d0: tau(T,P) coefficient
    :param d1: tau(T,P) coefficient
    :param d2: tau(T,P) coefficient
    :param h1: hysteresis correction coefficient
    :param h2: hysteresis correction coefficient
    :param h3: hysteresis correction coefficient
    """

    soc: float
    v_offset: float
    tau_20: float
    a: float
    b: float
    c: float
    e: float
    d0: float  # TODO: the sensor lists a d0 coefficient, but it doesn't seem to be used?
    d1: float
    d2: float
    h1: float
    h2: float
    h3: float

class PARCoefficients:
    """
    :param im: immersion coefficient
    :param a0: calibration slope
    :param a1: calibration offset
    :param multiplier: 1.0 for units of μEinsteins/m2 sec
    """

    im: float
    a0: float
    a1: float
    multiplier: float

class ECOCoefficients:
    """
    :param slope: units/count for digital, units/V for analog
    :param offset: dark counts for digital, dark voltage for analog
    """

    slope: float
    offset: float

# Selfmade
class Altimeter_Coefficients:
    scalefactor: float
    offset: float


class SPAR_Coefficients:
    conversionfactor: float
    ratiomultiplier: float