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


class Altimeter_Coefficients:
    scalefactor: float
    offset: float


class SPAR_Coefficients:
    conversionfactor: float
    ratiomultiplier: float



temp1_coefs = TemperatureCoefficients()
temp1_coefs.a0 = 4.33459406e-003
temp1_coefs.a1 = 6.28735680e-004
temp1_coefs.a2 = 2.05782221e-005
temp1_coefs.a3 = 1.68624541e-006

temp2_coefs = TemperatureCoefficients()
temp2_coefs.a0 = 4.34471154e-003
temp2_coefs.a1 = 6.38146193e-004
temp2_coefs.a2 = 2.04525401e-005
temp2_coefs.a3 = 1.60292340e-006


pres_coefs = PressureCoefficients()
pres_coefs.C1 = -4.212464e+004
pres_coefs.C2 = -1.972709e-001
pres_coefs.C3 = 1.365930e-002

pres_coefs.D1 = 3.604100e-002
pres_coefs.D2 = 0.000000e+000

pres_coefs.T1 = 3.046128e+001
pres_coefs.T2 = -3.791930e-004
pres_coefs.T3 = 4.202260e-006
pres_coefs.T4 = 3.386190e-009
pres_coefs.T5 = 0.000000e+000

pres_coefs.Slope = 1.00000000
pres_coefs.Offset = -0.20000
pres_coefs.AD590M = 1.279410e-002
pres_coefs.AD590B = -9.215740e+000



cond1_coefs = ConductivityCoefficients()
cond1_coefs.g = -9.94900904e+000
cond1_coefs.h = 1.42164741e+000
cond1_coefs.i = 1.26008747e-004
cond1_coefs.j = 6.16820694e-005
cond1_coefs.cpcor = -9.57000000e-008
cond1_coefs.ctcor = 3.2500e-006
cond1_coefs.wbotc = 0.00000000e+000

cond2_coefs = ConductivityCoefficients()
cond2_coefs.g = -9.77824825e+000
cond2_coefs.h = 1.31877668e+000
cond2_coefs.i = 2.26556605e-004
cond2_coefs.j = 6.59127902e-005
cond2_coefs.cpcor = -9.57000000e-008
cond2_coefs.ctcor = 3.2500e-006
cond2_coefs.wbotc = 0.00000000e+000


oxy1_coefs = Oxygen43Coefficients()
oxy1_coefs.soc = 4.7612e-001
oxy1_coefs.v_offset = -0.4891
oxy1_coefs.tau_20 = 1.7400
oxy1_coefs.a = -3.2939e-003
oxy1_coefs.b = 1.4617e-004
oxy1_coefs.c = -2.1408e-006
oxy1_coefs.e = 3.6000e-002
oxy1_coefs.d0 = 2.5826e+000
oxy1_coefs.d1 = 1.92634e-004
oxy1_coefs.d2 = -4.64803e-002
oxy1_coefs.h1 = -3.3000e-002
oxy1_coefs.h2 = 5.0000e+003
oxy1_coefs.h3 = 1.4500e+003

oxy2_coefs = Oxygen43Coefficients()
oxy2_coefs.soc = 5.4771e-001
oxy2_coefs.v_offset = -0.5020
oxy2_coefs.tau_20 = 1.2900
oxy2_coefs.a = -3.9126e-003
oxy2_coefs.b = 1.5415e-004
oxy2_coefs.c = -2.2575e-006
oxy2_coefs.e = 3.6000e-002
oxy2_coefs.d0 = 2.5826e+000
oxy2_coefs.d1 = 1.92634e-004
oxy2_coefs.d2 = -4.64803e-002
oxy2_coefs.h1 = -3.3000e-002
oxy2_coefs.h2 = 5.0000e+003
oxy2_coefs.h3 = 1.4500e+003



par_coefs = PARCoefficients()
par_coefs.im = 1.00000000
par_coefs.a0 = 11273957159.00000000
par_coefs.a1 = -0.09022432
par_coefs.multiplier = 1.00000000

spar_coefs = SPAR_Coefficients()
spar_coefs.conversionfactor = 1.6067e+003
spar_coefs.ratiomultiplier = 1.00000000

fluroeco_coefs = ECOCoefficients()
fluroeco_coefs.slope = 6.00000000e+000
fluroeco_coefs.offset = 0.0690

turbidity_coefs = ECOCoefficients()
turbidity_coefs.slope = 2.000000
turbidity_coefs.offset = 0.065000

ali_coefs = Altimeter_Coefficients()
ali_coefs.scalefactor = 15.000
ali_coefs.offset = 0.000