# calibration_coefficients.py
from ctd_dataclasses import TemperatureCoefficients, PressureCoefficients, ConductivityCoefficients, \
                        Oxygen43Coefficients, PARCoefficients, ECOCoefficients, \
                        Altimeter_Coefficients, SPAR_Coefficients

# Temperature coefficients
temp1_coefs = TemperatureCoefficients(
    a0=4.33459406e-003,
    a1=6.28735680e-004,
    a2=2.05782221e-005,
    a3=1.68624541e-006
)

temp2_coefs = TemperatureCoefficients(
    a0=4.34471154e-003,
    a1=6.38146193e-004,
    a2=2.04525401e-005,
    a3=1.60292340e-006
)

# Pressure coefficients
pres_coefs = PressureCoefficients(
    C1=-4.212464e+004,
    C2=-1.972709e-001,
    C3=1.365930e-002,
    D1=3.604100e-002,
    D2=0.000000e+000,
    T1=3.046128e+001,
    T2=-3.791930e-004,
    T3=4.202260e-006,
    T4=3.386190e-009,
    T5=0.000000e+000,
    Slope=1.00000000,
    Offset=-0.20000,
    AD590M=1.279410e-002,
    AD590B=-9.215740e+000
)

# Conductivity coefficients
cond1_coefs = ConductivityCoefficients(
    g=-9.94900904e+000,
    h=1.42164741e+000,
    i=1.26008747e-004,
    j=6.16820694e-005,
    cpcor=-9.57000000e-008,
    ctcor=3.2500e-006,
    wbotc=0.00000000e+000
)

cond2_coefs = ConductivityCoefficients(
    g=-9.77824825e+000,
    h=1.31877668e+000,
    i=2.26556605e-004,
    j=6.59127902e-005,
    cpcor=-9.57000000e-008,
    ctcor=3.2500e-006,
    wbotc=0.00000000e+000
)

# Oxygen coefficients
oxy1_coefs = Oxygen43Coefficients(
    soc=4.7612e-001,
    v_offset=-0.4891,
    tau_20=1.7400,
    a=-3.2939e-003,
    b=1.4617e-004,
    c=-2.1408e-006,
    e=3.6000e-002,
    d0=2.5826e+000,
    d1=1.92634e-004,
    d2=-4.64803e-002,
    h1=-3.3000e-002,
    h2=5.0000e+003,
    h3=1.4500e+003
)

oxy2_coefs = Oxygen43Coefficients(
    soc=5.4771e-001,
    v_offset=-0.5020,
    tau_20=1.2900,
    a=-3.9126e-003,
    b=1.5415e-004,
    c=-2.2575e-006,
    e=3.6000e-002,
    d0=2.5826e+000,
    d1=1.92634e-004,
    d2=-4.64803e-002,
    h1=-3.3000e-002,
    h2=5.0000e+003,
    h3=1.4500e+003
)

# Other coefficients
par_coefs = PARCoefficients(
    im=1.00000000,
    a0=11273957159.00000000,
    a1=-0.09022432,
    multiplier=1.00000000
)

spar_coefs = SPAR_Coefficients(
    conversionfactor=1.6067e+003,
    ratiomultiplier=1.00000000
)

fluroeco_coefs = ECOCoefficients(
    slope=6.00000000e+000,
    offset=0.0690
)

turbidity_coefs = ECOCoefficients(
    slope=2.000000,
    offset=0.065000
)

ali_coefs = Altimeter_Coefficients(
    scalefactor=15.000,
    offset=0.000
)