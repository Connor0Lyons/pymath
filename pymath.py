import cmath
from fractions import Fraction
import math
from random import *
"""
import CoolProp.CoolProp as cp
import matplotlib.pyplot as plt
"""

"""
import pint
unit = pint.UnitRegistry()
u = unit
"""

"""
import scipy as sp
from scipy.constants import *
"""

try:
    import numpy as np
    from numpy import matrix
    from numpy.linalg import *
except:
    print("Warning: Error importing NumPy. Some features may not work as intended")

from cmath import *
from math import *

#print("Using pymath")

#Universal constants:   ------------------------------------------------------------------------------------------
pi = math.pi                    #Circle constant pi (unit-less)
pi2 = pi**2                     #pi^2
pi3 = pi**3                     #pi^3
c = 299_792_458                 #Speed of light (m/s) - exact
c2 = c**2                       #Speed of light squared c^2 (m^2 / s^2)
c3 = c**3                       #Speed of light cubed c^3 (m^3 / s^3)
c4 = c**4                       #Speed of light to the fourth power c^4 (m^4 / s^4)
me = 9.109_383_701_5e-31        #Mass of electron (kg)
mp = 1.672_621_923_69e-27       #Mass of proton (kg)
mn = 1.674_927_498_04e-27       #Mass of neutron (kg)
mHe = 6.644_657_335_7e-27       #Mass of alpha particle ^4He^2+ (kg)
e = 1.602_176_634e-19           #Elementary charge (C) - exact
Na = 6.022_140_76e23            #Avogadro constant (1 / mol)
h = 6.626_070_15e-34            #Planck's constant (m^2 kg/s) - exact
hbar = h / (2*pi)               #Reduced Planck constant (m^2 kg/s)
kB = 1.380_649e-23              #Boltzmann constant (J / K) - exact
u0 = 4 * pi * 10**-7            #Vacuum permeability AKA Magnetic constant mu_not (H/m)
e0 = 1 / (u0 * c2)              #Vacuum permittivity epsilon_not (F / m)
k = 1 / (4 * pi * e0)           #Coulomb constant (N * m^2 / C^2)
G = 6.674_30e-11                #Gravitational constant (m^3 / kg s^2)
g = 9.81                        #Acceleration due to gravity (m / s^2)
uB = e * hbar / (2 * me)        #Bohr magneton mu_B (J / T)
uN = e * hbar / (2 * mp)        #Nuclear magneton mu_N (J/T)
H0 = 2.333_361e-18              #Hubble constant approx = 72 [km/s] / Mpc  (1/s) - inexact
msun = 1.988_47e30              #Mass of sun (kg)
F = Na * e                      #Faraday constant (C / mol)
Z0 = u0 * c                     #Characteristic impedance of vacuum or Impedance of free space (Ohms)
R = Na * kB                     #Molar gas constant AKA Universal gas constant (J / K mol)
alpha = e**2 / (2 * e0 * h * c)                     #Fine-structure constant alpha (unit-less)
a0 = 4 * pi * e0 * hbar**2 / (me * e**2)            #Bohr radius a_not (m)
Rinf = alpha**2 * me * c / (2 * h)                  #Rydberg constant R_infinity (1 / m)
sigma = 2 * pi**5 * kB**4 / (15 * h**3 * c2)        #Stefan-Boltzmann constant (W / m^2 K^4)

Rair = 287.05                   #Individual Gas Constant of Air [J / K kg]
Rwater = 461.52                 #Individual Gas Constant of Water Vapor [J / K kg]
Roxygen = 259.84                #Individual Gas Constant of Oxygen O2 [J / K kg]
Rnitrogen = 296.80              #Individual Gas Constant of Nitrogen N2 [J / K kg]
RCO2 = 118.92                   #Individual Gas Constant of Carbon Dioxide CO2 [J / K kg] 


#Conversion Factors:    ------------------------------------------------------------------------------------------
zepto = 1e-21                   #SI Small Prefixes
atto = 1e-18
femto = 1e-15
pico = 1e-12
nano = 1e-9
micro = 1e-6
milli = 1e-3
centi = 1e-2                    #--------------------
hecto = 1e2                     #SI Big Prefixes
kilo = 1e3
mega = 1e6
giga = 1e9
tera = 1e12
peta = 1e15
exa = 1e18
zetta = 1e21                    #--------------------

kibi = 2**10                    #Binary prefixes
mebi = 2**20
gibi = 2**30
tebi = 2**40
pebi = 2**50
exbi = 2**60
zebi = 2**70
yobi = 2**80                    #--------------------

# deg = pi / 180                #Degrees to radians factor (rad / degree)       #DEPRECIATED in favor of degrees function
MeV = e * 10**6                 #MeV conversion factor (J / MeV)
amu = 1.660_539_066_6e-27       #Atomic mass unit AMU conversion factor (kg / AMU) 
mol = Na                        #Mol to molecules conversion factor (molecules / mol)
min = 60                        #Time conversion (s / min)
hr = 60 * min                   #Time conversion (s / hr)
day = 24 * hr                   #Time conversion (s / day)
week = 7 * day                  #Time conversion (s / week)
month = 365 / 12 * day          #Time conversion (s / month)
yr = 365 * day                  #Time conversion (s / yr)
IN = 0.0254                     #Imperial length conversion (m / in)      #'in' is a reserved word in python so using IN instead
in2 = IN**2                     #Imperial area conversion (m^2 / in^2)
in3 = IN**3                     #Imperial volume conversion (m^3 / in^3)
in4 = IN**4                     #Imperial hyper volume conversion (m^4 / in^4)
ft = 0.3048                     #Imperial length conversion (m / ft)
ft2 = ft**2                     #Imperial area conversion (m^2 / ft^2)
ft3 = ft**3                     #Imperial volume conversion (m^3 / ft^3)
ft4 = ft**4                     #Imperial hyper volume conversion (m^4 / ft^4)
yd = 0.9144                     #Imperial length conversion (m / yd)
yd2 = yd**2                     #Imperial area conversion (m^2 / yd^2)
yd3 = yd**3                     #Imperial volume conversion (m^3 / yd^3)
yd4 = yd**4                     #Imperial hyper volume conversion (m^4 / yd^4)
mi = 1_609.344                  #Imperial length conversion (m / mi)
mi2 = mi**2                     #Imperial area conversion (m^2 / mi^2)
mi3 = mi**3                     #Imperial volume conversion (m^3 / mi^3)
mi4 = mi**4                     #Imperial hyper volume conversion (m^4 / mi^4)
kmh = 1e3 / hr                  #Km / hour -> m/s   ([m/s] / kmh)
mph = mi / hr                   #Miles per hour to m/s ([m/s] / mph)
au = 149_597_870_700            #Astronomical unit to meter (m / au)
pc =  3.085_677_581_28e16       #Parsec to meter   (m / Pc)
Mpc = 1e6 * pc                  #Mega parsecs to meter (m / MPc)
ly = c * 365.25 * day           #Lightyear to meter (m / ly)
mach = 340.5                    #one Mach (approx., at 15 C, 1 atm) in meters per second ([m/s] / mach)
angstrom = 1e-10                #Angstrom to meters (m / A)
lbf = 4.448_221_6               #Pound force to Newtons conversion (N / lbf)
kip = kilo * lbf                #KiloPound force to Newtons conversion (N / kip)
liter = 0.001                   #liter to m^3 conversion (m^3 / L)
mL = milli * liter              #milliliter to m^3 conversion (m^3 / mL)
bar = 100_000                   #bar to pascal conversion (pa / bar) = ([N / m^2] / bar)
atm = 101_325                   #standard atmosphere to pascal conversion ([N / m^2] / atm)
psi = lbf / in2                 #Pounds per square inch to pascal (pa / psi) = ([N / m^2] / [lbf / in^2])
ksi = kilo * psi                #Kips per square inch to pascal (pa / ksi) = ([N / m^2] / [kilolbf / in^2])
btu = 1_055.055_852_6           #British thermal unit to Joules conversion (J / btu)
therm = btu * 100_000           #therm to Joules conversion (J / therm)
lbm = 0.453_592_37              #Pound to kilogram conversion (kg / lbm)
ton = 2000 * lbf                #Ton to newton conversion (N/ ton)
hp = 745.75                     #Horsepower to Watt conversion (W / hp)
acre = 4046.856                 #acre to m^2 conversion (m^2 / acre)
ha = hecto * acre               #hectare to m^2 conversion (m^2 / ha)
floz = 28.413 * mL              #fluid ounce to m^3 conversion (m^3 / fl oz)
pint = 568.261 * mL             #pint to m^3 conversion (m^3 / pint)
quart = 1236.523 * mL           #quart to m^3 conversion (m^3 / quart) 
gallon = 4546.09 * mL           #gallon to m^3 conversion (m^3 / gal)
oz = 0.028_350                  #ounce to kg conversion (kg / oz)
rankine = 5/9                   #Rankine t conversion (K / R = C / F)
arcmin = radians(1/60)          #Arc Minute to Radians conversion (rad / arcmin)
arcsec = arcmin / 60            #Arc Second to Radians conversion (rad / arcsec)
mhr = milli * hr                #milli hours to seconds (s / mhr )
khr = kilo * hr                 #kilo hours to seconds (s / khr )
Mhr = mega * hr                 #Mega hours to seconds (s / Mhr )
cm2 = centi**2                  #area conversion (m^2 / cm^2)
cm3 = centi**3                  #volume conversion (m^3 / cm^3)
cm4 = centi**4                  #hyper volume conversion (m^4 / cm^4)
mm2 = milli**2                  #area conversion (m^2 / mm^2)
mm3 = milli**3                  #volume conversion (m^3 / mm^3)
mm4 = milli**4                  #hyper volume conversion (m^4 / mm^4)
celsius = 273.15                #Kelvin-celsius offset ([Kelvin] + celsius = [celsius] ; [celsius] - celsius = [Kelvin] )
calorie = 4.184                 #calorie to Joule conversion (J / calorie)
kcal = 4184                     #Calorie to Joule conversion (J / kcal)
rps = 2 * pi                    #Rotations per second to radians per second ([rad]/s / rot/s)
rpm = 2 * pi / 60               #Rotations per minute to radians per second ([rad]/s / rot/min)
lbft = lbf * ft                 #Pound force * feet to Newtons * meter ([N * m] / [lbf * ft])
lbin = lbf * IN                 #Pound force * inch to Newtons * meter ([N * m] / [lbf * in])
kipft = kip * ft                #Kilo pound force * feet to Newtons * meter ([N * m] / [kip * ft])
kipin = kip * IN                #Kilo pound force * inch to Newtons * meter ([N * m] / [kip * in])
barn = 1e-28                    #Barn to m^2 conversion (m^2 / b)
fb = femto * barn               #Femto Barn to m^2 conversion (m^2 / fb)

#Constant Aliases:    --------------------------------------------------------------------------------------------
m_e = me
m_p = mp
m_n = mn
ma, m_a, malpha, m_alpha, mhe, m_He, m_he = mHe , mHe, mHe, mHe, mHe, mHe, mHe          #Warning: potential conflict with mA = milli-Amp != ma
qe, q_e, Qe, Q_e, eV, ev = e, e, e, e, e, e                 #Warning: potential conflict with exponential growth constant e^x
N_a, NA, N_A, Avogadro, avogadro = Na, Na, Na, Na, Na
k_B, kb = kB, kB                                            #Warning: potential conflict with kibi = kilobyte 
mev, Mev = MeV, MeV
AMU, da, dalton = amu, amu, amu
u_0, mu0, mu_0 = u0, u0, u0
e_0, epsilon_not, epsilon_0, eps_0, eps0 = e0, e0, e0, e0, e0
u_B, muB, mu_B = uB, uB, uB
u_N, muN, mu_N = uN, uN, uN
mSun, m_sun, m_Sun = msun, msun, msun
Z_0, z0, z_0 = Z0, Z0, Z0
Rbar, rbar, gasconst, gasconstant = R, R, R, R
a_0 = a0
R_inf = Rinf

RAir = Rair
Rh2o, RH2o, RH2O, RH20, Rh20 = Rwater, Rwater, Rwater, Rwater, Rwater
RO2, Ro2, Roxygen = Roxygen, Roxygen, Roxygen
RCarbonDioxide, RCarbon_Dioxide, RcarbonDioxide, Rcarbon_Dioxide, Rcarbondioxide, Rcarbon_dioxide, Rco2, RCo2, RC02, rc02 = RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2

Zepto, zm, zs = zepto, zepto, zepto
Atto, am = atto, atto                           #'as' is reserved word in python
Femto, fm, fs, fempto, Fempto = (femto,)*5
Pico, pm, ps = pico, pico, pico
Nano, nm, ns, ppb, PPB, ppB = nano, nano, nano, nano, nano, nano
Micro, mum, mus, mu_m, mu_s, muC, muc, muF, muf, muA, um, us, uC, uc, uF, uf, uA, ppm, PPM, ppM = micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro
Milli, mm, ms, mA, mJ, mF, mf, mC, mv, mV, gram = milli, milli, milli, milli, milli, milli, milli, milli, milli, milli, milli            #Warning: potential conflict with ma = mass of alpha particle != mA
Centi, cm, cs = centi, centi, centi
Hecto, hm, hs = hecto, hecto, hecto
Kilo, km, ks, kPa, kpa, kW, kw, kJ, kj, Mg, kN, kn, kohm, thousand = kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo
Mega, Mm, Ms, MPa, MPA, Mpa, MW, MJ, Mj, Gg, km2, Mohm, million, mil = mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega
Giga, Gm, Gs, GPa, GPA, Gpa, gpa, GW, Gw, gw, GJ, Gj, gj, km3, billion, bil = giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga
Tera, Tm, Ts, km4, trillion = tera, tera, tera, tera, tera
Peta, Pm, Ps, quadrillion = peta, peta, peta, peta
Exa, Em, Es, quintillion = exa, exa, exa, exa
Zetta, Zm, Zs, sextillion = zetta, zetta, zetta, zetta

Ki, KiB, Kib, kib, KB, kilobyte, kiloB   = kibi, kibi, kibi, kibi, kibi, kibi, kibi             #Warning: potential conflict with boltzman constant kB
Mi, MiB, Mib, mib, MB, megabyte, megaB   = mebi, mebi, mebi, mebi, mebi, mebi, mebi
Gi, GiB, Gib, gib, GB, gigabyte, gigaB   = gibi, gibi, gibi, gibi, gibi, gibi, gibi
Ti, TiB, Tib, tib, TB, terabyte, teraB   = tebi, tebi, tebi, tebi, tebi, tebi, tebi
Pi, PiB, Pib, pib, PB, petabyte, petaB   = pebi, pebi, pebi, pebi, pebi, pebi, pebi
Ei, EiB, Eib, eib, EB, exabyte, exaB     = exbi, exbi, exbi, exbi, exbi, exbi, exbi
Zi, ZiB, Zib, zib, ZB, zettabyte, zettaB = zebi, zebi, zebi, zebi, zebi, zebi, zebi
Yi, YiB, Yib, yib, YB, yottabyte, yottaB = yobi, yobi, yobi, yobi, yobi, yobi, yobi

hrs, hour, hours, Whr, Wh, Ahr, Ah = hr, hr, hr, hr, hr, hr, hr
days = day
weeks, Week, Weeks = week, week, week
mon, Mon, Month, months, Months = month, month, month, month, month
yrs, year, years = yr, yr, yr
In, inch, inches, iN, inn, in1 = IN, IN, IN, IN, IN, IN
In2, inch2, inches2, iN2, inn2, IN2 = in2, in2, in2, in2, in2, in2
In3, inch3, inches3, iN3, inn3, IN3 = in3, in3, in3, in3, in3, in3
In4, inch4, inches4, iN4, inn4, IN4 = in4, in4, in4, in4, in4, in4
feet, Ft = ft, ft
feet2, Ft2 = ft2, ft2
feet3, Ft3 = ft3, ft3
yard, yards =  yd, yd
Mi, miles, mile = mi, mi, mi
Mi2, sqmi, SqMi, sqMi = mi2, mi2, mi2, mi2
parsec, Pc = pc, pc
megaparsec, MegaParsec, mpc, MPc, megaParsec = Mpc, Mpc, Mpc, Mpc, Mpc
lightyear, Lightyear, Ly = ly, ly, ly
A, AA, Angstrom = angstrom, angstrom, angstrom
kips, Kips, Kip, klbf, KIP, KIPS, = kip, kip, kip, kip, kip, kip
l, L, Liter, lit, Lit = liter, liter, liter, liter, liter
ml = mL
PSI = psi
BTU, Btu = btu, btu
HP, Hp, horsepower = hp, hp, hp
FlOz, Floz = floz, floz
pt = pint
qt = quart
gal, gallons, Gal = gallon, gallon, gallon
ounce, Oz, ounces = oz, oz, oz
degR, degr, degF, degf, Fahrenheit, fahrenheit, fahr, fahren  = rankine, rankine, rankine, rankine, rankine, rankine, rankine, rankine
mh, mWhr, mwhr, mWh, mwh, mAhr, mAhr, mAh, mAh, mWH, mwH, mAHr, mAHr, mAH, mAH = mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr
kh, kWhr, kwhr, kWh, kwh, kAhr, kAhr, kAh, kAh, kWH, kwH, kAHr, kAHr, kAH, kAH = khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr
Mh, MWhr, Mwhr, MWh, Mwh, MAhr, MAhr, MAh, MAh, MWH, MwH, MAHr, MAHr, MAH, MAH = Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr
Cel, cel, Celsius, Celcius, celcius, = celsius, celsius, celsius, celsius, celsius
cal = calorie
Cal, Calorie = kcal, kcal
Rps, RPS = rps, rps
Rpm, RPM = rpm, rpm
lbfft, ftlb, ftlbf = lbft, lbft, lbft
lbfin, lbinn, lbfinn, lbIn, lbfIn, lbIN, lbfIN, lbinch, lbfinch, inlb, innlb, inlbf, innlbf = lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin
ftkip, Kipft = kipft, kipft
Kipin, Kipinn, kipinn, kipIn, kipIN, kipinch, inKip, inkip, innkip, Inkip = kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin
Barn, b = barn, barn
fB = fb


#Intra Program Variables
ans = [0]
ans1 = 0
ans2 = 0
degreesFlag = False                         #Default in Radians Mode
floatDelta = nano                           #Accepatable percent error for floating point number comparisions

def quad(a1, b1, c1):
    global ans, ans1, ans2

    discriminant = b1 ** 2 - 4 * a1 * c1
    
    if discriminant == 0:
        ans1 = -b1 / (2 * a1)
        -b1 / (2 * a1)
        ans.insert(0, ans1)
        print(ans1)
        return
        
    ans1 = (-b1 + sqrt(discriminant) ) / (2 * a1)
    ans2 = (-b1 - sqrt(discriminant) ) / (2 * a1)
    ans.insert(0, ans1)
    ans.insert(0, ans2)
    print(ans1)
    print(ans2)
    return

def sqrt(x):
    try:
        return math.sqrt(x)
    except:
        return cmath.sqrt(x)
    
#Converts a parameter from radians to degrees, or sets the default behavior of normal trig functions to be in degrees if no paramaters passed
def degrees(val = None):
    if(val != None):
        return math.degrees(val)
    
    global degreesFlag 
    degreesFlag = True
    print("Now in Degrees Mode")

#Converts a parameter from degrees to radians, or sets the default behavior of normal trig functions to be in radians if no paramaters passed
def radians(val = None):
    if(val != None):
        return math.radians(val)
    
    global degreesFlag
    degreesFlag = False
    print("Now in Radians Mode")

#Chase did this, don't entirely know how it works
# sind, cosd, tand = map(lambda f: lambda x: f(radians(x)), (math.sin, math.cos, math.tan))

#Degree based trig functions
def sind(x):
    try:
        return math.sin(radians(x))
    except:
        return cmath.sin(radians(x))

def cosd(x):
    try:
        return math.cos(radians(x))
    except:
        return cmath.cos(radians(x))

def tand(x):
    try:
        return math.tan(radians(x))
    except:
        return cmath.tan(radians(x))

#Radian based Trig Funcitons
def sinr(x):
    try:
        return math.sin(x)
    except:
        return cmath.sin(x)

def cosr(x):
    try:
        return math.cos(x)
    except:
        return cmath.cos(x)

def tanr(x):
    try:
        return math.tan(x)
    except:
        return cmath.tan(x)

#Trig Functions that check whether in degree mode or radian mode, and give warnings 
def sin(x):
    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDelta < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (sind(x) calculated anyway)")
        return sind(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180 
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDelta < 1 ) ):
        print("Warning: In Radian Mode but possible degree number detected. (sinr(x) calculated anyway)")
    return sinr(x)

def cos(x):
    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDelta < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (cosd(x) calculated anyway)")
        return cosd(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDelta < 1 ) ):
        print("Warning: In Radian Mode but possible degree number detected. (cosr(x) calculated anyway)")
    return cosr(x)

def tan(x):
    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDelta < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (tand(x) calculated anyway)")
        return tand(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDelta < 1 ) ):
        print("Warning: In Radian Mode but possible degree number detected. (tanr(x) calculated anyway)")
    return tanr(x)

#returns arcsin(x), returns radian value if in radians mode or if complex number
def asin(x):
        ans = cmath.asin(x)
        if(ans.imag != 0):
            return ans
        
        return degrees(ans.real) if degreesFlag else ans.real
    
#returns arccos(x), returns radian value if in radians mode or if complex number
def acos(x):
        ans = cmath.acos(x)
        if(ans.imag != 0):
            return ans
        
        return degrees(ans.real) if degreesFlag else ans.real
    
#returns arctan(x), returns radian value if in radians mode or if complex number
def atan(x):
        ans = cmath.atan(x)
        if(ans.imag != 0):
            return ans
        
        return degrees(ans.real) if degreesFlag else ans.real

sec = lambda x: 1/cos(x)
csc = lambda x: 1/sin(x)
cot = lambda x: 1/tan(x)

sin2 = lambda x: sin(x)**2
cos2 = lambda x: cos(x)**2
tan2 = lambda x: sin(x)**2
sec2 = lambda x: cos(x)**-2
csc2 = lambda x: sin(x)**-2
cot2 = lambda x: tan(x)**-2

def ln(x):
    try:
        return math.log(x)
    except:
        return cmath.log(x)
    
def exp(x):
    try:
        return math.exp(x)
    except:
        return cmath.exp(x)

def lg(x):
    try:
        return math.log2(x)
    except:
        return cmath.log(x) / cmath.log(2)

E = lambda x : 10**x

sci = lambda x : f"{x : .4e}"       #Formats input in scientiific notation with 4 decimal places

#Temperature conversion functions

k2c = lambda x : x - 273.15                     #Kelvin to Celsius conversion
k2f = lambda x : (x - 273.15) * 9/5 + 32        #Kelvin to Fahrenheit conversion
k2r = lambda x : x * 1.8                        #Kelvin to Rankine conversion

c2k = lambda x : x + 273.15                     #Celsius to Kelvin conversion
c2f = lambda x : (x * 9/5) + 32                 #Celsius to Fahrenheit conversion
c2r = lambda x : (x * 9/5) + 491.67             #Celsius to Rankine conversion

f2k = lambda x : (x - 32) * 5/9 + 273.15        #Fahrenheit to Kelvin conversion
f2c = lambda x : (x - 32) * 5/9                 #Fahrenheit to Celcius conversion
f2r = lambda x : x + 459.67                     #Fahrenheit to Rankine conversion

r2k = lambda x : x * 5/9                        #Rankine to Kelvin conversion
r2c = lambda x : (x - 491.67) * 5/9             #Rankine to Celcius conversion
r2f = lambda x : x - 459.67                     #Rankine to Fahrenheit conversion

#area of circle given diameter
darea = lambda d : pi * (d/2)**2 

def factor(num: int):
    for i in range(1, math.floor(sqrt(num)+1) ):
        if( (num / i).is_integer() ):
            print(f"{i} * {(num / i) : .0f}")

def frac(num: float):
    return Fraction(num).limit_denominator()

det = lambda mat : np.linalg.det(mat)

def innerprod(array1:tuple, array2:tuple):
    sum = 0
    for i,j in zip(array1, array2):
        sum += np.conj(i) * j    
    return sum

def norm(array1:tuple):
    sum = 0
    for i in array1:
        sum += np.conj(i) * i    
    return math.sqrt(sum.real)

#Integrates function/lamda expression f in interval (a, b) using equal step sizes and trapezoid summation
#Example use: int(lambda x : 5 * x**2 + 10, 0, 5)
def int(f, a, b):
    steps = 1_000_000
    deltaX = (b-a)/steps
    sum = 0
    for i in range(steps):
        sum += deltaX * (f(a) + f(a+deltaX) ) / 2
        a += deltaX
    return sum

""" 
def integral_f(f, initial_step_size):
    \"""
    This algorithm calculates the definite integral of a function
    from 0 to 1, adaptively, by choosing smaller steps near
    problematic points.
    \"""
    x = 0.0
    h = initial_step_size
    accumulator = 0.0
    while x < 1.0:
        if x + h > 1.0:
            h = 1.0 - x  # At end of unit interval, adjust last step to end at 1.
        if error_too_big_in_quadrature_of_f_over_range(f, [x, x + h]):
            h = make_h_smaller(h)
        else:
            accumulator += quadrature_of_f_over_range(f, [x, x + h])
            x += h
            if error_too_small_in_quadrature_of_over_range(f, [x, x + h]):
                h = make_h_larger(h)  # Avoid wasting time on tiny steps.
    return accumulator
"""

#Returns sample standard deviation by default
def std(arr, population=False):
    if(population):
        return np.std(arr,ddof=1)
    return np.std(arr)

def popstd(arr):
    return np.std(arr)

#Returns sample variance by default
def var(arr, population=False):
    if(population):
        return np.var(arr,ddof=1)
    return np.var(arr)

def popvar(arr):
    return np.var(arr)



#Given an arbritrary number of days, returns the number of people required for Probability[two people share the same birthday] >= 0.5   
def birthdayProblem(days=365):
    return ceil( sqrt(2*days*ln(2)) + (3 - 2 * ln(2))/6 + (9 - 4*ln(2)**2 ) / (72 * sqrt(2*days*ln(2))) - 2 * ln(2)**2 / (135*days) )


#Function and Lambda Aliases
deg, Deg, Degrees, degreesMode, degreesmode, DegreesMode, DegreesMode, degMode, degmode = degrees, degrees, degrees, degrees, degrees, degrees, degrees, degrees, degrees
rad, Rad, Radians, radiansMode, radiansmode, Radiansmode, RadiansMode, radMode, radmode = radians, radians, radians, radians, radians, radians, radians, radians, radians
arcsin, sininverse, sininv = asin, asin, asin
arccos, cosinverse, cosinv = acos, acos, acos
arctan, taninverse, taninv = atan, atan, atan
K2c, K2C, k2C = k2c, k2c, k2c
K2f, K2F, k2F = k2f, k2f, k2f
K2r, K2R, k2R = k2r, k2r, k2r
C2k, C2K, c2K = c2k, c2k, c2k
C2f, C2F, c2F = c2f, c2f, c2f
C2r, C2R, c2R = c2r, c2r, c2r
F2k, F2K, f2K = f2k, f2k, f2k
F2c, F2C, f2C = f2c, f2c, f2c
F2r, F2R, f2R = f2r, f2r, f2r
R2k, R2K, r2K = r2k, r2k, r2k
R2c, R2C, r2C = r2c, r2c, r2c
R2f, R2F, r2F = r2f, r2f, r2f
InnerProduct, innerproduct, Innerproduct, inprod, inproduct, braket, BraKet, Braket = innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod
integral, integrate = int, int
stdDev, stddev, StdDev, standardDeviation = std, std, std, std
popstdDev, popstddev, popStdDev, populationStddev, populationStdDev, populationStd, populationstd, populationstandardDeviation, populationStandardDeviation = popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd
variance = var
populationVariance, popvariance, popVariance = popvar, popvar, popvar
birthdayproblem, birthdayparadox, birthdayParadox, generalizedBirthdayProblem = birthdayProblem, birthdayProblem, birthdayProblem, birthdayProblem