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

try:
    import scipy as sp
    #from scipy.constants import *
    from scipy.optimize import fsolve
except:
    print("Warning: Error importing SciPy. Some features may not work as intended")

try:
    import numpy as np
    from numpy import matrix
    from numpy.linalg import *
except:
    print("Warning: Error importing NumPy. Some features may not work as intended")

from cmath import *
from math import *

#print("Using pymath")

#Unless noted otherwise, units are in terms of official SI units: second [s], meter [m], kilogram [kg], ampere [A], kelvin [K], mole [mol], and candela [cd]

#Universal constants:   ------------------------------------------------------------------------------------------
pi = math.pi                    #Circle constant pi [unit-less]
pi2 = pi**2                     #pi^2
pi3 = pi**3                     #pi^3
c = 299_792_458                 #Speed of light [m/s] - exact
c2 = c**2                       #Speed of light squared c^2 [m^2 / s^2]
c3 = c**3                       #Speed of light cubed c^3 [m^3 / s^3]
c4 = c**4                       #Speed of light to the fourth power c^4 [m^4 / s^4]
me = 9.109_383_701_5e-31        #Mass of electron [kg]
mp = 1.672_621_923_69e-27       #Mass of proton [kg]
mn = 1.674_927_498_04e-27       #Mass of neutron [kg]
mHe = 6.644_657_335_7e-27       #Mass of alpha particle ^4He^2+ [kg]
e = 1.602_176_634e-19           #Elementary charge [C] - exact
Na = 6.022_140_76e23            #Avogadro constant [1 / mol]
h = 6.626_070_15e-34            #Planck's constant [m^2 kg/s] - exact
h2 = h**2                       #Planck's constant squared [m^2 kg/s]^2
h3 = h**3                       #Planck's constant cubed [m^2 kg/s]^3
h4 = h**4                       #Planck's constant to the fourth power [m^2 kg/s]^4
hbar = h / (2*pi)               #Reduced Planck constant [m^2 kg/s]
hbar2 = hbar**2                 #Reduced Planck constant squared [m^2 kg/s]^2
hbar3 = hbar**3                 #Reduced Planck constant cubed [m^2 kg/s]^3
hbar4 = hbar**4                 #Reduced Planck constant to the fourth power [m^2 kg/s]^4
kB = 1.380_649e-23              #Boltzmann constant [J / K] - exact
u0 = 4 * pi * 10**-7            #Vacuum permeability AKA Magnetic constant mu_not [H/m]
e0 = 1 / (u0 * c2)              #Vacuum permittivity epsilon_not [F / m]
k = 1 / (4 * pi * e0)           #Coulomb constant [N * m^2 / C^2)]
G = 6.674_30e-11                #Gravitational constant [m^3 / kg s^2]
g = 9.806_65                    #Acceleration due to gravity [m / s^2]
uB = e * hbar / (2 * me)        #Bohr magneton mu_B [J / T]
uN = e * hbar / (2 * mp)        #Nuclear magneton mu_N [J/T]
H0 = 2.333_361e-18              #Hubble constant approx = 72 [km / s * Mpc]  = _ [1 / s] - inexact
msun = 1.988_47e30              #Mass of sun [kg]
F = Na * e                      #Faraday constant [C / mol]
Z0 = u0 * c                     #Characteristic impedance of vacuum or Impedance of free space [Ohms]
R = Na * kB                     #Molar gas constant AKA Universal gas constant [J / K mol]
ln2 = math.log(2)               #Natural logarithm of 2 [unit-less]
phi = (1+sqrt(5))/2             #Golden ratio [unit-less]
alpha = e**2 / (2 * e0 * h * c)                     #Fine-structure constant alpha [unit-less]
a0 = 4 * pi * e0 * hbar**2 / (me * e**2)            #Bohr radius a_not [m]
Rinf = alpha**2 * me * c / (2 * h)                  #Rydberg constant R_infinity [1 / m]
sigma = 2 * pi**5 * kB**4 / (15 * h**3 * c2)        #Stefan-Boltzmann constant [W / m^2 K^4]

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

# deg = pi / 180                #Degrees to radians factor [rad / degree)       #DEPRECIATED in favor of degrees function
MeV = e * 10**6                 #MeV conversion factor [J / MeV)
amu = 1.660_539_066_6e-27       #Atomic mass unit AMU conversion factor [kg / AMU) 
mol = Na                        #Mol to molecules conversion factor [molecules / mol)
min = 60                        #Time conversion [s / min)
hr = 60 * min                   #Time conversion [s / hr)
day = 24 * hr                   #Time conversion [s / day)
week = 7 * day                  #Time conversion [s / week)
month = 365 / 12 * day          #Time conversion [s / month)
yr = 365 * day                  #Time conversion [s / yr)
IN = 0.0254                     #Imperial length conversion [m / in)      #'in' is a reserved word in python so using IN instead
in2 = IN**2                     #Imperial area conversion [m^2 / in^2)
in3 = IN**3                     #Imperial volume conversion [m^3 / in^3)
in4 = IN**4                     #Imperial hyper volume conversion [m^4 / in^4)
ft = 0.3048                     #Imperial length conversion [m / ft)
ft2 = ft**2                     #Imperial area conversion [m^2 / ft^2)
ft3 = ft**3                     #Imperial volume conversion [m^3 / ft^3)
ft4 = ft**4                     #Imperial hyper volume conversion [m^4 / ft^4)
yd = 0.9144                     #Imperial length conversion [m / yd)
yd2 = yd**2                     #Imperial area conversion [m^2 / yd^2)
yd3 = yd**3                     #Imperial volume conversion [m^3 / yd^3)
yd4 = yd**4                     #Imperial hyper volume conversion [m^4 / yd^4)
mi = 1_609.344                  #Imperial length conversion [m / mi)
mi2 = mi**2                     #Imperial area conversion [m^2 / mi^2)
mi3 = mi**3                     #Imperial volume conversion [m^3 / mi^3)
mi4 = mi**4                     #Imperial hyper volume conversion [m^4 / mi^4)
kmh = 1e3 / hr                  #Km / hour -> m/s   [[m/s] / kmh)
mph = mi / hr                   #Miles per hour to m/s [[m/s] / mph)
au = 149_597_870_700            #Astronomical unit to meter [m / au)
pc =  3.085_677_581_28e16       #Parsec to meter   [m / Pc)
Mpc = 1e6 * pc                  #Mega parsecs to meter [m / MPc)
ly = c * 365.25 * day           #Lightyear to meter [m / ly)
mach = 340.5                    #one Mach (approx., at 15 C, 1 atm) in meters per second [[m/s] / mach)
angstrom = 1e-10                #Angstrom to meters [m / A)
lbf = 4.448_221_6               #Pound force to Newtons conversion [N / lbf)
kip = kilo * lbf                #KiloPound force to Newtons conversion [N / kip)
liter = 0.001                   #liter to m^3 conversion [m^3 / L)
mL = milli * liter              #milliliter to m^3 conversion [m^3 / mL)
bar = 100_000                   #bar to pascal conversion [pa / bar) = [[N / m^2] / bar)
atm = 101_325                   #standard atmosphere to pascal conversion [[N / m^2] / atm)
psi = lbf / in2                 #Pounds per square inch to pascal [pa / psi) = [[N / m^2] / [lbf / in^2])
ksi = kilo * psi                #Kips per square inch to pascal [pa / ksi) = [[N / m^2] / [kilolbf / in^2])
btu = 1_055.055_852_6           #British thermal unit to Joules conversion [J / btu)
therm = btu * 100_000           #therm to Joules conversion [J / therm)
lbm = 0.453_592_37              #Pound to kilogram conversion [kg / lbm)
ton = 2000 * lbf                #Ton to newton conversion [N/ ton)
hp = 745.75                     #Horsepower to Watt conversion [W / hp)
acre = 4046.856                 #acre to m^2 conversion [m^2 / acre)
ha = hecto * acre               #hectare to m^2 conversion [m^2 / ha)
floz = 28.413 * mL              #fluid ounce to m^3 conversion [m^3 / fl oz)
pint = 568.261 * mL             #pint to m^3 conversion [m^3 / pint)
quart = 1236.523 * mL           #quart to m^3 conversion [m^3 / quart) 
gallon = 4546.09 * mL           #gallon to m^3 conversion [m^3 / gal)
oz = 0.028_350                  #ounce to kg conversion [kg / oz)
rankine = 5/9                   #Rankine t conversion [K / R = C / F)
arcmin = radians(1/60)          #Arc Minute to Radians conversion [rad / arcmin)
arcsec = arcmin / 60            #Arc Second to Radians conversion [rad / arcsec)
mhr = milli * hr                #milli hours to seconds [s / mhr )
khr = kilo * hr                 #kilo hours to seconds [s / khr )
Mhr = mega * hr                 #Mega hours to seconds [s / Mhr )
cm2 = centi**2                  #area conversion [m^2 / cm^2)
cm3 = centi**3                  #volume conversion [m^3 / cm^3)
cm4 = centi**4                  #hyper volume conversion [m^4 / cm^4)
mm2 = milli**2                  #area conversion [m^2 / mm^2)
mm3 = milli**3                  #volume conversion [m^3 / mm^3)
mm4 = milli**4                  #hyper volume conversion [m^4 / mm^4)
celsius = 273.15                #Kelvin-celsius offset [[Kelvin] + celsius = [celsius] ; [celsius] - celsius = [Kelvin] )   #Semi-depreciated, use temp conversion functions instead
calorie = 4.184                 #calorie to Joule conversion [J / calorie)
kcal = 4184                     #Calorie to Joule conversion [J / kcal)
rps = 2 * pi                    #Rotations per second to radians per second [[rad]/s / rot/s)
rpm = 2 * pi / 60               #Rotations per minute to radians per second [[rad]/s / rot/min)
lbft = lbf * ft                 #Pound force * feet to Newtons * meter [[N * m] / [lbf * ft])
lbin = lbf * IN                 #Pound force * inch to Newtons * meter [[N * m] / [lbf * in])
kipft = kip * ft                #Kilo pound force * feet to Newtons * meter [[N * m] / [kip * ft])
kipin = kip * IN                #Kilo pound force * inch to Newtons * meter [[N * m] / [kip * in])
barn = 1e-28                    #Barn to m^2 conversion [m^2 / b)
fb = femto * barn               #Femto Barn to m^2 conversion [m^2 / fb)
Torr = atm / 760                #Torr to pascal conversion [Pa / torr)
mmHg = 133.322_387_415          #Millimeters of mercury to pascal conversion [Pa / mmHg)  - exact
kgf = 9.806_65                  #kilogram force to newtons conversion [N / kgf)
at = 98066.5                    #Technical atmosphere [= kgf / cm^2) to pascal conversion [Pa / at)
dyne = 1e-5                     #Dyne to Newton conversion [N / dyn)
statcoulomb = 3.335_64e-10      #Statcoulomb to coulomb conversion [C / statC)
gft = 9.806_65 / ft             #Acceleration due to gravity [ft / s^2)     - Non SI
gin = 9.806_65 / IN             #Acceleration due to gravity [in / s^2)     - Non SI
slug = 14.593902937             #Slug to kilogram conversion [kg / slug)
story = 3.3                     #Story to meter conversion [m / story)

#Material constants:    -------------------------------------------------------------------------------------------
Rair = 287.05                   #Individual Gas Constant of Air [J / K kg]
Rwater = 461.52                 #Individual Gas Constant of Water Vapor [J / K kg]
Roxygen = 259.84                #Individual Gas Constant of Oxygen O2 [J / K kg]
Rnitrogen = 296.80              #Individual Gas Constant of Nitrogen N2 [J / K kg]
RCO2 = 118.92                   #Individual Gas Constant of Carbon Dioxide CO2 [J / K kg] 

pair = 1.204                    #Density of air at 20C 1atm [kg / m^3]
pwater = 998.19                 #Density of water at 20C 1atm [kg / m^3]
poxygen = 1.314                 #Density of oxygen at 20C 1atm [kg / m^3]
pnitrogen = 1.16                #Density of nitrogen at 20C 1atm [kg / m^3]
pCO2 = 1.815                    #Density of CO2 at 20C 1atm [kg / m^3]

EAluminum = 69 * giga           #Elastic modulus of Aluminum [Pa] 
EAl2014 = 73.1 * giga           #Elastic modulus of Aluminum 2014-T6 [Pa] 
EAl6061 = 68.9 * giga           #Elastic modulus of Aluminum 6061-T6 [Pa] 
EGrayCastIron = 67.0 * giga     #Elastic modulus of cast iron alloy Gray ASTM 20 [Pa] 
EMalleableCastIron = 172 * giga #Elastic modulus of cast iron alloy Malleable ASTM A-197 [Pa] 
ECopper = 117 * giga            #Elastic modulus of Copper [Pa] 
EBrass = 101 * giga             #Elastic modulus of Red Brass C83400 [Pa] 
EBronze = 103 * giga            #Elastic modulus of Bronze C86100 [Pa] 
EIron = 210 * giga              #Elastic modulus of Iron [Pa] 
ESteel = 200 * giga             #Elastic modulus of many Steel alloys (Structural A-36, Structural A992, Tool L2) [Pa] 
EStainlessSteel = 193 * giga    #Elastic modulus of steel alloy Stainless 304  [Pa] 
ENickel = 170 * giga            #Elastic modulus of Nickel [Pa]
ETitaniumAlloy = 120 * giga     #Elastic modulus of titanium alloy [Ti-6Al-4V] [Pa] 
ELowStrengthConcrete = 22.1 * giga      #Elastic modulus of Low Strength Concrete [Pa] 
EHighStrengthConcrete = 29.0 * giga     #Elastic modulus of High Strength Concrete [Pa] 
EDouglasFir = 13.1 * giga       #Elastic modulus of wood structural grade Douglas Fir [Pa] 
EWhiteSpruce = 9.65 * giga      #Elastic modulus of wood structural grade White Spruce [Pa] 

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
liters, l, L, Liter, lit, Lit = liter, liter, liter, liter, liter, liter
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
torr = Torr
mmhg = mmHg
dyn, Dyne = dyne, dyne
statC, statc, Fr, fr, franklin, esu = statcoulomb, statcoulomb, statcoulomb, statcoulomb, statcoulomb, statcoulomb
gimp, gimperial, gcust, gcus, gcustomary = gft, gft, gft, gft, gft
ginn, ginch, ginches, gIN, gIn = gin, gin, gin, gin, gin
stories, floors = story, story

RAir = Rair
RWater, Rh2o, RH2o, RH2O, RH20, Rh20 = Rwater, Rwater, Rwater, Rwater, Rwater, Rwater
RO2, Ro2, ROxygen = Roxygen, Roxygen, Roxygen
RCarbonDioxide, RCarbon_Dioxide, RcarbonDioxide, Rcarbon_Dioxide, Rcarbondioxide, Rcarbon_dioxide, Rco2, RCo2, RC02, rc02 = RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2

PAir, rhoair, rhoAir = pair, pair, pair
pWater, ph2o, pH2o, pH2O, pH20, ph20, rhoWater, rhowater, rhoH2O, rhoH2o, rhoh2o = pwater, pwater, pwater, pwater, pwater, pwater, pwater, pwater, pwater, pwater, pwater
pO2, po2, pOxygen, rhoO2, rhooxygen, rhoOxygen = poxygen, poxygen, poxygen, poxygen, poxygen, poxygen
pCarbonDioxide, pCarbon_Dioxide, pcarbonDioxide, pcarbon_Dioxide, pcarbondioxide, pcarbon_dioxide, pco2, pCo2, pC02, rhoCO2, rhoCarbonDioxide, rhocarbondioxide, rhoCarbondioxide = pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2, pCO2


EAl = EAluminum
EAluminum2014 = EAl2014 
EAluminum6061 = EAl6061
ECastIronGray, EIronGray, EGrayIron = EGrayCastIron, EGrayCastIron, EGrayCastIron
ECastIronMalleable, EIronMalleable, EMalleableIron = EMalleableCastIron, EMalleableCastIron, EMalleableCastIron
ECu = ECopper
EBrass
EBronze
EFe = EIron
ESt = ESteel
ESS, ESs, Ess = EStainlessSteel, EStainlessSteel, EStainlessSteel
ENi = ENickel
ETiAlloy, ETi6Al4V, ETiAlV = ETitaniumAlloy, ETitaniumAlloy, ETitaniumAlloy
ELSConcrete = ELowStrengthConcrete
EHSConcrete = EHighStrengthConcrete
EDouglasFir
EWhiteSpruce

#Intra Program Variables
ans = [0]
ans1 = 0
ans2 = 0
degreesFlag = False                         #Default in Radians Mode
floatDeltaPercent = nano                           #Accepatable percent error for floating point number comparisions
floatDeltaAbs = femto

# def floatComparision(float1, float2):
#     """returns true if percent error between 2 floats is less than acceptable percent error (floatDelta)

#     Args:
#         float1 (_type_): _description_
#         float2 (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     global floatDeltaAbs, floatDeltaPercent
#     absDiff = abs(float1 - float2) 
#     if(absDiff < floatDeltaPercent * float1 and absDiff < floatDeltaPercent * float2):
#         return True
#     #edge case: float1 ~= float2 ~= 0 -> return true if absolute val of both less than floatDeltaAbs
#     return (abs(float1) < floatDeltaAbs and abs(float2) < floatDeltaAbs)

def floatComparision(float1, float2, rel_tol=floatDeltaPercent):
    global floatDeltaAbs
    return math.isclose(float1, float2, rel_tol=rel_tol) or (abs(float1) < floatDeltaAbs and abs(float2) < floatDeltaAbs) 
    
    

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
    if(isinstance(x, complex)):
        return cmath.sqrt(x)
    return math.sqrt(x)
    
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
    if(isinstance(x, complex)):
        return cmath.sin(radians(x))
    return math.sin(radians(x))

def cosd(x):
    if(isinstance(x, complex)):
        return cmath.cos(radians(x))
    return math.cos(radians(x))

def tand(x):
    if(isinstance(x, complex)):
        return cmath.tan(radians(x))
    return math.tan(radians(x))

#Radian based Trig Funcitons
def sinr(x):
    if(isinstance(x, complex)):
        return cmath.sin(x)
    return math.sin(x)

def cosr(x):
    if(isinstance(x, complex)):
        return cmath.cos(x)
    return math.cos(x)

def tanr(x):
    if(isinstance(x, complex)):
        return cmath.tan(x)
    return math.tan(x)

#Trig Functions that check whether in degree mode or radian mode, and give warnings 
def sin(x):
    #if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            print("Warning: Complex number detected but in degree mode, sinr(z) calculated anyway. Use sind(z) if that is what you intended")
        return sinr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (sind(x) calculated anyway)")
        return sind(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180 
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
        print("Warning: In Radian Mode but possible degree number detected. (sinr(x) calculated anyway)")
    return sinr(x)

def cos(x):
    #if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            print("Warning: Complex number detected but in degree mode, cosr(z) calculated anyway. Use cosd(z) if that is what you intended")
        return cosr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (cosd(x) calculated anyway)")
        return cosd(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
        print("Warning: In Radian Mode but possible degree number detected. (cosr(x) calculated anyway)")
    return cosr(x)

def tan(x):
    #if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            print("Warning: Complex number detected but in degree mode, tanr(z) calculated anyway. Use tand(z) if that is what you intended")
        return tanr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            print("Warning: In Degree Mode but possible radian number detected. (tand(x) calculated anyway)")
        return tand(x)
    
    #Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
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
    if(isinstance(x, complex)):
        return cmath.log(x)
    return math.log(x)
    
def exp(x):
    if(isinstance(x, complex)):
        return cmath.exp(x)
    return math.exp(x)

def lg(x):
    if(isinstance(x, complex)):
        return cmath.log2(x)
    return math.log(x) / cmath.log(2)

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

#Volume of sphere given diameter
dvolume = lambda d : 4/3 * pi * (d/2)**3

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
    """Returns std dev of sample"""
    if(population):
        return np.std(arr,ddof=1)
    return np.std(arr)

def popstd(arr):
    """Returns population std dev"""
    return np.std(arr)

#Returns sample variance by default
def var(arr, population=False):
    if(population):
        return np.var(arr,ddof=1)
    return np.var(arr)

def popvar(arr):
    return np.var(arr)

#n choose r . Number of combinations (order doesn't matter) of r items chosen from n options
def nCr(n, r):
    return factorial(n) / (factorial(n-r) * factorial(r))

#n permute r . Number of permutations (order matters) of r items from list of n options
def nPr(n, r):
    return factorial(n) / factorial(n-r)

#Given an arbritrary number of days, returns the number of people required for Probability[two people share the same birthday] >= 0.5   
def birthdayProblem(days=365):
    return ceil( (3 - 2 * ln(2))/6 + sqrt(2*days*ln(2)) + (9 - 4*ln(2)**2 ) / (72 * sqrt(2*days*ln(2))) - 2 * ln(2)**2 / (135*days) )

# def singularity(magnitude, start, degree):
#     def fun(x):
#         if(x < start or degree < 0):
#             return 0
#         return magnitude * (x-start)**degree
#     return fun

def singularity3(magnitude, input, degree):
    """Singularity helper function : magnitude * <x - start>**degree

    Args:
        magnitude (float): 
        input (float): x - start
        degree (int): 

    Returns:
        _type_: float
    """
    if (input < 0 or degree < 0):
        return 0
    return magnitude * input**degree

def singularity(*args):
    """Singularity function : magnitude * <x - start>**degree

    Args:
        args[0]     = magnitude (float): optional argument
        args[0 / 1] = input (float): x - start
        args[1 / 2] = degree (int): 

    Returns:
        _type_: float
        
    Raises:
        Exception: Invalid number of parameters
    """
    if(len(args) == 2):
        input = args[0]
        degree = args[1]
        return singularity3(1, input, degree)
    elif(len(args) == 3):
        magnitude = args[0]
        input = args[1]
        degree = args[2]
        return singularity3(magnitude, input, degree)
    else:
        raise Exception('Invalid number of parameters')
    
def printfunction(f, start, end, stepSize, outputResolution = 4):
    """Prints values of f with range [start : step : end]

    Args:
        f (function): _description_
        start (float): start of range
        end (float): end of range
        stepSize (float): distance between successive x values
        outputResolution (int, optional): Number of decimal places to print output. Defaults to 4.
    """
    x = start
    inputResolution = max(0, ceil(-log10(stepSize)))     #resoultion = number of decimal places needed to format x, minimum of 0
    inputSize = (inputResolution > 0) + 1 + floor(max(log10(abs(start) + (start == 0))+ (start < 0), log10(abs(end) + (end == 0)) + (end < 0)))
    format = f"{inputSize+inputResolution}.{inputResolution}f"
    while (x < end or floatComparision(x, end)):
        print(f"{x :{format}}  :  {f(x) :9.{outputResolution}f} ")
        x += stepSize

def eig(mat):
    """np.linalg.eig redefinition - prints eigenvalues and eigenvectors in better format

    Args:
        mat (np.matrix): 
    """
    eigs = np.linalg.eig(mat)
    eigval = eigs[0]
    for i in range(len(eigval)):
        print(f"λ{i+1} = {eigval[i]}")
    print()
    print("Eigenvector Matrix (ans):")
    global ans
    print(ans := eigs[1])
    
def represent(x, deltaPercent = 0.000_1):
    exactPercent = 1e-18
    candidateList = []  #List of pair tuples with near but not exact representations, and the absolute distance to x
    constList = {1 : "1", pi: "pi", sqrt(2) : "sqrt(2)", sqrt(3) : "sqrt(3)"}
    for const in constList.keys():
        val = x / const
        checkedRatios = {0}     # Used to remove non reduced fractions, ie 2/4 when 1/2 has already been checked
        # for loop searches for fraction with denominator [1,200]
        for i in range(1, 401):
            if(round(val * i) / i in checkedRatios):
                continue
            # Checks if i is a viable denominator
            if( isclose(x, const * round(val * i) / i, rel_tol=deltaPercent ) ):
                # if exact match, print it and stop searching
                if( isclose(x, const * round(val * i) / i, exactPercent) ):
                    print(f"Found match: \n{constList[const]} * { round(val * i) } / {i}")
                    return
                else:
                    checkedRatios.add(round(val * i) / i)
                    candidateList.append((f"{constList[const]} * { round(val * i) } / {i}" , const * round(val * i) / i , abs(x - const * round(val * i) / i )))
    
    print("Did not find exact representation, here are closest matches:")
    candidateList.sort(key = lambda pair: pair[2])    #sort candidate list by smallest to largest distance
    print(f"{'True input':<30} : {f'{x:.8f}':^16}") 
    for i in candidateList:
        print(f"{i[0]:<30} : {f'{i[1]:.8f}':^16} : {f'{i[2]:12.9}':>12} ")

#Function and Lambda Aliases
root, roots = sp.optimize.fsolve, sp.optimize.fsolve
isclose, isClose, floatcomparsion, floatcomp = floatComparision, floatComparision, floatComparision, floatComparision
fact = factorial
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
InnerProduct, innerproduct, Innerproduct, inprod, inproduct, braket, BraKet, Braket, dotproduct, dotProduct, dotprod, dotProd, dot = innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod
magnitude, mag = norm, norm
integral, integrate = int, int
stdDev, stddev, StdDev, standardDeviation = std, std, std, std
popstdDev, popstddev, popStdDev, populationStddev, populationStdDev, populationStd, populationstd, populationstandardDeviation, populationStandardDeviation = popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd
variance = var
populationVariance, popvariance, popVariance = popvar, popvar, popvar
birthdayproblem, birthdayparadox, birthdayParadox, generalizedBirthdayProblem = birthdayProblem, birthdayProblem, birthdayProblem, birthdayProblem
sing = singularity
printFunction, printfun, printFun, plotFun, plotfun = printfunction, printfunction, printfunction, printfunction, printfunction
