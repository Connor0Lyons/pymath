import time
from time import sleep

import builtins

import cmath

import fractions
from fractions import Fraction

import math
from math import isclose
from math import factorial

import random
from random import random as rand

# import matplotlib.pyplot as plt

# import pint
# unit = pint.UnitRegistry()
# u = unit

import threading
from threading import Timer

import decimal

import sys

import os

import re

try:
    import CoolProp.CoolProp as CP
except(ModuleNotFoundError):
    print("Warning: Error importing CoolProp. Some features may not work as intended")

try:
    import scipy as sp
    # from scipy.constants import *
    # from scipy.optimize import fsolve
    import scipy.special
    import scipy.optimize as opt
    from scipy.optimize import fsolve
    from scipy.special import erf
    from scipy.special import erfc
    from scipy.special import j0
    from scipy.special import j1
    from scipy.special import jv
    from scipy.integrate import dblquad as iint
    from scipy.special import ndtri as invnorm
except(ModuleNotFoundError):
    print("Warning: Error importing SciPy. Some features may not work as intended")
    fsolve = lambda *args, **kwargs: print("SciPy not imported properly, fsolve function not found")
    j0 = lambda *args, **kwargs: print("SciPy not imported properly, j0 function not found")
    j1 = lambda *args, **kwargs: print("SciPy not imported properly, j1 function not found")
    jv = lambda *args, **kwargs: print("SciPy not imported properly, jv function not found")
    
try:
    import numpy as np
    from numpy import matrix
    from numpy import linalg
    from numpy.linalg import *
    from numpy import outer
    from numpy import cross
    from numpy import kron
    from numpy import sort
    from numpy import sort_complex
    from numpy import std as popstd
    # TODO Fix popvar
    from numpy import var as popvar
    from numpy import median
    from numpy import identity
except(ModuleNotFoundError):
    print("Warning: Error importing NumPy. Some features may not work as intended")
    outer = lambda *args, **kwargs: print("NumPy not imported properly, outer function not found")
    cross = lambda *args, **kwargs: print("NumPy not imported properly, cross function not found")
    sort = lambda *args, **kwargs: print("NumPy not imported properly, sort function not found. Maybe try list(arr).sort() instead.")
    popstd = lambda *args, **kwargs: print("NumPy not imported properly, popstd function not found.") 
    popvar = lambda *args, **kwargs: print("NumPy not imported properly, popvar function not found.") 
    def matrix(*args, **kwargs): 
        print("NumPy not imported properly, matrix constructor not found")
        return nan
    def kron(*args, **kwargs): 
        print("NumPy not imported properly, kron function not found")
        return nan
    def identity(*args, **kwargs): 
        print("NumPy not imported properly, identity function not found")
        return nan

try:
    import pandas as pd
except(ModuleNotFoundError):
    print("Warning: Error importing pandas. Some features may not work as intended")

from cmath import *

from math import *

import datetime as dt
from datetime import datetime


# Unless noted otherwise, units are in terms of official SI units: second [s], meter [m], kilogram [kg], ampere [A], kelvin [K], mole [mol], and candela [cd]

## Universal constants:   ------------------------------------------------------------------------------------------
pi = math.pi                    # Circle constant pi [unit-less]
pi2 = pi**2                     # pi^2
pi3 = pi**3                     # pi^3
ipi = pi * 1j                   # pi * sqrt(-1) 
c = 299_792_458                 # Speed of light [m/s] - exact
c2 = c**2                       # Speed of light squared c^2 [m^2 / s^2]
c3 = c**3                       # Speed of light cubed c^3 [m^3 / s^3]
c4 = c**4                       # Speed of light to the fourth power c^4 [m^4 / s^4]
me = 9.109_383_701_5e-31        # Mass of electron [kg]
mp = 1.672_621_923_69e-27       # Mass of proton [kg]
mn = 1.674_927_498_04e-27       # Mass of neutron [kg]
mHe = 6.644_657_335_7e-27       # Mass of alpha particle ^4He^2+ [kg]
e = 1.602_176_634e-19           # Elementary charge [C] - exact
Na = 6.022_140_76e23            # Avogadro constant [1 / mol] - exact
h = 6.626_070_15e-34            # Planck's constant [m^2 kg/s] - exact
h2 = h**2                       # Planck's constant squared [m^2 kg/s]^2
h3 = h**3                       # Planck's constant cubed [m^2 kg/s]^3
h4 = h**4                       # Planck's constant to the fourth power [m^2 kg/s]^4
hbar = h / (2*pi)               # Reduced Planck constant [m^2 kg/s]
hbar2 = hbar**2                 # Reduced Planck constant squared [m^2 kg/s]^2
hbar3 = hbar**3                 # Reduced Planck constant cubed [m^2 kg/s]^3
hbar4 = hbar**4                 # Reduced Planck constant to the fourth power [m^2 kg/s]^4
kB = 1.380_649e-23              # Boltzmann constant [J / K] - exact
u0 = 4 * pi * 10**-7            # Vacuum permeability AKA Magnetic constant mu_not [H/m]
e0 = 1 / (u0 * c2)              # Vacuum permittivity epsilon_not [F / m] = [C^2 / N m^2] = [A^2 s^4 / kg m^3]
k = 1 / (4 * pi * e0)           # Coulomb constant [N m^2 / C^2)]
G = 6.674_30e-11                # Gravitational constant [m^3 / kg s^2]
g = 9.806_65                    # Acceleration due to gravity [m / s^2]
uB = e * hbar / (2 * me)        # Bohr magneton mu_B [J / T]
uN = e * hbar / (2 * mp)        # Nuclear magneton mu_N [J/T]
H0 = 2.333_361e-18              # Hubble constant approx = 72 [km / s Mpc]  = _ [1 / s] - inexact
mSun = 1.988_47e30              # Mass of sun [kg]
rSun = 6.95700e8                # Radius of sun [m]
mEarth = 5.972_19e30            # Mass of Earth [kg]
rEarth =  6.378_1e6             # Radius of Earth [m]
F = Na * e                      # Faraday constant [C / mol]
Z0 = u0 * c                     # Characteristic impedance of vacuum or Impedance of free space [Ohms]
R = Na * kB                     # Molar gas constant AKA Universal gas constant [J / K mol]
ln2 = math.log(2)               # Natural logarithm of 2 [unit-less]
sqrt2 = math.sqrt(2)            # Square root of 2 [unit-less]
sqrt3 = math.sqrt(3)            # Square root of 3 [unit-less]
phi = (1 + math.sqrt(5)) / 2    # Golden ratio [unit-less]
alpha = e**2 / (2 * e0 * h * c)                     # Fine-structure constant alpha [unit-less]
a0 = 4 * pi * e0 * hbar**2 / (me * e**2)            # Bohr radius a_not [m]
Roo = alpha**2 * me * c / (2 * h)                   # Rydberg constant for large atoms R_infinity [1 / m]
RH = Roo * mp / (me + mp)                           # Rydberg constant for Hydrogen [1 / m]         # Warning: potential conflict with individual molar gas constant for hydrogen (not implemented)
Ry = h * c * Roo                                    # Rydberg Unit of energy [J]
sigma = 2 * pi**5 * kB**4 / (15 * h**3 * c2)        # Stefan-Boltzmann constant [W / m^2 K^4]

true = True
false = False


## Matrices / Quantum Gates
if "NumPy" in sys.modules:
    sigmax = matrix("0 1; 1 0")             # Pauli spin matrix x [unit-less]
    sigmay = matrix("0 -1j; 1j 0")          # Pauli spin matrix y [unit-less]
    sigmaz = matrix("1 0; 0 -1")            # Pauli spin matrix z [unit-less]
    H = matrix("1. 1.; 1. -1.") / sqrt(2)   # First Order Hadamard Matrix
    H2 = kron(H, H)                         # Second Order Hadamard Matrix
    H3 = kron(H, H2)                        # Third Order Hadamard Matrix
    CX = matrix([[1., 0., 0., 0.],          # Controlled X Gate
                [0., 1., 0., 0.],
                [0., 0., 0., 1.],
                [0., 0., 1., 0.]])
    CY = matrix([[1., 0.,  0., 0.],         # Controlled Y Gate
                [0., 1.,  0., 0.],
                [0., 0.,  0, -1j],
                [0., 0., +1j, 0.]])
    CZ = matrix([[1., 0., 0., 0.],          # Controlled Z Gate
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0, -1.]])
    zp = matrix("1. ; 0.")                  # Pauli +Z Basis Representation of Spin 1/2 +Z ket 
    zn = matrix("0. ; 1.")                  # Pauli +Z Basis Representation of Spin 1/2 -Z ket
    xp = matrix("1. ; 1.") / sqrt(2)        # Pauli +Z Basis Representation of Spin 1/2 +X ket
    xn = matrix("1. ; -1.") / sqrt(2)       # Pauli +Z Basis Representation of Spin 1/2 -X ket
    yp = matrix("1. ; 1j") / sqrt(2)        # Pauli +Z Basis Representation of Spin 1/2 +Y ket
    yn = matrix("1. ; -1j") / sqrt(2)       # Pauli +Z Basis Representation of Spin 1/2 -Y ket
    I2 = identity(2)                        # 2x2 Identity Matrix  - Note:Type = numpy.ndarray, which is compatible with but distinct from numpy.matrix
    I3 = identity(3)                        # 3x3 Identity Matrix  - Note:Type = numpy.ndarray, which is compatible with but distinct from numpy.matrix
    I4 = identity(4)                        # 4x4 Identity Matrix  - Note:Type = numpy.ndarray, which is compatible with but distinct from numpy.matrix
    SWAP = matrix([[1., 0., 0., 0.],        # Matrix representation of SWAP gate
                [0., 0., 1., 0.],
                [0., 1., 0., 0.],
                [0., 0., 0., 1.]])
    CSWAP = identity(8); CSWAP[5,5] = 0; CSWAP[5,6] = 1; CSWAP[6,6] = 0; CSWAP[6,5] = 1      # Matrix representation of Controlled SWAP gate  - Note:Type = numpy.ndarray, which is compatible with but distinct from numpy.matrix
    CCNOT = identity(8); CCNOT[6,6] = 0; CCNOT[6,7] = 1; CCNOT[7,7] = 0; CCNOT[7,6] = 1      # Matrix representation of Controlled Controlled NOT gate (AKA Toffoli gate) - Note:Type = numpy.ndarray, which is compatible with but distinct from numpy.matrix

## Conversion Factors:    ------------------------------------------------------------------------------------------
zepto = 1e-21                   # SI Small Prefixes
atto = 1e-18
femto = 1e-15
pico = 1e-12
nano = 1e-9
micro = 1e-6
milli = 1e-3
centi = 1e-2
deci = 1e-1                     #--------------------
deka = 1e1                      # SI Big Prefixes
hecto = 1e2                     
kilo = 1e3
mega = 1e6
giga = 1e9
tera = 1e12
peta = 1e15
exa = 1e18
zetta = 1e21                    #--------------------

kibi = 2**10                    # Binary prefixes
mebi = 2**20
gibi = 2**30
tebi = 2**40
pebi = 2**50
exbi = 2**60
zebi = 2**70
yobi = 2**80                    #--------------------

# deg = pi / 180                # Degrees to radians factor [rad / degree]       # DEPRECIATED in favor of degrees function
MeV = e * 10**6                 # MeV conversion factor [J / MeV]
amu = 1.660_539_066_6e-27       # Atomic mass unit AMU conversion factor [kg / AMU] 
mol = Na                        # Mol to molecules conversion factor [molecules / mol]
minute = 60                     # Time conversion [s / min]                    
hr = 60 * minute                # Time conversion [s / hr]
day = 24 * hr                   # Time conversion [s / day]
week = 7 * day                  # Time conversion [s / week]
month = 365 / 12 * day          # Time conversion [s / month]
yr = 365 * day                  # Time conversion [s / yr]
IN = 0.0254                     # Imperial length conversion [m / in]      # 'in' is a reserved word in python so using IN instead
in2 = IN**2                     # Imperial area conversion [m^2 / in^2]
in3 = IN**3                     # Imperial volume conversion [m^3 / in^3]
in4 = IN**4                     # Imperial hyper volume conversion [m^4 / in^4]
ft = 0.3048                     # Imperial length conversion [m / ft]
ft2 = ft**2                     # Imperial area conversion [m^2 / ft^2]
ft3 = ft**3                     # Imperial volume conversion [m^3 / ft^3]
ft4 = ft**4                     # Imperial hyper volume conversion [m^4 / ft^4]
yd = 0.9144                     # Imperial length conversion [m / yd]
yd2 = yd**2                     # Imperial area conversion [m^2 / yd^2]
yd3 = yd**3                     # Imperial volume conversion [m^3 / yd^3]
yd4 = yd**4                     # Imperial hyper volume conversion [m^4 / yd^4]
mi = 1_609.344                  # Imperial length conversion [m / mi]
mi2 = mi**2                     # Imperial area conversion [m^2 / mi^2]
mi3 = mi**3                     # Imperial volume conversion [m^3 / mi^3]
mi4 = mi**4                     # Imperial hyper volume conversion [m^4 / mi^4]
kmh = 1e3 / hr                  # Km / hour -> m/s   [[m/s] / kmh]
mph = mi / hr                   # Miles per hour to m/s [[m/s] / mph]
au = 149_597_870_700            # Astronomical unit to meter [m / au]
pc =  3.085_677_581_28e16       # Parsec to meter   [m / Pc]
Mpc = 1e6 * pc                  # Mega parsecs to meter [m / MPc]
ly = c * 365.25 * day           # Lightyear to meter [m / ly]
lightday = c * day              # Lightday to meter [m / ly]
lighthour = c * hr              # Lighthour to meter [m / ly]
mach = 340.5                    # one Mach (approx., at 15 C, 1 atm] in meters per second [[m/s] / mach]
angstrom = 1e-10                # Angstrom to meters [m / A]
lbf = 4.448_221_6               # Pound force to Newtons conversion [N / lbf]
kip = kilo * lbf                # KiloPound force to Newtons conversion [N / kip]
liter = 0.001                   # liter to m^3 conversion [m^3 / L]
mL = milli * liter              # milliliter to m^3 conversion [m^3 / mL]
bar = 100_000                   # bar to pascal conversion [pa / bar] = [[N / m^2] / bar]
atm = 101_325                   # standard atmosphere to pascal conversion [[N / m^2] / atm]
psi = lbf / in2                 # Pounds per square inch to pascal [pa / psi] = [[N / m^2] / [lbf / in^2]]
ksi = kilo * psi                # Kips per square inch to pascal [pa / ksi] = [[N / m^2] / [kilolbf / in^2]]
Msi = mega * psi                # Million pounds per square inch to pascal [pa / Msi] = [[N / m^2] / [kilolbf / in^2]]
btu = 1_055.055_852_6           # British thermal unit to Joules conversion [J / btu]
therm = btu * 100_000           # therm to Joules conversion [J / therm]
ton = 2000 * lbf                # Ton to newton conversion [N/ ton]
hp = 745.75                     # Horsepower to Watt conversion [W / hp]
acre = 4046.856                 # acre to m^2 conversion [m^2 / acre]
ha = hecto * acre               # hectare to m^2 conversion [m^2 / ha]
tsp = 4.92892159375 * mL        # US customary fluid teaspoon
tbsp = 14.78676478125 * mL      # US customary fluid tablespoon
floz = 29.5735295625 * mL       # US customary fluid ounce to m^3 conversion [m^3 / fl oz]
cup = 236.5882365 * mL          # US customary fluid cup to m^3 conversion (= 8 * floz) [m^3 / fl oz]
pint = 	473.176473 * mL         # US pint to m^3 conversion (= 16 * floz) [m^3 / pint]                  # Warning: potential conflict with point != pt
quart = 0.000946352946          # US quart to m^3 conversion (= 32 * floz) [m^3 / quart] 
gallon = 0.003785411784         # US gallon to m^3 conversion (= 128 * floz) [m^3 / gal]
lbm = 0.453_592_37              # Pound to kilogram conversion [kg / lbm]
oz = lbm / 16                   # International avoirdupois ounce to kg conversion [kg / oz]
rankine = 5/9                   # Rankine t conversion [K / R = C / F]
arcmin = pi / (180 * 60)        # Arc Minute to Radians conversion [rad / arcmin]
arcsec = arcmin / 60            # Arc Second to Radians conversion [rad / arcsec]
mhr = milli * hr                # milli hours to seconds [s / mhr ]
khr = kilo * hr                 # kilo hours to seconds [s / khr ]
Mhr = mega * hr                 # Mega hours to seconds [s / Mhr ]
cm2 = centi**2                  # area conversion [m^2 / cm^2]
cm3 = centi**3                  # volume conversion [m^3 / cm^3]
cm4 = centi**4                  # hyper volume conversion [m^4 / cm^4]
mm2 = milli**2                  # area conversion [m^2 / mm^2]
mm3 = milli**3                  # volume conversion [m^3 / mm^3]
mm4 = milli**4                  # hyper volume conversion [m^4 / mm^4]
celsius = 273.15                # Kelvin-celsius offset ([Kelvin] + celsius = [celsius] ; [celsius] - celsius = [Kelvin] )   # Semi-depreciated, use temp conversion functions instead
calorie = 4.184                 # calorie to Joule conversion [J / calorie]
kcal = 4184                     # Calorie to Joule conversion [J / kcal]
rps = 2 * pi                    # Rotations per second to radians per second [[rad]/s / rot/s]
rpm = 2 * pi / 60               # Rotations per minute to radians per second [[rad]/s / rot/min]
lbft = lbf * ft                 # Pound force * feet to Newtons * meter [[N * m] / [lbf * ft]]
lbin = lbf * IN                 # Pound force * inch to Newtons * meter [[N * m] / [lbf * in]]
kipft = kip * ft                # Kilo pound force * feet to Newtons * meter [[N * m] / [kip * ft]]
kipin = kip * IN                # Kilo pound force * inch to Newtons * meter [[N * m] / [kip * in]]
barn = 1e-28                    # Barn to m^2 conversion [m^2 / b]
fb = femto * barn               # Femto Barn to m^2 conversion [m^2 / fb]
Torr = atm / 760                # Torr to pascal conversion [Pa / torr]
mmHg = 133.322_387_415          # Millimeters of mercury to pascal conversion [Pa / mmHg]  - exact
inHg = 3386.39                  # Inches of mercury to pascal conversion [Pa / mmHg]  - exact
kgf = 9.806_65                  # kilogram force to newtons conversion [N / kgf]
at = 98066.5                    # Technical atmosphere [= kgf / cm^2] to pascal conversion [Pa / at]
dyne = 1e-5                     # Dyne to Newton conversion [N / dyn]
statcoulomb = 3.335_64e-10      # Statcoulomb to coulomb conversion [C / statC]
gft = 9.806_65 / ft             # Acceleration due to gravity [ft / s^2]     - Non SI
gin = 9.806_65 / IN             # Acceleration due to gravity [in / s^2]     - Non SI
slug = 14.593902937             # Slug to kilogram conversion [kg / slug]
story = 3.3                     # Story to meter conversion [m / story]
stoke = 1e-4                    # Stoke to pascal * second conversion [[Pa * s] / St ]
point =  127 / 36 * 1e-4        # Point to meter conversion [m / point]                         # Warning: potential conflict with pt = pint != point
mmH2O = 9.8066501               # Millimeters of water to pascal conversion (using rho_water = 1 kg / L)  [Pa / mmH2O]
inH2O = 249.08891               # Inches of water to pascal conversion (using rho_water = 1 kg / L)  [Pa / inH2O]
nauticalMile = 1852             # Nautical miles to meters conversion [m / NM]
knot = nauticalMile / hr        # Knot to meters per second conversion [m/s / knot]
thou = IN / 1000                # Thousandth of an inch to meters conversion [m / thou] == [m / mil]
ipm = IN / minute               # Inch per minute [ [m / s] / [in / min]]
sfm = ft / minute               # Surface feet per minute [ [m / s] / [ft / min]]

## Imperial Machine Screw Sizes [ANSI Units]
#! ANSI UNITS
screw0  = 0.0600
screw1  = 0.0730
screw2  = 0.0860
screw3  = 0.0990
screw4  = 0.1120
screw5  = 0.1250
screw6  = 0.1380
screw8  = 0.1640
screw10 = 0.1900
screw12 = 0.2160

## Material constants:    -------------------------------------------------------------------------------------------
RAir = 287.05                   # Individual Gas Constant of Air [J / K kg]
RWater = 461.52                 # Individual Gas Constant of Water Vapor [J / K kg]
ROxygen = 259.84                # Individual Gas Constant of Oxygen O2 [J / K kg]
RNitrogen = 296.80              # Individual Gas Constant of Nitrogen N2 [J / K kg]
RCO2 = 118.92                   # Individual Gas Constant of Carbon Dioxide CO2 [J / K kg] 
                                                                                                # Warning: Potential conflict with RH = Hydrogen Rydberg Constant != Individual gas constant of Hydrogen
rhoAir = 1.204                  # Density of air at 20C 1atm [kg / m^3]
rhoWater = 998.19               # Density of water at 20C 1atm [kg / m^3]
rhoOxygen = 1.314               # Density of oxygen at 20C 1atm [kg / m^3]
rhoNitrogen = 1.16              # Density of nitrogen at 20C 1atm [kg / m^3]
rhoCO2 = 1.815                  # Density of CO2 at 20C 1atm [kg / m^3]
rhoAl6061 = 2700                # Density of Aluminum 6061 [kg / m^3]
rhoAl3003 = 2730                # Density of Aluminum 3003 [kg / m^3]
rhoSteel = 7900                 # Density of Steel (general) [kg / m^3]         
rhoBP1 = 1860                   # Maximum density of lunar regolith simulant BP-1 [kg / m^3]         

cpAir = 1006.14                 # Constant pressure specific heat of air at 20C 1atm [J / K kg]
cpWater = 4184.05               # Constant pressure specific heat of water at 20C 1atm [J / K kg]
cpOxygen = 918.95               # Constant pressure specific heat of oxygen at 20C 1atm [J / K kg]
cpNitrogen = 1041.34            # Constant pressure specific heat of nitrogen at 20C 1atm [J / K kg]
cpCO2 = 846.05                  # Constant pressure specific heat of CO2 at 20C 1atm [J / K kg]

cvAir = 717.67                  # Constant volume specific heat of air at 20C 1atm [J / K kg]
cvWater = 4156.68               # Constant volume specific heat of water at 20C 1atm [J / K kg]
cvOxygen = 657.72               # Constant volume specific heat of oxygen at 20C 1atm [J / K kg]
cvNitrogen = 743.07             # Constant volume specific heat of nitrogen at 20C 1atm [J / K kg]
cvCO2 = 652.45                  # Constant volume specific heat of CO2 at 20C 1atm [J / K kg]

muAir = 18.13e-6                # Dynamic Viscosity of air at 20C [Pa s]
muWater = 0.0010005             # Dynamic Viscosity of water at 20C [Pa s]

nuAir = 15.06e-6                # Kinematic Viscosity of air at 20C [m^2 / s]
nuWater = 1.0023e-6             # Kinematic Viscosity of water at 20C [m^2 / s]

kAir = 25.87e-3                 # Thermal conductivity of air at 20C 1 bar [W / m K]
kWater = 0.59803                # Thermal conductivity of water at 20C 1 bar [W / m K]

PrAir = 0.7309                  # Prandtl number of air at 20C 1atm []

EAluminum = 69 * giga           # Elastic modulus of Aluminum [Pa] 
EAl2014 = 73.1 * giga           # Elastic modulus of Aluminum 2014-T6 [Pa] 
EAl6061 = 68.9 * giga           # Elastic modulus of Aluminum 6061-T6 [Pa] 
EAl3003 = 68.9 * giga           # Elastic modulus of Aluminum 3003-H14 [Pa] 
EGrayCastIron = 67.0 * giga     # Elastic modulus of cast iron alloy Gray ASTM 20 [Pa] 
EMalleableCastIron = 172 * giga # Elastic modulus of cast iron alloy Malleable ASTM A-197 [Pa] 
ECopper = 117 * giga            # Elastic modulus of Copper [Pa] 
EBrass = 101 * giga             # Elastic modulus of Red Brass C83400 [Pa] 
EBronze = 103 * giga            # Elastic modulus of Bronze C86100 [Pa] 
EIron = 210 * giga              # Elastic modulus of Iron [Pa] 
ESteel = 200 * giga             # Elastic modulus of many Steel alloys (Structural A-36, Structural A992, Tool L2) [Pa] 
EStainlessSteel = 193 * giga    # Elastic modulus of steel alloy Stainless 304  [Pa] 
ENickel = 170 * giga            # Elastic modulus of Nickel [Pa]
ETitaniumAlloy = 120 * giga     # Elastic modulus of titanium alloy [Ti-6Al-4V] [Pa] 
ELowStrengthConcrete = 22.1 * giga      # Elastic modulus of Low Strength Concrete [Pa] 
EHighStrengthConcrete = 29.0 * giga     # Elastic modulus of High Strength Concrete [Pa] 
EDouglasFir = 13.1 * giga       # Elastic modulus of wood structural grade Douglas Fir [Pa] 
EWhiteSpruce = 9.65 * giga      # Elastic modulus of wood structural grade White Spruce [Pa] 

GAl6061 = 26 * giga             # Shear modulus of Aluminum 6061 [Pa]
GAl3003 = 25 * giga             # Shear modulus of Aluminum 3003 [Pa]
GSteel = 80 * giga              # Shear modulus of Steel (general) [Pa]            

SyAl6061 = 276 * mega           # Tensile yield strength of Aluminum 6061-T6 [Pa]
SyAl3003 = 145 * mega           # Tensile yield strength of Aluminum 3003-H14 [Pa]

vAl6061 =  0.33                 # Poisson's ratio of Aluminum 6061 [unitless]
vAl3003 =  0.33                 # Poisson's ratio of Aluminum 3003 [unitless]
vSteel = 0.25                   # Poisson's ratio of Steel (general) [unitless]

## Dictionary Definitions:    ----------------------------------------------------------------------------------------------------------------------------------------
## Gauge Size Dictionaries
# Wire gauge and letter gauge drill size to inch diameter conversions [in] [ANSI Units]
drillGaugeDict = { 107 : 0.0019, 106 : 0.0023, 105 : 0.0027, 104 : 0.0031, 103 : 0.0035, 102 : 0.0039, 101 : 0.0043, 100 : 0.0047, 99 : 0.0051, 98 : 0.0055, 97 : 0.0059, 96 : 0.0063, 95 : 0.0067, 94 : 0.0071, 93 : 0.0075, 92 : 0.0079, 91 : 0.0083, 90 : 0.0087, 89 : 0.0091, 88 : 0.0095, 87 : 0.0100, 86 : 0.0105, 85 : 0.0110, 84 : 0.0115, 83 : 0.0120, 82 : 0.0125, 81 : 0.0130, 80 : 0.0135, 79 : 0.0145, 78 : 0.0160, 77 : 0.0180, 76 : 0.0200, 75 : 0.0210, 74 : 0.0225, 73 : 0.0240, 72 : 0.0250, 71 : 0.0260, 70 : 0.0280, 69 : 0.0292, 68 : 0.0310, 67 : 0.0320, 66 : 0.0330, 65 : 0.0350, 64 : 0.0360, 63 : 0.0370, 62 : 0.0380, 61 : 0.0390, 60 : 0.0400, 59 : 0.0410, 58 : 0.0420, 57 : 0.0430, 56 : 0.0465, 55 : 0.0520, 54 : 0.0550, 53 : 0.0595, 52 : 0.0635, 51 : 0.0670, 50 : 0.0700, 49 : 0.0730, 48 : 0.0760, 47 : 0.0785, 46 : 0.0810, 45 : 0.0820, 44 : 0.0860, 43 : 0.0890, 42 : 0.0935, 41 : 0.0960, 40 : 0.0980, 39 : 0.0995, 38 : 0.1015, 37 : 0.1040, 36 : 0.1065, 35 : 0.1100, 34 : 0.1110, 33 : 0.1130, 32 : 0.1160, 31 : 0.1200, 30 : 0.1285, 29 : 0.1360, 28 : 0.1405, 27 : 0.1440, 26 : 0.1470, 25 : 0.1495, 24 : 0.1520, 23 : 0.1540, 22 : 0.1570, 21 : 0.1590, 20 : 0.1610, 19 : 0.1660, 18 : 0.1695, 17 : 0.1730, 16 : 0.1770, 15 : 0.1800, 14 : 0.1820, 13 : 0.1850, 12 : 0.1890, 11 : 0.1910, 10 : 0.1935, 9 : 0.1960, 8 : 0.1990, 7 : 0.2010, 6 : 0.2040, 5 : 0.2055, 4 : 0.2090, 3 : 0.2130, 2 : 0.2210, 1 : 0.2280, "A" : 0.2340, "B" : 0.2380, "C" : 0.2420, "D" : 0.2460, "E" : 0.2500, "F" : 0.2570, "G" : 0.2610, "H" : 0.2660, "I" : 0.2720, "J" : 0.2770, "K" : 0.2810, "L" : 0.2900, "M" : 0.2950, "N" : 0.3020, "O" : 0.3160, "P" : 0.3230, "Q" : 0.3320, "R" : 0.3390, "S" : 0.3480, "T" : 0.3580, "U" : 0.3680, "V" : 0.3770, "W" : 0.3860, "X" : 0.3970, "Y" : 0.4040, "Z" : 0.4130 }

# American wire gauge to meter diameter conversion [m]
AWGDict = {"4/0" : 0.011684, "0000" : 0.011684, "3/0" : 0.010404, "000" : 0.010404, "2/0" : 0.0092658, "00" : 0.0092658, "1/0" : 0.0082515, "0" : 0.0082515, 0 : 0.0082515, 1 : 0.0073481, 2 : 0.0065437, 3 : 0.0058273, 4 : 0.0051894, 5 : 0.0046213, 6 : 0.0041154, 7 : 0.0036649, 8 : 0.0032636, 9 : 0.0029064, 10 : 0.0025882, 11 : 0.0023048, 12 : 0.0020525, 13 : 0.0018278, 14 : 0.0016277, 15 : 0.0014495, 16 : 0.0012908, 17 : 0.0011495, 18 : 0.0010237, 19 : 0.0009116, 20 : 0.0008118, 21 : 0.0007229, 22 : 0.0006438, 23 : 0.0005733, 24 : 0.0005106, 25 : 0.0004547, 26 : 0.0004049, 27 : 0.0003606, 28 : 0.0003211, 29 : 0.0002859, 30 : 0.0002546, 31 : 0.0002268, 32 : 0.0002019, 33 : 0.0001798, 34 : 0.0001601, 35 : 0.0001426, 36 : 0.000127, 37 : 0.0001131, 38 : 0.0001007, 39 : 0.0000897, 40 : 0.0000799 }

# Standard wire gauge to meter diameter conversion [m]. AKA British Standard Wire Gauge,  Imperial Wire Gauge,  British Standard Gauge
SWGDict = {"0000000": 0.012700, "7/0": 0.012700, "000000": 0.011786, "6/0": 0.011786, "00000": 0.010973, "5/0": 0.010973, "0000": 0.010160, "4/0": 0.010160, "000": 0.009449, "3/0": 0.009449, "00": 0.008839, "2/0": 0.008839, "1/0": 0.008230, "0": 0.008230, 0: 0.008230, 1: 0.007620, 2: 0.007010, 3: 0.006401, 4: 0.005893, 5: 0.005385, 6: 0.004877, 7: 0.004470, 8: 0.004064, 9: 0.003658, 10: 0.003251, 11: 0.002946, 12: 0.002642, 13: 0.002337, 14: 0.002032, 15: 0.001829, 16: 0.001626, 17: 0.001422, 18: 0.001219, 19: 0.001016, 20: 0.000914, 21: 0.000813, 22: 0.000711, 23: 0.000610, 24: 0.000559, 25: 0.0005080, 26: 0.0004572, 27: 0.0004166, 28: 0.0003759, 29: 0.0003454, 30: 0.0003150, 31: 0.0002946, 32: 0.0002743, 33: 0.0002540, 34: 0.0002337, 35: 0.0002134, 36: 0.0001930, 37: 0.0001727, 38: 0.0001524, 39: 0.0001321, 40: 0.0001219, 41: 0.0001118, 42: 0.0001016, 43: 0.0000914, 44: 0.0000813, 45: 0.0000711, 46: 0.0000610, 47: 0.0000508, 48: 0.0000406, 49: 0.0000305, 50: 0.0000254}

# U.S. Standard Gauge for Sheet and Plate Iron and Steel to meter conversion [m] 
sheetMetalGaugeDict = {"7/0": 0.0127, "6/0": 0.0119126, "5/0": 0.0110744, "4/0": 0.0103124, "3/0": 0.009525, "2/0": 0.0087376, "1/0": 0.0079502, "0000000": 0.0127, "000000": 0.0119126, "00000": 0.0110744, "0000": 0.0103124, "000": 0.009525, "00": 0.0087376, 0: 0.0079502, 1: 0.0071374, 2: 0.0067564, 3: 0.00635, 4: 0.0059436, 5: 0.0055626, 6: 0.0051562, 7: 0.0047752, 8: 0.0045466, 9: 0.0039624, 10: 0.0035814, 11: 0.003175, 12: 0.0027686, 13: 0.0023876, 14: 0.0019812, 15: 0.00178562, 16: 0.0015875, 17: 0.00143002, 18: 0.00127, 19: 0.00111252, 20: 0.0009525, 21: 0.00087376, 22: 0.00079502, 23: 0.00071374, 24: 0.000635, 25: 0.00055626, 26: 0.00047752, 27: 0.00043688, 28: 0.00039624, 29: 0.00035814, 30: 0.0003175, 31: 0.00027686, 32: 0.00026924, 33: 0.00023876, 34: 0.00021844, 35: 0.00019812, 36: 0.0001778, 37: 0.00016764, 38: 0.00016002, 39: 0.00014986, 40: 0.0001397}

# Stubs' Iron Wire Gauge to meters conversion [m]
stubsIronGaugeDict = {"4/0" : 0.0115316, "3/0" : 0.010795, "2/0" : 0.009652, "1/0" : 0.008636, "0000" : 0.0115316, "000" : 0.010795, "00" : 0.009652, "0" : 0.008636, 0 : 0.008636, 1 : 0.00762, 2 : 0.0072136, 3 : 0.0065786, 4 : 0.0060452, 5 : 0.005588, 6 : 0.0051562, 7 : 0.004572, 8 : 0.004191, 9 : 0.0037592, 10 : 0.0034036, 11 : 0.003048, 12 : 0.0027686, 13 : 0.002413, 14 : 0.0021082, 15 : 0.0018288, 16 : 0.001651, 17 : 0.0014732, 18 : 0.0012446, 19 : 0.0010668, 20 : 0.000889, 21 : 0.0008128, 22 : 0.0007112, 23 : 0.000635, 24 : 0.0005588, 25 : 0.000508, 26 : 0.0004572, 27 : 0.0004064, 28 : 0.0003556, 29 : 0.0003302, 30 : 0.0003048, 31 : 0.000254, 32 : 0.0002286, 33 : 0.0002032, 34 : 0.0001778, 35 : 0.000127, 36 : 0.0001016 }

# Stubs' Steel Wire Guage to meters conversion [m]
stubsSteelGaugeDict = {1 : 0.0057658, 2 : 0.0055626, 3 : 0.0053848, 4 : 0.0052578, 5 : 0.0051816, 6 : 0.0051054, 7 : 0.0050546, 8 : 0.0050038, 9 : 0.0049276, 10 : 0.0048514, 11 : 0.0047752, 12 : 0.004699, 13 : 0.0046228, 14 : 0.004572, 15 : 0.0045212, 16 : 0.004445, 17 : 0.0043688, 18 : 0.0042672, 19 : 0.0041656, 20 : 0.0040894, 21 : 0.0039878, 22 : 0.003937, 23 : 0.0038862, 24 : 0.0038354, 25 : 0.0037592, 26 : 0.0037084, 27 : 0.0036322, 28 : 0.0035306, 29 : 0.0034036, 30 : 0.0032258, 31 : 0.003048, 32 : 0.002921, 33 : 0.0028448, 34 : 0.002794, 35 : 0.0027432, 36 : 0.0026924, 37 : 0.0026162, 38 : 0.0025654, 39 : 0.0025146, 40 : 0.0024638 }

# Drill Fractions, wire gauge, and letter guage drill sizes to inch diameter conversions [in] [ANSI Units]
drillSizeDict = {str(Fraction(i / 64)) : i / 64 for i in range(1,65)}
drillSizeDict.update(drillGaugeDict)
# sort dictionary by values for use in drillSize function
drillSizeDict = {key: val for key, val in sorted(drillSizeDict.items(), key=lambda item: item[1]) }

## Constant Aliases:    ----------------------------------------------------------------------------------------------------------------------------------------
jpi, pii, pij = ipi, ipi, ipi 
lightsec = c
m_e = me
m_p = mp
m_n = mn
ma, m_a, malpha, m_alpha, mhe, m_He, m_he = mHe , mHe, mHe, mHe, mHe, mHe, mHe          # Warning: potential conflict with mA = milli-Amp != ma
qe, q_e, Qe, Q_e, eV, ev = e, e, e, e, e, e                 # Warning: potential conflict with exponential growth constant e^x
N_a, NA, N_A, Avogadro, avogadro = Na, Na, Na, Na, Na
k_B, kb = kB, kB                                            # Warning: potential conflict with kibi = kilobyte 
mev, Mev = MeV, MeV
AMU, da, dalton = amu, amu, amu
u_0, mu0, mu_0 = u0, u0, u0
e_0, epsilon_not, epsilon_0, eps_0, eps0 = e0, e0, e0, e0, e0
u_B, muB, mu_B = uB, uB, uB
u_N, muN, mu_N = uN, uN, uN
msun, m_sun, m_Sun = mSun, mSun, mSun
rsun, r_sun, r_Sun, Rsun, RSun, R_sun, R_Sun = rSun, rSun, rSun, rSun, rSun, rSun, rSun
mearth, m_earth, m_Earth = mEarth, mEarth, mEarth
rearth, r_earth, r_Earth, Rearth, REarth, R_earth, R_Earth = rEarth, rEarth, rEarth, rEarth, rEarth, rEarth, rEarth
Z_0, z0, z_0 = Z0, Z0, Z0
Rbar, rbar, gasconst, gasconstant = R, R, R, R
rt2 = sqrt2
rt3 = sqrt3
a_0 = a0
R_inf, Rinf, R_oo, R_infinity, Rinfinity = Roo, Roo, Roo, Roo, Roo
R_H = RH                                                # Warning: potential conflict with individual molar gas constant for hydrogen (not implemented)

if "NumPy" in sys.modules:
    sigx, sigma1, sig1 = sigmax, sigmax, sigmax
    sigy, sigma2, sig2 = sigmay, sigmay, sigmay
    sigz, sigmaz, sigz = sigmaz, sigmaz, sigmaz
    H1 = H
    Cx, cx, cX, CNOT, Cnot, cnot = CX, CX, CX, CX, CX, CX
    Cy, cy, cY = CY, CY, CY
    Cz, cz, cZ = CZ, CZ, CZ
    ket0 = zp
    ket1 = zn
    sig0 = I2
    Swap, swap = SWAP, SWAP
    CSwap, Cswap, cswap, Fredkin, fredkin = CSWAP, CSWAP, CSWAP, CSWAP, CSWAP
    Ccnot, ccnot, CCnot, CCNot, Toffoli, toffoli = CCNOT, CCNOT, CCNOT, CCNOT, CCNOT, CCNOT 

Zepto, zm, zs = zepto, zepto, zepto
Atto, am = atto, atto                           # 'as' is reserved word in python
Femto, fm, fs, fempto, Fempto = (femto,)*5
Pico, pm, ps, pF = pico, pico, pico, pico
Nano, nm, ns, ppb, PPB, ppB, nF = nano, nano, nano, nano, nano, nano, nano
Micro, mum, mus, mu_m, mu_s, muC, muc, muF, muf, muA, um, us, uC, uc, uF, uf, uA, ppm, PPM, ppM = micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro, micro
Milli, mm, ms, mA, mJ, mF, mf, mC, mv, mV, gram, cP, cp, centipoise = milli, milli, milli, milli, milli, milli, milli, milli, milli, milli, milli, milli, milli, milli            # Warning: potential conflict with ma = mass of alpha particle != mA
Centi, cm, cs = centi, centi, centi
poise, Poise = deci, deci
deka
Hecto, hm, hs = hecto, hecto, hecto
Kilo, km, ks, kPa, kpa, kW, kw, kJ, kj, Mg, kN, kn, kohm, kHz, khz, thousand = kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo, kilo
Mega, Mm, Ms, MPa, MPA, Mpa, MW, MJ, Mj, Gg, km2, Mohm, million, MHz, Mhz = mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega, mega
Giga, Gm, Gs, GPa, GPA, Gpa, gpa, GW, Gw, gw, GJ, Gj, gj, km3, billion, bil, GHz, Ghz, ghz, gHz = giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga, giga
Tera, Tm, Ts, km4, trillion = tera, tera, tera, tera, tera
Peta, Pm, Ps, quadrillion = peta, peta, peta, peta
Exa, Em, Es, quintillion = exa, exa, exa, exa
Zetta, Zm, Zs, sextillion = zetta, zetta, zetta, zetta

Ki, KiB, Kib, kib, KB, kilobyte, kiloB   = kibi, kibi, kibi, kibi, kibi, kibi, kibi             # Warning: potential conflict with boltzman constant kB
Mi, MiB, Mib, mib, MB, megabyte, megaB   = mebi, mebi, mebi, mebi, mebi, mebi, mebi
Gi, GiB, Gib, gib, GB, gigabyte, gigaB   = gibi, gibi, gibi, gibi, gibi, gibi, gibi
Ti, TiB, Tib, tib, TB, terabyte, teraB   = tebi, tebi, tebi, tebi, tebi, tebi, tebi
Pi, PiB, Pib, pib, PB, petabyte, petaB   = pebi, pebi, pebi, pebi, pebi, pebi, pebi
Ei, EiB, Eib, eib, EB, exabyte, exaB     = exbi, exbi, exbi, exbi, exbi, exbi, exbi
Zi, ZiB, Zib, zib, ZB, zettabyte, zettaB = zebi, zebi, zebi, zebi, zebi, zebi, zebi
Yi, YiB, Yib, yib, YB, yottabyte, yottaB = yobi, yobi, yobi, yobi, yobi, yobi, yobi

minutes = minute            # Note: min is also aliased to minute as a callableInt so it can also be used as the builtin minimum function.
hrs, hour, hours, Whr, Wh, Ahr, Ah = hr, hr, hr, hr, hr, hr, hr
days = day
weeks, Week, Weeks = week, week, week
mon, Mon, Month, months, Months = month, month, month, month, month
yrs, year, years = yr, yr, yr
In, inch, inches, iN, inn, in1 = IN, IN, IN, IN, IN, IN                             # 'in' is a reserved word in python
In2, inch2, inches2, iN2, inn2, IN2 = in2, in2, in2, in2, in2, in2
In3, inch3, inches3, iN3, inn3, IN3 = in3, in3, in3, in3, in3, in3
In4, inch4, inches4, iN4, inn4, IN4 = in4, in4, in4, in4, in4, in4
feet, Ft = ft, ft
feet2, Ft2 = ft2, ft2
feet3, Ft3, cfs = ft3, ft3, ft3
yds, yard, yards =  yd, yd, yd
Mi, miles, mile = mi, mi, mi
Mi2, sqmi, SqMi, sqMi = mi2, mi2, mi2, mi2
parsec, Pc = pc, pc
megaparsec, MegaParsec, mpc, MPc, megaParsec = Mpc, Mpc, Mpc, Mpc, Mpc
lightyear, Lightyear, Ly = ly, ly, ly
lightd = lightday
lighthr, lighthrs = lighthour, lighthour
A, AA, Angstrom, Ao = angstrom, angstrom, angstrom, angstrom
kips, Kips, Kip, klbf, KIP, KIPS, = kip, kip, kip, kip, kip, kip
liters, l, L, Liter, lit, Lit = liter, liter, liter, liter, liter, liter
ml = mL
patm, Patm = atm, atm
PSI = psi
kpsi = ksi
Mpsi = Msi
BTU, Btu = btu, btu
HP, Hp, horsepower = hp, hp, hp
FlOz, Floz = floz, floz
pt = pint                           # Warning: potential conflict with point
qt = quart
gal, gallons, Gal = gallon, gallon, gallon
ounce, Oz, ounces = oz, oz, oz
degR, degr, degF, degf, Fahrenheit, fahrenheit, fahr, fahren  = rankine, rankine, rankine, rankine, rankine, rankine, rankine, rankine
mh, mWhr, mwhr, mWh, mwh, mAhr, mAhr, mAh, mAh, mWH, mwH, mAHr, mAHr, mAH, mAH = mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr, mhr
kh, kWhr, kwhr, kWh, kwh, kAhr, kAhr, kAh, kAh, kWH, kwH, kAHr, kAHr, kAH, kAH = khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr, khr
Mh, MWhr, Mwhr, MWh, Mwh, MAhr, MAhr, MAh, MAh, MWH, MwH, MAHr, MAHr, MAH, MAH = Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr, Mhr
cc = cm3
cSt, cst, centistoke = mm2, mm2, mm2
Cel, cel, Celsius, Celcius, celcius, = celsius, celsius, celsius, celsius, celsius
cal = calorie
Cal, Calorie = kcal, kcal
Rps, RPS = rps, rps
Rpm, RPM = rpm, rpm
lbfft, ftlb, ftlbf, ftlbs, lbsft = lbft, lbft, lbft, lbft, lbft
lbfin, lbinn, lbfinn, lbIn, lbfIn, lbIN, lbfIN, lbinch, lbfinch, inlb, innlb, inlbf, innlbf, inlbs, lbsin = lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin, lbin
ftkip, Kipft = kipft, kipft
Kipin, Kipinn, kipinn, kipIn, kipIN, kipinch, inKip, inkip, innkip, Inkip = kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin, kipin
Barn, b = barn, barn
fB = fb
torr = Torr
mmhg, mmMercury, inmercury = mmHg, mmHg, mmHg
inhg, inMercury, inmercury = inHg, inHg, inHg
mmh2o, mmH2o, mmH20, mmh20, mmwater, mmWater = mmH2O, mmH2O, mmH2O, mmH2O, mmH2O, mmH2O
inh2o, inH2o, inH20, inh20, inwater, inWater, innh2o, innH2o, innH20, innh20, innwater, innWater = inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O, inH2O
dyn, Dyne = dyne, dyne
statC, statc, Fr, fr, franklin, esu = statcoulomb, statcoulomb, statcoulomb, statcoulomb, statcoulomb, statcoulomb
gimp, gimperial, gcust, gcus, gcustomary = gft, gft, gft, gft, gft
ginn, ginch, ginches, gIN, gIn = gin, gin, gin, gin, gin
stories, floors = story, story              # floor is also aliased as a callabeInt
St, st = stoke, stoke
NM, NMi, nmi = nauticalMile, nauticalMile, nauticalMile
Knot, knots = knot, knot
Thou, thous, mil, mils, Mils, Mil, MIL = thou, thou, thou, thou, thou, thou, thou
IPM = ipm
SFM = sfm

bolt0  = screw0 
bolt1  = screw1 
bolt2  = screw2 
bolt3  = screw3 
bolt4  = screw4 
bolt5  = screw5 
bolt6  = screw6 
bolt8  = screw8 
bolt10 = screw10
bolt12 = screw12

Rair = RAir
Rwater, Rh2o, RH2o, RH2O, RH20, Rh20 = RWater, RWater, RWater, RWater, RWater, RWater
RO2, Ro2, Roxygen = ROxygen, ROxygen, ROxygen
RN2, Rnitrogen = RNitrogen, RNitrogen
RCarbonDioxide, RCarbon_Dioxide, RcarbonDioxide, Rcarbon_Dioxide, Rcarbondioxide, Rcarbon_dioxide, Rco2, RCo2, RC02, rc02 = RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2, RCO2

densityAir, densAir, AirDensity, airDensity, airDens, PAir, pair, pAir, rhoAir = rhoAir, rhoAir, rhoAir, rhoAir, rhoAir, rhoAir, rhoAir, rhoAir, rhoAir
densityWater, densWater, WaterDensity, pwater, pWater, ph2o, pH2o, pH2O, pH20, ph20, rhowater, rhoH2O, rhoH2o, rhoh2o = rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater, rhoWater
densityO2, densO2, O2Density, pO2, po2, pOxygen, poxygen, rhoO2, rhooxygen = rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen, rhoOxygen
densityN2, densN2, N2Density, pN2, pNitrogen, rhoN2, pnitrogen, rhonitrogen = rhoNitrogen, rhoNitrogen, rhoNitrogen, rhoNitrogen, rhoNitrogen, rhoNitrogen, rhoNitrogen, rhoNitrogen
densityCarbonDioxide, densCarbonDioxide, CarbonDioxideDensity, pCarbonDioxide, pCarbon_Dioxide, pcarbonDioxide, pcarbon_Dioxide, pcarbondioxide, pcarbon_dioxide, pco2, pCo2, pC02, pCO2, rhoCarbonDioxide, rhocarbondioxide, rhoCarbondioxide = rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2, rhoCO2
densityAl6061, densal6061, al6061Density, density6061, dens6061, Al6061Density, pal6061, pAl6061, rhoal6061, pAl6, pAl60, p6061, rho6061 = rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061, rhoAl6061
densityAl3003, densal3003, al3003Density, density3003, dens3003, Al3003Density, pal3003, pAl3003, rhoal3003, p3003, rho3003 = rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003, rhoAl3003
densitySteel, densSteel, steelDensity, pSt, pSteel, rhoSt, rhosteel, psteel =  rhoSteel, rhoSteel, rhoSteel, rhoSteel, rhoSteel, rhoSteel, rhoSteel, rhoSteel

cpair = cpAir
cpwater, cph2o, cpH2o, cpH2O, cpH20, cph20 = cpWater, cpWater, cpWater, cpWater, cpWater, cpWater
cpO2, cpo2, cpoxygen = cpOxygen, cpOxygen, cpOxygen
cpN2, cpnitrogen = cpNitrogen, cpNitrogen
cpCarbonDioxide, cpCarbon_Dioxide, cpcarbonDioxide, cpcarbon_Dioxide, cpcarbondioxide, cpcarbon_dioxide, cpco2, cpCo2, cpC02, = cpCO2, cpCO2, cpCO2, cpCO2, cpCO2, cpCO2, cpCO2, cpCO2, cpCO2

cvair = cvAir
cvwater, cvh2o, cvH2o, cvH2O, cvH20, cvh20 = cvWater, cvWater, cvWater, cvWater, cvWater, cvWater
cvO2, cvo2, cvoxygen = cvOxygen, cvOxygen, cvOxygen
cvN2, cvnitrogen = cvNitrogen, cvNitrogen
cvCarbonDioxide, cvCarbon_Dioxide, cvcarbonDioxide, cvcarbon_Dioxide, cvcarbondioxide, cvcarbon_dioxide, cvco2, cvCo2, cvC02= cvCO2, cvCO2, cvCO2, cvCO2, cvCO2, cvCO2, cvCO2, cvCO2, cvCO2

muair, uair, uAir = muAir, muAir, muAir
muwater, uwater, uWater = muWater, muWater, muWater

nuair, vair, vAir = nuAir, nuAir, nuAir
nuwater, vwater, vWater = nuWater, nuWater, nuWater

kair = kAir
kwater = kWater

Prair = PrAir

EAl = EAluminum
EAluminum2014 = EAl2014 
EAluminum6061, E6061 = EAl6061, EAl6061
EAluminum3003, E3003 = EAl3003, EAl3003
ECastIronGray, EIronGray, EGrayIron = EGrayCastIron, EGrayCastIron, EGrayCastIron
ECastIronMalleable, EIronMalleable, EMalleableIron = EMalleableCastIron, EMalleableCastIron, EMalleableCastIron
ECu = ECopper
EBrass
EBronze
EFe = EIron
ESt, Est, Esteel = ESteel, ESteel, ESteel
ESS, ESs, Ess = EStainlessSteel, EStainlessSteel, EStainlessSteel
ENi = ENickel
ETiAlloy, ETi6Al4V, ETiAlV = ETitaniumAlloy, ETitaniumAlloy, ETitaniumAlloy
ELSConcrete = ELowStrengthConcrete
EHSConcrete = EHighStrengthConcrete
EDouglasFir
EWhiteSpruce

G6061 = GAl6061 
G3003 = GAl3003 
GSt = GSteel

Sy6061 = SyAl6061 
Sy3003 = SyAl3003 

v6061 = vAl6061
v3003 = vAl3003
vSt, nuSteel, nuSt = vSteel, vSteel, vSteel

## Intra Program Global Variables:    ----------------------------------------------------------------------------------------------------------------------------------------
ans = [0]
ans1 = 0
ans2 = 0
degreesFlag = False                             # Default trig functions in Radians Mode
expDegreesFlag = False                          # Default e^x in Radians Mode
floatDeltaPercent = nano                        # Accepatable percent error for floating point number comparisions
floatDeltaAbs = femto                           # Accepatable absolute error for floating point number comparisions
suppressWarningsFlag = False                    # Global boolean designed to toggle whether functions will print warnings or not.
warningTimer = Timer(0, lambda : "")            # Timer thread used for implementing temporary suppression of warnings.    
decimalContext = decimal.Context(prec=15)

## Class Definitions:    ----------------------------------------------------------------------------------------------------------------------------------------

class callableInt(builtins.int):
    """ The callableInt class is designed to be able to be used as both an int and a function. 
    I recognize that this is kinda stupid and terrible practice in the "real world", but for this script specifically it has some utility.
    
    Example: min = callableInt(minute, builtins.min)        # min can be used as both an alias for minute and to find the minimum element in an arraylike
    """
    def __new__(cls, value, function):
        return super(callableInt, cls).__new__(cls, value)
    
    def __init__(self, value, function):
        self.function = function
    
    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)

class callableFloat(builtins.float):
    """ The callableFloat class is designed to be able to be used as both an float and a function. 
    I recognize that this is kinda stupid and terrible practice in the "real world", but for this script specifically it has some utility.
    
    Example: min = callableInt(minute, builtins.min)        # min can be used as both an alias for minute and to find the minimum element in an arraylike
    """
    def __new__(cls, value, function):
        return super(callableFloat, cls).__new__(cls, value)
    
    def __init__(self, value, function):
        self.function = function
    
    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)
    
## Class Aliases:    ----------------------------------------------------------------------------------------------------------------------------------------
Int, integer, Integer, _int_, __int__ =  builtins.int, builtins.int,  builtins.int,  builtins.int,  builtins.int

## Object Definitions:    ----------------------------------------------------------------------------------------------------------------------------------------

min = callableInt(minute, builtins.min)         # min can be used as both an alias for minute and to find the minimum element in an arraylike
floor = callableFloat(story, math.floor)        # floor can be used as both an alias for story  

## Function Definitions:    ----------------------------------------------------------------------------------------------------------------------------------------

def warn(warning):
    """ Designed primarily to be a helper function for printing warning messages to the user in intra program functions.
    Will print warning if the global boolean variable suppressWarningsFlag = False, else will do nothing.
    
    I decided not to use the python warnings module because of the way it printed in the terminal. 

    Args:
        warning (str): warning to be printed
    """
    if (suppressWarningsFlag):
        return
    print("Warning: ", warning)
    print("\t[To suppress all warnings like this, use suppressWarnings()]")


def suppressWarnings(val="", time = math.inf):
    """ Sets the value of the global suppressWarningsFlag. If no value is provided, the current state will be toggled.
    If suppressWarningsFlag = True, no warnings will be printed.
    If a time is provided, then warnings will be suppressed for the specified duration (in seconds) and then toggled back on.
    
    Args:
        val (bool, optional): Value to set suppressWarningsFlag to. Default behavior is to toggle current value.
        time (float, optional): Duration used to temporarily turn off warning messages for a specified amount of time.
    """
    global suppressWarningsFlag, warningTimer
    
    if (0 <= time < math.inf):
        warningTimer.cancel()
        warningTimer = Timer(time, suppressWarnings, args=(False,))
        warningTimer.start()
        val = True
    
    if (val == ""):
        suppressWarningsFlag = not suppressWarningsFlag
    else:
        suppressWarningsFlag = val
    if (suppressWarningsFlag):
        print(f'Suppressing all warning messages {("for " + str(time) + " seconds") if time < inf else "indefinitely"}.')
    else:
        print("No longer suppressing all warning messages.")
        
silence = lambda time=math.inf : suppresswarning(True, time)
verbose = lambda time=math.inf : suppresswarning(False, time)

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
#     # edge case: float1 ~= float2 ~= 0 -> return true if absolute val of both less than floatDeltaAbs
#     return (abs(float1) < floatDeltaAbs and abs(float2) < floatDeltaAbs)


# Using the same function as math.isclose, but overriding the defualt parameters.
def isclose(a, b, *, rel_tol=floatDeltaPercent, abs_tol=floatDeltaAbs):
    return math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol) 

isclose.__doc__ = math.isclose.__doc__
    
def quad(a, b, c):
    """Quadratic formula find solutions to: a x^2 + b x + c = 0

    Args:
        a (float): 
        b (float): 
        c (float): 
    """
    global ans, ans1, ans2

    discriminant = b ** 2 - 4 * a * c
    
    if discriminant == 0:
        ans1 = -b / (2 * a)
        -b / (2 * a)
        ans.insert(0, ans1)
        print(ans1)
        return
        
    ans1 = (-b + sqrt(discriminant) ) / (2 * a)
    ans2 = (-b - sqrt(discriminant) ) / (2 * a)
    ans.insert(0, ans1)
    ans.insert(0, ans2)
    print(ans1)
    print(ans2)
    return

def sqrt(x):
    if(isinstance(x, complex) or x < 0):
        return cmath.sqrt(x)
    return math.sqrt(x)
    
# Converts a parameter from radians to degrees, or sets the default behavior of normal trig functions to be in degrees if no paramaters passed
def degrees(val = None):
    if(val != None):
        return math.degrees(val)
    
    global degreesFlag 
    degreesFlag = True
    print("Now in Degrees Mode")
    
    if(expDegreesFlag):
        print("Note, exp(1j * x) still will be computed as if x is measured in radians. To toggle this, use expDegrees()")

# Converts a parameter from degrees to radians, or sets the default behavior of normal trig functions to be in radians if no paramaters passed
def radians(val = None):
    if(val != None):
        return math.radians(val)
    
    global degreesFlag
    degreesFlag = False
    print("Now in Radians Mode")
    
    if(not expDegreesFlag):
        print("Note, exp(1j * x) still will be computed as if x is measured in degrees. To toggle this, use expRadians()")

def expDegrees():
    global expDegreesFlag
    expDegreesFlag = True
    print("Now in Exponential Degrees Mode. ( exp(180 * 1j) == -1 )")

def expRadians():
    global expDegreesFlag
    expDegreesFlag = False
    print("Now in Exponential Radians Mode. (Default;  exp(pi * 1j) == -1 )")
    

# Degree based trig functions
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

# Radian based Trig Funcitons
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

# Trig Functions that check whether in degree mode or radian mode, and give warnings 
def sin(x):
    # if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            warn("Complex number detected but in degree mode, sinr(z) calculated anyway. Use sind(z) if that is what you intended")
        return sinr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            warn("In Degree Mode but possible radian number detected. (sind(x) calculated anyway)")
        return sind(x)
    
    # Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180 
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
        warn("In Radian Mode but possible degree number detected. (sinr(x) calculated anyway)")
    return sinr(x)

def cos(x):
    # if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            warn("Complex number detected but in degree mode, cosr(z) calculated anyway. Use cosd(z) if that is what you intended")
        return cosr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            warn("In Degree Mode but possible radian number detected. (cosd(x) calculated anyway)")
        return cosd(x)
    
    # Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
        warn("In Radian Mode but possible degree number detected. (cosr(x) calculated anyway)")
    return cosr(x)

def tan(x):
    # if x is a complex number, will calculate result in radians regardless
    if(isinstance(x, complex)):
        if(degreesFlag):
            warn("Complex number detected but in degree mode, tanr(z) calculated anyway. Use tand(z) if that is what you intended")
        return tanr(x)

    if(degreesFlag):
        # Throws warning if in degree mode and x is near a multiple of pi/180
        if( (abs(x* 180/pi - round(x * 180/pi) ) / floatDeltaPercent < 1 ) ):
            warn("In Degree Mode but possible radian number detected. (tand(x) calculated anyway)")
        return tand(x)
    
    # Throws warning if in radian mode and x > 2pi and is not near a multiple of pi/180
    if( (x > 2 * pi) and not (abs(( x * 180/pi) - round( x * 180 / pi ) )/ floatDeltaPercent < 1 ) ):
        warn("In Radian Mode but possible degree number detected. (tanr(x) calculated anyway)")
    return tanr(x)

# returns arcsin(x), returns radian value if in radians mode or if complex number
def asin(x):
        ans = cmath.asin(x)
        if(ans.imag != 0):
            return ans
        
        return degrees(ans.real) if degreesFlag else ans.real
    
# returns arccos(x), returns radian value if in radians mode or if complex number
def acos(x):
        ans = cmath.acos(x)
        if(ans.imag != 0):
            return ans
        
        return degrees(ans.real) if degreesFlag else ans.real
    
# returns arctan(x), returns radian value if in radians mode or if complex number
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
    
def exp(x=1):
    """Returns e**x
    If x = 1j*theta, theta will be converted from degrees to radians if degreesFlag = true
    exp() will return exp(1) = 2.718281828459045

    Args:
        x (_type_, optional): Defaults to 1.

    Returns:
        float or complex: e**x
    """
    if(isinstance(x, complex)):
        if(x.real == 0):
            if(expDegreesFlag):
                warn("Exponential Degrees flag on, converting x from degrees to radians")
                return cos(x.imag) + 1j * sin(x.imag)
        return cmath.exp(x)
    return math.exp(x)

def lg(x):
    if(isinstance(x, complex)):
        return cmath.log2(x)
    return math.log(x) / cmath.log(2)

sci = lambda x : f"{x:.4e}"       # Formats input in scientific notation with 4 decimal places

# Temperature conversion functions

k2c = lambda x : x - 273.15                     # Kelvin to Celsius conversion
k2f = lambda x : (x - 273.15) * 9/5 + 32        # Kelvin to Fahrenheit conversion
k2r = lambda x : x * 1.8                        # Kelvin to Rankine conversion

c2k = lambda x : x + 273.15                     # Celsius to Kelvin conversion
c2f = lambda x : (x * 9/5) + 32                 # Celsius to Fahrenheit conversion
c2r = lambda x : (x * 9/5) + 491.67             # Celsius to Rankine conversion

f2k = lambda x : (x - 32) * 5/9 + 273.15        # Fahrenheit to Kelvin conversion
f2c = lambda x : (x - 32) * 5/9                 # Fahrenheit to Celcius conversion
f2r = lambda x : x + 459.67                     # Fahrenheit to Rankine conversion

r2k = lambda x : x * 5/9                        # Rankine to Kelvin conversion
r2c = lambda x : (x - 491.67) * 5/9             # Rankine to Celcius conversion
r2f = lambda x : x - 459.67                     # Rankine to Fahrenheit conversion

# area of circle given diameter
darea = lambda d : pi * (d/2)**2 

# Volume of sphere given diameter
dvolume = lambda d : 4/3 * pi * (d/2)**3

def factor(num: int):
    for i in range(1, math.floor(sqrt(num)+1) ):
        if( num % i == 0 ):
            print(f"{i} * {(num / i):.0f}")

def frac(num: float):
    return Fraction(num).limit_denominator()

# det = lambda mat : np.linalg.det(mat)         

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

def integral(f, start, end, steps = 1_000_000):
    """ Integrates function/lamda expression f in interval (start, end) using equal step sizes and trapezoid summation
        
        Example use: integral(lambda x : 5 * x**2 + 10, 0, 5)
    """
    deltaX = (end-start)/steps
    sum = 0
    x = start
    for i in range(int(steps)):
        sum += deltaX * (f(x) + f(x+deltaX) ) / 2
        x += deltaX
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

def mean(arr, *args, **kwargs):
    """Compute mean / averarge of arraylike. 
    If weights parameter is provided, then the return value is returned value of np.average(), else the return value is the returned value of np.mean()
    """
    if("weights" in kwargs):
        return np.average(arr, *args, **kwargs)
    return np.mean(arr, *args, **kwargs)


# TODO: Correct this to use *args, **kwargs
# Returns sample standard deviation by default
def std(arr, weights=None, population=False):
    """Returns standard deviation of array. Weights array is optional, but cumsum should be normalized to 1.
    If weight array is given, the standard deviation is calculated using the formula sig^2 = sum(x[i]**2 * p[i]) - mean(x, p)**2
    Else will use np.std.

    Args:
        arr (arraylike): Data values
        weights (arraylike, optional): Weights / probabilties for each value in arr. Must be normalized such that sum(weights) = 1. Defaults to None.
        population (bool, optional): Boolean to specify whether to calculate population (True) or sample (False) standard deviation. Defaults to False.

    Returns:
        float: Standard deviation of arr
    """
    if(weights):
        return sqrt( var(arr, weights))
    if(population):
        return np.std(arr)
    return np.std(arr,ddof=1)

# TODO: Correct this to use *args, **kwargs
# Returns sample variance by default
def var(arr, weights=None, population=False):
    """Returns variance of array. By default, the sample variance is returned unless population parameter is set to True.
    Weights array is optional, but weights sum should be normalized to 1.
    If weight array is given, the variance is calculated using the formula sig^2 = sum(x[i]**2 * p[i]) - mean(x, p)**2
    Else will use np.var.

    Args:
        arr (arraylike): Data values
        weights (arraylike, optional): Weights / probabilties for each value in arr. Must be normalized such that sum(weights) = 1. Defaults to None.
        population (bool, optional): Boolean to specify whether to calculate population (True) or sample (False) variance. Defaults to False.

    Returns:
        float: Variance of arr
    """
    if (weights):
        
        var = 0
        for x, p in zip(arr, weights):
            var += x**2 * p
        
        return var - mean(arr, weights)**2

        
    if(population):
        return np.var(arr)
    return np.var(arr,ddof=1)



def birthdayProblem(days=365):
    """ Given an arbritrary number of days, returns the number of people required for Probability[two people share the same birthday] >= 0.5   

    Args:
        days (int, optional): _description_. Defaults to 365.

    Returns:
        int: smallest number of people required for Probability[two people share the same birthday] >= 0.5   
    """
    return ceil( (3 - 2 * ln(2))/6 + sqrt(2*days*ln(2)) + (9 - 4*ln(2)**2 ) / (72 * sqrt(2*days*ln(2))) - 2 * ln(2)**2 / (135*days) )


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
    magnitude argument is optional, and will default to 1 if only 2 arguments are given.
    
    Examples:
        40 * singularity(x-2, 1)
        singularity(40, x-2, 1)

    Args:
        args[0]     = magnitude (float): optional argument
        args[0 or 1] = input (float): x - start
        args[1 or 2] = degree (int): 

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


def step(x):
    """ Step function

    Args:
        x (float): input

    Returns:
        float: 1 if x >= 0, else 0
    """
    return singularity3(1, x, 0)

def printfunction(f, start, end, stepSize = nan, steps = 50, outputResolution = 4, inputResolution = None, smallScientific=True):
    """Prints values of function f with range [start : stepSize : end]

    Args:
        f (callable): Single variable function
        start (float): start of range
        end (float): end of range
        stepSize (float, optional): distance between successive x values 
        steps (int, optional): number of steps of function to print (50 by default if stepSize not also specified)
        outputResolution (int, optional): Number of decimal places to print output. Defaults to 4.
        inputResolution (int, optional): Number of decimal places to print input. Defaults to automatically calculated based on start, end, and the magnitude of stepSize.
        smallScientific (bool, optional): Determines whether small nonzero output numbers will be automatically formatted in scientific notation instead of the default format if they would have otherwise be printed as 0.0000 . Defaults to True.
    """
    
    if(isnan(stepSize)):
        stepSize = (end - start) / steps
    
    if (inputResolution == None):
        # inputResolution = max(0, ceil(-log10(abs(stepSize))))     # resoultion = number of decimal places needed to format x, minimum of 0      // TODO: stepSize = 0.125 -> inputResolution = 1 (should be 3)
        if ((stepSize - floor(stepSize) == 0) and (start - floor(start) == 0) and (end - floor(end) == 0)):
            inputResolution = 0
        else:
            inputResolution = -2 + max( len(floatToStr(stepSize - floor(stepSize) ) ), len(floatToStr(start - floor(start) ) ), len(floatToStr(end - floor(end) ) ) )     # Maybe fixed
        if (abs(stepSize) >= 1):
            inputResolution = min(3, inputResolution)
        elif (abs(stepSize) >= 0.01 ):
            inputResolution = min(5, inputResolution)
        else:
            inputResolution = min(8, inputResolution)
        
    inputSize = (inputResolution > 0) + 1 + floor(max(log10(abs(start) + (start == 0)) + 2 * (start < 0), log10(abs(end) + (end == 0)) + 2 * (end < 0))) # maybe fixed?
    inputFormat = f"{inputSize+inputResolution}.{inputResolution}f"
    
    # print(f"{inputResolution = }, {inputSize = }, {inputFormat = }, {stepSize = }")
    threshold = 0.5 * 10**(-outputResolution)
    
    x = start
    yArr = []
    while (x < end or floatComparison(x, end)):
        try:
            yArr.append(f(x))
            # if(smallScientific and abs(temp) < threshold and temp != 0):
            #     print(f"{x :{inputFormat}}  : {' ' if temp > 0 else ''}{ temp :9.{outputResolution}e} ")
            # else:
            #     print(f"{x :{inputFormat}}  : {' ' if temp > 0 else ''}{ temp :9.{outputResolution}f} ")
        except (ZeroDivisionError):
            yArr.append(nan)

        x += stepSize
    
    # minSize / maxSize are the calculated lengths of the min and max elements, not including the trailing decimal values
    # If the magnitude of min and/or max << 1, cap the size to the desired number of decimals
    minSize = 0 if abs(np.nanmin(yArr)) == 0 else abs(max(-outputResolution, floor(log10(abs(np.nanmin(yArr)))))) + (np.nanmin(yArr) < 0)
    maxSize = 0 if abs(np.nanmax(yArr)) == 0 else abs(max(-outputResolution, floor(log10(abs(np.nanmax(yArr)))))) + (np.nanmax(yArr) < 0)
        
    # maxDigits = 1 + (decimal point) + (# of decimals) + (# of digits before decimal including negative sign)
    maxDigits = 1 + (outputResolution > 0) + outputResolution + max(0, maxSize, minSize)        
    sciFormat = ""
    threshold = 0.5 * 10**(-outputResolution)
    if (smallScientific):
        smallestElement = min( [abs(yArr[i]) for i in np.nonzero(yArr)[0]] )
        if (smallestElement < threshold):
            # maxDigits = 1 + (decimal point) + (# of decimals) + len("e-05") + (negative sign)
            maxDigits = max(maxDigits, 1 + 1 + outputResolution + 4 + (-smallestElement in yArr) )
    
    x = start
    for y in yArr:
        if(smallScientific and abs(y) < threshold and y != 0):
            print(f"{x :{inputFormat}} : {y :{maxDigits}.{outputResolution}e} ")
        else:
            print(f"{x :{inputFormat}} : {y :{maxDigits}.{outputResolution}f} ")
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
    
    
#TODO: small numbers still don't work - 1e-3 = no match
#TODO: implement other functions - sin, cos, tan, log10, lg, exp, 1/x
#TODO: function too good at finding representations - some are ridiculous
def represent(input, deltaPercent = 0.000_1, denom = 200, output=30, paramFunction=nan, paramCoef=nan):
    """Finds the math representation of input  (ie. determines if input can be represented as a fraction of pi) 

    Args:
        input (_type_): _description_
        deltaPercent (_type_, optional): Minium percent error for a number to be flagged as a near match. Defaults to 0.000_1.
        denom (_type_, optional): Max denominator that function will check. Defaults to 200.
        output (int, optional): Max number of candidates to output. Defaults to 30
        paramFunction (string, optional): #TODO
    """
    exactPercent = 1e-18
    candidateList = []  # List of canditate tuples with near but not exact representations. candidate = (str representation, candidateValue, abs(inputValue - candiateValue) )
    constList = {1 : "1", pi : "pi", sqrt(2) : "sqrt(2)", sqrt(3) : "sqrt(3)", sqrt(5) : "sqrt(5)", exp(1) : "exp(1)", 1 - exp(-1): "(1 - exp(-1))"}
    functionList = {lambda x: x : "{candidate}", lambda x : x**2 : "sqrt({candidate})", lambda x : exp(x) : "ln({candidate})"}
    functionInverseList = {"{candidate}" : lambda x: x, "sqrt({candidate})": lambda x: sqrt(x), "ln({candidate})" : lambda x: ln(x) }
    maxRepLen = 0
    candidate = ""
    for function in functionList.keys():
        val1 = function(input)
        functionInverse = functionInverseList[functionList[function]]
        for const in constList.keys():
            val2 = val1 / const
            checkedRatios = {0}     # Used to remove non reduced fractions, ie. 2/4 when 1/2 has already been checked
            
            # for loop searches for fraction with denominator [1,denom] , [1,200] by default
            for i in range(1, denom+1):
                if( (round(val2 * i) / i) in checkedRatios):
                    continue
            
                # Checks if i is a viable denominator
                if( isclose(input, functionInverse(const * round(val2 * i) / i ), rel_tol=deltaPercent ) ):
                    # if exact match, print it and stop searching
                    if( isclose(input, functionInverse(const * round(val2 * i) / i ), rel_tol=exactPercent) ):
                        if(const == 1):
                            candidate = f"{round(val2 * i)}{f' / {i}' if (i != 1) else ''}"
                        else:
                            candidate = f"{constList[const]}{f' * {round(val2 * i)}' if (round(val2 * i) != 1) else ''}{f' / {i}' if (i != 1) else ''}"
                        print("Found match: \n" + eval(f'f\"{functionList[function]}\"'))
                        return
            
                    else:
                        checkedRatios.add(round(val2 * i) / i)
                        if(const == 1):
                            candidate = f"{round(val2 * i)}{f' / {i}' if (i != 1) else ''}"
                        else:
                            candidate = f"{constList[const]}{f' * {round(val2 * i)}' if (round(val2 * i) != 1) else ''}{f' / {i}' if (i != 1) else ''}"
                        candidateList.append((eval(f"f'{functionList[function]}'") , functionInverse(const * round(val2 * i) / i) , abs(input - functionInverse(const * round(val2 * i) / i ))))    
                        maxRepLen = max(maxRepLen, len(eval(f"f'{functionList[function]}'")))
                        
    print("Did not find exact representation, here are closest matches:")
    candidateList.sort(key = lambda candidate: candidate[2])    # sort candidate list by smallest to largest distance
    inputSize = len(str(input))
    inputResolution = len(str(input)) - 1 - max(0, floor(log10(input)))
    print(f"{'True input':<{maxRepLen+1}} : {f'{input:{inputSize+2}.{inputResolution}f}':^{inputSize+2}}") 
    count = 0
    for candidate in candidateList:
        print(f"{candidate[0]:<{maxRepLen+1}} : {f'{candidate[1]:{inputSize+2}.{inputResolution}f}':^{inputSize+2}}  :  {f'{candidate[2]:g}':<12} ")
        count += 1
        if (count >= output):
            break 
        

def props(prop1: str, val1, prop2: str, val2, fluid: str, molFlag=False):
    """ Wrapper function for extracting data from CoolProp.
    Adds functionality enabling user to input specific volume as a fluid property.
    
    Patch 1/24/23: fluid can now be listed as first or last parameter
        i.e. props("T", 300, "P", 100_000, "water") === props("water", "T", 300, "P", 100_000)
    
    # Examples:
    # props("T", 300, "P", 100_000, "water")
    # props("T", 300, "P", 100_000, "water", molFlag=True)  # Outputs molar specific properties instead of mass specific properties
    # props("T", 400, "D", 996.556340388, "water")
    # props("P", 100000, "T", 400, "air")
    # props("T", 300, "v", 0.001003455559, "CO2")
    # props("", , "", , "")     #blank example

    Args:
        prop1 (str): _description_
        val1 (_type_): _description_
        prop2 (str): _description_
        val2 (_type_): _description_
        fluid (str): _description_
        molFlag (bool, optional): Toggles whether to output mass specific properties (False, default) or molar specific properties (True) 
        
    Valid CP inputs: 
    'P' = Pressure [Pa]
    'T' = Temperature [K]
    'Q' = Vapor quality [kg/kg]
    'D' = Density [kg / m^3]
    'S', 'Smass', 'SMASS' = Mass specific entropy [J/kg/K]
    'U', 'Umass', 'UMASS' = Mass specific internal energy [J/kg]
    'H', 'Hmass', 'HMASS' = Mass specific enthalpy [J/kg]
    Valid inputs with props function only:
    'V', 'v' = specific volume [m^3 / kg]
    
    Full Fluid List:
    1-Butene, Acetone, Air, Ammonia, Argon, Benzene, CarbonDioxide, CarbonMonoxide, CarbonylSulfide, cis-2-Butene, CycloHexane, Cyclopentane, CycloPropane, D4, D5, D6, Deuterium, Dichloroethane, DiethylEther, DimethylCarbonate, DimethylEther, Ethane, Ethanol, EthylBenzene, Ethylene, EthyleneOxide, Fluorine, HeavyWater, Helium, HFE143m, Hydrogen, HydrogenChloride, HydrogenSulfide, IsoButane, IsoButene, Isohexane, Isopentane, Krypton, m-Xylene, MD2M, MD3M, MD4M, MDM, Methane, Methanol, MethylLinoleate, MethylLinolenate, MethylOleate, MethylPalmitate, MethylStearate, MM, n-Butane, n-Decane, n-Dodecane, n-Heptane, n-Hexane, n-Nonane, n-Octane, n-Pentane, n-Propane, n-Undecane, Neon, Neopentane, Nitrogen, NitrousOxide, Novec649, o-Xylene, OrthoDeuterium, OrthoHydrogen, Oxygen, p-Xylene, ParaDeuterium, ParaHydrogen, Propylene, Propyne, R11, R113, R114, R115, R116, R12, R123, R1233zd(E), R1234yf, R1234ze(E), R1234ze(Z), R124, R1243zf, R125, R13, R134a, R13I1, R14, R141b, R142b, R143a, R152A, R161, R21, R218, R22, R227EA, R23, R236EA, R236FA, R245ca, R245fa, R32, R365MFC, R40, R404A, R407C, R41, R410A, R507A, RC318, SES36, SulfurDioxide, SulfurHexafluoride, Toluene, trans-2-Butene, Water, Xenon  
    """
    # Checks if user put fluid before properties, and corrects the variables, ie. props("water", "T", 300, "P", 100_000) 
    if(isinstance(val1, str) and not isinstance(fluid, str)):
        # prop1, val1, prop2, val2, fluid
        tempFluid =  prop1
        prop1 = val1
        val1 = prop2
        prop2 = val2
        val2 = fluid
        fluid = tempFluid
    
    print(f"\t Input:  {fluid} :  {prop1} = {val1}  ;  {prop2} = {val2} \n")
    
    proplst = [prop1, prop2]
    vallst = [val1, val2]
    # for loop to add functionality to CP such as ability to input specific volume
    for i in [0,1]:
        # Specific Volume
        if(proplst[i] in ["v", "V", "Vol", "vol"]):
            proplst[i] = 'D'
            vallst[i] = 1 / vallst[i]
        
    prop1 = proplst[0]
    prop2 = proplst[1]
    val1 = vallst[0]
    val2 = vallst[1]
        
    print("Temperature T [K] =  ", CP.PropsSI('T', prop1, val1, prop2, val2, fluid) )
    print("Pressure P [Pa] =  ", CP.PropsSI('P', prop1, val1, prop2, val2, fluid) )
    print("Density D [kg / m^3] =  ", CP.PropsSI('D', prop1, val1, prop2, val2, fluid) )
    print("Specific Volume ν [m^3 / kg] =  ", 1 / CP.PropsSI('D', prop1, val1, prop2, val2, fluid) )
    print("Vapor quality Q or χ [kg/kg] =  ", CP.PropsSI('Q', prop1, val1, prop2, val2, fluid) )
    print()
    
    if(molFlag):
        print("Molar Mass M [kg/mol] =  ", CP.PropsSI('molarmass', prop1, val1, prop2, val2, fluid) )
        print("Molar Enthalpy [J/mol] =  ", CP.PropsSI('Hmolar', prop1, val1, prop2, val2, fluid) )
        print("Molar Entropy S [J/mol/K] =  ", CP.PropsSI('Smolar', prop1, val1, prop2, val2, fluid) )
        print("Molar Internal Energy U [J/mol] =  ", CP.PropsSI('Umolar', prop1, val1, prop2, val2, fluid) )
        print("Molar Gibbs Free Energy g [J/mol] =  ", CP.PropsSI('Gmolar', prop1, val1, prop2, val2, fluid) )
    else:
        print("Enthalpy h [J/kg] =  ", CP.PropsSI('Hmass', prop1, val1, prop2, val2, fluid) )
        print("Entropy s [J/kg/K] =  ", CP.PropsSI('Smass', prop1, val1, prop2, val2, fluid) )
        print("Internal Energy u [J/kg] =  ", CP.PropsSI('Umass', prop1, val1, prop2, val2, fluid) )
        print("Gibbs Free Energy g [J/kg] =  ", CP.PropsSI('Gmass', prop1, val1, prop2, val2, fluid) )
    print()
    
    print("Constant-pressure specific heat c_p [J/kg/K] =  ", CP.PropsSI('Cpmass', prop1, val1, prop2, val2, fluid) )
    print("Constant-volume specific heat c_v [J/kg/K] =  ", CP.PropsSI('Cvmass', prop1, val1, prop2, val2, fluid) )
    print("Specific Heat Ratio c_p / c_v [] =  ", CP.PropsSI('Cpmass', prop1, val1, prop2, val2, fluid)  / CP.PropsSI('Cvmass', prop1, val1, prop2, val2, fluid) )
    print("Specific Molar Gas Consant R [J/kg/K] =  ", CP.PropsSI('gas_constant', prop1, val1, prop2, val2, fluid) / CP.PropsSI('molarmass', prop1, val1, prop2, val2, fluid) )
    print()
    
    print("Isentropic expansion coefficent [] =  ", CP.PropsSI('isentropic_expansion_coefficient', prop1, val1, prop2, val2, fluid) )
    print("Isobaric expansion coefficent [1/K] =  ", CP.PropsSI('isobaric_expansion_coefficient', prop1, val1, prop2, val2, fluid) )
    print("Isothermal compressibility [1/Pa] =  ", CP.PropsSI('isothermal_compressibility', prop1, val1, prop2, val2, fluid) )
    print()
    
    print("Speed of sound [m/s] =  ", CP.PropsSI('speed_of_sound', prop1, val1, prop2, val2, fluid) )
    print()
    
    # print(" =  ", CP.PropsSI('', prop1, val1, prop2, val2, fluid) )


def thermoPlot(xAxis: str, xInitial: float, xFinal: float, yAxis: str, constParam: str, constParamVal: float, fluid: str, title: str = ""):
    """ Wrapper function for creating plots with Matplotlib and CoolProp.
    
    Will attempt to look at xAxis and yAxis labels and use non base-SI units if present within [] brackets
        - ie: xAxis="Density [kg/liter]"
        - only works on units with conversion variables predefined
        - Warning - uses python eval() on unsanitized string inputs
    
    # Examples:
    thermoPlot("Temp", 300, 400, "P", "D", 900, "CO2")      # Plots pressure vs temperature of CO2 at constant density 900 kg/m^3 from 300K to 400K
    thermoPlot("T", 300, 400, "P", "D", 900, "CO2", "Pressure vs temperature of CO2 at 900 kg/m^3")      # Plots pressure vs temperature of CO2 at constant density 900 kg/m^3 from 300K to 400K, and adds title to graph
    thermoPlot("Pressure [MPa]", 2*Mpa, 80*Mpa, "Density [kg / liter]", "T", 300, "CO2")    # Plots density (in units kg / liter) vs (pressure in units of MPa) of CO2 at constant temperatrue 300 K from 2 MPa to 80 MPa

    Args:
        xAxis (str): Label for x-axis - can use standard CoolProp labels if plotting in SI (for instance, xAxis="P" will default to x axis label "Pressure [Pa]")
                     For nonstandard units, type out full xAxis label followed by units in []
                     Non base SI unit example: xAxis="Density [kg/liter]"
        xInitial (float): starting point of x axis to be plotted
        xFinal (float): end point of x axis to be plotted
        yAxis (str): Label for y-axis - can use standard CoolProp labels if plotting in SI (for instance, yAxis="P" will default to x axis label "Pressure [Pa]")
                     For nonstandard units, type out full yAxis label followed by units in []
                     Non base SI unit example: yAxis="Pressure [MPa]"
        constParam (str): thermodynamic quantity to remain constant in this plot
        constParamVal (float): value of thermodynamic quanity to remain constant
        fluid (str): fluid to be plotted
        title (str, optional): Plot title. Defaults to "".
    """
    
    print("Importing Matplotlib - may take a while")
    import matplotlib.pyplot as plt
    print("Matplotlib successfully imported")
    
    steps = 100
    deltaX = (xFinal - xInitial) / steps
    xscale = 1
    yscale = 1
    xplt = []
    yplt = []
    
    plt.figure()  # create a new figure
    
    # Give the x axis a descriptive label
    
    Txunit = "K"
    
    if('[' in xAxis):
        units = xAxis[xAxis.find('[')+1 : xAxis.find(']')]
        J, kg, K, m, m2, m3, Pa, C = 1, 1, 1, 1, 1, 1, 1, 1
        # if(['C', 'R', 'F'] not in units):
        print('xUnits = ', units)
        if(xAxis == 'T' or 'temp' in xAxis.lower()  ):
            Txunit = units
        else:
            xscale = 1/eval(units)
        print('xscale = ', xscale)
    
    xvolflag = False
    
    if(xAxis == 'T' or 'temp' in xAxis.lower()  ):
        if(Txunit == "K"):
            plt.xlabel('Temperature [K]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'T'
    elif(xAxis == 'P' or 'pres' in xAxis.lower() ):
        if(xscale == 1):
            plt.xlabel('Pressure [Pa]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'P'
    elif(xAxis == 'Q'):
        if(xscale == 1):
            plt.xlabel('Vapor Quality [kg/kg]')
        else:
            plt.xlabel(xAxis)
    elif(xAxis == 'D' or 'dens' in xAxis.lower() ):
        if(xscale == 1):
            plt.xlabel('Density [kg/m^3]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'D'
    elif(xAxis == 'S' or xAxis == 'Smass' or 'entropy' in xAxis.lower()):
        if(xscale == 1):
            plt.xlabel('Mass Specific Entropy [J/kg/K]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'S'
    elif(xAxis == 'U' or xAxis == 'Umass' or 'internal energy' in xAxis.lower() ):
        if(xscale == 1):
            plt.xlabel('Mass Specific Internal Energy [J/kg]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'U'
    elif(xAxis == 'H' or xAxis == 'Hmass' or 'enthalpy' in xAxis.lower() ):
        if(xscale == 1):
            plt.xlabel('Mass Specific Enthalpy [J/kg]')
        else:
            plt.xlabel(xAxis)
        xAxis = 'U'
    elif(xAxis.lower() == 'v' or 'vol' in xAxis.lower() ):
        if(xscale == 1):
            plt.xlabel('Specific Volume m^3 / kg')
        else:
            plt.xlabel(xAxis)
        xvolflag = True
        xAxis = 'D'
        xInitial = xInitial**-1
        xFinal = xFinal**-1
        deltaX = (xFinal - xInitial) / steps
    else:
        plt.xlabel(xAxis)
    
    # Give y-axis a descriptive label
    
    Tyunit = "K"
    
    if('[' in yAxis):
        units = yAxis[yAxis.find('[')+1 : yAxis.find(']')]
        J, kg, K, m, m2, m3, Pa = 1, 1, 1, 1, 1, 1, 1
        # if(['C', 'R', 'F'] not in units):
        print('yUnits = ', units)
        if(yAxis == 'T' or 'temp' in yAxis.lower()  ):
            Tyunit = units
        else:
            yscale = 1/eval(units)
        print('yscale = ', yscale)
    
    yvolflag = False
    
    if(yAxis == 'T' or 'temp' in yAxis.lower()  ):
        if(Tyunit == "K"):
            plt.ylabel('Temperature [K]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'T'
    elif(yAxis == 'P' or 'pres' in yAxis.lower() ):
        if(yscale == 1):
            plt.ylabel('Pressure [Pa]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'P'
    elif(yAxis == 'Q'):
        if(yscale == 1):
            plt.ylabel('Vapor Quality [kg/kg]')
        else:
            plt.ylabel(yAxis)
    elif(yAxis == 'D' or 'dens' in yAxis.lower() ):
        if(yscale == 1):
            plt.ylabel('Density [kg/m^3]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'D'
    elif(yAxis == 'S' or yAxis == 'Smass' or 'entropy' in yAxis.lower()):
        if(yscale == 1):
            plt.ylabel('Mass Specific Entropy [J/kg/K]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'S'
    elif(yAxis == 'U' or yAxis == 'Umass' or 'internal energy' in yAxis.lower() ):
        if(yscale == 1):
            plt.ylabel('Mass Specific Internal Energy [J/kg]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'U'
    elif(yAxis == 'H' or yAxis == 'Hmass' or 'enthalpy' in yAxis.lower() ):
        if(yscale == 1):
            plt.ylabel('Mass Specific Enthalpy [J/kg]')
        else:
            plt.ylabel(yAxis)
        yAxis = 'H'
    elif(yAxis.lower() == 'v' or 'vol' in yAxis.lower() ):
        if(yscale == 1):
            plt.ylabel('Specific Volume m^3 / kg')
        else:
            plt.ylabel(yAxis)
        yvolflag = True
        yAxis = 'D'
    else:
        plt.ylabel(yAxis)
    
    for i in range(steps):
        xStep = i * deltaX + xInitial
        xplt.append(xStep * xscale)
        yplt.append(CP.PropsSI(yAxis, xAxis, xStep, constParam, constParamVal, fluid) * yscale)
    
    if(xvolflag):
        for i in range(steps):
            xplt[i] = xplt[i]**-1

    if(yvolflag):
        for i in range(steps):
            yplt[i] = yplt[i]**-1
    
    if(Txunit != "K"):
        for i in range(steps):
            if("c" in Txunit.lower()):
                xplt[i] = k2c(xplt[i])
            if("r" in Txunit.lower()):
                xplt[i] = k2r(xplt[i])
            if("f" in Txunit.lower()):
                xplt[i] = k2f(xplt[i])
    
    if(Tyunit != "K"):
        for i in range(steps):
            if("c" in Tyunit.lower()):
                yplt[i] = k2c(yplt[i])
            if("r" in Tyunit.lower()):
                yplt[i] = k2r(yplt[i])
            if("f" in Tyunit.lower()):
                yplt[i] = k2f(yplt[i])
                
    # print("xvolflag = ", xvolflag)
    # print("yvolflag = ", yvolflag)
    
    plt.title(title)
    
    plt.plot(xplt, yplt)
    plt.show()


def int(*args, **kwargs):
    """Casts to int if given single input, computes integral if given multiple inputs. 
    
    See integral(...) function for more details. 
    """
    if(len(args) == 1 and len(kwargs) == 0):
        x = args[0]
        if(isinstance(x, bool)):
            return 1 if x else 0
        if(isinstance(x, builtins.int)):
           return x
        if(isinstance(x, float)):
            return x.__int__()
        if(isinstance(x, str)):
            return float(x).__int__()
        if(isinstance(x, complex)):
            return (x.real).__int__()
        # Check if unknown object has builtin cast to int function
        try: 
            return x.__int__()
        except:
            # Per the python standard, __index__() is used to cast to int when __int__() is not defined.
            return x.__index__()
        
    return integral(*args, **kwargs)

# Angular Frequency to Frequency conversions
w2f = lambda w: w / tau 
f2w = lambda f: f * tau 


def roots(func, leftBound=-10_000, rightBound=10_000, numberOfRoots = nan, iterations = 1_000, repetitions = 2, hardCodedFlag = True):
    """Finds all roots of a function within specified range using fsolve

    Args:
        func (): single parameter function
        leftBound (_type_, optional): Defaults to -10_000.
        rightBound (_type_, optional): Defaults to 10_000.
        numberOfRoots (_type_, optional): Defaults to nan.
        iterations (int, optional): Iterations per recursion. Defaults to 1_000.
        repetitions (int, optional): Number of times algorithm repeats. Defaults to 1.
        hardCodedFlag (bool, optional): Flag for whether or not the algorithm will try hard coded guesses. Defaults to True.

    Returns:
        _type_: _description_
    """
    # Don't look too closely at this code, I am not proud of it
    
    if (numberOfRoots == nan):
        numberOfRoots = iterations + 1
    
    # For reasons, it makes sense to temporarily include the leftBound and rightBound in the rootsFound array even if they are not true roots. They may not be defined for the function, so putting in try/except blocks.
    try:
        if(not isClose(func(leftBound), 0)):
            numberOfRoots += 1
    except:
        pass
    try:
        if(not isClose(func(leftBound), 0)):
            numberOfRoots += 1
    except:
        pass
    
    rootsFound = [leftBound, rightBound]
    r = 0
    
    # In order to increase effectiveness of the algorithm, I'm hard coding some numbers in to be used as guesses
    if(hardCodedFlag):
       hardCodedGuesses = [-100, -50, -10, -5, -2, -1, 0, 1, 2, 5, 10, 50, 100]
    else:
        hardCodedGuesses = []
        
    # To keep track of whether the program has iterated through the hard coded values yet or not
    guess = 0
    
    # Adding a few iterations to solve edge case when all hard coded guesses are out of range and specified iterations number is less than 
    for i in range(iterations + len(hardCodedGuesses)):
        if(len(rootsFound) >= numberOfRoots):
                break
        
        if (hardCodedFlag):
            # First try hard coded guesses
            # Ensure that hard coded guess is within specified range
            while (guess < len(hardCodedGuesses) and (hardCodedGuesses[guess] < leftBound or hardCodedGuesses[guess] > rightBound)):
                guess += 1
            if (guess >= len(hardCodedGuesses)):
                hardCodedFlag = false
                continue
            sol, info, ier, msg = fsolve(func, hardCodedGuesses[guess], full_output=True)
            r = sol[0]
            guess += 1
        else:
            # After trying all hard coded guesses, try uniform random numbers within specified range
            randTemp = random.uniform(leftBound, rightBound)
            sol, info, ier, msg = fsolve(func, randTemp, full_output=True)
            r = sol[0]

        # when fsolve can't find a root, it returns a number anyway. Ensure that r is an actual root. (ier == 1)
        if (ier == 1):
            if ((r > rightBound) or (r < leftBound)):
                continue
            
            # Ensure no duplicates in rootsFound
            flag = false
            for x in rootsFound:
                if (isClose(x,r)):
                    flag = true
                    # Check if r is a better estimate than x
                    if(abs(func(r)) < abs(func(x)) ):
                        rootsFound.remove(x)
                        rootsFound.append(r)
                    break
            
            if(not flag):
                rootsFound.append(r)
                
    rootsFound.sort() 
    # Repeat algorithm with upper/lower bounds set to found roots to find roots between roots
    for i in range(repetitions):
        # No need to repeat, no roots found between upper and lower bound
        if(len(rootsFound) == 2):
            return rootsFound
        # Divide iterations up evenly among each range so algorithm doesn't take too long
        secondIterations = ceil(iterations / len(rootsFound))
        j = 0
        while (j < len(rootsFound) - 1):
            # Successive reptitions of algorithm sometimes improve the guess of the root value and therefore may generate duplicates 
            # Remove duplicate by removing both values and generating a better guess from their average
            if(isClose(rootsFound[j], rootsFound[j+1], rel_tol=1e-7)):
                r = (rootsFound[j] + rootsFound[j+1])/2
                # print(r)
                # print(rootsFound)
                rootsFound.remove(rootsFound[j])
                rootsFound[j] = fsolve(func, r)[0]
                rootsFound.sort()
                continue
            moreRoots = roots(func, leftBound=rootsFound[j], rightBound=rootsFound[j+1], numberOfRoots = numberOfRoots - len(rootsFound), iterations= secondIterations, repetitions=0, hardCodedFlag=False)
            for k in range(1, len(moreRoots) - 1):
                rootsFound.append(moreRoots[k])
            rootsFound.sort()
            j+=1
        
        try:
            func(rootsFound[0])
            if(not isClose(fsolve(func, rootsFound[0])[0], rootsFound[0])):
                rootsFound.remove(rootsFound[0])
        except:
            rootsFound.remove(rootsFound[0])
        try:
            func(rootsFound[-1])
            if(not isClose(fsolve(func, rootsFound[-1])[0], rootsFound[-1])):
                rootsFound.remove(rootsFound[-1])
        except:
            rootsFound.remove(rootsFound[-1])
            
    # Final pass through to improve guesses
    for i in range(len(rootsFound)):
        rootsFound[i] = fsolve(func, rootsFound[i])[0]

    # Remove any last duplicates generated in last step
    i = 0
    while (i < len(rootsFound) - 1):
        if(isClose(rootsFound[i], rootsFound[i+1], rel_tol=1e-7)):
            rootsFound.remove(rootsFound[i])
            continue
        i += 1

    return rootsFound


def boltChart():
    print(f"screw0 \t :  {screw0:.4f}") 
    print(f"screw1 \t :  {screw1:.4f}") 
    print(f"screw2 \t :  {screw2:.4f}") 
    print(f"screw3 \t :  {screw3:.4f}") 
    print(f"screw4 \t :  {screw4:.4f}") 
    print(f"screw5 \t :  {screw5:.4f}") 
    print(f"screw6 \t :  {screw6:.4f}") 
    print(f"screw8 \t :  {screw8:.4f}") 
    print(f"screw10\t :  {screw10:.4f}")
    print(f"screw12\t :  {screw12:.4f}")
    i = 1/4
    while (i < 1):
        print(f"{frac(i)}\t :  {i:.4f}")
        i += 1/32
    
# TODO: Clean up this code. Possibly create helper function for pretty printing numbers.
def printList(arr, decimals = 4, format=">{maxDigits}.{decimals}f", smallScientific=True, scientificDecimals = -2):
    """Print arraylike of items and attempt to automatically handle formatting. 
    
    Currently implemented: 
        - any 1D arraylike consisting of only ints 
        - any 1D arraylike consisting of floats and/or ints
        - any 2D np.ndarray consisting of only numbers
    
    If arr is not an implemented format, then every eleement of arr will be printed with default python print formatting.
    
    WARNING: uses python eval() on unsanitized string input for format parameter.

    Args:
        arr (arraylike): 
        decimals (int, optional): _description_. Defaults to 4.
        format (str, optional): Format string for printing numbers if arr contains only ints and floats. Defaults to ">{maxDigits}.{decimals}f" - maxDigits is automatically computed and is equal to the maximum string length of numbers in arr. 
        smallScientific (bool, optional): Determines whether small numbers will be automatically formatted in scientific notation instead of the default format if they would have otherwise be printed as 0.0000 . Defaults to True.
        scientificDecimals(int, optional): Determines number of decimals used to display small values in scientific when smallScientific=True. If a negative integer is given, then the number of decimals is calculated based on the number  
    """
    
    # Test if arr is 1 dimensional and all elements of arr are numbers
    if(all(isinstance(n, builtins.int) or isinstance(n, float) for n in arr)):
        # If all the element values are integers (type independent), don't print any decimals
        if (all(isinstance(n, builtins.int) for n in arr)  or  all(float(n).is_integer() for n in arr)):
            decimals = 0
            smallScientific = False
        # minSize / maxSize are the calculated lengths of the min and max elements, not including the trailing decimal values
        # If the magnitude of min and/or max << 1, cap the size to the desired number of decimals
        minSize = 0 if abs(np.nanmin(arr)) == 0 else abs(max(-decimals, floor(log10(abs(np.nanmin(arr)))))) + (np.nanmin(arr) < 0)
        maxSize = 0 if abs(np.nanmax(arr)) == 0 else abs(max(-decimals, floor(log10(abs(np.nanmax(arr)))))) + (np.nanmax(arr) < 0)
        
        # maxDigits = 1 + (decimal point) + (# of decimals) + (# of digits before decimal including negative sign)
        maxDigits = 1 + (decimals > 0) + decimals + max(0, maxSize, minSize)        
        sciFormat = ""
        threshold = 0.5 * 10**(-decimals)
        if (smallScientific):
            smallestElement = min( [abs(arr[i]) for i in np.nonzero(arr)[0]] )
            if (smallestElement < threshold):
                # maxDigits = 1 + (decimal point) + (# of decimals) + len("e-05") + (negative sign)
                maxDigits = max(maxDigits, 1 + 1 + decimals + 4 + (-smallestElement in arr) )
            
        format = eval('f\"' + format + '\"')
        sciFormat = format[:-1] + "e"
        
        for i in arr:
            if(smallScientific and abs(i) < threshold and i != 0):
                print(f"{i:{sciFormat}}")   
                continue
            
            print(f"{i:{format}}")   
            
        return
    
    # Test if arr is a two dimensional np.ndarray of numbers
    if (isinstance(arr, np.ndarray)  and  np.issubdtype(arr.dtype, np.number)  and  len(arr.shape) == 2 ):
        if (np.issubdtype(arr.dtype, np.integer)  or  all(float(n).is_integer() for n in np.nditer(arr)) ):
            decimals = 0
            smallScientific = False
        minSize = 0 if abs(np.min(arr)) == 0 else floor(log10(abs(np.min(arr)))) + (np.min(arr) < 0)
        maxSize = 0 if abs(np.max(arr)) == 0 else floor(log10(abs(np.max(arr)))) + (np.max(arr) < 0)
        
        # maxDigits = 1 + (decimal point) + (# of decimals) + (# of digits before decimal including negative sign)
        maxDigits = 1 + (decimals > 0) + decimals + max(0, maxSize, minSize)        
        threshold = 0.5 * 10**(-decimals)
        if (smallScientific):
            smallestElement = min( [abs(arr[i,j]) for i,j in zip(np.nonzero(arr)[0], np.nonzero(arr)[1])] )
            if (smallestElement < threshold):
                # maxDigits = 1 + (decimal point) + (# of decimals) + (negative sign)
                maxDigits = max(maxDigits, 1 + 1 + decimals + 4 + (-smallestElement in arr) )
            
        format = eval('f\"' + format + '\"')
        sciFormat = format[:-1] + "e"
        
        for row in arr:
            print("[", end="") 
            for i in row:
                if(smallScientific and abs(i) < threshold and i != 0):
                    print(f"{i:{sciFormat}}",  end=", ")   
                    continue
            
                print(f"{i:{format}}", end=", ")  
            print("\b\b]")      # Using backspace character instead of coming up with a better and more elegant solution. \shrug
        return
                
            
    # If arr is not one of the hardcoded formats above, default to simply print each element in arr
    else:
        for i in arr:
            print(i)


def expandWeightedList(data, weights):
    """Returns expanded data set given a set of data points and the frequency of each data point

    Args:
        data (arraylike): list of unique data points
        weights (arraylike): frequency list in same order as data list

    Returns:
        list: expanded data set
    """
    if (len(data) != len(weights)):
        raise ValueError("Data length does not match weight length")
    arr = []
    for i,j in zip(data, weights):
        arr += [i for _ in range(j)]
        
    return arr


def stats(arr):
    return pd.DataFrame(arr).describe()
    
    
def derivative(f, x, step = 1e-6):
    return (f(x+step) - f(x-step)) / (2 * step)     


def binompdf(n, p, X):
    """ Binomial probability mass function.

    Args:
        n (int): Number of trials
        p (float): Probability of success
        X (int): Number of successes

    Raises:
        ValueError: n and X should be integers such that n >= X
        ValueError: Probabilty should be a decimal [0, 1]

    Returns:
        float: probability of having X success given n trials
    """
    if (X > n or not isinstance(X, Integer) or not isinstance(n, Integer)):
        raise ValueError("n and X should be integers such that n >= X")
    if (p > 1 or p < 0):
        raise ValueError("Probabilty should be a decimal [0, 1]")
    return nCx(n, X) * p**X * (1 - p)**(n-X)


def binomcdf(n, p, X):
    """ Binomial cumulative density function.

    Args:
        n (int): Number of trials
        p (float): Probability of success
        X (int): Number of successes

    Raises:
        ValueError: n and X should be integers such that n >= X
        ValueError: Probabilty should be a decimal [0, 1]

    Returns:
        float: probability of having X successes or less given n trials
    """
    if (X > n or not isinstance(X, Integer) or not isinstance(n, Integer)):
        raise ValueError("n and X should be integers such that n >= X")
    if (p > 1 or p < 0):
        raise ValueError("Probabilty should be a decimal [0, 1]")
    sum = 0
    for x in range(X + 1):
        sum += nCx(n, x) * p**x * (1 - p)**(n-x)
    return sum


def poissonpdf(lam, t, x):
    """ Poisson probability mass function.

    Args:
        lam (int): Number of events in unit time/space
        t (float): Amount of sampled time/space 
        X (int): Number of successes
    
    Returns:
        float: probability of having X successes in t space
    """
    return exp(-lam * t) * (lam * t)**x / factorial(x)

def poissoncdf(lam, t, X):
    """ Poisson cumulative density function.

    Args:
        lam (int): Number of events in unit time/space
        t (float): Amount of sampled time/space 
        X (int): Number of successes
    
    Returns:
        float: probability of having X or less successes in t space
    """
    sum = 0
    for x in range(X + 1):
        sum += exp(-lam * t) * (lam * t)**x / factorial(x)
    return sum


def normalpdf(x, xbar = 0, sigx = 1):
    if (xbar != 0 or sigx != 1):
        x = (x - xbar) / sigx

    return exp(-x**2 / 2) / sqrt(2 * pi)

def normalcdf(x, xbar = 0, sigx = 1):
    if (xbar != 0 or sigx != 1):
        x = (x - xbar) / sigx

    return integral(lambda t: exp(-t**2 / 2) / sqrt(2 * pi), -100, x)

def exponentialcdf(lam, t): 
    return 1 - exp(-lam * t)

def isInteger(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


def fun2arr(f, start, end, stepSize=nan):
    if(isnan(stepSize)):
        stepSize = (end - start) / 10_000

    x = start
    arr = np.zeros( 1 + ceil((end - start) / stepSize ))
    i = 0
    while (x <= end):
        arr[i] = (f(x))
        x += stepSize
        i += 1
    return arr

def arr2fun(arr, start, end, stepSize=nan):
    if(isnan(stepSize)):
        stepSize = (end - start) / 10_000
        
    return lambda x: arr[round((x - start) / stepSize )]

def floatToStr(f):
    """Convert the given float to a string, without resorting to scientific notation
    """
    d1 = decimalContext.create_decimal(repr(f))
    return format(d1, 'f')

# TODO
def drillSize(input):
    if (isinstance(input, float) and 0 < input < 1):
        i = 0
        while (i < len(drillSizeDict.values())):
            pass
    else:
        if (isinstance(input, str)):
            input = input.upper()
        return drillSizeDict[input]    


def date(date:str, format = ""):
    adjustedStr = ""
    if (format == ""):
        # Parsing date
        date = date.replace("/","-").replace("\\","-")
        if (match := re.search(r"\d+-\d+-?\d*", date)):
            format += "%m-%d"
            adjustedStr += match.group()

            # Check if year is in format YY, else assume YYYY*
            if (re.search(r"\d+-\d+-\d\d(\D|\Z)", date)):
                format += "-%y"
            elif (re.search(r"\d+-\d+-\d{4}", date)):
                format += "-%Y"

        # Parsing Time
        if (match := re.search(r"\d+:\d+:?\d*", date)):
            adjustedStr += " " + match.group()
            format += " %H:%M"
            if (re.search(r"\d+:\d+:\d+", date)):
                format += ":%S"
            # Check for am/pm signature, else assume 24 hr time
            if (match := re.search(r"(A|P)\.?M", date.upper())):
                format = format.replace("%H", "%I")
                format += " %p"
                adjustedStr += " " + match.group().replace('.', '')
        
    if (format == ""):
        raise ValueError("Unknown date format")

    print(f"Parsed: {adjustedStr}")
    return datetime.strptime(adjustedStr, format)

def dateDifference(date1, date2 = "", format=""):
    if (isinstance(date1, dt.timedelta)):
        diff = abs(date1)
        # The duration in months can't be calulated without the original dates 
        monthFlag = False
    else:
        if (not isinstance(date1, datetime)):
            date1 = date(date1, format)
        if (not isinstance(date2, datetime)):
            date2 = date(date2, format)
        
        if (date1 < date2):
            date1, date2 = date2, date1
        
        leapYearFlag = (not (date1.year % 400)) or ((date1.year % 100) and not (date1.year % 4))
        
        daysInMonths = [31, 28 + leapYearFlag, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        monthFlag = True
        diff = date1 - date2
        monthDiff = date1.month - date2.month 
        dayDiff = date1.day - date2.day
        
        if (dayDiff < 0):
            # Add days from previous month, and reduce month count by 1
            dayDiff += daysInMonths[date1.month - 2] 
            monthDiff -= 1
        if (monthDiff < 0):
            monthDiff += 12
    
    secondStr = ""
    if (diff.seconds != 0):
        secondStr = f", {diff.seconds // 3600} hours, {diff.seconds % 3600 // 60} minutes, {diff.seconds % 60} seconds"
    
    output = ""
    if (diff.days > 0):
        # TODO is // 365 correct?
        output += f"{diff.days // 365} years, {monthDiff} months, {dayDiff // 7} weeks, {dayDiff % 7} days{secondStr}\n" if (monthFlag and (diff.days > 364)) else ""
        output += f"{diff.days // 365} years, {monthDiff} months, {dayDiff} days{secondStr}\n" if (monthFlag and (diff.days > 364)) else ""
        output += f"{diff.days // 365} years, {(diff.days % 365) // 7} weeks, {(diff.days % 365) % 7 } days{secondStr}\n" if (diff.days > 364) else ""
        output += f"{diff.days // 365} years, {(diff.days % 365) } days{secondStr}\n" if (diff.days > 364) else ""
        output += f"{diff.days // 365 * 12 + monthDiff} months, {dayDiff} days{secondStr}\n" if (monthFlag and (monthDiff or diff.days > 31)) else ""
        output += f"{diff.days // 7} weeks, {diff.days % 7} days{secondStr}\n" if (diff.days > 7) else ""
        output += f"{diff.days} days{secondStr}\n" if (diff.days) else ""
    
    output += f"{diff.days * 24 + diff.seconds // 3600} hours, {diff.seconds % 3600 // 60} minutes, {diff.seconds % 60} seconds\n" if (diff.days > 0 or diff.seconds > 3600) else ""
    output += f"{diff.days * 24 * 60 + diff.seconds // 60} minutes, {diff.seconds % 60} seconds\n" if (diff.days > 0 or diff.seconds > 60) else ""
    output += f"{diff.total_seconds()} seconds"

    print(output)

    return diff.total_seconds()
    
    

# Function Aliases:    ----------------------------------------------------------------------------------------------------------------------------------------
minimum = builtins.min
wait = sleep
floatComparison, isClose, floatcomparison, floatcomp = isclose, isclose, isclose, isclose
fact = factorial
J0 = j0
J1 = j1
Jv, jn, Jn, bessel, bessel1, firstBessel = jv, jv, jv, jv, jv, jv
outerProd, outerprod, outerProduct, outerproduct = outer, outer, outer, outer
crossprod, crossProd, crossProduct, crossproduct = cross, cross, cross, cross
stdDev, stddev, StdDev, standardDeviation = std, std, std, std
popstdDev, popstddev, popStdDev, populationStddev, populationStdDev, populationStd, populationstd, populationstandardDeviation, populationStandardDeviation = popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd, popstd
# root, roots = fsolve, fsolve
supressWarnings, suppresswarnings, supresswarnings, supressWarning, suppresswarning, supresswarning = suppressWarnings, suppressWarnings, suppressWarnings, suppressWarnings, suppressWarnings, suppressWarnings
silentWarnings, silenceWarnings, silentwarnings, silencewarnings, silentWarning, silenceWarning, silentwarning, silencewarning, silent = silence, silence, silence, silence, silence, silence, silence, silence, silence
verboseWarnings, verbosewarnings, verboseWarning, verbosewarning = verbose, verbose, verbose, verbose
nCr, nCk, nCx, ncr, nck, ncx, combinations = math.comb, math.comb, math.comb, math.comb, math.comb, math.comb, math.comb
nPr, nPk, nPx, npr, npk, npx, permutations = math.perm, math.perm, math.perm, math.perm, math.perm, math.perm, math.perm
deg, Deg, Degrees, degreesMode, degreesmode, DegreesMode, DegreesMode, degMode, degmode = degrees, degrees, degrees, degrees, degrees, degrees, degrees, degrees, degrees
rad, Rad, Radians, radiansMode, radiansmode, Radiansmode, RadiansMode, radMode, radmode = radians, radians, radians, radians, radians, radians, radians, radians, radians
expdeg, expDeg, expdegrees, expdegreesMode, expdegreesmode, expDegreesMode, expDegreesMode, expdegMode, expdegmode = expDegrees, expDegrees, expDegrees, expDegrees, expDegrees, expDegrees, expDegrees, expDegrees, expDegrees
exprad, expRad, expradians, expradiansMode, expradiansmode, expRadiansmode, expRadiansMode, expradMode, expradmode = expRadians, expRadians, expRadians, expRadians, expRadians, expRadians, expRadians, expRadians, expRadians
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
InnerProduct, innerproduct, Innerproduct, inprod, inproduct, braket, BraKet, Braket, dotproduct, dotProduct, dotprod, dotProd, dot, inner = innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod, innerprod
magnitude, mag = norm, norm
integrate = integral            # int is also (sort of) aliased to integral - see int documentation
average, avg = mean, mean
variance = var
populationVariance, popvariance, popVariance = popvar, popvar, popvar
birthdayproblem, birthdayparadox, birthdayParadox, generalizedBirthdayProblem = birthdayProblem, birthdayProblem, birthdayProblem, birthdayProblem
sing = singularity
printFunction, printfun, printFun, plotFun, plotfun = printfunction, printfunction, printfunction, printfunction, printfunction
eigen = eig
representation = represent
coolpropPlot, coolPropPlot, CoolPropPlot, CoolpropPlot, cpPlot, CPPlot, plot = thermoPlot, thermoPlot, thermoPlot, thermoPlot, thermoPlot, thermoPlot, thermoPlot
printlist, printTable, printtable, printMatrix, printMat, printArr, printarr, printArray, printarray = printList, printList, printList, printList, printList, printList, printList, printList, printList
boltchart, screwChart, screwchart, boltsizes, screwsizes, boltSizes, screwSizes = boltChart, boltChart, boltChart, boltChart, boltChart, boltChart, boltChart
expandList, expandlist, expandArr, expandarr, weightedList, weightedlist, weightList, weightlist = expandWeightedList, expandWeightedList, expandWeightedList, expandWeightedList, expandWeightedList, expandWeightedList, expandWeightedList, expandWeightedList
ddx, ddt, deriv, dydx, dydt, dfdx, dfdt = derivative, derivative, derivative, derivative, derivative, derivative, derivative
binompmf, binpdf, binpmf = binompdf, binompdf,  binompdf
binomcmf, bincdf, bincmf = binomcdf, binomcdf,  binomcdf
normcdf = normalcdf
normpdf = normalpdf
isInt = isInteger
Date, DATE, dateParser, dateParse, parseDate, Time, TIME, timeParser, timeParse, parseTime = date, date, date, date, date, date, date, date, date, date
dateDiff, timeDifference, timeDiff, timediff, datediff = dateDifference, dateDifference, dateDifference, dateDifference, dateDifference

# print(f"total time = {total}")