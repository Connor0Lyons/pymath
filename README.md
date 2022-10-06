# pymath
### Intended to be a helper script for using an interactive python terminal as a calculator.

## Features:
- Variables for universal constants
- Variables for unit conversions to base SI units
- Variables for some material constants
- Redefined functions (such as sqrt, sin, etc) to have more functionality and practicality
- Misc other functions I find useful

## Use:
With python3 installed, your terminal/shell of choice to execute the following command(s)
```bash
cd [Directory where pymath.py is located]
python3 -i pymath.py
```
OR
```bash
python3 -i [file address of directory where pymath.py is located]/pymath.py
```

Within the python terminal, write any valid python expression and the result will output on the following newline. '_' is a variable that gets updated after every command and contains the output of the previous line, similar to the Ans button on a traditional calculator.
Variables and functions were given aliases so you don't have to remember what the exact variable name is and you can usually hopefully guess correctly (ie. yr, yrs, year, years are all equal to the number of seconds in a year). 

## Warnings
- I tried to comment everything clearly within pymath.py, so hopefully everything is at least clear in its purpose.
- The variable 'e' refers to the elementary charge. For the exponential growth constant e^x, use exp(x)
- Unit variables are designed to convert non-SI untis to SI. Unit variables are defined such that multiplying a number by the unit variable specifies that number is measured in that unit (ie 10 yards converted to meters would be written as 10 * yards = 9.144)
- All variables are in SI units unless specified otherwise. Divide by a unit variable to convert to that unit (ie to convert 9.144 meter to yards, 9.144/yards = 1.093613)
- 'in' is a reserved word in python, and therefore cannot be used as variable name for inches. The following variables can be used instead: In, inch, inches, iN, inn, in1

## Disclaimer
This script was designed for my personal use, and as such all features were implemented as I considered them convenient and useful. This is public on github due to requests from friends. Feature requests and bug reports may be considered if I consider them personally useful.
