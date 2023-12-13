# EasyEppler
Scripts and tools to make using the Eppler Code easier and less error prone, plus enabling automization
The Eppler Code is NOT part of this, since it does not belong to me in any way

## Content
### EpplerFoil.py 
a wrapper class that makes creating the input files for the Eppler Code easy and can execute the Code
(!This requires proper placement of the files!)

### BasicEasyFoil.py
is a skript that generates a basic laminar airfoil from a small amount of inputs.
Its based on the NACA design phylosophy and additionaly features a function to reduce suction peaks
This can be used as a very baseline starting point for an airfoil design

## Installation
If you just wana create the input file nothing special is necessary just put all these files in a folder and work in that.
Currently only EpplerFoil.py is really of value, FoilClassShowcase.py is a tutorial on basic functions and could be deleted.
If you want to use the ExecuteEppler command and automatically run the code with the created file, your python musst run in
the folder the Eppler Executable is in. I mostly do this by starting Visual Studio Code in the folder the Eppler Exes are in.

## Requirements (that must be manually installed)
EpplerFoil.py -> Numpy
BasicEasyFoil.py -> Numpy, Aerosandbox