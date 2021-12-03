# calibrate-hpge
Small program that generates quadratic calibration parameters using analysis trees. Original code by Stephen Gillespie.

# Installation:
Compile the program using the command
```
g++ calib-ge.cxx -std=c++0x  -o Calibrate
```

# Usage: 
To generate quadratic calibration parameters for GRIFFIN dataset first pass the known source calibration files to the program preceded by the source name, e.g. 
```
./Calibrate 60Co /path/to/analysistree.root 152Eu /path/to/analysistree.root 133Ba /path/to/analysistree.root
```
The program requires a 60Co source, but the others are optional and you can include and many as you want. You may need to define energies and errors in the ```quad_energy.C``` file but the most common calibration sources are already defined. **Note:** you must have a calibration file in the same directory as the code titled ```CalibrationFile.cal```; this is necessary to get the proper channel addessses for the sorting of the analysis trees charge histograms. 

# Outputs
The program outputs multiple files you can use for diagnostics and the final parameters are written to ```quad_energy_coeff.txt``` in the form of three arrays. 

# Further Notes:
More information about this program and the calibration method for GRIFFIN's HPGe detectors is located at: https://grsi.wiki.triumf.ca/index.php/Stephens_page_of_things
