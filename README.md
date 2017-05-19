# simulation-related 
A few examples of codes used for the data analysis.

### rename.py --> for outputs (e.g. simulations with Gasoline) which havethe format sim_name.XXXXX
=========

This program renames simulation files. If you output every Xth run of your simulation, i.e. the snapshot numbers follow a sequence with a step X (start, start+X, start+2X, ...), with this program you can change the ordering to be an array of integers with a step X=1. 

### APEX --> X-ray luminosity calculator. Uses APEC output models.
=========

Calculation of the X-ray luminosity of your gasoline file.
Requires packages: astropy, pynbody, numpy, pickle, gc
Tested for atomdb_v3.0.2 (download here: http://www.atomdb.org/download.php).

Goal: compute total X-ray luminosity in the 0.5-2keV band within a virial radius of a halo, as a sum of the contribution from lines and a continuum.
