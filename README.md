# FILL IT!

A Python interface to compute pit filling using priority-flood on structured and unstructured meshes.

## Dependencies

### FillIT

+ [numpy](http://numpy.org)
+ [meshplex](https://github.com/nschloe/meshplex)
+ fortran compiler, preferably [gfortran](https://gcc.gnu.org/wiki/GFortran)

### Notebooks

To run the Examples (jupyter notebooks) it is recommended to use the `gscape-docker` image which has all the libraries already installed.

[https://hub.docker.com/u/geodels/dashboard/](https://hub.docker.com/u/geodels/dashboard/)

## Installation 

```bash
python setup.py install
```

## Usage

Two classes are included as part of the **FillIT** package:

+ `depressionFillingZhou`
+ `depressionFillingBarnes`

These classes share similar methods for both structured and unstructured mesh and can be easily interchanged.

For regular grids, the following `classes` are available:
+ `depressionFillingZhou`
+ `depressionFillingBarnes`

To call one of these classes, you will typically do as follows:

``` python
import fillit as pitfill

# Class initialisation
pitClass = pitfill.depressionFillingBarnes(dem)

# Performing pit filling
fillZ = pitClass.performPitFillingSruct()

```

## References

Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
Vol 62, Jan 2014, pp 117–127.

Zhou, Sun, Fu. "An efficient variant of the Priority-Flood algorithm for filling
depressions in raster digital elevation models". Computers & Geosciences.
Vol 90, Feb 2016, pp 87–96.
