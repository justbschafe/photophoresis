# Overview

This is the repository for analytical and numerical code used to calculate the photophoretic forces in "Photophoretic flight of perforated structures in near-space conditions" by Benjamin C. Schafer, Jong-hyoung Kim, Felix Sharipov, Gyeong-Seok Hwang, Joost J. Vlassak, and David W. Keith.

# Analytical model

Mathematica 14.1 was used to write and compile the analytical model. For system requirements of the latest version of Mathematica, see https://www.wolfram.com/mathematica/system-requirements/index.php.en. For an introduction to using Mathematica, see https://www.wolfram.com/wolfram-u/courses/wolfram-language/hands-on-start-to-mathematica-wl005/. The analytical code in this repository was run successfully on MacOS 15.1 and Windows 11. 

The Mathematica notebook has the following structure: 
* Model functions
  * Constants and functions
  * DSMC data for a permeable membrane
  * Hybrid force models 
* Lofting force calculations
* 2D temperature models
  * Microscale
  * Macroscale   

**Model functions.** 
In this section, the analytical model is defined.

**Constants and functions.** 
Physical constants for the relevant environment (e.g. Earth's atmosphere) are defined. Mathematica's Standard Atmosphere package (https://reference.wolfram.com/language/StandardAtmosphere/tutorial/StandardAtmosphere.html.en) is used to define atmospheric properties, and those in all other environments are given in the paper's Supplementary Information.

**DSMC data for a permeable membrane.** 
DSMC data from the numerical model is input. This data can be found in the .csv file in this repository. Sample data for a porosity of 0.5 is imported.

**Hybrid force models.** 
The analytical model is defined as stated in the methods of the paper's main text. 

**Lofting force calculations.** 
In this section, the model is used to calculate photophoretic forces of interest. Descriptions of the model's inputs are given, and sample calculations of the lofting force on a structure with porosity 0.5 are shown. Calculating the lofting force for a single set of input parameters typically takes a few ms. 

**2D Temperature calculations.** 
In this section, the 2D model described in the paper's Supplementary Information is defined on two scales: microscale in the plane of the sample and macroscale around a uniform thin disk.

**Microscale.** 
The temperature difference between the two membranes is calculated on the unit cell that defines the hexagonal hole pattern.

**Microscale.** 
The local temperature around a thin disk is calculated by assuming a known heat flux on the surface of the disk and a fixed temperature at the edge of the domain.


# Numerical model

The files are used to calculate the dimensionless radiomentric force and heat fluxes on a thin disk.

**The minimum requirements to run the code are:**
* A cluster with 60 cores
* Memory of 150 MB per node
* MPI fortran compiler

**List of files:**
* ```disc.for``` contains Fortran code for the main calculations.
* ```p.inc``` is an auxiliary file containing input data for ```disc.for```.
* ```Res.for``` contains secondary Fortran code to form a results file.
* ```d.sh``` is a script file that submits the main code to the queue.
* ```Res.dat``` is a file created after the secondary code (```Res.for```) runs.
* The ```MATRIX``` directory contains databases of viscosities and deflection angles.

**Steps to perform the calculations:**
1. Compile the main code:
```mpif77 -w -fallow-argument-mismatch -O3 -o d.exe disc.for```
2. Compile the secondary code:
```gfortran -o R.exe Res.for```
3. Run the main codes via the queue. If allowed by the cluster, run it directly:
```mpirun -np 60 d.exe```
4. Run the secondary code (```Res.for```):
```./R.exe```
5. The last step generates the files ```Res.dat``` and ```F.dat```. The former contains the total force, total energy fluxes, local pressures and local energy fluxes. The latter contains the flow-field.


