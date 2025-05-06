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

**Model functions.** In this section, the analytical model is defined.

**Constants and functions** In this subsection, physical constants for the relevant environment (e.g. Earth's atmosphere) are defined. Mathematica's Standard Atmosphere package (https://reference.wolfram.com/language/StandardAtmosphere/tutorial/StandardAtmosphere.html.en) is used to define atmospheric properties, and those in all other environments are given in the paper's Supplementary Information.

**DSMC data for a permeable membrane** In this subsection, DSMC data from the numerical model is input. This data can be found in the .csv file in this repository. Sample data for a porosity of 0.5 is imported.

**Hybrid force models** In this subsection, the analytical model is defined as stated in the methods of the paper's main text. 

**Lofting force calculations** In this section, the model is used to calculate photophoretic forces of interest. Descriptions of the model's inputs are given, and sample calculations of the lofting force on a structure with porosity 0.5 are shown. Calculating the lofting force for a single set of input parameters typically takes a few ms. 

**2D Temperature calculations** 

# Numerical model

Fortran 77 was used to write and compile the numerical model.
