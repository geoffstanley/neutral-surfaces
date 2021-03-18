[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4615747.svg)](https://doi.org/10.5281/zenodo.4615747)

# Neutral Surfaces
## Software for approximately neutral surfaces and geostrophic streamfunctions

### Approximately Neutral Surfaces

In the absence of irreversible mixing, fluid parcels in the ocean move such that they are always at their level of neutral buoyancy with their environment.  This constrains the direction of flow lie in a plane, called the neutral tangent plane. This plane is well-defined at every point in the ocean.  A neutral surface is an extensive 2D surface that is everywhere aligned with the neutral tangent plane (McDougall 1987).  However, neutral surfaces are not well-defined 2D surfaces -- a consequence of the non-linear equation of state for the density of seawater, and the resulting non-zero neutral helicity in the ocean (McDougall and Jackett 1988).  Hence, physical oceanographers craft approximately neutral surfaces, which are well-defined extensive 2D surfaces that are everywhere approximately aligned with the neutral tangent plane.  The following approximately neutral surfaces may be calculated with this software package. 

#### Potential density surfaces

The potential density (Wust 1935) is the density a seawater parcel would have if adiabatically and isentropically moved to a given reference pressure.  A potential density surface is an isosurface of the 3D potential density field. 

#### Specific volume anomaly surfaces

The specific volume anomaly (Montgomery 1937) is the difference between the  in-situ specific volume and the specific volume of a seawater parcel at the local pressure but having a reference practical / Absolute salinity and a reference potential / Conservative temperature.   A specific volume anomaly surface is an isosurface of the 3D specific volume anomaly field.  In the Boussinesq approximation, it is more useful to consider the in-situ density anomaly (defined similarly) and its isosurfaces. 

#### Omega surfaces

Omega surfaces (Klocker et al., 2009;  Stanley et al. 2021) are highly accurate approximately neutral surfaces, formed from an iterative procedure that solves a global least squares problem to minimize the neutrality error. 

#### Topobaric surfaces

Topobaric surfaces (Stanley 2019a) are approximately neutral surfaces that are highly accurate, fast to compute, and possess an exact geostrophic streamfunction (though it is ill-defined, as it is for truly neutral surfaces). Topobaric surfaces are built from a multivalued functional relationship between the in-situ density and the pressure that exists on a truly neutral surface, or on a topobaric surface. The single-valued branches of this function are valid on regions that are determined by the Reeb graph of the pressure on the surface. Topobaric surfaces are constructed from an initial approximately neutral surface of any quality using an iterative procedure.  Each iteration calculates the Reeb graph, empirically fits the density to the pressure using simple functions, then updates the surface such that the density on the surface matches that given by the these simple functions. 

Modified topobaric surfaces are similar to topobaric surfaces, but their empirical fits have extra constraints that ensure these surfaces possess an exact geostrophic streamfunction that is well-defined. 

#### Orthobaric surfaces

Topobaric surfaces are the topologically correct extension of isosurfaces of orthobaric density (de Szoeke et al. 2000) to possess geographical dependence of the functional relationship between pressure and in-situ density.  This software can compute orthobaric surfaces as a special case of topobaric surfaces, though it does not compute isosurfaces of the 3D orthobaric density defined by de Szoeke et al. (2000). 

### Geostrophic streamfunctions

Some approximately neutral surfaces possess an exact geostrophic streamfunction (GSF), while others do not.  Those that do include specific volume anomaly surfaces, orthobaric density surfaces, and modified topobaric surfaces.  (Truly neutral surfaces and topobaric surfaces both possess exact GSFs but they are ill-defined.) 

The exact GSF on a specific volume anomaly surface is the Montgomery (1937) potential.  Using this on other surfaces introduces an error in the geostrophic velocity estimation.  Zhang and Hogg (1992) noted this and derived a simply modification to the Montgomery (1937) potential to minimize this error.  McDougall and Klocker (2010) extend the Zhang and Hogg (1992) GSF for approximately neutral surfaces, while Cunningham (2000) present a GSF useful for inverse estimates.

The "orthobaric Montgomery potential" defined by Stanley (2019b) generalizes the Zhang and Hogg (1992) extension of the Montgomery (1937) potential, by approximating the pressure as a single-valued function of density on the surface.  This is a simple and highly accurate GSF. 

McDougall (1989) proved that an exact GSF does exist, at least locally, on a truly neutral surface.  Stanley (2019b) derived the form of this exact GSF, based on the multivalued functional relation between density and pressure on a neutral surface, but showed that this exact GSF is ill-defined (does not exist globally).  Nonetheless, Stanley (2019b) defined the "topobaric GSF" that is well-defined and approximates the exact GSF for neutral surfaces.  The topobaric GSF can be used on any approximately neutral surface with extremely good accuracy, and is exact on modified topobaric surfaces. 


## Contents:
- `./fex/                             `- software from the MATLAB file exchange
- `./lib/                             `- libraries (e.g. surface analysis, equations of state, geostrophic streamfunctions)
- `./run/                             `- example scripts
- `./omega-surface/                   `- create omega surfaces
- `./potential-density-surface/       `- create potential density surfaces
- `./specific-volume-anomaly-surface/ `- create specific volume anomaly surfaces
- `./topobaric-surface/               `- create topobaric surfaces, modified topobaric surfaces, and orthobaric surfaces
- `./LICENSE                          `- license
- `./ns_add_to_path.m                 `- function to add relevant subfolders to MATLAB's path
- `./ns_install.m                     `- function to install this package
- `./README.md                        `- this file

In addition, the following functions calculate geostrophic streamfunctions:
- `./topobaric-surface/topobaric_geostrf.m `- create topobaric GSF (Stanley 2019b)
- `./lib/gsf/zhanghogg92.m           		   `- create Zhang and Hogg (1992) GSF
- `./lib/gsf/cunningham00.m         		   `- create Cunningham (2000) GSF
- `./lib/gsf/mcdougallklocker10.m    		   `- create McDougall and Klocker (2010) GSF
- `./lib/gsf/orthobaric_montgomery.m 	     `- create orthobaric Montgomery potential (Stanley 2019b)


## Requirements:
MATLAB 2016b or higher (tested on 2017b, 2018b, and 2020a) with the Optimization Toolbox


## Installation:
Run the following command in MATLAB, replacing `~/work/neutral-surfaces/` with the path to this README.md file:
```
>> run('~/work/neutral-surfaces/ns_install.m')
```
If this does not succeed, see manual_install.md. 


## Usage:
The `./run` folder contains scripts to produce the results shown in various papers.

Note: The scripts for old papers are not necessarily updated as the codebase evolves; to guarantee reproducibility one must roll back the codebase to the time of the appropriate script.

To reproduce the figures of Stanley (2019a,b) exactly, use the v1.0 Topobaric-Surface code available at https://github.com/geoffstanley/Topobaric-Surface


## ReCon
Topobaric Surface uses ReCon to compute the Reeb Graph.

ReCon is available at
http://vgl.serc.iisc.ernet.in/software/software.php?pid=003

For more information on ReCon, see
Doraiswamy, H. & Natarajan, V. Computing Reeb Graphs as a Union of Contour Trees. IEEE Transactions on Visualization and Computer Graphics 19, 249–262 (2013).

The ReCon code included with Topobaric Surface has been modified to work with double precision and to allow the input and output to be from memory rather than files on the hard disk.  This portion of the Neutral Surfaces toolbox is licensed under the GNU Lesser General Public License (LGPL).

Note that ReCon requires Java 1.5 or higher, though it must be built for the version of the Java Runtime Environment that is packaged with MATLAB.


## References:
Cunningham, S.A., 2000. Circulation and volume flux of the North Atlantic using synoptic hydrographic data in a Bernoulli inverse. Journal of marine research 58, 1–35. https://doi.org/10.1357/002224000321511188

de Szoeke, R.A., Springer, S.R., Oxilia, D.M., 2000. Orthobaric density: A thermodynamic variable for ocean circulation studies. Journal of physical oceanography 30, 2830–2852.

Doraiswamy, H. & Natarajan, V. Computing Reeb Graphs as a Union of Contour Trees. IEEE Transactions on Visualization and Computer Graphics 19, 249–262 (2013).

Klocker, A., McDougall, T.J., Jackett, D.R., 2009. A new method for forming approximately neutral surfaces. Ocean Science 5, 155–172. https://doi.org/10.5194/os-5-155-2009

McDougall, T.J., 1989. Streamfunctions for the lateral velocity vector in a compressible ocean. Journal of Marine Research 47, 267–284. https://doi.org/10.1357/002224089785076271

McDougall, T.J., 1987. Neutral Surfaces. Journal of Physical Oceanography. https://doi.org/10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

McDougall, T.J., Jackett, D.R., 1988. On the helical nature of neutral trajectories in the ocean. Progress in Oceanography 20, 153–183. https://doi.org/10.1016/0079-6611(88)90001-8

Montgomery, R., 1937. A suggested method for representing gradient flow in isentropic surfaces. Bull. Amer. Meteor. Soc 18, 210–212.

Stanley, G.J., 2019a. Neutral surface topology. Ocean Modelling 138, 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

Stanley, G.J., 2019b. The exact geostrophic streamfunction for neutral surfaces. Ocean Modelling 138, 107–121. https://doi.org/10.1016/j.ocemod.2019.04.002

Stanley, G.J., McDougall, T.J., Barker, P.M., 2021. Algorithmic improvements to finding approximately neutral surfaces. J. Adv. Model. Earth Syst. submitted.

Wüst, G., 1935. The stratosphere of the Atlantic ocean. Scientific Results of the German Atlantic Expedition of the Research Vessel “Meteor” 1925–27 6.

Zhang, H.-M., Hogg, N.G., 1992. Circulation and water mass balance in the Brazil Basin. Journal of marine research 50, 385–420. https://doi.org/10.1357/002224092784797629


## Copyright:
MIT License

Copyright (c) 2021 Geoff Stanley

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  

Author(s) : Geoff Stanley

Email     : g.stanley@unsw.edu.au

Email     : geoffstanley@gmail.com
