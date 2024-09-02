Table of contents
=================

* [Model description](#model-description)
* [References](#references)
* [Model parameters](#model-parameters)
* [Required files](#required-files)
* [Installation and usage in XSPEC](#installation-and-usage-in-xspec)


Model description
=================

This model computes the emission and its polarisation properties from
an X-ray source of power-law emission of arbitrary incident polarisation 
that is reprocessed in axially symmetric structures. The reprocessing
is precomputed in a form of FITS fits files that are required for
usage. These files were computed for three different states
of incident polarisation and this model interpolates for any primary
polarisation state. In addition the model possesses other properties of the
reprocessing tables that were precomputed for nearly neutral disc reflection.

In this case of xsstokes_disc, it is reflection from distant regions
of a geometrically thin optically thick accretion disc residing in the
equatorial plane around a central black hole. The source of emission
is a hot X-ray corona illuminating isotropically between
0 <= cos_incl_incident <= M_i, to represent distant disc reflection.
The local disc reflection was computed using the STOKES code and more
details are given in Podgorný J. et al. (2022). More details on the
integration and implementation for arbitrary polarisation state in XSPEC
are given in Podgorný J. et al. (2024). Apart from the incident
polarisation degree and angle, the model then also depends on
observer's inclination angle, the primary power-law index Gamma,
the integration upper limit M (i.e. representing the corona size),
the position angle on the sky, and the overall Dopler shift.
No relativistic effects inside the system are taken into account.
All components are static.

For any issues regarding the use of xsstokes_disc, please, contact J. Podgorný at 
[jakub.podgorny@asu.cas.cz](mailto:jakub.podgorny@asu.cas.cz) or M. Dovčiak
[michal.dovciak@asu.cas.cz](mailto:michal.dovciak@asu.cas.cz).


References
==========

Podgorný J, Dovčiak M, Marin F, Goosmann RW and Różańska A (2022)
_Spectral and polarization properties of reflected X-ray emission from black hole accretion discs_
[MNRAS, 510, pp.4723-4735](https://doi.org/10.1093/mnras/stab3714)
[[arXiv:2201.07494](https://arxiv.org/abs/2201.07494)]

Podgorný J, Dovčiak M, Marin F (2024)
_Simple numerical X-ray polarization models of reflecting axially symmetric structures around accreting compact objects_  
[MNRAS, 530, pp.2608-2626](https://doi.org/10.1093/mnras/stae1009)
[[arXiv:2310.15647](https://arxiv.org/abs/2310.15647)]

Model parameters
================

* **par1 ... Size**
  - the upper limit M_i in cos_incl_incident integration,
  - representing the corona size
* **par2 ... PhoIndex** 
  - photon index of the primary power-law X-ray flux
* **par3 ... cos_incl** 
  - cosine of the observer inclination (1.-pole, 0.-disc)
* **par4 ... poldeg**  
  - intrinsic polarisation degree of primary radiation between 0. and 1.
* **par5 ... chi**  
  - intrinsic polarisation angle (in degrees, -90 < chi < 90) of primary radiation
  - the orientation is degenarate by 180 degrees
* **par6 ... pos_ang**  
  - orientation of the system (-90 < pos_ang < 90), 
  - the position angle (in degrees) of the system 
  - rotation axis with direction up,
  - the orientation is degenarate by 180 degrees
* **par7 ... zshift**
  - overall Doppler shift
* **par8 ... Stokes**
  - defines the output of the model:
    - -1: the output is defined according to the XFLT0001 keyword of the 
          SPECTRUM extension of the data file, where "Stokes:0" means photon 
          number density flux, "Stokes:1" means Stokes parameter Q devided by 
          energy and "Stokes:2" means Stokes parameter U devided by energy
    -  0: array of photon number density flux per bin (array of Stokes parameter 
          I devided by energy) with the polarisation computations switched off
    -  1: array of photon number density flux per bin (array of Stokes parameter 
          I devided by energy), with the polarisation computations switched on
    -  2: array of Stokes parameter Q devided by energy
    -  3: array of Stokes parameter U devided by energy
    -  4: array of Stokes parameter V devided by energy
    -  5: array of degree of polarisation
    -  6: array of polarisation angle &psi;=0.5*atan(U/Q)
    -  7: array of "Stokes" angle 
          &beta; = 0.5*asin(V/sqrt(Q<sup>2</sup>+U<sup>2</sup>+V<sup>2</sup>))
    -  8: array of Stokes parameter Q devided by I
    -  9: array of Stokes parameter U devided by I
    - 10: array of Stokes parameter V devided by I


Further output of the model
===========================

When the fit is finished the XSPEC command xset can be used to see the following 
additional information of the model:
 
* **inc_degrees**
  - inclination in degrees measured from the disc normal, i.e. "acos(cos_incl)/PI*180."


Required files
==============

* **Source code files**
  - xsstokes_disc.c
  - lmodel-stokesni.dat  
* **reprocessing tables**
  (for reprocessing of photons in distant discs, more details in Podgorný J. et al. (2024))
  - [23576508.zip](https://doi.org/10.6084/m9.figshare.23576508)  


Installation and usage in XSPEC
===============================

1. **Download the source code files**
   into a directory where you want to install the model, e.g. '/path/to/xsstokes_disc-master/'

2. **Download the FITS files** 
   [23576508.zip](https://doi.org/10.6084/m9.figshare.23576508) 
   with the tables into the directory with xsstokes_disc, i.e. '/path/to/xsstokes_disc-master'.

3. **Unzip the reflection tables** with polarisation information, e.g. by the command:

   `unzip 23576508.zip`
   
4. **Compile the code** inside XSPEC (needed to be done only once):

   The code is compiled inside XSPEC with the following command (assuming all 
   the source files and FITS tables are in the directory 
   '/path/to/xsstokes_disc-master'):

   `initpackage stokesnidisc lmodel-stokesni.dat /path/to/xsstokes_disc-master`

   **Note**:
   Your XSPEC installation must have been originally installed from the source 
   code distribution. Local models, like xsstokes_disc, cannot be used if the XSPEC 
   was originally installed from the pre-compiled binary distribution.

5. **Load the xsstokes_disc model** into XSPEC:
   
   To use the xsstokes_disc model inside XSPEC, first the model package needs to be 
   loaded and you can also setup a directory containing the xsstokes_disc set:

   `lmod stokesnidisc /path/to/xsstokes_disc-master`
   `xset XSDIR /path/to/xsstokes_disc-master`

6. Then the **xsstokes_disc model may be used** in the usual way, e.g.:

   `mo stokesdisc`

   **Note**:
   In case of segmentation fault, one may need to increase the stack size before 
   launching XSPEC, e.g. with the command:

   `ulimit -s unlimited`

   or 

   `ulimit -s 65532`
