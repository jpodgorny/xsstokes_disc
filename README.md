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

   `initpackage stokesdisc lmodel-stokesdisc.dat /path/to/xsstokes_disc-master`

   **Note**:
   Your XSPEC installation must have been originally installed from the source 
   code distribution. Local models, like xsstokes_disc, cannot be used if the XSPEC 
   was originally installed from the pre-compiled binary distribution.

5. **Load the xsstokes_disc model** into XSPEC:
   
   To use the xsstokes_disc model inside XSPEC, first the model package needs to be 
   loaded and you can also setup a directory containing the xsstokes_disc set:

   `lmod stokesdisc /path/to/xsstokes_disc-master`

   `xset XSDIR /path/to/xsstokes_disc-master`

6. Then the **xsstokes_disc model may be used** in the usual way, e.g.:

   `mo stokesdisc`

   **Note**:
   In case of segmentation fault, one may need to increase the stack size before 
   launching XSPEC, e.g. with the command:

   `ulimit -s unlimited`

   or 

   `ulimit -s 65532`


Viewing the model in XSPEC
==========================

One usually needs to have polarisation data sets loaded for all three stokes
parameters to view the polarisation properties predicted by a model inside XSPEC.
The model is then shown in the energy range covered by, and with the energy binning
of, these data sets. To overcome these disadvantages, the xsstokes_disc model provides
a parameter (par8) that defines the output of the model, e.g. to show the model
prediction for the polarisation degree, one needs to set par8 to 5, see Section
[Model parameters](#model-parameters). In this case, dummy response may be used
without any data loaded, e.g. `'dummyrsp 1. 100. 200 log'`. Then the model may be
viewed (after it has been loaded in XSPEC) in a usual way using `'plot model'`. Note
that when showing polarisation degree, angle and normalised Stokes parameters, i.e.
when par8 is set to 5-10, the normalisation of the model needs to be set to unity!

To use the original way of showing results of polarisation model inside XSPEC, we
also provide a fake null data sets for all three Stokes parameters, i, q and u,
binned in 300 channels, together with corresponding unit response matrices, rmf, arf
and mrf, defined in 0.1 to 100 keV in 300 channels, see
[fake_null_iqu_300ch.tar.gz](https://owncloud.asu.cas.cz/index.php/s/Flk6cwYLISmw0D5).
One can see the model predicted polarisation properties in the following way:

1. **Download the package containing the fake null data and unit repsonses** -
   [fake_null_iqu_300ch.tar.gz](https://owncloud.asu.cas.cz/index.php/s/Flk6cwYLISmw0D5)

2. **Uncompress the package**, e.g. by the command:

   `tar -xzf fake_null_iqu_300ch.tar.gz`

   This will uncompress the following files:

   - `fake_null_i_300ch.pha`
   - `fake_null_q_300ch.pha`
   - `fake_null_u_300ch.pha`
   - `fake_unit_response_0.1-100keV_300ch.rmf`
   - `fake_unit_response_0.1-100keV_300ch.arf`
   - `fake_unit_response_0.1-100keV_300ch.mrf`

3. **Load the fake null data into XSPEC**:

   `data 1 fake_null_i_300ch.pha`
   `data 2 fake_null_q_300ch.pha`
   `data 3 fake_null_u_300ch.pha`

   the fake unit responses will be automatically loaded. For convenience, change also the energy range of the data

   `ignore *:**-1.`

4. **Set xsstokes_disc par8 to -1**

   To be able to use the traditional way of using the polarisation data with the model in XSPEC, one needs to set the xsstokes_disc parameter par8 to -1:

   `newpar 8 -1 -1 -1 -1 10 10`

   Here, we have also redefined the boundaries of this parameter in case -1 value is not allowed by the default settings.

5. **View the polarisation degree and angle** in XSPEC:

   `plot polfrac`
   `plot polangle`

   Note that the model has to be loaded into XSPEC first, as described in Section
   [Installation and usage in XSPEC](#installation-and-usage-in-xspec).
