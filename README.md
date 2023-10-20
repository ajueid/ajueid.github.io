# Toolbox for QCD uncertainties on particle spectra from dark-matter annihilation

## Introduction

We provide the spectra of stable particles in dark-matter (DM) annihilation or decay in a tabulated form using PYTHIA version 8306. In addition to the central prediction, we estimate the QCD uncertainties both due to hadronization as well as to the parton-shower variations. The DM masses considered here vary between 5 GeV and 100 TeV. We consider 14 primary annihilation channels:

```console
DM DM -> e+e-, mu+ mu-, tau tau, uu, dd, ss, cc, bb, tt, WW, ZZ, gg, and hh     (1)
```

The data is provided in such a way to faciliate including the uncertainties in DM fits. In other words, for each file, we provide the central prediction for the flux of the particle specie in addition to +- uncertainties on the flux (three independent uncertainties). We also provide the Tables that contains only the central values of the fluxes and not the uncertainties. The Tables are provided for six final states: photons, positrons, electron antineutrinos, muon antineutrinos, tau antineutrinos, and antiprotons.

## Structure of the Tables

For each stable final state, there is one dedicated file in ascii format. The tables are organised into two folders: 
```console
$ data/woUncertainty
$ data/wUncertainty
```

#### 1) data/woUncertainty

In this folder, we provide the spectra without QCD uncertainties, i.e. just the nominal predictions. Each file contains sixteen columns: The DM mass in GeV, the fraction x -- defined as the kinetic energy divided by the DM mass. The rest of the columns are organised as 

```console
dN/dx [ee]  dN/dx [mumu]   dN/dx [tautau]    dN/dx [uu]    dN/dx [dd]    dN/dx [ss]     dN/dx [cc]      dN/dx [bb]     dN/dx [tt]     dN/dx [gaga]     dN/dx [zz]         dN/dx [ww]      dN/dx [gg]      dN/dx [hh]     
```

#### 2) data/wUncertainty

In this folder, we provide the particle spectra along with the QCD uncertainties on the central value. Each file contains one hundred columns: The DM mass in GeV, the fraction x -- defined as the kinetic energy of the particle divided by the DM mass --. The rest of the columns are organised as 

```console
dN/dx [XX]    +DHad [XX]   -DHad [XX]   +DScale [XX]   -DScale [XX]    +DcNS [XX]    -DcNS [XX]  

dN/dx [XX]: is the flux at the best-fit point of the PYTHIA8 parameters.
+-DHad [XX]: are the uncertainties on the flux from the variations of the hadronisation model parameters.
+-DScale [XX]: are the shower uncertainties on the flux from the variation of the shower evolution variable by a factor of 2.
+-DcNS [XX]: are the uncertainties from the variations of the non-singular terms of the DGLAP splitting kernels. 
```

In the two folders, the files are named as:
```console
AtProduction-FS.dat; FS = Ga, Positrons, Nuel, Numu, Nuta, AntiP.
```

## Interpolation

The file 'Interpolate.py' containts five functions which can provide either an interpolated annihilation spectrum, a simulated diffusion spectrum if the final-state particles are antiprotons, or one of three types of DM uncertainties. The functions can be called by simply importing the 'Interpolate.py' file.

## Antiproton propagation

In order to compute the propagation of cosmic antiprotons originating from DM annihilation we use the DRAGON 2 code. We use the following diffusion parameters:
D_0 = 4.2
delta = 0.45
v_A = 11.8
eta = -1.9
alpha = 2.29

All other parameters have been left at their default values. The pre-computed diffusion tables can be provided by the user for all functions that utilize this functionality if different diffusion parameters are desired.

## Annihilation and Diffused spectra

The annihilation spectra including QCD uncertainties can be obtained by calling the 'annihilation_spectrum' function. The input parameters are the DM mass, annihilation channel, and final state particle, and optionally the path to the computed annihilation spectra if the location file structure has been altered. 
The diffued spectra for antiprotons can be obtained with the function 'diffused spectrum'. This function takes the same arguments as 'annihilation_spectrum' with an additional optional argument of the path to the pre-computed diffusion table. All interpolation is done linearly.

## DM uncertainties

To compute DM uncertainties there are three available functions: 'sigmav_uncertainty', 'mass_uncertainty', and 'spectrum uncertainty'. These functions take the DM mass, annihilation channel, and final state particles as mandatory input parameters. The path to the computed annihilation spectra can be manually provided, as can the path to the pre-computed diffusion tables. Finally, the minimum and maximum energy of the spectra considered in the fit can be provided.

## Sampling code

The Sampler.py code allows for the easy sampling of non-standard neutral particle spectra from Dark Matter annihilation or decay. The code can be run with Run.py; either directly, which requires manual input, or by passing input.txt as an argument. The format for an input file is exactly the same as for the manual input.  
  
The required input for the code is:

• How the spectrum is produced: DM decay (1) or DM annihilation (2)  
• The particle type: νe, νμ, ντ , or γ.  
• Which process: DMDM → ν/γX (1) or DMDM → XX (2).  
• The DM mass MDM in GeV.  
• The X mass MX in GeV.  
• The decay modes of X with the format being ‘BR daughter1 daughter2’. The total branching ratio must sum to 1.  
• The number of samples from the spectrum.  
• The path to the csv file where the points are saved. When no input is given, no save is made. The keyword ‘plot’ or ’logplot’ can be entered to plot the sampled data linearly or logarithmically. 
  
## Citations

If you use these Tables please cite the following references:

- [S. Amoroso, S. Caron, A. Jueid, R. Ruiz de Austri, P. Skands, JCAP 05 (2019) 007](https://arxiv.org/abs/1812.07424)
- [A. Jueid, J. Kip, R. Ruiz de Austri, P. Skands, JCAP 04 (2023) 068](https://arxiv.org/abs/2202.11546)
- [A. Jueid, J. Kip, R. Ruiz de Austri, P. Skands, e-Print:2303.11363](https://arxiv.org/2302.11363)

If you use the sampler code please cite:

- [W. Beenakker, S. Caron, J. Kip, R. Ruiz de Austri, Zhongyi Zhang, e-Print: 2306.16523](https://arxiv.org/abs/2306.16523)

And optionally
- [S. Mrenna, P. Skands, Phys.Rev.D 94 (2016) 7, 074005](https://arxiv.org/abs/1605.08352)
- [S. Amoroso, S. Caron, A. Jueid, R. Ruiz de Austri, P. Skands,  PoS TOOLS2020 (2021) 028](https://arxiv.org/abs/2012.08901)
- [A. Jueid, J. Phys.: Conf. Ser. 2156 012057](https://arxiv.org/abs/2110.09747)

## Authors
- [Adil Jueid](adil.hep@gmail.com)
- [Jochem Kip](jochem.kip@ru.nl)
- [Roberto Ruiz de Austri](rruiz@ific.uv.es)
- [Peter Skands](peter.skands@monash.edu)
