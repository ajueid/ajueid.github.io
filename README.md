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

## Citations

If you use these Tables please cite the following references:

- [S. Amoroso, S. Caron, A. Jueid, R. Ruiz de Austri, P. Skands, JCAP 05 (2019) 007](https://arxiv.org/abs/1812.07424)
- [A. Jueid, J. Kip, R. Ruiz de Austri, P. Skands, e-Print: 2202.11546](https://arxiv.org/abs/2202.11546)

And optionally
- [S. Mrenna, P. Skands, Phys.Rev.D 94 (2016) 7, 074005](https://arxiv.org/abs/1605.08352)
- [S. Amoroso, S. Caron, A. Jueid, R. Ruiz de Austri, P. Skands,  PoS TOOLS2020 (2021) 028](https://arxiv.org/abs/2012.08901)
- [A. Jueid, J. Phys.: Conf. Ser. 2156 012057](https://arxiv.org/abs/2110.09747)

## Authors
- [Adil Jueid](adil.hep@gmail.com)
- [Jochem Kip](jochem.kip@ru.nl)
- [Roberto Ruiz de Austri](rruiz@ific.uv.es)
- [Peter Skands](peter.skands@monash.edu)
