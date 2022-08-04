[![DOI](https://zenodo.org/badge/397937474.svg)](https://zenodo.org/badge/latestdoi/397937474)

# Osmotic transport at the aqueous graphene and hBN interfaces:  scaling laws from a unified, first principles description.

Laurent Joly, Rober Meissner, Marcella Iannuzzi, Gabriele Tocci

## Abstract

Osmotic transport in nanoconfined aqueous electrolytes provides new venues for water desalination and blue energy harvesting; the osmotic response of nanofluidic systems is controlled by the interfacial structure of water and electrolyte solutions in the so-called electrical double layer (EDL), but a molecular-level picture of the EDL is to a large extent still lacking.  Particularly, the role of the electronic  structure has not been considered in the description of electrolyte/surface interactions. Here, we report enhanced sampling simulations
based on ab initio molecular dynamics,  aiming at unravelling the free energy of prototypical  ions adsorbed at the aqueous graphene and hBN interfaces, and  its consequences on nanofluidic osmotic transport.  Specifically, we predicted the zeta potential, the diffusio-osmotic mobility and the diffusio-osmotic conductivity for a wide range of salt concentrations from the ab initio water and ion spatial distributions through an analytical
framework based on Stokes equation and a modified Poisson-Boltzmann equation. We observed concentration-dependent scaling laws,  together with  dramatic differences in osmotic transport  between the two interfaces, including diffusio-osmotic flow and current reversal on hBN, but not on graphene. We could rationalize the results for the three osmotic responses with a simple model based on characteristic length scales for ion and water adsorption at the surface, which are quite different on graphene and on hBN. Our work provides first principles insights into the structure and osmotic transport of aqueous electrolytes on two-dimensional materials and explores new pathways for efficient water desalination and osmotic energy conversion.

### Description

We have included the notebooks, python scripts and CP2K input files to reproduce the main results of the following paper:

* L. Joly, R. H. Meißner, M. Iannuzzi, G. Tocci, *Osmotic Transport at the Aqueous Graphene and hBN Interfaces: Scaling Laws from a Unified, First-Principles Description* ACS Nano **15**, 9 15249-15258 (2021); DOI: [10.1021/acsnano.1c05931](https://doi.org/10.1021/acsnano.1c05931)

### Dependencies

* bash, python3, Wolfram Mathematica 12, numpy, scipy, pandas, jupyter

### Executing programs

* In the directory ```calculate_free_energies``` execute ```./extract_free_energies.sh``` to get the free energy profiles from the collective variable (CV) files for each umbrella sampling simulation;
* Run the mathematica notebook ```solv_osm_mpb.nb``` to solve the mPB equation with the finite elements method (FEM) and to calculate the integrals of the osmotic transport coefficients numerically taking as input the water density profiles and the ions' free energy profiles;
* Run the jupyter notebook ```osmotic_transport.ipynb``` to produce figures 2 and 3 in the main text by taking as input the data that has been pre-calculated from the previous two steps

