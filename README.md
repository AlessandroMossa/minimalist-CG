# minimalist-CG
From the output of an atomistic simulation to a minimalist coarse-grained model ready to be simulated in LAMMPS

The minimalist coarse-grained (CG) model for protein dynamics that I have in mind is the one described in F. Delfino, Y. Porozov, E. Stepanov, G. Tamazian and V. Tozzini, "Evolutionary switches structural transitions via coarse-grained models" J Comput Biol 27(2), 189-199 (2020) doi: [10.1089/cmb.2019.0338](https://dx.doi.org/10.1089/cmb.2019.0338) and in "Structural transition states explored with minimalist coarse grained models: Applications to calmodulin" Front Mol Biosci 6, 104 (2019) doi: [10.3389/fmolb.2019.00104](https://dx.doi.org/10.3389/fmolb.2019.00104) by the same authors. The [Julia](https://julialang.org) script in this page takes as an input a representative structure in the file md_0.gro and builds a data file ready to be used for a CG molecular dynamics with [LAMMPS](https://www.lammps.org/).

## Minimalist Tutorial

Let's say we want to build a minimalist CG model of a "regular" protein (i.e., not significantly disordered, no exotic atoms). Let's extract from an atomistic simulation a representative snapshot (here is the file `md_0.gro`, which contains a frame from an atomistic simulation of the protein ACE2 ([PDB:1R42](https://www.rcsb.org/structure/1R42))). If we have Julia already installed (if not, follow the instructions [here](https://julialang.org/install/)), that takes no time at all: turn on the REPL
