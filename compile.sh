#!/bin/bash
lj_param_chff=def_params/lj_param_chff.dat
bond_file=def_params/bond_data.dat
one_four_params=def_params/one_four_lj.dat
psf_file=def_params/ab_0.psf
crd_path_file=crd/crd_names.dat
ref_file=data/dft_ref_energies_avg.dat
opt_file=def_params/cons_species.dat
one_four_species=def_params/one_four_species.dat
nfiles=198
num_pep_atoms=34
n_onefour=1
ngenerations=20000

gfortran -g -ffpe-trap=zero,invalid,overflow,underflow differential_evolution.f90 graphf.f90 dijkstra.f90 utils_lj_fit.f90 lj_fit.f90 -o lj_fit.x
./lj_fit.x  $lj_param_chff $one_four_params $one_four_species $psf_file $bond_file $crd_path_file $opt_file $ref_file $nfiles $num_pep_atoms $n_onefour $ngenerations
rm lj_fit.x
