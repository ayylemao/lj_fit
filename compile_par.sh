#!/bin/bash
lj_param_chff=def_params/lj_param_tip4.dat
bond_file=def_params/bond_data.dat
one_four_params=def_params/one_four_tip4.dat
psf_file=def_params/ab_0.psf
crd_path_file=crd/crd_names_no_ra.dat
ref_file=data/4conf_iamo/4conf_iamo_comb.dat
opt_file=def_params/iamo_species.dat
one_four_species=def_params/one_four_species.dat
nfiles=341
num_pep_atoms=34
n_onefour=5
ngenerations=100
l_bound=0.1
u_bound=2.0
rmse_penalty=2.0


gfortran -fopenmp -O2 differential_evolution.f90 graphf.f90 dijkstra.f90 utils_lj_fit.f90 lj_fit.f90 -o lj_fit.x
./lj_fit.x  $lj_param_chff $one_four_params $one_four_species $psf_file $bond_file $crd_path_file $opt_file $ref_file $nfiles $num_pep_atoms $n_onefour $ngenerations $l_bound $u_bound $rmse_penalty


