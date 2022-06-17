#!/bin/bash
lj_param_chff=def_params/lj_param_chff.dat
bond_file=def_params/bond_data.dat
one_four_params=def_params/one_four_lj.dat
psf_file=def_params/ab_0.psf
crd_path_file=crd/crd_names.dat
ref_file=data/4conf_tip3.dat
opt_file=def_params/opt_species.dat
one_four_species=def_params/one_four_species.dat
nfiles=366
num_pep_atoms=34
n_onefour=5
ngenerations=200000
l_bound=0.1
u_bound=2.0
rmse_penalty=2.0
test_x=def_params/test_x.dat



gfortran -fopenmp -O2 differential_evolution.f90 graphf.f90 dijkstra.f90 utils_lj_fit.f90 lj_fit.f90 -o lj_fit.x
#./lj_fit.x  $lj_param_chff $one_four_params $one_four_species $psf_file $bond_file $crd_path_file $opt_file $ref_file $nfiles $num_pep_atoms $n_onefour $ngenerations $l_bound $u_bound $rmse_penalty $test_x


