
program lj_fit
use differential_evolution
use utils
implicit none
character(len=100) :: base_name, lj_param_file, onefour_file, opt_file, junk
character(len=100) :: psf_file, bond_file, crd_dir, onefour_species_file
character(len=100) :: ref_name
character(len=100), allocatable, dimension(:) :: crd_names
real*8, dimension(:,:), allocatable :: search_range
real*8, dimension(36) :: init_val_search, x
real*8 :: energy
integer :: i, j
logical :: verbose = .true.








! ============ START MAIN =================================================


! File name declarations
base_name = "data/test"
lj_param_file = "def_params/lj_param_chff.dat"
bond_file = "def_params/bond_data.dat"
onefour_file = "def_params/one_four_lj.dat"
psf_file = "def_params/ab_0.psf"
crd_dir = "crd/"
ref_name = "data/dft_ref_energies.dat"
opt_file = "def_params/opt_species.dat"
onefour_species_file = "def_params/one_four_species.dat"
! initialize parameters and other things
call init_params(num_files          =        198,&
                 cut_on             =   10.0d0,&
                 cut_off            =   12.0d0,&
                 num_pep_atoms      =   34)

call init_system(psf_file, lj_param_file, onefour_file, bond_file,&
opt_file, onefour_species_file)
call load_data(crd_dir, ref_name)


init_val_search=0
open(69, file="def_params/test_x.dat", status = 'old')

do i = 1, nopt
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i) 
end do
do i = nopt+1, nopt + 5 
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i)
end do
close(69) 

call init_search_range(search_range)


call DE_init(set_range               = search_range,     &
             set_popSize             = 100,              &
             set_maxGens             = 50,               &
             set_maxChilds           = 1,                &
             set_forceRange          = .true.,         &
             set_mutationStrategy    = DErand1,  &
             set_crossProb           = 1.d0,             &
             set_verbose             = verbose,          &
             set_Nprint              = 2)
           

call DE_optimize(opt_func, feasible, sumconstr, x, init_pop=init_pop)
do i = 1, nopt
    write(*,*) x(i)
end do
 

!write(*,*) opt_func(init_val_search) 
! =============== END MAIN ================================================

contains
    
real*8 function opt_func(y)
    real*8, dimension(:) :: y
    real*8, dimension(nfiles) :: energies
    real*8 :: curr_energy, min_val, rmse
    integer i
    do i = 1, nfiles       
    call get_lj_energy(crd_data(:,:,i), curr_energy, y) 
        energies(i) = curr_energy 
    !    write(*,*) curr_energy
    end do
    min_val = MINVAL(energies)
    energies = energies - min_val
    rmse = 0.0d0
    do i = 1, nfiles
        rmse = rmse + (ref_energies(i) - energies(i))**2
    end do
    write(*,*) "RMSE:", sqrt(rmse), sqrt(rmse/real(nfiles))  
    opt_func = sqrt(rmse) 
end function opt_func 

subroutine init_pop(pop)
    real*8, dimension(:,:), intent(out) :: pop
    real*8 :: ran
    integer :: popSize, i, j

    popSize = size(pop, dim=2)
    do i = 1, popSize
        do j = 1, 2*(nopt+5)
            call random_number(ran)
            pop(j,i) = (1.2*init_val_search(j)-0.8*init_val_search(j))*ran
            pop(j,i) = pop(j,i) + init_val_search(j)*0.8
        end do
    end do
end subroutine init_pop


subroutine init_system(psf_file, lj_param_file, onefour_file, bond_file,&
opt_file, onefour_species_file)
    character(len=100) :: psf_file, lj_param_file, onefour_file, bond_file
    character(len=100) :: opt_file, onefour_species_file
    call load_psf(psf_file)
    call load_lj_params(lj_param_file)
    call load_one_four_params(onefour_file, onefour_species_file)
    call get_excl_array(bond_file)
    call load_opt_spec(opt_file)
end subroutine init_system

subroutine init_search_range(search_range)
    implicit none
    real*8, allocatable, dimension(:,:) :: search_range
    real*8 :: eps_min, eps_max, rmin_min, rmin_max
    if (.not. allocated(search_range)) allocate(search_range(2, 2*(nopt+5)))
    do i = 1, 2*(nopt+5)
        search_range(1, i) = 1.5*init_val_search(i)
        search_range(2, i) = 0.5*init_val_search(i)
    end do
end subroutine init_search_range


logical function feasible(y)
    implicit none    
    real*8, dimension(:) :: y 
    do i = 1, 2*(nopt+5)
        if (y(i) .ge. 1.5*init_val_search(i)) then
            feasible = .false.
            return
        else if (y(i) .le. 0.5*init_val_search(i)) then
            feasible = .false.
            return
        else
            feasible = .true.
        end if
    end do
end function feasible

real*8 function sumconstr(y)
    implicit none   
    real*8, dimension(:) :: y
    sumconstr = 0.0
end function sumconstr

end program lj_fit

    


