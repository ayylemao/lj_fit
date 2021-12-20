
program lj_fit
use differential_evolution
use utils
implicit none
character(len=100) :: base_name, lj_param_file, onefour_file, opt_file, junk
character(len=100) :: psf_file, bond_file, crd_dir, onefour_species_file
character(len=100) :: ref_name
real*8, dimension(:,:), allocatable :: search_range
real*8, dimension(36) :: init_val_search, x
character(len=3), dimension(36) :: print_helper
real*8 :: energy, rmse
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
call init_params(num_files          =      198,&
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
print_helper(1:13) = opt_species
print_helper(14) = "CT1"
print_helper(15) = "CT3"
print_helper(16) = "NH1"
print_helper(17) = "O  "
print_helper(18) = "OB "


call init_search_range(search_range)


call DE_init(set_range               = search_range,     &
             set_popSize             = 100,              &
             set_maxGens             = 5000,               &
             set_maxChilds           = 5,                &
             set_forceRange          = .false.,         &
             set_mutationStrategy    = DErand1,  &
             set_crossProb           = 0.9d0,             &
             set_verbose             = verbose,          &
             set_Nprint              = 10)

call DE_optimize(opt_func, feasible, sumconstr, x, init_pop=init_pop)
write(*,*) "BEST SOLUTION:"
do i = 1, nopt+5
    write(*,*) x(2*i-1), x(2*i)
end do
write(*,*) "CALC ENER", "FIT ENERGY"
do i = 1,nfiles 
    call get_lj_energy(i, energy, x)
    write(*,*) energy, ref_energies(i), sqrt((energy-ref_energies(i))**2)
end do
write(*,*) "FINAL RMSE", opt_func(x)
write(*,*) "FINAL LJ PARAMS:"
do i = 1, nopt+5
    write(*,'(A4,F10.5,2F10.5)') print_helper(i), 0.0, x(2*i-1), x(2*i)
end do

! =============== END MAIN ================================================

contains
    
real*8 function opt_func(y)
    real*8, dimension(:) :: y
    real*8, dimension(nfiles) :: energies
    real*8 :: curr_energy, ref_val, rmse
    integer i, j
    real*8 :: t1, t2
    energies = 0
    call calc_full_lj(energies, y) 
    rmse = 0.0d0   
    do i = 1,nfiles 
        rmse = rmse + (ref_energies(i) - energies(i))**2
    end do 
    opt_func = sqrt(rmse) 
end function opt_func 

subroutine calc_full_lj(energies, y)
    implicit none
    real*8, dimension(:) :: y
    real*8, dimension(:) :: energies
    integer :: i
    do i = 1, nfiles
        call get_lj_energy(i, energies(i), y)
    end do
end subroutine calc_full_lj

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
    call load_opt_spec(opt_file)
    call get_excl_array(bond_file)
    
end subroutine init_system

subroutine init_search_range(search_range)
    implicit none
    real*8, allocatable, dimension(:,:) :: search_range
    real*8 :: eps_min, eps_max, rmin_min, rmin_max
    integer i,j
    if (.not. allocated(search_range)) allocate(search_range(2, 2*(nopt+5)))
    do i = 1, nopt+5
        !eps parameters search range
        search_range(1, 2*i-1) = 1.5*init_val_search(2*i-1)
        search_range(2, 2*i-1) = 0.5*init_val_search(2*i-1)
        !rmin parameters search range
        search_range(1, 2*i) = 0.5*init_val_search(2*i)
        search_range(2, 2*i) = 1.5*init_val_search(2*i)
    end do
end subroutine init_search_range


logical function feasible(y)
    implicit none    
    real*8, dimension(:) :: y 
    integer :: i,j
     
    do i = 1, 2*(nopt+5)
        if (abs(y(i)) .ge. 2.0*abs(init_val_search(i))) then
            feasible = .false.
            return
        else if (abs(y(i)) .le. 0.25*(init_val_search(i))) then
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
    sumconstr = 0.0d0 
end function sumconstr

end program lj_fit

    


