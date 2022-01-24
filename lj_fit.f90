
program lj_fit
use differential_evolution
use utils
implicit none
character(len=100) :: base_name, lj_param_file, onefour_file, opt_file, junk
character(len=100) :: psf_file, bond_file, crd_dir, onefour_species_file
character(len=100) :: ref_name, crd_name_file
real*8, dimension(:,:), allocatable :: search_range
real*8, dimension(:), allocatable :: init_val_search, x
character(len=3), dimension(:), allocatable :: print_helper
real*8 :: energy, rmse
integer :: i, j
logical :: verbose = .true.

! ============ START MAIN =================================================


! File name declarations
base_name = "data/test"
lj_param_file = "wat_dim/lj_param_chff.dat"
bond_file = "wat_dim/bond_data.dat"
onefour_file = "wat_dim/one_four_lj.dat"
psf_file = "wat_dim/wat_dim_0.psf"
crd_name_file = "crd_wat_dim/crd_names.dat"
crd_dir = "crd_wat_dim/"
ref_name = "wat_dim_data/wat_dim_ref.dat"
opt_file = "wat_dim/opt_species.dat"
onefour_species_file = "wat_dim/one_four_species.dat"

! initialize parameters and other things
call init_params(num_files          =      50,&
                 cut_on             =   10.0d0,&
                 cut_off            =   12.0d0,&
                 num_pep_atoms      =   6)

call init_system(psf_file, lj_param_file, onefour_file, bond_file,&
opt_file, onefour_species_file)
call load_data(crd_name_file, crd_dir, ref_name)

if (.not. allocated(init_val_search)) allocate(init_val_search(2*(nopt+nonefour)))
if (.not. allocated(x)) allocate(x(2*(nopt+nonefour)))
if (.not. allocated(print_helper)) allocate(print_helper(nopt+nonefour))

init_val_search=0
open(69, file="wat_dim/test_x.dat", status = 'old')

do i = 1, nopt
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i) 
end do
do i = nopt+1, nopt + nonefour 
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i)
end do
close(69) 



print_helper(1:2) = opt_species
!print_helper(14) = "CT1"
!print_helper(15) = "CT3"
!print_helper(16) = "NH1"
!print_helper(17) = "O  "
!print_helper(18) = "OB "
call calc_look_ups()






call init_search_range(search_range)


call DE_init(set_range               = search_range,     &
             set_popSize             = 100,              &
             set_maxGens             = 10000,               &
             set_maxChilds           = 1,                &
             set_forceRange          = .false.,         &
             set_mutationStrategy    = DErand1,  &
             set_crossProb           = 0.9d0,             &
             set_verbose             = verbose,          &
             set_Nprint              = 1)

call DE_optimize(opt_func, feasible, sumconstr, x, init_pop=init_pop) 

!call DE_optimize(opt_func, feasible, sumconstr, x, guess=init_val_search) 
write(*,*) "BEST SOLUTION:"
do i = 1, nopt+nonefour
    write(*,*) x(2*i-1), x(2*i)
end do
write(*,*) "CALC ENER", "FIT ENERGY"
do i = 1,nfiles 
    call get_lj_energy(i, energy, x)
    write(*,*) energy, ref_energies(i), sqrt((energy-ref_energies(i))**2)
end do
write(*,*) "FINAL RMSE", opt_func(x)/SQRT(float(nfiles))
write(*,*) "FINAL LJ PARAMS:"
do i = 1, nopt+nonefour
    write(*,'(A4,F10.5,2F10.5)') print_helper(i), 0.0, -abs(x(2*i-1)), abs(x(2*i))
end do
write(*,*)
do i = 1, nopt
    write(*,'(A14,2F10.5)') "   ", -abs(init_val_search(2*i-1)), abs(init_val_search(2*i))
end do

write(*,*)
!do i = 1, nopt
!    write(*,'(A14,2F10.5)') "   ", -abs(init_val_search(2*i-1))/-abs(x(2*i-1)), &
!                                   abs(init_val_search(2*i))/abs(x(2*i))
!end do


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
    ref_val = sum(energies)/size(energies)
!    energies = energies - ref_val

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
        do j = 1, 2*(nopt+nonefour)
            call random_number(ran)
            pop(j,i) = (10*init_val_search(j)-0.01*init_val_search(j))*ran
            pop(j,i) = pop(j,i) + init_val_search(j)*0.01
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
    if (.not. allocated(search_range)) allocate(search_range(2, 2*(nopt+nonefour)))
    do i = 1, nopt+nonefour
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
     
    do i = 1, 2*(nopt+nonefour)
        if (abs(y(i)) .ge. 10*abs(init_val_search(i))) then
            feasible = .false.
            return
        else if (abs(y(i)) .le. 0.01*abs(init_val_search(i))) then
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

    


