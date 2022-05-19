
program lj_fit
use differential_evolution
use utils
implicit none
character(len=100) :: lj_param_chff, one_four_params, opt_file, junk
character(len=100) :: psf_file, bond_file, crd_path_file, one_four_species
character(len=100) :: ref_file, up_bound, low_bound, inp_multi
character(len=100) :: string
real*8, dimension(:,:), allocatable :: search_range
real*8, allocatable, dimension(:) :: init_val_search, x
real*8 :: energy, rmse, l_bound, u_bound, cons_multi
integer :: i, j
integer :: num_files, num_pep_atoms, n_onefour, ngenerations
character(len=100) :: snum_files, snum_pep_atoms, snone_four, sngenerations
logical :: verbose = .true.

! ============ START MAIN =================================================

!COMMAND LINE INPUTS
if (COMMAND_ARGUMENT_COUNT().ne.15) then
    write(*,*) "ERROR: NEED 15 ARGUMENTS"
    CALL EXIT(1)
end if
call GET_COMMAND_ARGUMENT(1, lj_param_chff) ! CHFF LJ parameters (normal int)
call GET_COMMAND_ARGUMENT(2, one_four_params) ! CHFF LJ one four params
call GET_COMMAND_ARGUMENT(3, one_four_species) ! CHFF one four species
call GET_COMMAND_ARGUMENT(4, psf_file) ! psf file of system
call GET_COMMAND_ARGUMENT(5, bond_file) ! bond data of system
call GET_COMMAND_ARGUMENT(6, crd_path_file) ! path to crd path file
call GET_COMMAND_ARGUMENT(7, opt_file) ! Atom species to optimize file
call GET_COMMAND_ARGUMENT(8, ref_file) ! reference energies
call GET_COMMAND_ARGUMENT(9, snum_files)
call GET_COMMAND_ARGUMENT(10, snum_pep_atoms)
call GET_COMMAND_ARGUMENT(11, snone_four)
call GET_COMMAND_ARGUMENT(12, sngenerations)
call GET_COMMAND_ARGUMENT(13, low_bound)
call GET_COMMAND_ARGUMENT(14, up_bound)
call GET_COMMAND_ARGUMENT(15, inp_multi)


read(snum_files, *) num_files
read(snum_pep_atoms, *) num_pep_atoms
read(snone_four, *) n_onefour
read(sngenerations, *) ngenerations
read(up_bound, *) u_bound
read(low_bound, *) l_bound
read(inp_multi, *) cons_multi
! initialize parameters and other things
call init_params(num_files          =       num_files,&
                 cut_on             =       10.0d0,&
                 cut_off            =       12.0d0,&
                 num_pep_atoms      =       num_pep_atoms,&
                 n_onefour          =       n_onefour,&
                 l_bound            =       l_bound,&
                 u_bound            =       u_bound,&
                 con_multi         =       cons_multi)

call init_system(psf_file, lj_param_chff, one_four_params, bond_file,&
opt_file, one_four_species)

call load_data(crd_path_file, ref_file)

allocate(init_val_search(2*(nopt+num_one_four)))
allocate(x(2*(nopt+num_one_four)))


init_val_search=0

! THIS IS STILL HARDCODED. IT IS A FILE WITH THE ORIGINAL LJ PARAMETERS WE ARE
! OPTIMIZING. IT WAS ORIGINALY ONLY FOR DEBUGGING PURPOSES BUT HAS TO BE THERE
! NOW SINCE WE ARE GENERATING OUR INITIAL GUESS FROM THIS
open(69, file="def_params/test_tip4.dat", status = 'old')

do i = 1, nopt
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i) 
end do
do i = nopt+1, nopt + num_one_four 
    read(69,*) junk, init_val_search(2*i-1), init_val_search(2*i)
end do
close(69) 

call init_print_helper()

call calc_look_ups()

call init_search_range(search_range)



write(*,*) "FOLLOWING ATOM TYPES ARE BEING OPTIMIZED:"
do i = 1,nopt
    write(*,*) print_helper(i)
end do
write(*,*) "RMSE TO CHFF PARAMS", opt_func(init_val_search)


call DE_init(set_range               = search_range,     &
             set_popSize             = 100,              &
             set_maxGens             = ngenerations,               &
             set_maxChilds           = 1,                &
             set_forceRange          = .false.,         &
             set_mutationStrategy    = DErand1,  &
             set_crossProb           = 0.9d0,             &
             set_verbose             = verbose,          &
             set_Nprint              = 10)


call DE_optimize(opt_func, feasible, sumconstr, x, guess=init_val_search)

write(*,*) "BEST SOLUTION:"
do i = 1, nopt+num_one_four
    write(*,*) x(2*i-1), x(2*i)
end do
write(*,*) "CALC ENER", "FIT ENERGY"
do i = 1,nfiles 
    call get_lj_energy(i, energy, x)
    write(*,*) energy, ref_energies(i), sqrt((energy-ref_energies(i))**2)
end do
write(*,*) "INITIAL RMSE: ", opt_func(init_val_search)
write(*,*) "FINAL RMSE WITH PENALTY: ", opt_func(x)
write(*,*) "RMSE PENALTY: ", rmse_penalty(x, init_val_search, rmse_multi)
write(*,*) "RMSE W/O PENALTY: ", opt_func(x) - rmse_penalty(x, init_val_search,&
                                                                rmse_multi)
write(*,*) "FINAL LJ PARAMS:"
do i = 1, nopt+num_one_four
    write(*,'(A4,F10.5,2F10.5)') print_helper(i), 0.0, -abs(x(2*i-1)), abs(x(2*i))
end do

! =============== END MAIN ================================================

contains
! objective function 
real*8 function opt_func(y)
    real*8, dimension(:) :: y
    real*8, dimension(nfiles) :: energies
    real*8 :: curr_energy, ref_val, rmse
    integer i, j
    real*8 :: t1, t2
    energies = 0
    call calc_full_lj(energies, y) 
    ref_val = sum(energies)/size(energies)
    energies = energies - ref_val
    rmse = 0.0d0   
    do i = 1,nfiles 
        rmse = rmse + (ref_energies(i) - energies(i))**2
    end do  
    opt_func = sqrt(rmse/nfiles)
    opt_func = opt_func + rmse_penalty(y, init_val_search, rmse_multi)
end function opt_func 

! RMSE penalty function to constrain the DE search space to only modify
! parameters that impact the RMSE on a larger scale
real*8 function rmse_penalty(y, ref, multi)
    real*8, dimension(:) :: y, ref
    real*8 :: multi
    rmse_penalty = 0.d0
    do i = 1, 2*(nopt+num_one_four) 
        rmse_penalty = rmse_penalty + ((y(i) - ref(i))*multi)**2
    end do 
    rmse_penalty = sqrt(rmse_penalty/(2*(nopt+num_one_four)))
end function rmse_penalty  

! calculates the lj energies for all files
subroutine calc_full_lj(energies, y)
    implicit none
    real*8, dimension(:) :: y
    real*8, dimension(:) :: energies
    integer :: i
    do i = 1, nfiles 
        call get_lj_energy(i, energies(i), y)
    end do
end subroutine calc_full_lj

! population initialization for DE
subroutine init_pop(pop)
    real*8, dimension(:,:), intent(out) :: pop
    real*8 :: ran
    integer :: popsize, i, j
    popsize = size(pop, dim=2)
    do i = 1, popsize
        do j = 1, 2*(nopt+num_one_four)
            call random_number(ran)
            pop(j,i) = (1.2*init_val_search(j)-0.8*init_val_search(j))*ran
            pop(j,i) = pop(j,i) + init_val_search(j)*0.8
        end do
    end do
end subroutine init_pop

! calls all functions and subroutines needed for initialization of the system
subroutine init_system(psf_file, lj_param_chff, one_four_params, bond_file,&
opt_file, one_four_species)
    character(len=100) :: psf_file, lj_param_chff, one_four_params, bond_file
    character(len=100) :: opt_file, one_four_species
    call load_psf(psf_file)
    call load_lj_params(lj_param_chff)
    call load_one_four_params(one_four_params, one_four_species)
    call load_opt_spec(opt_file)
    call get_excl_array(bond_file)
    
end subroutine init_system

! initializes DE search range
subroutine init_search_range(search_range)
    implicit none
    real*8, allocatable, dimension(:,:) :: search_range
    real*8 :: eps_min, eps_max, rmin_min, rmin_max
    integer i,j
    if (.not. allocated(search_range)) allocate(search_range(2,&
                                        2*(nopt+num_one_four)))
    do i = 1, nopt+num_one_four
        !eps parameters search range
        search_range(1, 2*i-1) = 1.5*init_val_search(2*i-1)
        search_range(2, 2*i-1) = 0.5*init_val_search(2*i-1)
        !rmin parameters search range
        search_range(1, 2*i) = 0.5*init_val_search(2*i)
        search_range(2, 2*i) = 1.5*init_val_search(2*i)
    end do
end subroutine init_search_range

! checks feasability of solution, the upper and lower bound variable are passed
! to the program from the command line
logical function feasible(y)
    implicit none    
    real*8, dimension(:) :: y 
    integer :: i,j
     
    do i = 1, 2*(nopt+num_one_four)
        if (abs(y(i)) .ge. upper_bound*abs(init_val_search(i))) then
            feasible = .false.
            return
        else if (abs(y(i)) .le. lower_bound*abs(init_val_search(i))) then
            feasible = .false.
            return
        else
            feasible = .true.
        end if
    end do

end function feasible

! sum constraint, currently not in use
real*8 function sumconstr(y)
    implicit none   
    real*8, dimension(:) :: y
    sumconstr = 0.0d0 
end function sumconstr


end program lj_fit

    


