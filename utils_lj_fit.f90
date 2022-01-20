!utilities and look up tables for lj fitting trialanine
MODULE utils 
use graphf
implicit none

character(len=4), allocatable, dimension(:) :: ordering_array, lj_species
character(len=4), allocatable, dimension(:) :: opt_species, o_f_species
real*8, allocatable, dimension(:,:) :: chff_lj_params 
real*8, allocatable, dimension(:,:) :: o_f_array
real*8, allocatable, dimension(:,:,:) :: crd_data, dist_array
real*8, allocatable, dimension(:) :: ref_energies
character(len=100), allocatable, dimension(:) :: crd_names
integer :: natoms, nfiles, nspecies, nbonds, nonefour, nopt, npep_atoms
integer, allocatable, dimension(:,:) :: bond_array
integer, allocatable, dimension(:,:) :: excl_array
integer, allocatable, dimension(:) :: stan_lj_index, spec_lj_index
logical, allocatable, dimension(:) :: is_opt_arr
real*8 :: r_off, r_on
public ordering_array, lj_species, bond_array, nbonds, o_f_species, dist_array
public natoms, nfiles, nspecies, chff_lj_params, nopt, crd_data, npep_atoms
public excl_array, o_f_array, r_on, r_off, opt_species, ref_energies, crd_names
public stan_lj_index, spec_lj_index, is_opt_arr
public nonefour


contains
! initializes parameters of system

subroutine init_params(num_files, cut_on, cut_off, num_pep_atoms)
    implicit none
    integer :: num_atoms, num_files, num_pep_atoms, num_species 
    real*8 :: cut_on, cut_off 
    nfiles = num_files
    r_on = cut_on
    r_off = cut_off 
    npep_atoms = num_pep_atoms
   
end subroutine init_params

subroutine read_crd_file(file_name, crd_conf)
    implicit none
    character(len=100) :: file_name, junk
    real*8, allocatable, dimension(:,:) :: crd_conf    
    integer :: i, num_at
    if(.NOT. allocated(crd_conf)) allocate(crd_conf(natoms, 3)) 
    open(69, file=file_name, status = 'old')
    read(69,*)
    read(69,*)
    read(69,*)
    read(69,*) num_at, junk 
    do i = 1, num_at
        read(69,*) junk, junk, junk, junk, crd_conf(i,1), crd_conf(i,2), crd_conf(i,3) 
    end do
    
    close(69)
end subroutine read_crd_file


subroutine load_data(crd_name_file, base_name, ref_file)
    implicit none
    character(len=100) :: base_name, enum, file_name, ref_file, crd_name_file
    real*8, allocatable, dimension(:,:) :: crd_conf
    integer :: i, iatom, idir
    
    allocate(crd_data(natoms, 3, nfiles))
    allocate(crd_names(nfiles))
    open(69, file=crd_name_file, status = 'old')
    do i = 1, nfiles
        read(69,*) crd_names(i)
    end do
    crd_data = 0
    do i = 1, nfiles
        
        file_name = trim(base_name) // trim(crd_names(i))
        
        call read_crd_file(file_name, crd_conf)

        do iatom = 1, natoms
            do idir = 1, 3 
                crd_data(iatom, idir, i) = crd_conf(iatom, idir)
            end do
        end do 
    end do
    close(69)
    open(69, file=ref_file, status = 'old')
    allocate(ref_energies(nfiles))
    do i = 1, nfiles
        read(69, *) ref_energies(i)
    end do
    close(69)
    call crt_dist_array()
end subroutine load_data

! loads relevant parameters like atom ordering and lj species
subroutine load_lj_params(file_name)
    implicit none
    character(len=100) :: file_name
    integer :: i
    real*8 :: junk
    open(69, file=file_name, status = 'old')
    read(69, *) nspecies
    if (.not. allocated(lj_species)) allocate(lj_species(nspecies))
    if (.not. allocated(chff_lj_params)) allocate(chff_lj_params(nspecies, 2))
    do i = 1, nspecies
       read(69,*) lj_species(i), junk, chff_lj_params(i,1), chff_lj_params(i,2)
    end do
    close(69)
end subroutine load_lj_params

subroutine load_one_four_params(file_name, name_file)
    implicit none
    character(len=100) :: file_name, junk, name_file
    integer :: i
    open(69, file=file_name, status = 'old')
    read(69, *) junk 
    if (.not. allocated(o_f_array)) allocate(o_f_array(nspecies, 2))
    
    do i = 1, nspecies
        read(69, *) junk, o_f_array(i, 1), o_f_array(i, 2)
    end do
    close(69) 
    open(69, file=name_file, status = 'old')
    read(69,*) nonefour
    if (.not. allocated(o_f_species)) allocate(o_f_species(nonefour))
    do i = 1, nonefour
        read(69,*) o_f_species(i)
    end do 
    close(69) 
end subroutine load_one_four_params

subroutine load_psf(file_name)
    implicit none
    character(len=100) :: file_name, junk
    integer :: i
    open(69, file=file_name, status = 'old')
    do i = 1, 6
        read(69,*)
    end do
    read(69, *) natoms
    if (.not. allocated(ordering_array)) allocate(ordering_array(natoms)) 
    do i = 1, natoms
        read(69,*) junk, junk, junk, junk, junk, ordering_array(i)
    end do
    close(69)
end subroutine load_psf     

subroutine load_opt_spec(file_name)
    implicit none
    character(len=100) file_name
    integer :: i
    open(69, file=file_name, status = 'old')
    read(69, *) nopt
    allocate(opt_species(nopt))
    do i = 1, nopt
        read(69, *) opt_species(i)
    end do
end subroutine load_opt_spec

! converts atom index in psf file to its LJ label
character(len=3) function at_to_label(iatom)
    integer :: iatom
    at_to_label = ordering_array(iatom) 
end function at_to_label

! converts atom LJ label to index in LJ array
integer function lab_to_lj_i(label)
    character(len=3) :: label
    integer :: i    
    do i = 1, nspecies
        if (trim(label) == lj_species(i)) then
            lab_to_lj_i = i
        end if
    end do
end function lab_to_lj_i


! calculates euclidian distance between crd1 and crd2
real*8 function get_distance(crd1, crd2)
    real*8, dimension(3) :: crd1, crd2
    get_distance = NORM2(crd1 - crd2)
end function  get_distance

! calculates lj interaction energy between a pair of atoms at distance dist
real*8 function calc_lj_pair(eps1, rmin1, eps2, rmin2, dist)
    real*8 :: eps1, rmin1, eps2, rmin2, epsij, rminij, dist 

    epsij = SQRT(eps1 * eps2)
    rminij = rmin1 + rmin2  
    calc_lj_pair = epsij*((rminij/dist)**12 - 2*((rminij/dist)**6)) 
end function calc_lj_pair


logical function is_opt(iatom)
    integer :: iatom, i
    do i=1, nopt 
        if (at_to_label(iatom) == opt_species(i)) then 
            is_opt = .true.
            return
        else
            is_opt = .false.
        end if
    end do
end function is_opt

logical function is_o_f(iatom)
    integer :: iatom, i
    if (nonefour == 0) then
        is_o_f = .false.
        return
    end if
    do i = 1, nonefour 
        if (at_to_label(iatom) == o_f_species(i)) then
            is_o_f = .true.
            return
        else
            is_o_f = .false.
        end if
    end do
end function is_o_f

integer function label_to_opt_index(label)
    character(len=3) :: label
    integer :: i
    do i = 1, nopt
        if(trim(label) == opt_species(i)) then
            label_to_opt_index = i
            return
        end if
    end do
end function label_to_opt_index

!THIS NEEDS TO BE MADE GENERAL
integer function label_to_x_vec_o_f(label)
    character(len=3) :: label
    integer :: i, index_x
    if (label == "CT1") then
        label_to_x_vec_o_f = 14
    else if (label == "CT3") then
        label_to_x_vec_o_f = 15
    else if (label == "NH1") then
        label_to_x_vec_o_f = 16
    else if (label == "O  ") then
        label_to_x_vec_o_f = 17
    else if (label == "OB ") then
        label_to_x_vec_o_f = 18
    end if
end function label_to_x_vec_o_f

subroutine calc_look_ups() 
    implicit none
    integer :: iatom, lj_index, nopt_atoms
    integer :: i, j
    if (.not. allocated(is_opt_arr)) allocate(is_opt_arr(iatom))
    do iatom = 1, natoms
       is_opt_arr(iatom) = is_opt(iatom) 
    end do   
    if (.not. allocated(stan_lj_index)) allocate(stan_lj_index(natoms))
    do iatom = 1, natoms
        if (is_opt_arr(iatom) .eqv. .true.) then
            stan_lj_index(iatom) = label_to_opt_index(at_to_label(iatom))
        else
            stan_lj_index(iatom) = lab_to_lj_i(at_to_label(iatom)) 
        end if
    end do
end subroutine calc_look_ups
    
    

real*8 function get_eps_stand(iatom, x)
    integer :: iatom, lj_index
    logical :: optimize
    real*8, dimension(36) :: x
    optimize = is_opt_arr(iatom) 
    if (optimize .eqv. .true.) then
        lj_index = stan_lj_index(iatom) 
        get_eps_stand = -1*abs(x(2*lj_index-1))
        return
    else
        lj_index = stan_lj_index(iatom) 
        get_eps_stand = chff_lj_params(lj_index, 1)
        return
    end if
end function get_eps_stand


real*8 function get_rmin_stand(iatom, x)
    integer :: iatom, lj_index
    logical :: optimize
    real*8, dimension(36) :: x
    optimize = is_opt_arr(iatom) 
    if (optimize .eqv. .true.) then
        lj_index = stan_lj_index(iatom) 
        get_rmin_stand = abs(x(2*lj_index))
        return
    else
        lj_index = stan_lj_index(iatom) 
        get_rmin_stand = chff_lj_params(lj_index, 2)
        return
    end if
end function get_rmin_stand


real*8 function get_eps_spec(iatom, x)
    integer :: iatom, lj_index
    logical :: optimize
    real*8, dimension(36) :: x
    optimize = is_opt_arr(iatom) 
    if ((optimize .eqv. .true.) .and. (is_o_f(iatom) .eqv. .true.)) then
        lj_index = label_to_x_vec_o_f(at_to_label(iatom))
        get_eps_spec = -1*abs(x(2*lj_index-1))
        return
    else if ((optimize .eqv. .true.) .and. (is_o_f(iatom) .eqv. .false.)) then
        lj_index = stan_lj_index(iatom) 
        get_eps_spec = -1*abs(x(2*lj_index-1))
        return
    else  
        lj_index = stan_lj_index(iatom) 
        get_eps_spec = o_f_array(lj_index, 1)
    end if
end function get_eps_spec


real*8 function get_rmin_spec(iatom, x)
    integer :: iatom, lj_index
    logical :: optimize
    real*8, dimension(36) :: x
    optimize = is_opt_arr(iatom) 
    if ((optimize .eqv. .true.) .and. (is_o_f(iatom) .eqv. .true.)) then
        lj_index = label_to_x_vec_o_f(at_to_label(iatom))
        get_rmin_spec = abs(x(2*lj_index))
        return
    else if ((optimize .eqv. .true.) .and. (is_o_f(iatom) .eqv. .false.)) then
        lj_index = stan_lj_index(iatom) 
        get_rmin_spec = abs(x(2*lj_index))
        return
    else 
        lj_index = stan_lj_index(iatom) 
        get_rmin_spec = o_f_array(lj_index, 2)
    end if
end function get_rmin_spec

subroutine get_lj_energy(ifile, energy, x)
    implicit none
    real*8 :: energy, dist_ij, curr_lj_energy, eps1, eps2, rmin1, rmin2 
    integer :: iatom, jatom, bond_dist, i, ifile
    real*8, dimension(nopt*2) :: x, curr_sol
    real*8 :: t1, t2 
    energy = 0.0 
    curr_lj_energy = 0.0
    do iatom = 1, npep_atoms
        do jatom = iatom + 1, natoms
            bond_dist = excl_array(iatom, jatom) 
            dist_ij = dist_array(iatom, jatom, ifile)
            if (bond_dist .ge. 4) then
                eps1 = get_eps_stand(iatom, x)
                eps2 = get_eps_stand(jatom, x)
                rmin1 = get_rmin_stand(iatom, x)
                rmin2 = get_rmin_stand(jatom, x)  
                if ((eps1 .ge. 0.0d0) .or. (eps2 .ge. 0.0d0)) then
                do i = 1, nopt+5
                    write(*,*)
                    write(*,*) x(2*i-1), x(2*i) 
                end do
                end if
                curr_lj_energy = calc_lj_pair(eps1, rmin1, eps2, rmin2, dist_ij)
                curr_lj_energy = curr_lj_energy*vswitch(dist_ij, r_on, r_off)
            else if (bond_dist == 3) then
                eps1 = get_eps_spec(iatom, x)
                eps2 = get_eps_spec(jatom, x)
                rmin1 = get_rmin_spec(iatom, x)
                rmin2 = get_rmin_spec(jatom, x)
                curr_lj_energy = calc_lj_pair(eps1, rmin1, eps2, rmin2, dist_ij)    
                curr_lj_energy = curr_lj_energy*vswitch(dist_ij, r_on, r_off) 
            end if
            energy = energy + curr_lj_energy  
            curr_lj_energy = 0.0                      
        end do    
    end do
end subroutine get_lj_energy

real*8 function vswitch(dist, r_on, r_off)
    real*8 :: dist, r_on, r_off
    
    if (dist .le. r_on) then
        vswitch = 1.0
        return
    else if ((dist .gt. r_on) .and. (dist .le. r_off)) then
        vswitch = (r_off**2 - dist**2)**2*(r_off**2 + 2*dist**2 - 3*r_on**2)
        vswitch = vswitch / (r_off**2 - r_on **2)**3
        return
    else if (dist .gt. r_off) then
        vswitch = 0.0
        return
    end if
end function vswitch


subroutine crt_dist_array()
    implicit none
    integer :: iatom, jatom, ifile
    if (.not. allocated(dist_array)) allocate(dist_array(natoms, natoms, nfiles))
      do ifile = 1, nfiles
        do iatom = 1, natoms
          do jatom = 1, natoms
            dist_array(iatom, jatom, ifile) = get_distance(crd_data(iatom, :, ifile),&
                                                           crd_data(jatom, :, ifile))
          end do
        end do
      end do
end subroutine crt_dist_array


subroutine get_excl_array(bond_file)
    implicit none
    type(graph) :: G
    real*8 :: distance
    character(len=100) :: bond_file
    integer :: i, j, start, sto 
    logical :: bonded
    
    call G%setOrder(natoms)  
    open(69, file=bond_file, status = 'old')
    read(69,*) nbonds 
    if(.not. allocated(bond_array)) allocate(bond_array(nbonds, 2))
    do i = 1, nbonds
        read(69, *) bond_array(i, 1), bond_array(i, 2)
    end do 
    do i = 1, nbonds
        call G%addEdge(bond_array(i, 1), bond_array(i, 2))
    end do
    if (.not. allocated(excl_array)) allocate(excl_array(natoms, natoms)) 
    do i = 1, natoms
        do j = 1, natoms
            call exampleDijkstra(G, i, j, distance, bonded)
            if (bonded .eqv. .TRUE.) then
                excl_array(i,j) = INT(distance)
            else
                excl_array(i,j) = 1000 
            end if
        end do
    end do
end subroutine get_excl_array
    
    
end module utils
