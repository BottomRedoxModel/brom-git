! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov,
!                     Elizaveta Protsenko, Phil Wallhead
!-----------------------------------------------------------------------


    module io_ascii

    use fabm, only: type_model, fabm_get_bulk_diagnostic_data
    use fabm_driver
    use yaml_types
    use yaml,yaml_parse=>parse,yaml_error_length=>error_length
    use fabm_types, only: attribute_length, rk

    implicit none
    private

    type brom_par ! Type for sparse matrix
        real(rk) :: realvalue
        character(len=64) :: initname
        character(len=64) :: outname
        type(brom_par), pointer :: next
    end type brom_par

    type brom_name
        character(len=64) :: outname
        character(len=64) :: initname
        type(brom_name), pointer :: next
    end type brom_name

    type (brom_par), pointer :: first
    type (brom_name), pointer :: first_name

    public find_index, porting_initial_state_variables, init_common, saving_state_variables, saving_state_variables_diag,&
           get_brom_par, get_brom_name, svan, make_vert_grid, input_primitive_physics, input_ascii_physics, &
           make_physics_bbl_sed, make_physics_bbl_sed_old

    contains


!=======================================================================================================================
    subroutine init_common()

    class (type_node), pointer       :: node
    character(len=yaml_error_length) :: yaml_error
    integer                          :: unit_eff
    character(len=256)               :: path_eff

    unit_eff = 1
    path_eff = 'brom.yaml'
    ! Parse YAML file.
    node => yaml_parse(trim(path_eff),unit_eff,yaml_error)
    if (yaml_error/='') call fatal_error('fabm_create_model_from_yaml_file',trim(yaml_error))
    if (.not.associated(node)) call fatal_error('fabm_create_model_from_yaml_file', &
         'No configuration information found in '//trim(path_eff)//'.')
    select type (node)
         class is (type_dictionary)
            call tree(node)
         class is (type_node)
            call fatal_error('brom', trim(path_eff)//' must contain a dictionary &
               & at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
    end select

    end subroutine init_common
!=======================================================================================================================





!=======================================================================================================================
    subroutine tree(mapping)

    class (type_dictionary), intent(in) :: mapping
    class (type_node), pointer          :: node
    type (type_key_value_pair), pointer :: pair
    character(len=64)                   :: instancename

    node => mapping%get('instances')
    if (.not.associated(node)) &
        call fatal_error('create_model_tree_from_dictionary', 'No "instances" dictionary found at root level.')
    select type (node)
        class is (type_dictionary)
            pair => node%first
        class is (type_node)
            nullify(pair)
            call fatal_error('create_model_tree_from_dictionary',trim(node%path)// &
               ' must be a dictionary with (model name : information) pairs, not a single value.')
    end select

    ! only brom acceptance, iterate throw other models with no action
    do while (associated(pair))
        instancename = trim(pair%key)
        if (instancename == "brom") then
            select type (dict=>pair%value)
                class is (type_dictionary)
                    call from_tree(instancename, dict)
                class is (type_null)
                    call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                        trim(instancename)//'" must be a dictionary, not a single value or nothing')
                class is (type_node)
                    call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                        trim(instancename)//'" must be a dictionary, not a single value.')
                end select
        end if
        pair => pair%next
    end do

    end subroutine tree
!=======================================================================================================================






!=======================================================================================================================
    subroutine from_tree(instancename, node)

    character(len=*),       intent(in)          :: instancename
    class (type_dictionary),intent(in)          :: node
    class (type_dictionary),pointer             :: childmap
    type (type_error),pointer                   :: config_error

    childmap => node%get_dictionary('initialization',required=.false.,error=config_error)
    if (associated(config_error)) call fatal_error('create_model_from_dictionary',config_error%message)
    if (associated(childmap)) call parse_initialization(childmap)

    end subroutine from_tree
!=======================================================================================================================






!=======================================================================================================================
    subroutine parse_initialization(node)

    class (type_dictionary),intent(in)  :: node

    type (type_key_value_pair), pointer :: pair
    type (brom_par), pointer :: current
    type (brom_name), pointer :: current_name

    logical             :: success
    real(rk)            :: realvalue
    character(len=64)   :: initname
    character(len=64)   :: outname

    ! Transfer user-specified initial state to the model.
    nullify(first)
    pair => node%first
    do while (associated(pair))
        select type (value=>pair%value)
            class is (type_scalar)
                initname = trim(pair%key)
                realvalue = value%to_real(default=real(0,real_kind),success=success)
                if (.not.success) then !Assume the value is a string
                    realvalue = huge(1.0_rk)
                    outname = value%string
                else
                    outname = repeat('z',64)
                end if
            class is (type_null)
                call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
            class is (type_dictionary)
                call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
            end select
        allocate (current)
        current = brom_par(realvalue, initname, outname, first)
        first => current
        pair => pair%next
    end do

    end subroutine parse_initialization
!=======================================================================================================================






!=======================================================================================================================
    function get_brom_par(initname,defvalue) result(realvalue)

    real(rk)                       :: realvalue
    character(len=*), intent(in)   :: initname
    real(rk), optional, intent(in) :: defvalue
    type (brom_par), pointer       :: current

    current => first
    do
        if (.not.associated(current)) exit
        if (trim(initname) == trim(current%initname)) then
            realvalue = current%realvalue
            return
        end if
        current => current%next
    end do

    if (present(defvalue)) then
        realvalue = defvalue
    else
        print *, "Check brom.yaml or name of the input variable for ", initname
        stop
    end if

    end function get_brom_par
!=======================================================================================================================






!=======================================================================================================================
    function get_brom_name(initname) result(outname)

    character(len=64) :: outname
    character(len=*), intent(in)   :: initname
    type (brom_par), pointer      :: current

    current => first
    do
        if (.not.associated(current)) exit
        if (trim(initname) == trim(current%initname)) then
            outname = trim(current%outname)
            return
        end if
        current => current%next
    end do

    print *, "Check brom.yaml or name of the input variable for ", initname
    stop

    end function get_brom_name
!=======================================================================================================================






!=======================================================================================================================
    function find_index(names, name) result(index)

    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: name
    integer :: index

    do index = 1, size(names)
        if (names(index) == name) return
    end do
    index = -1

    end function
!=======================================================================================================================






!=======================================================================================================================
    subroutine porting_initial_state_variables(icfile_name, model_year, julianday, i_max, k_max, par_max, &
                                                par_name, cc, vv)

    !Reads in initial conditions from an ascii file <icfilename>
    !Note: format should match that used in subroutine saving_state_variables

    !Input variables
    character(len=*), intent(in)            :: icfile_name
    integer, intent(in)                     :: model_year, julianday, i_max, k_max, par_max
    character(len=*), intent(in)            :: par_name(:)

    !Output variables
    real(rk), dimension(:,:,:), pointer, intent(out) :: cc, vv

    !Local variables
    integer                                 :: d_year_start, julianday_start, k_max_start, par_max_start
    integer, allocatable                    :: column2state(:)
    character(20000)                        :: labels
    character(200)                          :: comments
    integer                                 :: i, j, k, ip, istart, istop, foo, m
    real(rk)                                :: value

    open(9,file=icfile_name)
    read(9,'( 1x,i4,1x,i5,1x,i5,1x,i5 )') d_year_start, julianday_start, k_max_start, par_max_start
    if (k_max_start /= k_max) stop 'number of levels in start-up file does not match number of levels in model'

    allocate(column2state(par_max))
    read(9, '( 1x, a )') labels

    istart = 1
    do ip = 1, par_max_start
        istop = istart - 1 + index(labels(istart:), ' ')
        column2state(ip) = find_index(par_name, labels(istart:istop-1))
        istart = istop + 1
    end do

    read(9, *) comments
    do i=1,i_max
        do k=1,k_max
            read(9, '( i5 )', advance = 'no') foo
            read(9, '( i5 )', advance = 'no') foo
            read(9, '( 1x,f8.4 )', advance = 'no') value !we don't read again z()
            do m=1,4
                read(9, '( 1x,f15.9 )', advance = 'no') value !we don't read again t,s,Kz,hz
            end do
            read(9, '( 1x,f15.9 )', advance = 'no') vv(i,k,1)
            do ip=1,par_max
                read(9, '( 1x, f17.9 )', advance = 'no') value
                if (ip /= -1) then
                    if (value /= 0.) then
                        cc(i,k,ip) = value
                    else
                        cc(i,k,ip) = 0.000
                    end if
                end if
            end do
            read(9, *)
        end do
    end do
    close(9)

    end subroutine porting_initial_state_variables
!=======================================================================================================================





!=======================================================================================================================
    subroutine saving_state_variables_diag(outfile_name, model_year, julianday, i_max, k_max, par_max, par_name, &
                                        z, hz, cc, vv, t, s, kz, model, extend_out)

    !Saves final conditions in an ascii file <outfilename>

    !Input variables
    character (len = *), intent(in)            :: outfile_name
    integer, intent(in)                        :: model_year, julianday, i_max, k_max, par_max
    character(len=*), intent(in)               :: par_name(:)
    real(rk), dimension(:), intent(in)         :: z, hz
    real(rk), dimension(:,:,:), intent(in)     :: cc
    real(rk), dimension(:,:,:), intent(in)     :: vv
    real(rk), dimension (:,:,:), intent(in)    :: t, s, kz
    type (type_model), intent(in)              :: model
    logical, intent(in)                        :: extend_out

    !Local variables
    integer :: ip, k, i, ilast
    real    :: temp_matrix(i_max,k_max)

    open(10,file=outfile_name)

    ! First line with year, julian day, number of levels, number of state variables
    write(10,'( i4,1x,i5,1x,i5,1x,i5 )') model_year, julianday, k_max, par_max
    do ip=1, par_max
        write(10,'(1x,a)',advance='NO') trim(par_name(ip))
    end do
    write(10,*)

    write(10,'(3h i ,3h k ,6hDepth ,4h t  ,4h s   ,4h Kz ,4h hz ,4hVol )',advance='NO')
    do ip=1,par_max
        write(10,'(1x,a)',advance='NO') trim(par_name(ip)(15:))
    end do

    if (extend_out .eqv. .true.) then
        do ip = 1, size(model%diagnostic_variables)
            if (model%diagnostic_variables(ip)%save) then
                ilast = index(model%diagnostic_variables(ip)%path,'/',.true.)
                write(10,'(1x,a)',advance='NO') trim(model%diagnostic_variables(ip)%path(ilast+1:))
            end if
        end do
    end if

    write(10,*)
    !Subsequent lines: one for each depth level
    do i=1,i_max
        do k=1,k_max
            write(10,'(i5)',advance='NO') i
            write(10,'(i5)',advance='NO') k
            write(10,'(1x,f8.4)',advance='NO') z(k)
            write(10,'(1x,f15.9)',advance='NO') t(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') s(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') kz(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') hz(k)
            write(10,'(1x,f15.9)',advance='NO') vv(i,k,1)
            do ip=1,par_max
                write(10,'(1x,f15.9)',advance='NO') cc(i,k,ip)
            end do
            if (extend_out .eqv. .true.) then
                do ip = 1, size(model%diagnostic_variables)
                    if (model%diagnostic_variables(ip)%save) then
                        !temp_matrix = fabm_get_bulk_diagnostic_data(model, ip)
                        write(10,'( 1x,f15.9 )',advance='NO') fabm_get_bulk_diagnostic_data(model, ip)

                    end if
                end do
            end if
            write(10,*)
        end do
    end do

    close(10)

    end subroutine saving_state_variables_diag
!=======================================================================================================================






!=======================================================================================================================
    subroutine saving_state_variables(outfile_name, model_year, julianday, i_max, k_max, par_max, par_name, &
                                        z, hz, cc, vv, t, s, kz)

    !Saves final conditions in an ascii file <outfilename>

    !Input variables
    character (len = *), intent(in)            :: outfile_name
    integer, intent(in)                        :: model_year, julianday, i_max, k_max, par_max
    character(len=*), intent(in)               :: par_name(:)
    real(rk), dimension(:), intent(in)         :: z, hz
    real(rk), dimension(:,:,:), intent(in)     :: cc
    real(rk), dimension(:,:,:), intent(in)     :: vv
    real(rk), dimension (:,:,:), intent(in)    :: t, s, kz

    !Local variables
    integer :: ip, k, i


    open(10,file=outfile_name)

    ! First line with year, julian day, number of levels, number of state variables
    write(10,'( i4,1x,i5,1x,i5,1x,i5 )') model_year, julianday, k_max, par_max
    do ip=1, par_max
        write(10,'(1x,a)',advance='NO') trim(par_name(ip))
    end do
    write(10,*)

    write(10,'(3h i ,3h k ,6hDepth ,4h t  ,4h s   ,4h Kz ,4h hz ,4hVol )',advance='NO')
    do ip=1,par_max
        write(10,'(1x,a)',advance='NO') trim(par_name(ip)(15:))
    end do
    write(10,*)

    !Subsequent lines: one for each depth level
    do i=1,i_max
        do k=1,k_max
            write(10,'(i5)',advance='NO') i
            write(10,'(i5)',advance='NO') k
            write(10,'(1x,f8.4)',advance='NO') z(k)
            write(10,'(1x,f15.9)',advance='NO') t(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') s(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') kz(i,k,julianday)
            write(10,'(1x,f15.9)',advance='NO') hz(k)
            write(10,'(1x,f15.9)',advance='NO') vv(i,k,1)
            do ip=1,par_max
                write(10,'(1x,f17.9)',advance='NO') cc(i,k,ip)
            end do
            write(10,*)
        end do
    end do

    close(10)

    end subroutine saving_state_variables
!=======================================================================================================================







!=======================================================================================================================
    subroutine make_vert_grid(z, dz, hz, z_w, dz_w, hz_w, k_wat_bbl, k_max, k_bbl_sed)

    !Constructs full vertical grid (water + sediments) given water column parameters and parameters from brom.yaml

    !Input variables
    real(rk), dimension(:), intent(in)  :: z_w, dz_w, hz_w    !Water column layer midpoints (z_w), increments between them (dz_w), and layer thicknesses (hz_w)
    integer, intent(in)                 :: k_wat_bbl          !Number of layers above the water-BBL interface
    integer, intent(in)                 :: k_max              !Total number of layers (set in init_brom_transport using k_wat_bbl and k_points_below_water from brom.yaml)

    !Output variables
    real(rk), dimension(:), intent(out) :: z, dz, hz          !Full column layer midpoints (z), increments between them (dz), and layer thicknesses (hz)
    integer, intent(out)                :: k_bbl_sed          !Number of layers above the BBL-sediment interface (SWI)

    !Local variables
    real(rk)                            :: bbl_thickness      !Input BBL thickness [m]
    real(rk)                            :: hz_bbl_min         !Minimum hz in the BBL [m]
    real(rk)                            :: hz_sed_min         !Minimum hz in the sediments [m]
    real(rk)                            :: hz_sed_max         !Maximum hz in the sediments [m]
    real(rk)                            :: scale_fac, scale_fac_next !Scale factors to ensure that total BBL thickness = bbl_thickness
    integer                             :: k

    ! The grid structure is:
    !
    !======================= (air-sea interface)
    !           o            z(1), hz(1)
    !-----------------------
    !           o            z(2)=z(1)+dz(1), hz(2)
    !           :                 :
    !           :                 :
    !           o            z(k_wat_bbl), hz(k_wat_bbl)   NOTE: k_wat_bbl is set by brom.yaml or (with priority) by netcdf input
    !======================= (water-bbl interface)
    !           o            z(k_wat_bbl+1)=z(k_wat_bbl)+dz(k_wat_bbl), hz(k_wat_bbl+1)
    !-----------------------
    !           o            z(k_wat_bbl+2), hz(k_wat_bbl+2)
    !           :                 :
    !           :                 :
    !           o            z(k_bbl_sed), hz(k_bbl_sed)  NOTE: k_bbl_sed is set by this routine
    !======================= (bbl-sediment interface)
    !           o            z(k_bbl_sed+1)=z(k_bbl_sed)+dz(k_bbl_sed), hz(k_bbl_sed+1)
    !-----------------------
    !           o            z(k_bbl_sed+2), hz(k_bbl_sed+2)
    !           :                 :
    !           :                 :
    !           o            z(k_max)=z(k_max-1)+dz(k_max-1), hz(k_max)  NOTE: k_max is set by this routine using k_bbl_sed and input k_sed
    !======================= (bottom boundary)



    !Get parameters from brom.yaml
    bbl_thickness = get_brom_par("bbl_thickness")
    hz_bbl_min = get_brom_par("hz_bbl_min")
    hz_sed_min = get_brom_par("hz_sed_min")
    hz_sed_max = get_brom_par("hz_sed_max")


    !Define layer midpoint depths [m] and their thicknesses/increments in the water column
    z(1:k_wat_bbl) = z_w(1:k_wat_bbl)
    hz(1:k_wat_bbl) = hz_w(1:k_wat_bbl)
    dz(1:k_wat_bbl-1) = dz_w(1:k_wat_bbl-1)


    !Adjust the input grid to allow high-resolution BBL to be inserted (if not too thick)
    if (bbl_thickness.gt.hz(k_wat_bbl)) then
        write(*,*) "Requested BBL thickness exceeds input thickness of bottom pelagic layer: cannot insert BBL"
        stop
    else
        z(k_wat_bbl) = z_w(k_wat_bbl) - 0.5_rk*bbl_thickness
        hz(k_wat_bbl) = hz_w(k_wat_bbl) - bbl_thickness
        dz(k_wat_bbl-1) = dz_w(k_wat_bbl-1) - 0.5_rk*bbl_thickness
    end if


    !Define the grid in the BBL
    do k=k_wat_bbl+1,k_max
        hz(k)   = bbl_thickness*(0.5_rk**(k-k_wat_bbl))
        dz(k-1) = 0.5_rk*(hz(k-1)+hz(k))
        z(k)    = z(k-1) + dz(k-1)   !The prescribed BBL thickness is traversed with layers geometrically decreasing in thickness
        scale_fac_next = bbl_thickness / (0.5_rk*hz(k)+sum(hz(k_wat_bbl+1:k))) !Rescaling factor after the NEXT layer is included
        if (0.5_rk*hz(k)*scale_fac_next<hz_bbl_min) then !If the NEXT layer will be too thin (after rescaling), complete the BBL by rescaling all layers so that total thickness = bbl_thickness
            scale_fac = bbl_thickness/sum(hz(k_wat_bbl+1:k))
            hz(k_wat_bbl+1:k) = scale_fac * hz(k_wat_bbl+1:k)
            hz(k)=hz(k)-hz_sed_min  !here we add one thin (hz_sed_min) layer above the SWI
            hz(k+1)=hz_sed_min
            k_bbl_sed = k+1
!            k_bbl_sed = k
            goto 1
        end if
    end do
1   continue


    !Recalculate dz and z using the rescaled hz
    do k=k_wat_bbl+1,k_bbl_sed
        dz(k-1) = 0.5_rk*(hz(k-1)+hz(k))
        z(k)    = z(k-1) + dz(k-1)
    end do


    !Define the grid in the sediments
    do k=k_bbl_sed+1,k_max
        hz(k)   = min(hz_sed_min*(2.0_rk**(k-k_bbl_sed-1)), hz_sed_max)
        dz(k-1) = 0.5_rk*(hz(k-1)+hz(k))
        z(k)    = z(k-1) + dz(k-1)  !Layer thickness is increased geometrically from minimum value hz(k_bbl_sed+1) = hz_sed_min
    end do
    dz(k_max) = 0.0_rk


    !Record the vertical grid in an ascii output file
    open(10,FILE = 'Vertical_grid.dat')
    write(10,'(5h k   ,5hz[k]   ,6hdz[k]   ,6hhz[k] )')
    do k=1,k_max
        write(10,'(1x,i4, 1x,f8.4, 5h _o_ , 1x,f8.4)') k,z(k),hz(k)
        write(10,'(2x,i4, 5h ===  , 1x,f8.4,5h ===  )') k,dz(k)
    end do
    close(10)

    end subroutine make_vert_grid
!=======================================================================================================================







!=======================================================================================================================
    subroutine input_primitive_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_max, days_in_yr)

    !Constructs (temperature, salinity, diffusivity) forcings for the water column assuming a simple hydrophysical scenario
    !with hypothetical sinusoidal variations

    !Input variables
    integer, intent(in)                                :: k_wat_bbl, i_max, days_in_yr
    real(rk), intent(in)                               :: water_layer_thickness

    !Output variables
    real(rk), allocatable, dimension (:,:,:), intent(out)  :: t_w, s_w, kz_w
    real(rk), allocatable, dimension (:), intent(out)      :: z_w, dz_w, hz_w

    !Local variables
    real(rk), allocatable, dimension(:)                :: sal_sum, sal_win, tem_sum, tem_win, dens
    integer                                            :: i, k, iday


    !Allocate local variables
    allocate(sal_sum(k_wat_bbl))
    allocate(sal_win(k_wat_bbl))
    allocate(tem_sum(k_wat_bbl))
    allocate(tem_win(k_wat_bbl))
    allocate(dens(k_wat_bbl))
    allocate(z_w(k_wat_bbl))
    allocate(dz_w(k_wat_bbl))
    allocate(t_w(i_max,k_wat_bbl,days_in_yr))
    allocate(s_w(i_max,k_wat_bbl,days_in_yr))
    allocate(kz_w(i_max,k_wat_bbl,days_in_yr))


!!! HERE THE USER MUST ADAPT FOR HIS/HER DESIRED SCENARIO --------------------------------------------------------------
!
    !Make grid parameters assuming uniform grid defined by water_layer_thickness and k_wat_bbl in brom.yaml
    hz_w = water_layer_thickness/k_wat_bbl
    dz_w = hz_w
    z_w(1) = 0.5_rk*hz_w(1)
    do k=2,k_wat_bbl
        z_w(k) = z_w(k-1) + dz_w(k-1)
    end do


    !Set summer and winter vertical distributions
    do k=1,k_wat_bbl
        sal_sum(k) = 35.0_rk + 0.01_rk*k
        sal_win(k) = 35.0_rk + 0.001_rk*k
        tem_sum(k) = 25.0_rk - 0.1_rk*k
        tem_win(k) = 10.0_rk - 0.001_rk*k
    enddo

    !Seasonal variability of T and S
    do iday=1,days_in_yr
        do k=1,k_wat_bbl
            s_w(i_max,k,iday) = sal_win(k)+0.5*(1.+sin(2*3.14*(iday+240.)/365.))*(sal_sum(k)-sal_win(k))
            t_w(i_max,k,iday) = tem_win(k)+0.5*(1.+sin(2*3.14*(iday+240.)/365.))*(tem_sum(k)-tem_win(k))
        enddo
    enddo

    !Calculate  vertical turbulent coefficient Kz from density (Gargett)
    do iday=1,days_in_yr
        do k=1, k_wat_bbl
            call svan(s_w(i_max,k,iday),t_w(i_max,k,iday),z_w(k),dens(k)) !calculate density as f(p,t,s)
        end do
        do k=1, k_wat_bbl-1
            Kz_w(i_max,k,iday)=0.5E-6 /&
                ((9.81/(1000.+(dens(k)+dens(k+1))/2.)&
                *(abs(dens(k+1)-dens(k))/dz_w(k)) &
                )**0.5)
        end do
    end do
!!!---------------------------------------------------------------------------------------------------------------------

    !Write to output ascii file
    open(12,file = 'Hydrophysics.dat')
        do iday=1,days_in_yr
            do k=1,k_wat_bbl
                write(12,'(2(1x,i4))',advance='NO') iday, k
                do i=1,i_max
                    write (12,'(1x,i4, 2(1x,f8.4),1x,f15.11,2f7.3)',advance='NO') i, t_w(i_max,k,iday), s_w(i_max,k, iday), kz_w(i_max,k,iday), dz_w(k), z_w(k)
                end do
                write(12,*)
           end do
        end do
    close(12)

    end subroutine input_primitive_physics
!=======================================================================================================================







!=======================================================================================================================
    subroutine input_ascii_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_max, days_in_yr)

    !Constructs (temperature, salinity, diffusivity) forcings for water column using input physics in ascii format

    !Input variables
    integer, intent(in)                                     :: k_wat_bbl, i_max, days_in_yr
    real(rk), intent(in)                                    :: water_layer_thickness

    !Output variables
    real(rk), allocatable, dimension (:,:,:), intent(out)  :: t_w, s_w, kz_w
    real(rk), allocatable, dimension (:), intent(out)      :: z_w, dz_w, hz_w

    !Local variables
    real(rk), allocatable, dimension(:)                     :: sal_sum, sal_win, tem_sum, tem_win, dens
    integer                                                 :: i, j, k, iday, k_max_TEL
!    real(rk) :: tem_TEL(55,365), sal_TEL(55,365), kz_TEL(55,365) ! this is for TELEMARK model outputs for Berre
    real(rk) :: tem_TEL(100,365), sal_TEL(100,365), kz_TEL(100,365) ! this is for TELEMARK model outputs for Berre

    !Allocate local variables
    allocate(sal_sum(k_wat_bbl))
    allocate(sal_win(k_wat_bbl))
    allocate(tem_sum(k_wat_bbl))
    allocate(tem_win(k_wat_bbl))
    allocate(dens(k_wat_bbl))
    allocate(z_w(k_wat_bbl))
    allocate(dz_w(k_wat_bbl))
    allocate(hz_w(k_wat_bbl))
    allocate(t_w(i_max,k_wat_bbl,days_in_yr))
    allocate(s_w(i_max,k_wat_bbl,days_in_yr))
    allocate(kz_w(i_max,k_wat_bbl,days_in_yr))

!!! HERE THE USER MUST ADAPT FOR HIS/HER PARTICULAR APPLICATION --------------------------------------------------------
    goto 2
!
1   continue ! CASE BERRE:
!
! The code below is provided as an example, reading (T,S) data from an ascii file of interpolated output from a model
! (TELEMARC) and using the resulting seawater density to compute Kz for a particular application (Berre)
!____________________________________________________________________
!   TELEMARK produced file
!____________________________________________________________________
    !Make grid parameters assuming uniform grid defined by water_layer_thickness and k_wat_bbl in brom.yaml
    hz_w = water_layer_thickness/k_wat_bbl
    dz_w = hz_w
    z_w(1) = 0.5_rk*hz_w(1)
    do k=2,k_wat_bbl
        z_w(k) = z_w(k-1) + dz_w(k-1)
    end do

    !Read (T,S) data for Berre
    k_max_TEL= 35
    open(19, file='tem2_new.dat')
    do i=1,days_in_yr
        do j=1,k_max_TEL
            read(19, *) tem_TEL(j,i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close(19)
    open(20, file='sal2_new.dat')
    do i=1,days_in_yr
        do j=1,k_max_TEL
            read(20, *) sal_TEL(j,i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close(20)
    do iday=1,days_in_yr
        do k=1,k_wat_bbl
            s_w(i_max,k,iday) = sal_TEL(k,iday)
            t_w(i_max,k,iday) = tem_TEL(k,iday)
        enddo
    enddo

    !Calculate  vertical turbulent coefficient Kz from density (Gargett)
    do iday=1,days_in_yr
        do k=1, k_wat_bbl
            call svan(s_w(i_max,k,iday),t_w(i_max,k,iday),z_w(k),dens(k)) !calculate density as f(p,t,s)
        end do
        do k=1, k_wat_bbl-1
            Kz_w(i_max,k,iday)=0.5E-6 /&
                ((9.81/(1000.+(dens(k)+dens(k+1))/2.)&
                *(abs(dens(k+1)-dens(k))/dz_w(k)) &
                )**0.5)
        end do
            Kz_w(i_max,k_wat_bbl,iday)=Kz_w(i_max,k_wat_bbl-1,iday)
    end do
!---end  CASE BERRE-----------------------------------------------------------------------------------------------------
2 continue ! CASE SPANGEREID:
!
! The code below is provided as an example, reading (T,S) data from an ascii file of interpolated output from NODC 
! from the North Sea and using the resulting seawater density to compute Kz for a particular application (Spangereid)
!____________________________________________________________________
    !Make grid parameters assuming uniform grid defined by water_layer_thickness and k_wat_bbl in brom.yaml
    hz_w = water_layer_thickness/k_wat_bbl
    dz_w = hz_w
    z_w(1) = 0.5_rk*hz_w(1)
    do k=2,k_wat_bbl
        z_w(k) = z_w(k-1) + dz_w(k-1)
    end do

    !Read (T,S) data for Berre
    k_max_TEL= 100
    open(19, file='spa_t.dat')
    do j=1,k_max_TEL
        do i=1,days_in_yr
            read(19, *) k,k,tem_TEL(j,i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close(19)
    open(20, file='spa_s.dat')
    do j=1,k_max_TEL
        do i=1,days_in_yr
            read(20, *) k,k,sal_TEL(j,i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close(20)
    do iday=1,days_in_yr
        do k=1,k_wat_bbl
            s_w(i_max,k,iday) = sal_TEL(k,iday)
            t_w(i_max,k,iday) = tem_TEL(k,iday)
        enddo
    enddo

    !Calculate  vertical turbulent coefficient Kz from density (Gargett)
    do iday=1,days_in_yr
        do k=1, k_wat_bbl
            call svan(s_w(i_max,k,iday),t_w(i_max,k,iday),z_w(k),dens(k)) !calculate density as f(p,t,s)
        end do
        do k=1, k_wat_bbl-1
            Kz_w(i_max,k,iday)=1.0E-5 /& !0.5E-4
                ((9.81/(1000.+(dens(k)+dens(k+1))/2.)&
                *(abs(dens(k+1)-dens(k))/dz_w(k)) &
                )**0.4)
        end do
            Kz_w(i_max,k_wat_bbl,iday)=Kz_w(i_max,k_wat_bbl-1,iday)
    end do
!---end  CASE SPANGEREID--------------------------------------------------------------------------------------------------

    
    
    end subroutine input_ascii_physics
!=======================================================================================================================






!=======================================================================================================================
    subroutine make_physics_bbl_sed(t, s, kz, hmix_rate, cc_hmix, u_x, u_x_w, t_w, s_w, kz_w, cc_hmix_w, kz_mol, kz_bio, &
        z, dz, hz, i_min, i_max, k_wat_bbl, k_bbl_sed, k_max, par_max, steps_in_yr, alpha, is_solid, phi, phi1, phi_inv, &
        pF1, pF2, mu0_musw, tortuosity, w_b, u_b, rho, dt, freq_turb, par_name, diff_method, bioturb_across_SWI, &
        ip_sol, ip_par, dphidz_SWI, wat_content, pWC, input_step)
    
    !Constructs profiles for (temperature, salinity, vertical diffusivity etc.) for the full column (water + sediments)
    !given water column variability and parameters from brom.yaml

    !Input variables
    integer, intent(in)                       :: i_min, i_max, k_wat_bbl, k_bbl_sed, k_max, par_max, freq_turb
    integer, intent(in)                       :: diff_method, bioturb_across_SWI, ip_sol, ip_par
    integer                                   :: input_step, steps_in_yr
    integer, dimension(:), intent(in)         :: is_solid
    real(rk), intent(in)                      :: dt
    real(rk), dimension(:), intent(in)        :: z, dz, hz
    real(rk), dimension(:,:,:), intent(in)    :: t_w, s_w, kz_w, u_x_w
    real(rk), dimension(:,:,:,:), intent(in)  :: cc_hmix_w
    character(len=attribute_length), dimension(:), intent(in) :: par_name
    
    !Output variables
    real(rk)                                  :: mu0_musw, dphidz_SWI
    real(rk), dimension(:), intent(out)       :: rho
    real(rk), dimension(:,:), intent(out)     :: kz_bio, alpha, phi, phi_inv, tortuosity, w_b, u_b, wat_content
    real(rk), dimension(:,:,:), intent(out)   :: t, s, kz, u_x, hmix_rate, kz_mol, pF1, pF2, pWC
    real(rk), dimension(:,:,:,:), intent(out) :: cc_hmix

    !Local variables
    integer                                   :: i, k, ip, sel, iday, istep
    integer                                   :: k_sed(k_max-k_bbl_sed), k_sed1(k_max+1-k_bbl_sed), k_bbl1(k_bbl_sed-k_wat_bbl)
    real(rk)                                  :: z_wat_bbl, z_bbl_sed, kz_gr, z1(k_max+1), z_s1(k_max+1), phi1(i_max,k_max+1)
    integer                                   :: kz_bbl_type, dynamic_kz_bbl
    real(rk)                                  :: kz_bbl_max, dbl_thickness, kz_mol0
    real(rk)                                  :: a1_bioirr, a2_bioirr
    real(rk)                                  :: kz_bioturb_max, z_const_bioturb, z_decay_bioturb
    real(rk)                                  :: phi_0, phi_inf, z_decay_phi, w_binf, rho_def, wat_con_0, wat_con_inf
    real(rk)                                  :: kzCFL(k_bbl_sed-1,steps_in_yr), kz_molCFL(k_max-1,par_max)

    write(*,*) "Sed 1"

    !!Get diffusivity parameters from brom.yaml
    !Turbulence in the benthic boundary layer
    kz_bbl_type = get_brom_par("kz_bbl_type")
    kz_bbl_max = get_brom_par("kz_bbl_max")
    dbl_thickness = get_brom_par("dbl_thickness")
    dynamic_kz_bbl = get_brom_par("dynamic_kz_bbl")
    write(*,*) "Sed 2"
    !Molecular diffusivity of solutes (single constant value, infinite dilution)
    kz_mol0 = get_brom_par("kz_mol0")
    mu0_musw = get_brom_par("mu0_musw")
    write(*,*) "Sed 3"
    !Bioturbation
    kz_bioturb_max = get_brom_par("kz_bioturb_max")
    z_const_bioturb = get_brom_par("z_const_bioturb")
    z_decay_bioturb = get_brom_par("z_decay_bioturb")
    write(*,*) "Sed 4"
    !Bioirrigation
    a1_bioirr = get_brom_par("a1_bioirr")
    a2_bioirr = get_brom_par("a2_bioirr")
    write(*,*) "Sed 5"
    !Porosity
    phi_0 = get_brom_par("phi_0")
    phi_inf = get_brom_par("phi_inf")
    z_decay_phi = get_brom_par("z_decay_phi")
    wat_con_0 = get_brom_par("wat_con_0")
    wat_con_inf = get_brom_par("wat_con_inf")
    write(*,*) "Sed 6"
    !Vertical advection in the sediments
    w_binf = get_brom_par("w_binf")

    !Density of particles
    rho_def = get_brom_par("rho_def")
    rho = 0.0_rk
    do ip=1,par_max
        if (is_solid(ip).eq.1) then
            rho(ip) = get_brom_par('rho_' // trim(par_name(ip)),rho_def)
            write(*,*) "Assumed density of ", trim(par_name(ip)), " = ", rho(ip)
        end if
    end do

    write(*,*) "Sed 7"
    !!Set useful parameters for calculations
    !Useful depths
    z_wat_bbl = z(k_wat_bbl+1) - 0.5_rk*hz(k_wat_bbl+1) !Depth of water-BBL interface
    z_bbl_sed = z(k_bbl_sed+1) - 0.5_rk*hz(k_bbl_sed+1) !Depth of BBL-sediment interface (SWI)
    z1(1:k_max) = z(:)-0.5_rk*hz(:)         !Depth of layer interfaces
    z1(k_max+1) = z(k_max)+0.5_rk*hz(k_max)
    z_s1 = z1-z_bbl_sed                     !Depth of interfaces wrt SWI
    write(*,*) "Sed 8"
    !Useful index vectors
    k_sed = (/(k,k=k_bbl_sed+1,k_max)/) !Indices of layer midpoints in the sediments
    k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/) !Indices of layer interfaces in the sediments (including the SWI)
    k_bbl1 = (/(k,k=k_wat_bbl+1,k_bbl_sed)/) !Indices of layer interfaces in the BBL (including the top)

    write(*,*) "Sed 9"

    !!Calculate physical forcings
    !Assume (t, s, cc, hmix_rate) at layer midpoints; (kz, w_b) on interfaces (as in GOTM and ROMS grids):
    !Note: kz vertical index starts from 1, not 0 as in e.g. GOTM, ROMS
    !      Using index starting from 0 leads to array misalignment passing between subroutines (PWA, 11/03/2016)
    !
    !========= (air-sea interface) kz(1), w_b(1) (unused)
    !    o     t(1), s(1), cc(1), hmix_rate(1)
    !--------- kz(2), w_b(2) = 0
    !    o     t(2), s(2), cc(2), hmix_rate(2)
    !    :
    !    :
    !    o     t(k_wat_bbl), s(k_wat_bbl), cc(k_wat_bbl), hmix_rate(k_wat_bbl)
    !========= (water-bbl interface) kz(k_wat_bbl+1) = kz_bbl or kz_bbl_max
    !    o     t(k_wat_bbl+1), s(k_wat_bbl+1), cc(k_wat_bbl+1), hmix_rate(k_wat_bbl+1) = 0
    !--------- kz(k_wat_bbl+2) = kz_bbl  (if kz_bbl_type = 0)
    !    o     t(k_wat_bbl+2), s(k_wat_bbl+2), cc(k_wat_bbl+2), hmix_rate(k_wat_bbl+2) = 0
    !    :
    !    :
    !    o     t(k_bbl_sed), s(k_bbl_sed), cc(k_bbl_sed), hmix_rate(k_bbl_sed) = 0
    !========= (bbl-sediment interface) kz(k_bbl_sed+1) = 0, w_b(k_bbl_sed+1)
    !    o     t(k_bbl_sed+1), s(k_bbl_sed+1), cc(k_bbl_sed+1), hmix_rate(k_bbl_sed+1) = 0
    !--------- kz(k_bbl_sed+2) = 0, w_b(k_bbl_sed+2)
    !    o     t(k_bbl_sed+2), s(k_bbl_sed+2), cc(k_bbl_sed+2), hmix_rate(k_bbl_sed+2) = 0
    !    :
    !    :
    !    o     t(k_max), s(k_max), cc(k_max), hmix_rate(k_max) = 0
    !========= (bottom) kz(k_max+1), w_b(k_max+1)
    kz = 0.0_rk
    hmix_rate = 0.0_rk
    hmix_rate = 0.0_rk
    cc_hmix = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    iday = 1
    
    write(*,*) "Sed 10.1"
    do istep=1,steps_in_yr 

        if (mod(istep,int(86400/input_step)).eq.0) then
		    iday = iday + 1
		end if
		write(*,*) "Sed 10.2"
        !Salinity (s)
        s(i,1:k_wat_bbl,istep)       = s_w(i,1:k_wat_bbl,istep)
        s(i,k_wat_bbl+1:k_max,istep) = s_w(i,k_wat_bbl,istep) !Assume constant below z_wat_bbl
        write(*,*) "Sed 10.3"
        !Temperature (t)
        t(i,1:k_wat_bbl,istep)       = t_w(i,1:k_wat_bbl,istep)
        t(i,k_wat_bbl+1:k_max,istep) = t_w(i,k_wat_bbl,istep) !Assume constant below z_wat_bbl
        write(*,*) "Sed 10.4"
        !Horizontal advection (u_x)
        u_x(i,1:k_wat_bbl,istep)       = u_x_w(i,1:k_wat_bbl,istep)
        u_x(i,k_wat_bbl+1:k_max,istep) = 0.0_rk !Assume zero below z_wat_bbl

        write(*,*) "Sed 10.5"
        !Vertical diffusivity in water column (kz)
        kz(i,1:k_wat_bbl,istep)      = kz_w(i,1:k_wat_bbl,istep) !Use all values on upper layer interfaces in water column

        if (dynamic_kz_bbl.eq.0) then !Static kz_bbl
            if (kz_bbl_type.eq.0) then !Constant kz across BBL
                do k=k_wat_bbl+1,k_bbl_sed
                    if (z1(k) < (z_bbl_sed-dbl_thickness)) then !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
                        kz(i,k,istep) = kz_bbl_max
                    else
                        kz(i,k,istep) = 0.0_rk !eddy diffusivity kz is assumed to be zero within the diffusive boundary layer
                    end if
                end do
            end if
            if (kz_bbl_type.eq.1) then !Linear kz across BBL (~=> log-layer for velocity, Holtappels & Lorke, 2011)
                kz_gr = (0.0_rk-kz_bbl_max) / (z_bbl_sed-dbl_thickness-z_wat_bbl)
                !Note: kz is assumed to reach zero at height dbl_thickness above the SWI
                do k=k_wat_bbl+1,k_bbl_sed
                    kz(i,k,istep) = max(0.0_rk, kz_bbl_max + kz_gr * (z1(k)-z_wat_bbl)) !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
                end do
            end if
        end if
        write(*,*) "Sed 10.6"
        if (dynamic_kz_bbl.eq.1) then !Dynamic kz_bbl
            kz_gr = (0.0_rk-kz(i,k_wat_bbl,istep)) / (z_bbl_sed-dbl_thickness-z1(k_wat_bbl))
            do k=k_wat_bbl+1,k_bbl_sed
                kz(i,k,istep) = max(0.0_rk, kz(i,k_wat_bbl,istep) + kz_gr * (z1(k)-z1(k_wat_bbl))) !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
            end do
        end if
        write(*,*) "Sed 10.7"
        !Warning if the diffusive CFL condition is > 0.5
        kzCFL(:,istep) = (kz(i,2:k_bbl_sed,istep)*dt/freq_turb)/(dz(1:k_bbl_sed-1)**2)
        if (diff_method.eq.0.and.maxval(kzCFL(:,istep)).gt.0.5_rk) then
            write(*,*) "WARNING!!! CFL condition due to eddy diffusivity exceeds 0.5 on day", istep
            write(*,*) "z_L, kz, kzCFL = "
            do k=1,k_bbl_sed-1
                if (kzCFL(k,istep).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kz(i,k+1,istep), kzCFL(k,istep)
            end do
        end if
        write(*,*) "Sed 10.8"
        !Horizontal mixing rate (hmix_rate) and reservoir concentrations (cc_hmix)
        hmix_rate(i,1:k_wat_bbl,istep) = 0.0_rk 
        cc_hmix(i,1:par_max,1:k_wat_bbl,iday) = 0.0_rk !cc_hmix_w(i,1:par_max,1:k_wat_bbl,iday)
    end do
  enddo
  write(*,*) "Sed 11"
    !Porosity (phi) (assumed constant in time)
    phi = 1.0_rk
    phi_inv = 1.0_rk/phi
    wat_content = 1.0_rk
    phi1 = 1.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    phi(i,k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
    dphidz_SWI = -1.0_rk * (phi_0-phi_inf) / z_decay_phi
    !water content (wat_content) (assumed constant in time)
    wat_content(i,k_sed) = wat_con_inf + (wat_con_0-wat_con_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
    !Porosity on layer interfaces (phi1)
    phi1(i,k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-0.5_rk*hz(k_sed)-z_bbl_sed)/z_decay_phi)
    phi1(i,k_max+1) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_max)+0.5_rk*hz(k_max)-z_bbl_sed)/z_decay_phi)
  enddo
  write(*,*) "Sed 12"
    !Porosity factors used in diffusivity calculations (pF1, pF2)
    !(assumed constant in time but will vary between solutes vs. solids)
    !These allow us to use a single equation to model diffusivity updates in the water column and sediments, for both solutes and solids:
    ! dC/dt = d/dz(pF2*kzti*d/dz(pF1*C)) where C has units [mass per unit total volume (water+sediments)]
    pF1 = 1.0_rk
    pF2 = 1.0_rk
    pWC = 1.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    do ip=1,par_max
        if (is_solid(ip).eq.0) then !Factors for solutes
            pF1(i,k_sed,ip) = 1.0_rk/phi(i,k_sed)            !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments
            pF2(i,k_sed1,ip) = phi1(i,k_sed1)                !Porosity-related area restriction factor for fluxes across layer interfaces
            pWC(i,k_sed,ip) = 1.0_rk/wat_content(i,k_sed)    !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments (for analytical data)
            ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
            ! dC/dt = d/dz(phi*kzti*d/dz(C/phi))                     in the sediments
        end if
        if (is_solid(ip).eq.1) then !Factors for solids
            pF1(i,k_sed,ip) = 1.0_rk/(1.0_rk - phi(i,k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments
            pF2(i,k_sed1,ip) = 1.0_rk - phi1(i,k_sed1)       !Porosity-related area restriction factor for fluxes across layer interfaces
            pWC(i,k_sed,ip) = 1.0_rk/(1.0_rk - wat_content(i,k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments (for analytical data)
            ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
            ! dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi)))             in the sediments
        end if
    end do
    !Tortuosity on layer interfaces (following Boudreau 1996, equatoin 4.120)
    tortuosity(i,:) = sqrt(1.0_rk - 2.0_rk*log(phi1(i,:)))
  enddo
  
  write(*,*) "Sed 13"

    kz_mol = 0.0_rk
    kz_bio = 0.0_rk
    alpha = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    !Vertical diffusivity due to effective molecular diffusivity of solutes in the sediments (kz_mol)
    !(assumed constant in time but in general may vary between solutes)
    do ip=1,par_max
        write(*,*) "Sed 14.1"
        if (is_solid(ip).eq.0) then
            kz_mol(i,1:k_max+1,ip) = kz_mol0
            kz_mol(i,k_sed1,ip) = mu0_musw * kz_mol0/(tortuosity(i,k_sed1)**2)

            !Warning if the diffusive CFL condition is > 0.5
            kz_molCFL(:,ip) = (kz_mol(i,2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
            if (diff_method.eq.0.and.maxval(kz_molCFL(:,ip)).gt.0.5_rk) then
                write(*,*) "WARNING!!! CFL condition due to molecular diffusivity exceeds 0.5 for variable", trim(par_name(ip))
                write(*,*) "z_L, kz_mol, kz_molCFL = "
                do k=1,k_max-1
                    write(*,*) "Sed 14.2"
                    if (kz_molCFL(k,ip).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kz_mol(i,k+1,ip), kz_molCFL(k,ip)
                end do
            end if
        end if
    end do
    write(*,*) "Sed 15"
    !Vertical diffusivity due to bioturbation (kz_bio)
    !(assumed constant in time and between solutes/solids)
    do k=k_bbl_sed+(2-bioturb_across_SWI),k_max+1  !Note: bioturbation diffusivity is assumed to be non-zero on the SWI
        if (z_s1(k)<z_const_bioturb) then
            kz_bio(i,k) = kz_bioturb_max
        else
            kz_bio(i,k) = kz_bioturb_max*exp(-1.0_rk*(z_s1(k)-z_const_bioturb)/z_decay_bioturb)
        end if
    end do
    write(*,*) "Sed 15.1"
    !Bioirrigation rate (alpha)
    !(assumed constant in time and between solutes)
    alpha(i,k_sed) = a1_bioirr*exp(-1.0_rk*a2_bioirr*(z(k_sed)-z_bbl_sed)) !Schluter et al. (2000), Eqn. 2
  enddo
  write(*,*) "Sed 16"
    !Background vertical advective velocities of particulates and solutes on layer interfaces in the sediments (w_b, u_b)
    !(these assume steady state compaction and neglect reaction terms)
    w_b = 0.0_rk
    u_b = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    w_b(i,k_sed1) = ((1.0_rk - phi_inf)/(1.0_rk - phi1(i,k_sed1))) * w_binf !Boudreau (1997), Eqn 3.67; Holzbecher (2002) Eqn 3
    u_b(i,k_sed1) = (phi_inf/phi1(i,k_sed1)) * w_binf !Boudreau (1997), Eqn 3.68; Holzbecher (2002) Eqn 12
    !do k=k_bbl_sed-1,k_max
    !  w_b(i,:) =  w_binf
    !  u_b(i,:) =  w_binf
    !enddo
  enddo 
  write(*,*) "Sed 17"
    !!Write to output file
    open(12,FILE='Hydrophysics.dat')
    write(12,'(6hiday  ,5hk    ,5hi    ,11hhz[m]      ,11hz[m]       ,11hz_H[m]     ,9ht[degC]  ,10hs[psu]    ,16hKz_H[m2/s]      ,18hKz_mol1_H[m2/s]   , 18hKz_bio_H[m2/s]    ,15halpha[/s]      ,5hphi  ,14htortuosity_H  ,15hhmix_rate[/d]  ,17hw_bH [1E-10 m/s]  ,17hu_bH [1E-10 m/s]  )')
        iday=1
        do istep=1,steps_in_yr
            if (mod(istep,int(86400/input_step)).eq.0) then
		        iday = iday + 1
		    end if
            do k=1,k_max
                write (12,'(2(1x,i4))',advance='NO') istep, k
                do i=1,i_max 
                    write(12,'(1x,i4,f10.4,1x,f10.4,1x,f10.4,2(1x,f8.4),1x,f15.11,1x,f17.13,1x,f17.13,1x,f15.11,f7.4,1x,f7.4,12x,f7.4,12x,f7.4,12x,f7.4,12x,f7.4,12x,i5)',advance='NO') &
                        i, hz(k), z(k), (z(k)-0.5_rk*hz(k)),&
                        t(i,k,istep), s(i,k,istep), kz(i,k,istep), kz_mol(i,k,1), &
                        kz_bio(i,k), alpha(i,k), phi(i,k), tortuosity(i,k), &
                        hmix_rate(i,k,istep), 1.E10*w_b(i,k), 1.E10*u_b(i,k), &
                        u_x(i,k,istep),istep
                end do
                write(12,*)
            end do
        end do
    close(12)
    write(*,*) "Sed 18"

    end subroutine make_physics_bbl_sed
!=======================================================================================================================






!=======================================================================================================================
    subroutine make_physics_bbl_sed_old(t, s, kz, hmix_rate, cc_hmix, u_x, u_x_w, t_w, s_w, kz_w, cc_hmix_w, kz_mol, kz_bio, &
        z, dz, hz, i_min, i_max, k_wat_bbl, k_bbl_sed, k_max, par_max, steps_in_yr, alpha, is_solid, phi, phi1, phi_inv, &
        pF1, pF2, mu0_musw, tortuosity, w_b, u_b, rho, dt, freq_turb, par_name, diff_method, bioturb_across_SWI, &
        ip_sol, ip_par, dphidz_SWI, wat_content,pWC,input_step)

    !Constructs profiles for (temperature, salinity, vertical diffusivity etc.) for the full column (water + sediments)
    !given water column variability and parameters from brom.yaml

    !Input variables
    integer, intent(in)                       :: i_min, i_max, k_wat_bbl, k_bbl_sed, k_max, par_max, steps_in_yr, freq_turb
    integer, intent(in)                       :: diff_method, bioturb_across_SWI, ip_sol, ip_par, input_step
    integer, dimension(:), intent(in)         :: is_solid
    real(rk), intent(in)                      :: dt
    real(rk), dimension(:), intent(in)        :: z, dz, hz
    real(rk), dimension(:,:,:), intent(in)    :: t_w, s_w, kz_w, u_x_w
    real(rk), dimension(:,:,:,:), intent(in)  :: cc_hmix_w
    character(len=attribute_length), dimension(:), intent(in) :: par_name

    !Output variables
    real(rk)                                  :: mu0_musw, dphidz_SWI
    real(rk), dimension(:), intent(out)       :: rho
    real(rk), dimension(:,:), intent(out)     :: kz_bio, alpha, phi, phi_inv, tortuosity, w_b, u_b, wat_content
    real(rk), dimension(:,:,:), intent(out)   :: t, s, kz, u_x, hmix_rate, kz_mol, pF1, pF2, pWC
    real(rk), dimension(:,:,:,:), intent(out) :: cc_hmix

    !Local variables
    integer                                   :: i, k, ip, sel, iday, istep
    integer                                   :: k_sed(k_max-k_bbl_sed), k_sed1(k_max+1-k_bbl_sed), k_bbl1(k_bbl_sed-k_wat_bbl)
    real(rk)                                  :: z_wat_bbl, z_bbl_sed, kz_gr, z1(k_max+1), z_s1(k_max+1), phi1(i_max,k_max+1)
    integer                                   :: kz_bbl_type, dynamic_kz_bbl
    real(rk)                                  :: kz_bbl_max, dbl_thickness, kz_mol0
    real(rk)                                  :: a1_bioirr, a2_bioirr
    real(rk)                                  :: kz_bioturb_max, z_const_bioturb, z_decay_bioturb
    real(rk)                                  :: phi_0, phi_inf, z_decay_phi, w_binf, rho_def, wat_con_0, wat_con_inf
    real(rk)                                  :: kzCFL(k_bbl_sed-1,steps_in_yr), kz_molCFL(k_max-1,par_max)



    !!Get diffusivity parameters from brom.yaml
    !Turbulence in the benthic boundary layer
    kz_bbl_type = get_brom_par("kz_bbl_type")
    kz_bbl_max = get_brom_par("kz_bbl_max")
    dbl_thickness = get_brom_par("dbl_thickness")
    dynamic_kz_bbl = get_brom_par("dynamic_kz_bbl")

    !Molecular diffusivity of solutes (single constant value, infinite dilution)
    kz_mol0 = get_brom_par("kz_mol0")
    mu0_musw = get_brom_par("mu0_musw")

    !Bioturbation
    kz_bioturb_max = get_brom_par("kz_bioturb_max")
    z_const_bioturb = get_brom_par("z_const_bioturb")
    z_decay_bioturb = get_brom_par("z_decay_bioturb")

    !Bioirrigation
    a1_bioirr = get_brom_par("a1_bioirr")
    a2_bioirr = get_brom_par("a2_bioirr")

    !Porosity
    phi_0 = get_brom_par("phi_0")
    phi_inf = get_brom_par("phi_inf")
    z_decay_phi = get_brom_par("z_decay_phi")
    wat_con_0 = get_brom_par("wat_con_0")
    wat_con_inf = get_brom_par("wat_con_inf")

    !Vertical advection in the sediments
    w_binf = get_brom_par("w_binf")

    !Density of particles
    rho_def = get_brom_par("rho_def")
    rho = 0.0_rk
    do ip=1,par_max
        if (is_solid(ip).eq.1) then
            rho(ip) = get_brom_par('rho_' // trim(par_name(ip)),rho_def)
            write(*,*) "Assumed density of ", trim(par_name(ip)), " = ", rho(ip)
        end if
    end do


    !!Set useful parameters for calculations
    !Useful depths
    z_wat_bbl = z(k_wat_bbl+1) - 0.5_rk*hz(k_wat_bbl+1) !Depth of water-BBL interface
    z_bbl_sed = z(k_bbl_sed+1) - 0.5_rk*hz(k_bbl_sed+1) !Depth of BBL-sediment interface (SWI)
    z1(1:k_max) = z(:)-0.5_rk*hz(:)         !Depth of layer interfaces
    z1(k_max+1) = z(k_max)+0.5_rk*hz(k_max)
    z_s1 = z1-z_bbl_sed                     !Depth of interfaces wrt SWI

    !Useful index vectors
    k_sed = (/(k,k=k_bbl_sed+1,k_max)/) !Indices of layer midpoints in the sediments
    k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/) !Indices of layer interfaces in the sediments (including the SWI)
    k_bbl1 = (/(k,k=k_wat_bbl+1,k_bbl_sed)/) !Indices of layer interfaces in the BBL (including the top)



    !!Calculate physical forcings
    !Assume (t, s, cc, hmix_rate) at layer midpoints; (kz, w_b) on interfaces (as in GOTM and ROMS grids):
    !Note: kz vertical index starts from 1, not 0 as in e.g. GOTM, ROMS
    !      Using index starting from 0 leads to array misalignment passing between subroutines (PWA, 11/03/2016)
    !
    !========= (air-sea interface) kz(1), w_b(1) (unused)
    !    o     t(1), s(1), cc(1), hmix_rate(1)
    !--------- kz(2), w_b(2) = 0
    !    o     t(2), s(2), cc(2), hmix_rate(2)
    !    :
    !    :
    !    o     t(k_wat_bbl), s(k_wat_bbl), cc(k_wat_bbl), hmix_rate(k_wat_bbl)
    !========= (water-bbl interface) kz(k_wat_bbl+1) = kz_bbl or kz_bbl_max
    !    o     t(k_wat_bbl+1), s(k_wat_bbl+1), cc(k_wat_bbl+1), hmix_rate(k_wat_bbl+1) = 0
    !--------- kz(k_wat_bbl+2) = kz_bbl  (if kz_bbl_type = 0)
    !    o     t(k_wat_bbl+2), s(k_wat_bbl+2), cc(k_wat_bbl+2), hmix_rate(k_wat_bbl+2) = 0
    !    :
    !    :
    !    o     t(k_bbl_sed), s(k_bbl_sed), cc(k_bbl_sed), hmix_rate(k_bbl_sed) = 0
    !========= (bbl-sediment interface) kz(k_bbl_sed+1) = 0, w_b(k_bbl_sed+1)
    !    o     t(k_bbl_sed+1), s(k_bbl_sed+1), cc(k_bbl_sed+1), hmix_rate(k_bbl_sed+1) = 0
    !--------- kz(k_bbl_sed+2) = 0, w_b(k_bbl_sed+2)
    !    o     t(k_bbl_sed+2), s(k_bbl_sed+2), cc(k_bbl_sed+2), hmix_rate(k_bbl_sed+2) = 0
    !    :
    !    :
    !    o     t(k_max), s(k_max), cc(k_max), hmix_rate(k_max) = 0
    !========= (bottom) kz(k_max+1), w_b(k_max+1)
    kz = 0.0_rk
    hmix_rate = 0.0_rk
    cc_hmix = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    iday = 1
    do istep=1,steps_in_yr
		
		

        if (mod(istep,int(86400/input_step)).eq.0) then
!		if (mod(istep,24).eq.0) then
		    iday = iday + 1
		end if
		
        !Salinity (s)
        s(i,1:k_wat_bbl,istep)       = s_w(i,1:k_wat_bbl,istep)
        s(i,k_wat_bbl+1:k_max,istep) = s_w(i,k_wat_bbl,istep) !Assume constant below z_wat_bbl

        !Temperature (t)
        t(i,1:k_wat_bbl,istep)       = t_w(i,1:k_wat_bbl,istep)
        t(i,k_wat_bbl+1:k_max,istep) = t_w(i,k_wat_bbl,istep) !Assume constant below z_wat_bbl

        !Horizontal advection (u_x)
        u_x(i,1:k_wat_bbl,istep)       = u_x_w(i,1:k_wat_bbl,istep)
        u_x(i,k_wat_bbl+1:k_max,istep) = 0.0_rk !Assume zero below z_wat_bbl

        !Vertical diffusivity in water column (kz)
        kz(i,1:k_wat_bbl,istep)      = kz_w(i,1:k_wat_bbl,istep) !Use all values on upper layer interfaces in water column

        if (dynamic_kz_bbl.eq.0) then !Static kz_bbl
            if (kz_bbl_type.eq.0) then !Constant kz across BBL
                do k=k_wat_bbl+1,k_bbl_sed
                    if (z1(k) < (z_bbl_sed-dbl_thickness)) then !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
                        kz(i,k,istep) = kz_bbl_max
                    else
                        kz(i,k,istep) = 0.0_rk !eddy diffusivity kz is assumed to be zero within the diffusive boundary layer
                    end if
                end do
            end if
            if (kz_bbl_type.eq.1) then !Linear kz across BBL (~=> log-layer for velocity, Holtappels & Lorke, 2011)
                kz_gr = (0.0_rk-kz_bbl_max) / (z_bbl_sed-dbl_thickness-z_wat_bbl)
                !Note: kz is assumed to reach zero at height dbl_thickness above the SWI
                do k=k_wat_bbl+1,k_bbl_sed
                    kz(i,k,istep) = max(0.0_rk, kz_bbl_max + kz_gr * (z1(k)-z_wat_bbl)) !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
                end do
            end if
        end if

        if (dynamic_kz_bbl.eq.1) then !Dynamic kz_bbl
            kz_gr = (0.0_rk-kz(i,k_wat_bbl,istep)) / (z_bbl_sed-dbl_thickness-z1(k_wat_bbl))
            do k=k_wat_bbl+1,k_bbl_sed
                kz(i,k,istep) = max(0.0_rk, kz(i,k_wat_bbl,istep) + kz_gr * (z1(k)-z1(k_wat_bbl))) !Note that kz(k) is at depth z1(k) = z(k)-hz(k)/2
            end do
        end if

        !Warning if the diffusive CFL condition is > 0.5
        kzCFL(:,istep) = (kz(i,2:k_bbl_sed,istep)*dt/freq_turb)/(dz(1:k_bbl_sed-1)**2)
        if (diff_method.eq.0.and.maxval(kzCFL(:,istep)).gt.0.5_rk) then
            write(*,*) "WARNING!!! CFL condition due to eddy diffusivity exceeds 0.5 on day", istep
            write(*,*) "z_L, kz, kzCFL = "
            do k=1,k_bbl_sed-1
                if (kzCFL(k,istep).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kz(i,k+1,istep), kzCFL(k,istep)
            end do
        end if

        !Horizontal mixing rate (hmix_rate) and reservoir concentrations (cc_hmix)
        hmix_rate(i,1:k_wat_bbl,istep) = 0.0_rk 
        cc_hmix(i,1:par_max,1:k_wat_bbl,iday) = 0.0_rk !cc_hmix_w(i,1:par_max,1:k_wat_bbl,iday)
    end do
  enddo
    !Porosity (phi) (assumed constant in time)
    phi = 1.0_rk
    phi_inv = 1.0_rk/phi
    wat_content = 1.0_rk
    phi1 = 1.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    phi(i,k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
    dphidz_SWI = -1.0_rk * (phi_0-phi_inf) / z_decay_phi
    !water content (wat_content) (assumed constant in time)
    wat_content(i,k_sed) = wat_con_inf + (wat_con_0-wat_con_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
    !Porosity on layer interfaces (phi1)
    phi1(i,k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-0.5_rk*hz(k_sed)-z_bbl_sed)/z_decay_phi)
    phi1(i,k_max+1) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_max)+0.5_rk*hz(k_max)-z_bbl_sed)/z_decay_phi)
  enddo
    !Porosity factors used in diffusivity calculations (pF1, pF2)
    !(assumed constant in time but will vary between solutes vs. solids)
    !These allow us to use a single equation to model diffusivity updates in the water column and sediments, for both solutes and solids:
    ! dC/dt = d/dz(pF2*kzti*d/dz(pF1*C)) where C has units [mass per unit total volume (water+sediments)]
    pF1 = 1.0_rk
    pF2 = 1.0_rk
    pWC = 1.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    do ip=1,par_max
        if (is_solid(ip).eq.0) then !Factors for solutes
            pF1(i,k_sed,ip) = 1.0_rk/phi(i,k_sed)            !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments
            pF2(i,k_sed1,ip) = phi1(i,k_sed1)                !Porosity-related area restriction factor for fluxes across layer interfaces
            pWC(i,k_sed,ip) = 1.0_rk/wat_content(i,k_sed)    !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments (for analytical data)
            ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
            ! dC/dt = d/dz(phi*kzti*d/dz(C/phi))                     in the sediments
        end if
        if (is_solid(ip).eq.1) then !Factors for solids
            pF1(i,k_sed,ip) = 1.0_rk/(1.0_rk - phi(i,k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments
            pF2(i,k_sed1,ip) = 1.0_rk - phi1(i,k_sed1)       !Porosity-related area restriction factor for fluxes across layer interfaces
            pWC(i,k_sed,ip) = 1.0_rk/(1.0_rk - wat_content(i,k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments (for analytical data)
            ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
            ! dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi)))             in the sediments
        end if
    end do
    !Tortuosity on layer interfaces (following Boudreau 1996, equatoin 4.120)
    tortuosity(i,:) = sqrt(1.0_rk - 2.0_rk*log(phi1(i,:)))
  enddo

    kz_mol = 0.0_rk
    kz_bio = 0.0_rk
    alpha = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    !Vertical diffusivity due to effective molecular diffusivity of solutes in the sediments (kz_mol)
    !(assumed constant in time but in general may vary between solutes)
    do ip=1,par_max
        if (is_solid(ip).eq.0) then
            kz_mol(i,1:k_max+1,ip) = kz_mol0
            kz_mol(i,k_sed1,ip) = mu0_musw * kz_mol0/(tortuosity(i,k_sed1)**2)

            !Warning if the diffusive CFL condition is > 0.5
            kz_molCFL(:,ip) = (kz_mol(i,2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
            if (diff_method.eq.0.and.maxval(kz_molCFL(:,ip)).gt.0.5_rk) then
                write(*,*) "WARNING!!! CFL condition due to molecular diffusivity exceeds 0.5 for variable", trim(par_name(ip))
                write(*,*) "z_L, kz_mol, kz_molCFL = "
                do k=1,k_max-1
                    if (kz_molCFL(k,ip).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kz_mol(i,k+1,ip), kz_molCFL(k,ip)
                end do
            end if
        end if
    end do
    !Vertical diffusivity due to bioturbation (kz_bio)
    !(assumed constant in time and between solutes/solids)
    do k=k_bbl_sed+(2-bioturb_across_SWI),k_max+1  !Note: bioturbation diffusivity is assumed to be non-zero on the SWI
        if (z_s1(k)<z_const_bioturb) then
            kz_bio(i,k) = kz_bioturb_max
        else
            kz_bio(i,k) = kz_bioturb_max*exp(-1.0_rk*(z_s1(k)-z_const_bioturb)/z_decay_bioturb)
        end if
    end do

    !Bioirrigation rate (alpha)
    !(assumed constant in time and between solutes)
    alpha(i,k_sed) = a1_bioirr*exp(-1.0_rk*a2_bioirr*(z(k_sed)-z_bbl_sed)) !Schluter et al. (2000), Eqn. 2
  enddo
    !Background vertical advective velocities of particulates and solutes on layer interfaces in the sediments (w_b, u_b)
    !(these assume steady state compaction and neglect reaction terms)
    w_b = 0.0_rk
    u_b = 0.0_rk
    !Calculation below are for water column horizontal index
  do i = i_min, i_max
    w_b(i,k_sed1) = ((1.0_rk - phi_inf)/(1.0_rk - phi1(i,k_sed1))) * w_binf !Boudreau (1997), Eqn 3.67; Holzbecher (2002) Eqn 3
    u_b(i,k_sed1) = (phi_inf/phi1(i,k_sed1)) * w_binf !Boudreau (1997), Eqn 3.68; Holzbecher (2002) Eqn 12
    !do k=k_bbl_sed-1,k_max
    !  w_b(i,:) =  w_binf
    !  u_b(i,:) =  w_binf
    !enddo
  enddo 

    !!Write to output file
    open(12,FILE='Hydrophysics.dat')
    write(12,'(6hiday  ,5hk    ,5hi    ,11hhz[m]      ,11hz[m]       ,11hz_H[m]     ,9ht[degC]  ,10hs[psu]    ,16hKz_H[m2/s]      ,18hKz_mol1_H[m2/s]   , 18hKz_bio_H[m2/s]    ,15halpha[/s]      ,5hphi  ,14htortuosity_H  ,15hhmix_rate[/d]  ,17hw_bH [1E-10 m/s]  ,17hu_bH [1E-10 m/s]  )')
        iday=1
        do istep=1,steps_in_yr
            if (mod(istep,int(86400/input_step)).eq.0) then
!            if (mod(istep,24).eq.0) then
		        iday = iday + 1
		    end if
            do k=1,k_max
                write (12,'(2(1x,i4))',advance='NO') istep, k
                do i=1,i_max 
                    write(12,'(1x,i4,f10.4,1x,f10.4,1x,f10.4,2(1x,f8.4),1x,f15.11,1x,f17.13,1x,f17.13,1x,f15.11,f7.4,1x,f7.4,12x,f7.4,12x,f7.4,12x,f7.4,12x,f7.4,12x,i5)',advance='NO') &
                        i, hz(k), z(k), (z(k)-0.5_rk*hz(k)),&
                        t(i,k,istep), s(i,k,istep), kz(i,k,istep), kz_mol(i,k,1), &
                        kz_bio(i,k), alpha(i,k), phi(i,k), tortuosity(i,k), &
                        hmix_rate(i,k,istep), 1.E10*w_b(i,k), 1.E10*u_b(i,k), &
                        u_x(i,k,istep),istep
                end do
                write(12,*)
            end do
        end do
    close(12)


    end subroutine make_physics_bbl_sed_old
!=======================================================================================================================




!=======================================================================================================================
    subroutine svan(s, t, po, sigma)
!/*
!c------specific volume anomaly based on 1980 equation of state for
!c      seawater and 1978 practical salinity scale
!c      pressure          PO     decibars
!c      temperature        T     degree celsius (IPTS-68)
!c      salinitty          S     (PSS-78)
!c      spec.vol.anom.  SVAN     1.0E-8 m**3/Kg
!c      density anom.   SIGMA    Kg/m**3
!*/

    real(rk)  p,sig,sr,r1,r2,r3,s,t,po,sigma
    real(rk)  a,b,c,d,e,a1,b1,aw,bw,k,ko,kw,k35,v350p,sva,dk,gam,pk,dr35p,dvan

    real(rk) r3500, r4 ,dr350
    data  r3500 /1028.1063/, r4/4.8314E-4/,dr350/28.106331/

    p=po/10.
    sr=sqrt(abs(s))

    r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t &
           -9.095290E-3)*t+6.793952E-2)*t-28.263737
    r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t &
          +8.24493E-1
    r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3

    sig=(r4*s + r3*sr + r2)*s +r1

    v350p=1.0/r3500
    sva=-sig*v350p/(r3500+sig)
    sigma= sig + dr350

    if (p.eq.0.0) return

    e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
       bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
    b = bw + e*s

    d= 1.91075E-4
    c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
    aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
    a = (d*sr + c)*s + aw

    b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
    a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
    kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t &
           +148.4206)*t-1930.06
    ko = (b1*sr + a1)*s + kw

    dk = (b*p+a)*p+ko
    k35 = (5.03217E-5*p+3.359406)*p+21582.27
    gam=p/k35
    pk=1.0-gam
    sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))

    v350p= v350p*pk
    dr35p=gam/v350p
    dvan= sva/(v350p*(v350p+sva))
    sigma = dr350 + dr35p -dvan
    return
    end subroutine svan
!=======================================================================================================================


    end module io_ascii
