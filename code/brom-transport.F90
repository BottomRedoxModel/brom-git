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
!                     Anfisa Berezina 
!-----------------------------------------------------------------------

    module brom_transport

    ! AB check fabm subroutines
    use fabm
    use fabm_config
    use fabm_types, only: attribute_length, rk ! check attribute_length
    use fabm_standard_variables ! T, S, light, CO2 concentration etc
    use io_netcdf
    use io_ascii 
    use mtridiagonal, only: init_tridiagonal,clean_tridiagonal
    use ids         !Provides access to variable indices id_O2 etc.


    implicit none
    private
    public init_brom_transport, do_brom_transport, clear_brom_transport
	
	real(rk), parameter :: pi=3.141592653589793_rk

    !FABM model with all data and procedures related to biogeochemistry
    type (type_model) :: model
    type (type_horizontal_standard_variable),parameter :: &
        ice_thickness = type_horizontal_standard_variable(name=' ice_thickness',units='m') ! horizontal - 2D
    type (type_bulk_standard_variable),parameter :: & ! check if used!!
        volume_of_cell = type_bulk_standard_variable(name=' volume_of_cell')


    !Solution parameters to be read from brom.yaml (see brom.yaml for details)
    integer   :: i_min, i_max       !x-axis related 
    integer   :: k_min, k_wat_bbl, k_bbl_sed !z-axis related
    integer   :: k_points_below_water, k_max !z-axis related
    integer   :: par_max                     !no. BROM variables
    integer   :: i_day, year, days_in_yr, freq_turb, freq_sed, freq_float, last_day   !time related ! ?? freq_sed, freq_turb
    integer   :: diff_method, kz_bbl_type, bioturb_across_SWI  !vertical diffusivity related
    integer   :: h_adv, h_relax, h_turb  !horizontal transport  (advection and relaxation) switches
    integer   :: use_Eair, use_hice, use_gargett ! use input for light, ice, calculate Kz
    integer   :: input_type, port_initial_state, ncoutfile_type !I/O related
    integer   :: use_leak, i_leak, start_leak !leak related
    real(rk)  :: dt, water_layer_thickness, cc_leak, cc_leak2, w_leak_adv
    real(rk)  :: K_O2s, gargett_a0, gargett_q, mult_Kz
    
    ! Free timestep input and output
    integer   :: ist, i_step, input_step, steps_in_yr, output_step ! ist - steps in the course of the day, i_step - step in the course of the year

    character(len=64) :: icfile_name, outfile_name, ncoutfile_name
    character :: hmix_file

    !Forcings to be provided to FABM: These must have the POINTER attribute
    real(rk), pointer, dimension(:)            :: pco2atm, windspeed, hice
    real(rk), pointer, dimension(:,:)          :: surf_flux, bott_flux, bott_source, Izt, pressure, depth
    real(rk), pointer, dimension(:,:,:)        :: t, s, u_x
    real(rk), pointer, dimension(:,:,:)        :: vv, dvv, cc, cc_out, dcc, dcc_R, wbio ! add the description cc - all params, dcc - volumes of solids
    ! vv -> 2d

    !Surface and bottom forcings, used within brom-transport only
    real(rk), allocatable, dimension(:)      :: Eair
    real(rk), allocatable, dimension(:,:,:)    :: cc_top, cc_bottom

    !Horizontal mixing forcings, used within brom-transport only
    real(rk), allocatable, dimension(:,:,:)    :: kl
    real(rk), allocatable, dimension(:,:,:,:)  :: cc_hmix ! relaxation

    !Grid parameters and forcings for water column only
    real(rk), allocatable, dimension(:)        :: z_w, dz_w, hz_w
    real(rk), allocatable, dimension(:,:,:)    :: t_w, s_w, kz_w, u_x_w
    real(rk), allocatable, dimension(:,:,:,:)  :: cc_hmix_w

    !Grid parameters and forcings for full column including water and sediments
    real(rk), allocatable, dimension(:)        :: z, dz, hz, x, dx, hx, z1, z_s1
    real(rk), allocatable, dimension(:,:)      :: bc_top, bc_bottom, kz_bio, alpha, phi, phi1, phi_inv, tortuosity, kztCFL
    real(rk), allocatable, dimension(:,:)      :: wCFL, w_b, u_b, wat_content
    real(rk), allocatable, dimension(:,:,:)    :: kz, kzti, fick, fick_per_day, sink, sink_per_day, bcpar_top, bcpar_bottom
    real(rk), allocatable, dimension(:,:,:)    :: kz_mol, pF1, pF2, wti, w_bub, pWC
    integer, allocatable, dimension(:)         :: is_solid, is_gas, k_wat, k_sed, k_sed1, k_bbl1
    integer, allocatable, dimension(:,:)       :: bctype_top, bctype_bottom, hmixtype
    integer, allocatable, dimension(:,:)       :: kzCFL, kz_molCFL
    character(len=attribute_length), allocatable, dimension(:)    :: par_name

            ! variables for building the grid
    integer                                   :: sel, iday, istep
!    integer                                   :: k_sed(k_max-k_bbl_sed), k_sed1(k_max+1-k_bbl_sed), k_bbl1(k_bbl_sed-k_wat_bbl)
    real(rk)                                  :: z_wat_bbl, z_bbl_sed, kz_gr !, z1(k_max+1), z_s1(k_max+1), phi1(i_max,k_max+1)
    integer                                   :: dynamic_kz_bbl
    real(rk)                                  :: kz_bbl_max, dbl_thickness, kz_mol0
    real(rk)                                  :: a1_bioirr, a2_bioirr
    real(rk)                                  :: kz_bioturb_max, z_const_bioturb, z_decay_bioturb
    real(rk)                                  :: phi_0, phi_inf, z_decay_phi, w_binf, rho_def, wat_con_0, wat_con_inf
!    real(rk)                                  :: kzCFL(k_bbl_sed-1,steps_in_yr), kz_molCFL(k_max-1,par_max)
    
    
    !Constant forcings that can be read as parameters from brom.yaml
    real(rk) :: wind_speed, pco2_atm, mu0_musw, dphidz_SWI , dx_min, dx_max ! dx_adv
    real(rk), allocatable, dimension(:)        :: rho

    !Counters
    integer   :: i, k, ip, ip_sol, ip_par, kj, ij !, i_dummy AB
    real(rk)  :: i_dummy ! AB


    contains

!=======================================================================================================================
!  
!  name: brom_transport.init_brom_transport
!  @param
!  @return
!  
    subroutine init_brom_transport()

    !Initialises the offline vertical transport model BROM-transport

    use ids         !Provides access to variable indices id_O2 etc

    implicit none


    !Reading brom.yaml
    call init_common()

    !Get grid and numerical solution parameters from from brom.yaml
    dt = get_brom_par("dt")
    freq_turb = get_brom_par("freq_turb")
    freq_sed  = get_brom_par("freq_sed ")
    freq_float  = get_brom_par("freq_float ")
    last_day = get_brom_par("last_day")
    water_layer_thickness = get_brom_par("water_layer_thickness")
    k_min = get_brom_par("k_min")
    k_wat_bbl = get_brom_par("k_wat_bbl")
    k_points_below_water = get_brom_par("k_points_below_water")
    i_min = get_brom_par("i_min") 
    i_max = get_brom_par("i_max")
    year = get_brom_par("year")
    days_in_yr = get_brom_par("days_in_yr")
    
    ! for free length output (assumed to be a day fraction)
    input_step = get_brom_par("input_step")
    output_step = get_brom_par("output_step")
    
    diff_method = get_brom_par("diff_method")
    bioturb_across_SWI = get_brom_par("bioturb_across_SWI")
    input_type = get_brom_par("input_type")
    use_Eair = get_brom_par("use_Eair")
    use_hice = get_brom_par("use_hice")
    port_initial_state = get_brom_par("port_initial_state")
    icfile_name = get_brom_name("icfile_name")
    outfile_name = get_brom_name("outfile_name")
    ncoutfile_name = get_brom_name("ncoutfile_name")
    ncoutfile_type = get_brom_par("ncoutfile_type")
    K_O2s = get_brom_par("K_O2s")
    h_adv =  get_brom_par("h_adv")
    h_relax =  get_brom_par("h_relax")
    h_turb =  get_brom_par("h_turb")
!    dx_adv = get_brom_par("dx_adv")
	
	! vertical grid params
    
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

	! horizontal grid params
	dx_min = get_brom_par("dx_min")
	dx_max = get_brom_par("dx_max")
    use_gargett = get_brom_par("use_gargett")
    gargett_a0 = get_brom_par("gargett_a0")
    gargett_q = get_brom_par("gargett_q")
    mult_Kz = get_brom_par("mult_Kz")
    use_leak = get_brom_par("use_leak")
    i_leak = get_brom_par("i_leak")
    start_leak = get_brom_par("start_leak")
    w_leak_adv = get_brom_par("w_leak_adv")
    cc_leak = get_brom_par("cc_leak")
    cc_leak2 = get_brom_par("cc_leak2")
    !Initialize FABM model from fabm.yaml
    call fabm_create_model_from_yaml_file(model)
    par_max = size(model%state_variables)

    steps_in_yr = days_in_yr*24*3600/input_step ! determine how much timesteps in the course of the year

    !Allocate biological variables now that par_max is known
    allocate(surf_flux(i_max,par_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bott_flux(i_max,par_max))     !bottom flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bott_source(i_max,k_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bc_top(i_max,par_max))
    allocate(bc_bottom(i_max,par_max))
    allocate(bctype_top(i_max,par_max))
    allocate(bctype_bottom(i_max,par_max))
    allocate(bcpar_top(i_max,par_max,3))
    allocate(bcpar_bottom(i_max,par_max,3))
    allocate(par_name(par_max))
    allocate(pco2atm(i_max))
    allocate(windspeed(i_max))
    allocate(cc_top(i_max,par_max,days_in_yr))
    allocate(cc_bottom(i_max,par_max,days_in_yr))
    allocate(is_solid(par_max))
    allocate(is_gas(par_max))
    allocate(hmixtype(i_max,par_max))
    allocate(rho(par_max))


    !Retrieve the parameter names from the model structure
    do ip=1,par_max
        par_name(ip) = model%state_variables(ip)%name
    end do
    !Make the named parameter indices id_O2 etc.
    call get_ids(par_name)


    !Get boudary condition parameters from brom.yaml:
    !bctype = 0, 1, 2, 3 for no flux (default), Dirichlet constant, Dirichlet sinusoid, and Dirichlet netcdf input respectively
    do ip=1,par_max
        bctype_top(i_max,ip) = get_brom_par('bctype_top_' // trim(par_name(ip)),0.0_rk)
        if (bctype_top(i_max,ip).eq.1) then
            bc_top(i_max,ip) = get_brom_par('bc_top_' // trim(par_name(ip)))
            write(*,*) "Constant Dirichlet upper boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3)') " = ", bc_top(i_max,ip)
        else if (bctype_top(i_max,ip).eq.2) then     !Model: bc_top = a1top + a2top*sin(omega*(julianday-a3top))
            bcpar_top(i_max,ip,1) = get_brom_par('a1top_' // trim(par_name(ip)))
            bcpar_top(i_max,ip,2) = get_brom_par('a2top_' // trim(par_name(ip)))
            bcpar_top(i_max,ip,3) = get_brom_par('a3top_' // trim(par_name(ip)))
            write(*,*) "Sinusoidal Dirichlet upper boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_top(i_max,ip,1), " + ", &
                  bcpar_top(i_max,ip,2), "*sin(omega*(julianday -", bcpar_top(i_max,ip,3), "))"
        else if (bctype_top(i_max,ip).eq.3) then     !Read from netcdf
            write(*,*) "NetCDF specified Dirichlet upper boundary condition for " // trim(par_name(ip))
        else if (bctype_top(i_max,ip).eq.4) then     !Read from ascii or calc from calintiiy
            write(*,*) "Upper boundary condition from NODC " // trim(par_name(ip))
        end if        
        bctype_bottom(i_max,ip) = get_brom_par('bctype_bottom_' // trim(par_name(ip)),0.0_rk)
        if (bctype_bottom(i_max,ip).eq.1) then
            bc_bottom(i_max,ip) = get_brom_par('bc_bottom_' // trim(par_name(ip)))
            write(*,*) "Constant Dirichlet lower boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3)') " = ", bc_bottom(i_max,ip)
        else if (bctype_bottom(i_max,ip).eq.2) then  !Model: bc_bottom = a1bottom + a2bottom*sin(omega*julianday-a3bottom))
            bcpar_bottom(i_max,ip,1) = get_brom_par('a1bottom_' // trim(par_name(ip)))
            bcpar_bottom(i_max,ip,2) = get_brom_par('a2bottom_' // trim(par_name(ip)))
            bcpar_bottom(i_max,ip,3) = get_brom_par('a3bottom_' // trim(par_name(ip)))
            write(*,*) "Sinusoidal Dirichlet lower boundary condition for " // trim(par_name(ip))
            write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_bottom(i_max,ip,1), " + ", &
            bcpar_bottom(i_max,ip,2), "*sin(omega*(julianday -", bcpar_bottom(i_max,ip,3), "))"
        else if (bctype_bottom(i_max,ip).eq.3) then  !Read from netcdf
            write(*,*) "NetCDF specified Dirichlet lower boundary condition for " // trim(par_name(ip))
        end if
        bctype_top(:,ip) = bctype_top(i_max,ip)
        bc_top(:,ip) = bc_top(i_max,ip)
        bctype_bottom(:,ip) = bctype_bottom(i_max,ip)
        bc_bottom(:,ip) = bc_bottom(i_max,ip)
        do k=1,3
          bcpar_top(:,ip,k) = bcpar_top(i_max,ip,k)
          bcpar_bottom(:,ip,k) = bcpar_bottom(i_max,ip,k)
        enddo
    end do
    
    write(*,*) "All other boundary conditions use surface and bottom fluxes from FABM"

    
    !Input forcing data
    if (input_type.eq.0) then !Input sinusoidal seasonal changes (hypothetical)
        call input_primitive_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_max, days_in_yr)
        allocate(Eair(days_in_yr))
        allocate(hice(days_in_yr))
        Eair = 0.0_rk
        hice = 0.0_rk
        write(*,*) "Done sinusoidal input"
    end if
    if (input_type.eq.1) then !Input physics from ascii
        call input_ascii_physics(z_w, dz_w, hz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_max, days_in_yr)
        allocate(Eair(days_in_yr))
        allocate(hice(days_in_yr))
        allocate(cc_hmix_w(i_max,par_max,k_wat_bbl,days_in_yr))
        cc_hmix_w = 0.0_rk
        Eair = 0.0_rk
        hice = 0.0_rk
        write(*,*) "Done ascii input"
    end if
    if (input_type.eq.2) then !Input water column physics from netcdf
        call input_netcdf_2(z_w, dz_w, hz_w, t_w, s_w, kz_w, Eair, use_Eair, &
        hice, use_hice, gargett_a0, gargett_q, use_gargett, &
        year, i_max, steps_in_yr, k_wat_bbl, u_x_w)
        kz_w=kz_w*mult_Kz
        write(*,*) "Done netcdf input"
        !Note: This uses the netCDF file to set z_w = layer midpoints, dz_w = increments between layer midpoints, hz_w = layer thicknesses
    end if
 
    !Determine total number of vertical grid points (layers) now that k_wat_bbl is determined
    k_max = k_wat_bbl + k_points_below_water

    !Allocate full grid variables now that k_max is knownk
    allocate(z(k_max))
    allocate(dz(k_max))
    allocate(hz(k_max))
    allocate(x(i_max))
    allocate(dx(i_max))
    allocate(hx(i_max))
    allocate(t(i_max,k_max,steps_in_yr))
    allocate(s(i_max,k_max,steps_in_yr))
    allocate(kz(i_max,k_max+1,steps_in_yr))
    allocate(kl(i_max,k_max,steps_in_yr))
    allocate(u_x(i_max,k_max,steps_in_yr))
    allocate(cc_hmix(i_max,par_max,k_max,days_in_yr))
    allocate(kz_mol(i_max,k_max+1,par_max))
    allocate(kz_bio(i_max,k_max+1))
    allocate(pF1(i_max,k_max,par_max))
    allocate(pF2(i_max,k_max+1,par_max))
    allocate(pWC(i_max,k_max+1,par_max)) ! water content??
    allocate(alpha(i_max,k_max))
    allocate(phi(i_max,k_max))
    allocate(wat_content(i_max,k_max))
    allocate(phi1(i_max,k_max+1))
    allocate(phi_inv(i_max,k_max))
    allocate(tortuosity(i_max,k_max+1))
    allocate(w_b(i_max,k_max+1))
    allocate(u_b(i_max,k_max+1))
    allocate(wti(i_max,k_max+1,par_max)) ! vertical velocity
    allocate(w_bub(i_max,k_max+1,par_max))
    allocate(cc(i_max,k_max,par_max))
    allocate(cc_out(i_max,k_max,par_max))
    allocate(dcc(i_max,k_max,par_max))
    allocate(dcc_R(i_max,k_max,par_max))
    allocate(fick(i_max,k_max+1,par_max))
    allocate(fick_per_day(i_max,k_max+1,par_max))
    allocate(wbio(i_max,k_max,par_max))    !vertical velocity (m/s, negative for sinking)
    allocate(sink(i_max,k_max+1,par_max))  !sinking flux (mmol/m2/s, positive downward)
    allocate(sink_per_day(i_max,k_max+1,par_max))
    allocate(vv(i_max,k_max,1))
    allocate(dvv(i_max,k_max,1))
    allocate(Izt(i_max,k_max))
    allocate(pressure(i_max,k_max))
    allocate(depth(i_max,k_max))
    allocate(kzti(i_max,k_max+1,par_max))
    allocate(kztCFL(k_max-1,par_max))
    allocate(wCFL(k_max-1,par_max))
!    allocate(k_sed(k_max-k_bbl_sed))
!    allocate(k_sed1(k_max+1-k_bbl_sed))
    allocate(k_bbl1(k_bbl_sed-k_wat_bbl))
    allocate(z1(k_max+1))
    allocate(z_s1(k_max+1))
    allocate(kzCFL(k_bbl_sed-1,steps_in_yr))
    allocate(kz_molCFL(k_max-1,par_max))

    if (k_points_below_water==0) then  !This is to "unlock" BBL and sediments for a "classical" water column model
        k_max=k_wat_bbl
        z=z_w
        dz=dz_w
        hz=hz_w
        k_bbl_sed=k_wat_bbl !needed for Irradiance calculations
    else
        !Construct the full vertical grid
        call make_vert_grid(z, dz, hz, z_w, dz_w, hz_w, k_wat_bbl, k_max, k_bbl_sed)
        write(*,*) "Made vertical grid"
        allocate(k_wat(k_bbl_sed))
        allocate(k_sed(k_max-k_bbl_sed))
        allocate(k_sed1(k_max+1-k_bbl_sed))
        k_wat = (/(k,k=1,k_bbl_sed)/)       !Index vector for all points in the water column
        k_sed = (/(k,k=k_bbl_sed+1,k_max)/) !Index vector for all points in the sediments
        k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/) !Indices of layer interfaces in the sediments (including the SWI)
    endif

    !Specify horizontal transport
    !dx(:)= dx_adv !horizontal resolution in m
    write(*,*) dx_min, dx_max
    if (i_max.gt.1) then
        do i=1,i_max 
            dx(i) = dx_min + (dx_max-dx_min)*(cos((i)*2*pi/(i_max+1))+1)*0.5
            write(*,*) i, dx(i)
        enddo
    endif
    if (h_adv.eq.1.or.h_turb.eq.1) then
        x(1)=0.0_rk
        !u_x_w(:,:,:)=0.02  !horizontal velocity in m/s if not read from forcing file...
        do i=2,i_max            !works when we have more than 2 columns
            x(i)=x(i-1)+dx(i)
        enddo
    else
        !allocate u_x_w
        x=0.0_rk
        u_x_w=0.0_rk
    endif
    !Initialize tridiagonal matrix if necessary
    if (diff_method.gt.0) then
        call init_tridiagonal(k_max)
        write(*,*) "Initialized tridiagonal matrix"
    end if

    !Set model domain
    call fabm_set_domain(model,i_max,k_max)

    !Specify vertical index of surface and bottom
    call model%set_surface_index(k_min)
    call model%set_bottom_index(k_max)

    !Initial volumes of layers:
    vv(i:i_max,1:k_max,1) = 1.0_rk

    !Make they (full) pressure variable to pass to FABM
    !This is used by brom_eqconst.F90 to compute equilibrium constants for pH calculations
    !and by brom_carb.F90 to compute equilibrium constants for saturation states (subroutine CARFIN)
    do i=1,i_max
        pressure(i,:) = z(:) + 10.0_rk
    end do

    !Point FABM to array slices with biogeochemical state.
    do ip=1,par_max
        call fabm_link_bulk_state_data(model, ip, cc(:,:,ip))
    end do


    !Link temperature and salinity data to FABM (needs to be redone every time julianday is updated below)
    call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,1))
    call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,1))
    !Link other data needed by FABM
    call fabm_link_bulk_data(model, standard_variables%downwelling_photosynthetic_radiative_flux, Izt)  !W m-2
    call fabm_link_bulk_data(model, standard_variables%pressure, pressure)                              !dbar
    call fabm_link_bulk_data(model, standard_variables%depth, depth)                                    !m
    call fabm_link_horizontal_data(model, standard_variables%wind_speed, windspeed)                     !m s-1
    call fabm_link_horizontal_data(model, standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2atm)  !ppm
    if (use_hice.eq.1) call fabm_link_horizontal_data(model, ice_thickness, hice)
    call fabm_link_bulk_data(model, volume_of_cell, vv(:,:,1))


    !Check FABM is ready
    call fabm_check_ready(model)


    !Allow FABM models to use their default initialization (this sets cc)
    do k=1,k_max
        call fabm_initialize_state(model, 1, i_max, k)
    end do


    !Read initial values from ascii file if req'd
    if (port_initial_state.eq.1) call porting_initial_state_variables(trim(icfile_name), year, &
                                                i_day, i_max, k_max, par_max, par_name, cc, vv)
                                                
    if (port_initial_state.eq.2) then 
       call porting_initial_state_variables(trim(icfile_name), year, &
                                           i_day, 1, k_max, par_max, par_name, cc, vv)
       do i=2,i_max
          cc(i,:,:) = cc(1,:,:)
          vv(i,:,:) = vv(1,:,:) 
       enddo
    endif

    !Check biological parameters for non-zero values
        do ip=1,par_max
            if (ip.eq.id_Phy.or.ip.eq.id_Het.or.ip.eq.id_Baae.or.ip.eq.id_Baan.or.ip.eq.id_Bhae.or.ip.eq.id_Bhan) then 
                do i=1,i_max
                    do k=1,k_max
                        if(cc(i,k,ip).le.0.0_rk) cc(i,k,ip)= 1.0E-11
                    enddo
                enddo
            endif
        enddo 
    

    !Initialize output
    call init_netcdf(trim(ncoutfile_name), i_max, k_max, model, use_Eair, use_hice, year)

    !Establish which variables will be treated as solid phase in the sediments, based on the biological velocity (sinking/floating) from FABM.
    wbio = 0.0_rk
    call fabm_get_vertical_movement(model, i_max, i_max, k_wat_bbl, wbio(i_max:i_max,k_wat_bbl,:))
    do i=i_min,i_max-1
        wbio(i,:,:)=wbio(i_max,:,:)
    enddo
    wbio = -1.0_rk * wbio !FABM returns NEGATIVE wbio for sinking; sign change here means that wbio is POSITIVE for sinking
    is_solid = 0
    is_gas = 0
    ip_sol = 0
    ip_par = 0
    write(*,*) "The following variables are assumed to join the solid phase in the sediments"
    do ip=1,par_max
      do i=i_min,i_max
        if (wbio(i,k_wat_bbl,ip).gt.0.0_rk) then
            is_solid(ip) = 1 !Any variables that SINKS (wbio>0) in the bottom cell of the water column will become "solid" in the sediments
!            write(*,*) trim(par_name(ip))
            if (ip_par.eq.0) ip_par = ip !Set ip_par = index of first particulate variable
        elseif (wbio(i,k_wat_bbl,ip).lt.0.0_rk) then
            is_gas(ip) = 1 !Any variables that floats (wbio>1) in the bottom cell of the water column will become "solid" in the sediments
!            write(*,*) trim(par_name(ip))           
        else
            if (ip_sol.eq.0) ip_sol = ip !Set ip_par = index of first solute variable
        end if
      enddo
    end do


    !Density of particles
    rho_def = get_brom_par("rho_def")
    rho = 0.0_rk
    do ip=1,par_max
        if (is_solid(ip).eq.1) then
            rho(ip) = get_brom_par('rho_' // trim(par_name(ip)),rho_def)
            write(*,*) "Assumed density of ", trim(par_name(ip)), " = ", rho(ip)
        end if
    end do


    !Complete hydrophysical forcings
    cc_hmix=0.0_rk
    
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
!    hmix_rate = 0.0_rk
!    cc_hmix = 0.0_rk
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
        !Linear u_x across BBL (to be checked!)
                !kz_gr = (0.0_rk-kz_bbl_max) / (z_bbl_sed-dbl_thickness-z_wat_bbl)
                !Note: u_x is assumed to reach zero at height dbl_thickness above the SWI
                do k=k_wat_bbl+1,k_bbl_sed-1
                    !u_x(i,k,istep) = max(0.0_rk, u_x_w(i,k_wat_bbl,istep) + kz_gr * (z1(k)-z_wat_bbl)) 
                    ! velosity profile (Cushman-Roisin, Beckers, 2011) U(Z)=(U*/K)*ln(Z-Z0)
                    ! were Z-distance from bottom, Z0-roughness, K=0.4 von Karman const, U*- friction velosity.
                    ! 0.07 here is for (U*/K) normalization  
                    u_x(i,k,istep) = (1.0_rk/log((z(k_bbl_sed+1)-z(k_wat_bbl+1))/dbl_thickness)) &
                        *u_x_w(i,k_wat_bbl,istep)*log((z(k_bbl_sed+1)-z(k))/dbl_thickness)
                    !u_x(i,k,istep) =0.0_rk
                end do
!        u_x(i,k_wat_bbl+1:k_max,istep) = 0.0_rk !Assume zero below z_wat_bbl
        u_x(i,k_bbl_sed+1:k_max,istep) = 0.0_rk !Assume zero below z_wat_bbl

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
!        hmix_rate(i,1:k_wat_bbl,istep) = 0.0_rk 
!        cc_hmix(i,1:par_max,1:k_wat_bbl,iday) = 0.0_rk !cc_hmix_w(i,1:par_max,1:k_wat_bbl,iday)
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
    !Calculation below are for water column horizontal index
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
    write(12,'(6hiday  ,5hk    ,5hi    ,11hhz[m]      ,11hz[m]       ,11hz_H[m]     ,9ht[degC]  ,10hs[psu]    ,16hKz_H[m2/s]      ,18hKz_mol1_H[m2/s]   , 18hKz_bio_H[m2/s]    , 18hu[m/s]            , 18hstep              ,15halpha[/s]      ,5hphi  ,14htortuosity_H  ,17hw_bH [1E-10 m/s]  ,17hu_bH [1E-10 m/s]  )')
        iday=1
        do istep=1,steps_in_yr
            if (mod(istep,int(86400/input_step)).eq.0) then
!            if (mod(istep,24).eq.0) then
		        iday = iday + 1
		    end if
            do k=1,k_max
                write (12,'(2(1x,i4))',advance='NO') istep, k
                do i=1,i_max 
                    write(12,'(1x,i4,f10.4,1x,f10.4,1x,f10.4,2(1x,f8.4),1x,f15.11,1x,f17.13,1x,f17.13,12x,f7.4,12x,i5,1x,f15.11,f7.4,1x,f7.4,12x,f7.4,12x,f7.4)',advance='NO') &
                        i, hz(k), z(k), (z(k)-0.5_rk*hz(k)),&
                        t(i,k,istep), s(i,k,istep), kz(i,k,istep), kz_mol(i,k,1), &
                        kz_bio(i,k), u_x(i,k,istep), istep, &
                        alpha(i,k), phi(i,k), tortuosity(i,k), &
                        1.E10*w_b(i,k), 1.E10*u_b(i,k) 
                end do
                write(12,*)
            end do
        end do
    close(12)

    write(*,*) "Made physics of BBL, sediments"

    !Get horizontal relaxation parameters from brom.yaml:
    !hmixtype = 0, 1 or 2  for no horizontal relaxation (default), box model mixing respectively
        do ip=1,par_max
            hmixtype(i_max,ip) = get_brom_par('hmix_' // trim(par_name(ip)),0.0_rk)
            hmixtype(:,ip)=hmixtype(i_max,ip)   
                        
            if (hmixtype(i_max,ip).eq.1) then
                write(*,*) "Horizontal relaxation assumed for " // trim(par_name(ip))            
            end if
            if (hmixtype(i_max,ip).eq.2) then
    !            hmix_file = get_brom_name("hmix_filename_" // trim(par_name(ip)(13:))) !niva_oxydep_NUT
                write(*,*) "Horizontal relaxation (ASCII) assumed for " // trim(par_name(ip))
            !else if (hmixtype(i_max,ip).eq.2) then  ! read relaxation files in two cases, also for top bondary condition  #bctype_top(i_max,ip).eq.4.or.  
                open(20, file= get_brom_name("hmix_filename_" // trim(par_name(ip))))!'' // hmix_file
                write(*,*) ("hmix_filename_" // trim(par_name(ip)))
                do k=1,k_wat_bbl
                    do i_day=1,days_in_yr
                        read(20, *) i_dummy,i_dummy,cc_hmix(i_max,ip,k,i_day) ! NODC data (i_max,par_max,k_max,days_in_yr)
                    end do
                end do 
                close(20)
                do i=i_min,i_max-1
                    cc_hmix(i,:,:,:)=cc_hmix(i_max,:,:,:)
                enddo
               
            else if (bctype_top(i_max,ip).eq.4.) then    
                open(40, file = get_brom_name("bc_top_filename_" // trim(par_name(ip)))) !to boundary condition file  
                do i_day=1,days_in_yr
                    read(40, *) cc_top(i,ip,i_day)
                end do   
                close(40)
                do i=i_min,i_max-1
                    cc_top(i,:,:)=cc_top(i_max,:,:)
                enddo
            end if
        end do
    
    open(8,FILE = 'burying_rate.dat')
    !Set constant forcings
    wind_speed = get_brom_par("wind_speed")    ! 10m wind speed [m s-1]
    pco2_atm   = get_brom_par("pco2_atm")      ! CO2 partical pressure [ppm]
    pco2atm(:) = pco2_atm                      ! Same pCO2 values for all the surface gridpoints
    windspeed(:) = wind_speed                  ! Same wind speed values to all the surface gridpoints

    end subroutine init_brom_transport
!=======================================================================================================================






!=======================================================================================================================
    subroutine do_brom_transport()

    !Executes the offline vertical transport model BROM-transport

    use calculate, only: calculate_light, calculate_phys, calculate_sed, calculate_sed_eya, calculate_bubble

    implicit none

    integer      :: id, idt, idf                     !time related
    integer      :: surf_flux_with_diff              !1 to include surface fluxes in diffusion update, 0 to include in bgc update
    integer      :: bott_flux_with_diff              !1 to include bottom fluxes in diffusion update, 0 to include in bgc update
    integer      :: model_w_sed                      !1 to assume porosity effects for solutes and solids, 0 - sumplified approach
    integer      :: constant_w_sed                   !1 to assume constant burial (advection) velocities in the sediments
    integer      :: dynamic_w_sed                    !1 to assume dynamic burial (advection) velocities in the sediments depending on dVV(k_bbl_sed)
    integer      :: show_maxmin, show_kztCFL, show_wCFL, show_nan, show_nan_kztCFL, show_nan_wCFL     !options for runtime output to screen
    integer      :: bc_units_convert, sediments_units_convert !options for conversion of concentrations units in the sediment 
    integer      :: julianday, model_year
    
    integer      :: ist, istep ! AB
    integer      :: i_count_pr, i_count_s, i_sec_pr
    
    integer      :: k_inj,i_inj,inj_swith,inj_num,start_inj,stop_inj    !#number of layer and column to inject into, start day, stop day number 
    integer      :: start_inj2,stop_inj2,start_inj3,stop_inj3    !#number of layer and column to inject into, start day, stop day number 
    real(rk)     :: cnpar                            !"Implicitness" parameter for GOTM vertical diffusion (set in brom.yaml)
    real(rk)     :: cc0                              !Resilient concentration (same for all variables)
    real(rk)     :: omega                            !angular frequency of sinusoidal forcing = 2*pi/365 rads/day
    real(rk)     :: O2stat                           !oxygen status of sediments (factor modulating the bioirrigation rate)
    real(rk)     :: a1_bioirr                        !to detect whether or not bioirrigation is activated
    real(rk)     :: kl_unifrom                !uniform horizontal relaxation if prescribed
    real(rk)     :: tau_relax                        ! relaxation time (AB)
    real(rk)     :: fresh_PM_poros                   ! porosity of fresh precipitated PM (i.e. dVV) 
    real(rk)     :: w_binf                   ! baseline rate of burying
    real(rk)     :: w_bub0                   ! baseline rate of floating up of bubbles
    real(rk)     :: bu_co                   ! "Burial coeficient" for setting velosity exactly to the SWI proportional to the 
                                           !   settling velocity in the water column (0<bu_co<1), 0 - for no setting velosity, (nd)
    real(rk)     :: N_bubbles              ! Number of rising up bubbles
    real(rk)     :: injection_rate                   ! injection rate
    real(rk)     :: injection_rate_ini               ! injection rate initial constant or multiplier if injection rate is a function
    character(len=attribute_length), allocatable, dimension(:)    :: inj_var_name

    omega = 2.0_rk*pi/365.0_rk

    !Get parameters for the time-stepping and vertical diffusion / sedimentation
    cnpar = get_brom_par("cnpar")
    model_w_sed = get_brom_par("model_w_sed")
    dynamic_w_sed = get_brom_par("dynamic_w_sed")
    constant_w_sed = get_brom_par("constant_w_sed")
    fresh_PM_poros = get_brom_par("fresh_PM_poros")
    w_binf = get_brom_par("w_binf")
    N_bubbles = get_brom_par("N_bubbles")
    w_bub0 = get_brom_par("w_bub0")
    bu_co = get_brom_par("bu_co")
    cc0 = get_brom_par("cc0")
    a1_bioirr = get_brom_par("a1_bioirr")
    surf_flux_with_diff = get_brom_par("surf_flux_with_diff")
    bott_flux_with_diff = get_brom_par("bott_flux_with_diff")
    show_maxmin = get_brom_par("show_maxmin")
    show_kztCFL = get_brom_par("show_kztCFL")
    show_wCFL = get_brom_par("show_wCFL")
    show_nan = get_brom_par("show_nan")
    show_nan_kztCFL = get_brom_par("show_nan_kztCFL")
    show_nan_wCFL = get_brom_par("show_nan_wCFL")
    bc_units_convert = get_brom_par("bc_units_convert")
    sediments_units_convert = get_brom_par("sediments_units_convert")
    kl_unifrom = get_brom_par("kl_unifrom")
    tau_relax = get_brom_par("tau_relax")
    injection_rate_ini = get_brom_par("injection_rate_ini") 
    k_inj = get_brom_par("k_injection") 
    i_inj = get_brom_par("i_injection") 
    inj_swith = get_brom_par("injection_swith")
    start_inj = get_brom_par("start_inj")
    stop_inj = get_brom_par("stop_inj")    
    start_inj2 = get_brom_par("start_inj2")
    stop_inj2 = get_brom_par("stop_inj2")    
    start_inj3 = get_brom_par("start_inj3")
    stop_inj3 = get_brom_par("stop_inj3")
    !i_sec_pr = get_brom_par("i_sec_pr")    
    idt = int(86400._rk/dt)                                      !number of cycles per day
    model_year = 0
    kzti = 0.0_rk
    sink=0.0_rk
    wti=0.0_rk
    w_bub=0.0_rk
    dVV = 0.0_rk

        !convert bottom boundary values from 'mass/pore water ml' for dissolved and 'mass/mass' for solids into 'mass/total volume'
    if (bc_units_convert.eq.1) then
        do i=i_min,i_max
            do ip=1,par_max
                if (bctype_bottom(i,ip).eq.1) bc_bottom(i,ip) = bc_bottom(i,ip)/pF1(i,k_max,ip)
            enddo
        enddo    
    end if

    !Uniform horizontal relaxation if prescribed
    if(kl_unifrom>0.0_rk)    kl=kl_unifrom
    !Master time step loop over days
    write(*,*) "Starting time stepping"
    do i=1, i_max
        cc(i,:,:)=cc(1,:,:)
    enddo
    
 !_______BIG Cycle ("i_day"=0,...,last_day-1)________!
	ist=1
	istep=1
  do i_day=0,(last_day-1)

        julianday = i_day - int(i_day/days_in_yr)*days_in_yr + 1    !"julianday" (1,2,...,days_in_yr)
        if (julianday==1) model_year = model_year + 1

        write (*,'(a, i4, a, i4, a, f9.4)') " model year:", model_year, "; julianday:", julianday,"; w_sed (cm/yr):", wti(1,k_bbl_sed+2,1)*365.*8640000.
        !If including ice using horizontal coordinate, set ice volume in top cell of ice column
            !vv = 0.0_rk
            !if (i_max.eq.2) then
            !if (use_hice.eq.1) vv(1,k_min,1) = max(0.0_rk, hice(julianday))   ! i.e. vv() is an amount of solid matter in the cell
        !end if

!###########################################################################################    
    do id=1,idt !in the course of 1 day
        !Note: The numerical approach here is Operator Splitting with tracer transport processes assumed to be
        !numerically more demanding than the biogeochemistry (hence freq_turb, freq_sed >= 1) (Butenschon et al., 2012)


        !Set time-varying Dirichlet boundary conditions for current julianday
        do ip=1,par_max
            do i=i_min, i_max
                !Sinusoidal variations
                if (bctype_top(i,ip).eq.2) bc_top(i,ip) = bcpar_top(i,ip,1) + &
                    bcpar_top(i,ip,2)*sin(omega*(julianday-bcpar_top(i,ip,3)))
                if (bctype_bottom(i,ip).eq.2) bc_bottom(i,ip) = bcpar_bottom(i,ip,1) + &
                    bcpar_bottom(i,ip,2)*sin(omega*(julianday-bcpar_bottom(i,ip,3)))
                
                !Variations read from netcdf
                if (bctype_top(i,ip).eq.3) bc_top(i,ip) = cc_top(i,ip,julianday)
                if (bctype_bottom(i,ip).eq.3) bc_bottom(i,ip) = cc_bottom(i,ip,julianday)
                
                !Variations read from ascii file and/or calculated as a function of something
                if (bctype_top(i,ip).eq.4) then                                
                    bc_top(i,ip) = cc_hmix(i,ip,1,julianday)
                end if 
                
                !SO4 in mmol/m3, SO4/Salt from Morris, A.W. and Riley, J.P.(1966) quoted in Dickson et al.(2007)              
                if (bctype_top(i,ip).eq.5) bc_top(i,ip)=(0.1400_rk/96.062_rk)*(s(1,1,julianday)/1.80655_rk)*1.e6_rk !.and.ip.eq.id_SO4
                !if (bctype_top(i,ip).eq.4.and.ip.eq.id_Alk) bc_top(i,ip)=0.068*s(1,1,julianday)  
                !         !Alk in mmol/m3, Alk/Salt from Murray, 2014
            enddo
        enddo        
        
        
      !Calculate Izt = <PAR(z)>_24hr for this day
      do i=i_min, i_max
          call calculate_light(julianday, i, k_bbl_sed, k_max, par_max, hz, &
                               Eair, use_Eair, hice, use_hice, cc, is_solid, rho, Izt)
      enddo

      if (mod(int(id*dt),input_step).eq.0) then  
        ist = int((id*int(dt))/input_step)
        istep = int((julianday-1)*(24*3600/input_step)) + ist
        
        ! Reload changes in t and s (needed for FABM/biology)
        call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,istep))
        call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,istep))
      endif
      
	  i_count_pr=0

        !_______vertical diffusion________!
            do i=i_min, i_max
                call calculate_phys(i, k_max, par_max, model, cc, kzti, fick, &
                    dcc, bctype_top, bctype_bottom, bc_top, bc_bottom, &
                    surf_flux, bott_flux, bott_source, k_bbl_sed, dz, hz, kz, &
                    kz_mol, kz_bio, julianday, id_O2, K_O2s, dt, freq_turb, &
                    diff_method, cnpar, surf_flux_with_diff,bott_flux_with_diff, &
                    bioturb_across_SWI, pF1, pF2, phi_inv, is_solid, cc0)
            enddo

        !_______bioirrigation_____________!
            if (a1_bioirr.gt.0.0_rk) then
                do i=i_min, i_max
                    dcc = 0.0_rk
                    !Oxygen status of sediments set by O2 level just above sediment surface
 !                   O2stat = cc(i,k_bbl_sed,id_O2) / (cc(i,k_bbl_sed,id_O2) + K_O2s)   
                    O2stat = cc(i,k_bbl_sed,id_Oxy) / (cc(i,k_bbl_sed,id_Oxy) + K_O2s)   
                    do ip=1,par_max
                        if (is_solid(ip).eq.0) then
                            !Calculate tendencies dcc
                        
                            !Schluter et al. (2000), Meile et al. (2001)
                            !Note use of factor pF1 = 1/phi to convert cc from [mass per unit total volume]
                            !to [mass per unit volume pore water]                                                
                            dcc(i,k_sed,ip) = O2stat*alpha(i,k_sed)*phi(i,k_sed) * &
                                (cc(i,k_bbl_sed,ip) - pF1(i,k_sed,ip)*cc(i,k_sed,ip))

                            !Bottom cell of water column receives -1 * sum of all exchange fluxes (conservation of mass)                       
                            dcc(i,k_bbl_sed,ip) = -1.0_rk * sum(dcc(i,k_sed,ip)*hz(k_sed)) / hz(k_bbl_sed) 
                            !Bottom cell of water column receives -1 * sum of all exchange fluxes (conservation of mass)
                        
                            !Update concentrations
                            !Simple Euler time step for all k in the sediments (index vector k_sed)
                            cc(i,k_sed,ip) = cc(i,k_sed,ip) + dt*dcc(i,k_sed,ip)                       
                            cc(i,k_bbl_sed,ip) = cc(i,k_bbl_sed,ip) + dt*dcc(i,k_bbl_sed,ip)
                            if (bctype_bottom(i,ip).gt.0) cc(i,k_max,ip) = bc_bottom(i,ip) !Reassert Dirichlet BC if required
                            cc(i,k_sed,ip) = max(cc0, cc(i,k_sed,ip)) !Impose resilient concentration
                        end if
                    end do
                enddo
            end if

            !_____water_biogeochemistry_______!
            dcc = 0.0_rk
            do i=i_min, i_max
                do k=1,k_max
                call fabm_do(model, i, i, k, dcc(i:i,k,:))   !to ask Jorn
                !Note: We MUST pass the range "i:i" to fabm_do -- a single value "i" will produce compiler error
            end do
            !Add surface and bottom fluxes if treated here
            if (surf_flux_with_diff.eq.0) then
                surf_flux = 0.0_rk
                call fabm_do_surface(model, i, i, surf_flux(i:i,:))
                fick(i,k_min,:) = surf_flux(i,:)
                do ip=1,par_max
                    dcc(i,k_min,ip) = dcc(i,k_min,ip) + surf_flux(i,ip) / hz(k_min)
                end do
            end if
            
            if (bott_flux_with_diff.eq.0) then
                bott_flux = 0.0_rk
                bott_source = 0.0_rk                
                call fabm_do_bottom(model, i, i, bott_flux(i:i,:),bott_source(i:i,:))
                sink(i,k_max+1,:) = bott_flux(i,:)
                do ip=1,par_max
                    dcc(i,k_max,ip) = dcc(i,k_max,ip) + bott_flux(i,ip) / hz(k_max)
                end do
            end if
            if (dynamic_w_sed.eq.1) dcc_R(i,:,:) = dcc(i,:,:) !Record biological reaction terms for use in calculate_sed
            !Euler time step due to FABM biogeochemistry
            do k=1,k_max
                do ip=1,par_max
                    cc(i,k,ip) = cc(i,k,ip) + dt*dcc(i,k,ip)
                end do
            end do
            !Reassert Dirichlet BCs
            do ip=1,par_max
                if (bctype_top(i,ip).gt.0) then
                    cc(i,i_min,ip) = bc_top(i,ip)
                end if
                if (bctype_bottom(i,ip).gt.0) then
                    cc(i,k_max,ip) = bc_bottom(i,ip)
                end if
            end do
            cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration
            enddo

        !_______Particles sinking_________!
        do i=i_min, i_max
            if(model_w_sed.ge.1) then
                call calculate_sed(i, k_max, par_max, model, cc, wti, sink, &
                dcc, dcc_R, bctype_top, bctype_bottom, bc_top, bc_bottom, hz, &
                dz, k_bbl_sed, wbio, w_b, u_b, julianday, dt, freq_sed, &
                dynamic_w_sed, is_solid, rho, phi1, fick, k_sed1, K_O2s, &
                kz_bio, id_O2, dphidz_SWI, cc0, bott_flux, bott_source)
            else
                call calculate_sed_eya(i, k_max, par_max, model, cc, wti, &
                sink, dcc, dVV,bctype_top, bctype_bottom, bc_top, &
                bc_bottom, hz, dz, k_bbl_sed, wbio, w_b, u_b, julianday, &
                dt, freq_sed, dynamic_w_sed, constant_w_sed, is_solid, &
                rho, phi1, fick, k_sed1, K_O2s, kz_bio, fresh_PM_poros, &
                id_O2, dphidz_SWI, cc0, bott_flux, bott_source, w_binf, bu_co)
            endif
        enddo

       !_______Calculate changes of volumes of cells_________!
        dVV(:,:,1)=0.0
        do i=i_min, i_max
            k=k_bbl_sed
           do ip=1,par_max !Sum over contributions from each particulate variable
              if (is_solid(ip).eq.1) then
              !change of Volume of a cell as a function of biology, salt precipitation and sinking
                dVV(i,k,1)= dVV(i,k,1)+(sink(i,k-1,ip))/rho(ip)
              end if
           end do
        enddo  
          if (id.eq.1) write (8,'(a, i8,a, i4, a, i4, a, e, a, e,a, e, a, e)') " i_day:", &
              i_day, " model_year:", model_year, "; julianday:", julianday, &
              "; dVV(k_bbl_sed):", dVV(1,k_bbl_sed,1), "; w_sed_bl(cm/yr):", &
              wti(1,k_bbl_sed+2,1)*365.*8640000., "; w_sed_inj(cm/yr):", &
              wti(i_inj,k_bbl_sed+2,1)*365.*8640000., ";sink_of_POML:" , sink(i,k_bbl_sed-1,12)    
        
       !________Horizontal relaxation_________!
        if (h_relax.eq.1) then          
            dcc = 0.0_rk
            do i=i_min, i_max
                do ip=1,par_max
                    if  (hmixtype(i,ip).ge.1) then
                        do k=1,k_wat_bbl
                        !Calculate tendency dcc (water column only) 
                            dcc(i,k,ip) = (cc_hmix(i,ip,k,julianday)-cc(i,k,ip))/tau_relax
                        !Update concentration (water column only)
                            cc(i,k,ip) = cc(i,k,ip) + dt*dcc(i,k,ip)
                        end do
                        cc(i,:,ip) = max(cc0, cc(i,:,ip)) !Impose resilient concentration
                    end if
                end do
            enddo
        endif
       
       !________Horizontal turbulence_________!
        if (h_turb.eq.1.and.i_max.gt.1) then
            dcc = 0.0_rk
            do ip=1,par_max
                do i=i_min+1, i_max-1
                    dcc(i,:,ip) = kl(i,:,istep)*(cc(i+1,:,ip)+cc(i-1,:,ip)-2.0_rk*cc(i,:,ip))/dx(i)/dx(i)
                end do    
                    dcc(i_min,:,ip) = kl(i_min,:,istep)*(cc(i_min+1,:,ip)+cc(i_max,:,ip)-2.0_rk*cc(i_min,:,ip))/dx(i_min)/dx(i_min)
                    dcc(i_max,:,ip) = kl(i_max,:,istep)*(cc(i_min,:,ip)+cc(i_max-1,:,ip)-2.0_rk*cc(i_max,:,ip))/dx(i_max)/dx(i_max)
                do i=i_min, i_max
                    do k=1,k_wat_bbl
                        cc(i,k,ip) = cc(i,k,ip) + dt*dcc(i,k,ip) !Simple Euler time step
                    enddo
                    cc(i,:,ip) = max(cc0, cc(i,:,ip)) !Impose resilient concentration
                enddo        
            enddo          
        end if
       
       !________Horizontal advection_________!
        if (h_adv.eq.1.and.i_max.gt.1) then
            dcc = 0.0_rk
            do ip=1,par_max
                do i=i_min+1, i_max-1
                    do k=1,k_wat_bbl
                      dcc(i,k,ip) = max(0.0_rk, u_x(i,k,istep))*(cc(i-1,k,ip)-cc(i,k,ip))/dx(i)  &
                                  + max(0.0_rk,-u_x(i,k,istep))*(cc(i+1,k,ip)-cc(i,k,ip))/dx(i)  
!                      dcc(i,k,ip) = max(0.0_rk, u_x(i,k,istep))*(cc(i-1,k,ip)-cc(i+1,k,ip))/dx(i)/2._rk  &
!                                  + max(0.0_rk,-u_x(i,k,istep))*(cc(i+1,k,ip)-cc(i-1,k,ip))/dx(i)/2._rk  
                    enddo              
                end do
                    do k=1,k_wat_bbl
                dcc(i_min,k,ip) = max(0.0_rk,u_x(i_min,k,istep))*(cc(i_max,k,ip)-cc(i_min,k,ip))/dx(i_min)  &
                       + max(0.0_rk,-u_x(i_min,k,istep))*(cc(i_min+1,k,ip)-cc(i_min,k,ip))/dx(i_min)
                dcc(i_max,k,ip) = max(0.0_rk,u_x(i_max,k,istep))*(cc(i_max-1,k,ip)-cc(i_max,k,ip))/dx(i_max)  &
                       + max(0.0_rk,-u_x(i_max,k,istep))*(cc(i_min,k,ip)-cc(i_max,k,ip))/dx(i_max)
                    enddo
                do i=i_min, i_max
                    do k=1,k_wat_bbl
                        cc(i,k,ip) = cc(i,k,ip) + dt*dcc(i,k,ip) !Simple Euler time step
                    enddo
                    cc(i,:,ip) = max(cc0, cc(i,:,ip)) !Impose resilient concentration
                enddo        
            enddo          
        endif
        
        !________Vertical advection in the sediments due to leak_________
        if (use_leak.gt.0.and.i_day.gt.start_leak) then          
            dcc = 0.0_rk
            cc(i_leak,k_max,id_DIC)=cc_leak
            cc(i_leak,k_max,id_H2S)=cc_leak2
                do ip=1,par_max
                   if (is_solid(ip).ne.1) then
                        !Calculate tendency dcc (liquids in the sediment  only)
                        do k=k_wat_bbl,k_max-1
                          dcc(i_leak,k,ip) = -w_leak_adv*(cc(i_leak,k,ip)-cc(i_leak,k+1,ip))/hz(k)
                        enddo
                        do k=k_wat_bbl,k_max-1
                          cc(i_leak,k,ip) = cc(i_leak,k,ip) + dt*dcc(i_leak,k,ip)
                        end do
!                        cc(i_leak,:,ip) = max(cc0, cc(i_leak,:,ip)) !Impose resilient concentration
                    endif
                end do
        endif
        
        !________Injection____________________!
        !            !Source of "acetate" 1292 mmol/sec, should be devided to the volume of the grid cell, i.e. dz(k)*dx(i)*dx(i)
        if (i_day.gt.start_inj.and.i_day.le.stop_inj.and.i_max.gt.0) then
            injection_rate = 0.5*max(0., 90.2614*sin(0.0054*(i_day-start_inj))+6.3394*sin(2*0.0054*(i_day-start_inj)) &
                              +10.7102*sin(3*0.0054*(i_day-start_inj))-1.0403*sin(4*0.0054*(i_day-start_inj)) &
                              -22.2653*sin(5*0.0054*(i_day-start_inj))+8.0852*sin(6*0.0054*(i_day-start_inj)) ) &
                              *injection_rate_ini
            if (inj_swith.eq.1)  then
                do ip = 1, par_max
                    if (par_name(ip).eq.get_brom_name("inj_var_name")) exit 
                    inj_num = ip+1           
                end do 

!from 26.02.2020 we use  (after Wang et al., 2012) :			
!Inj O2 / Inj Waste-N [molO2/molN] = 34,1			
!Inj DIC / Inj Waste-N [molC/molN] = 27.3			
!Inj NH4 / Inj Waste-N  [molN/molN] = 2.50 		strict definition	
!Inj DOML (DON) / Inj Waste-N [molN/molN]= 0.80 			
!Inj PO4 / Inj Waste-N  [molN/molN] = 0.0998

              do kj = k_inj-1, k_inj
                do ij= i_inj-2,i_inj+2, 2
                  cc(ij,kj,inj_num)=cc(ij,kj,inj_num)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                  cc(ij,kj,id_O2)=cc(ij,kj,id_O2)    -dt*34.0_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                  cc(ij,kj,id_PO4)=cc(ij,kj,id_PO4)  +dt*0.10_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                  cc(ij,kj,id_NH4)=cc(ij,kj,id_NH4)  +dt*2.5_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                  cc(ij,kj,id_DIC)=cc(ij,kj,id_DIC)  +dt*27.3_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                  cc(ij,kj,id_DOML)=cc(ij,kj,id_DOML)+dt*0.8_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
                enddo
              enddo
              
            end if 
        end if
        
!        if (i_day.gt.start_inj2.and.i_day.le.stop_inj2.and.i_max.gt.0) then
!            injection_rate = 0.5*max(0., 90.2614*sin(0.0054*(i_day-start_inj2))+6.3394*sin(2*0.0054*(i_day-start_inj2)) &
!                              +10.7102*sin(3*0.0054*(i_day-start_inj2))-1.0403*sin(4*0.0054*(i_day-start_inj2)) &
!                              -22.2653*sin(5*0.0054*(i_day-start_inj2))+8.0852*sin(6*0.0054*(i_day-start_inj2)) )&
!                              *injection_rate_ini
!            if (inj_swith.eq.1)  then
!                do ip = 1, par_max
!                    if (par_name(ip).eq.get_brom_name("inj_var_name")) exit 
!                    inj_num = ip+1           
!                end do 
!
!              do kj = k_inj-1, k_inj
!                do ij= i_inj,i_inj !-2,i_inj+2, 2
!                  cc(ij,kj,inj_num)=cc(ij,kj,inj_num)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!!                  cc(ij,kj,id_DOML)=cc(ij,kj,id_DOML)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_O2)=cc(ij,kj,id_O2)  -dt*2._rk*6.625_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_PO4)=cc(ij,kj,id_PO4)+dt*2._rk/16.0_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_NH4)=cc(ij,kj,id_NH4)+dt*2._rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_DIC)=cc(ij,kj,id_DIC)+dt*2._rk*6.625_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                enddo
!              enddo
!
!            end if 
!        end if

!        if (i_day.gt.start_inj3.and.i_day.le.stop_inj3.and.i_max.gt.0) then
!            injection_rate = 0.5*max(0., 90.2614*sin(0.0054*(i_day-start_inj3))+6.3394*sin(2*0.0054*(i_day-start_inj3)) &
!                              +10.7102*sin(3*0.0054*(i_day-start_inj3))-1.0403*sin(4*0.0054*(i_day-start_inj3)) &
!                              -22.2653*sin(5*0.0054*(i_day-start_inj3))+8.0852*sin(6*0.0054*(i_day-start_inj3)) )&
!                              *injection_rate_ini
!            if (inj_swith.eq.1)  then
!                do ip = 1, par_max
!                    if (par_name(ip).eq.get_brom_name("inj_var_name")) exit 
!                    inj_num = ip+1           
!                end do 

!              do kj = k_inj-1, k_inj
!                do ij= i_inj-2,i_inj+2, 2
!                  cc(ij,kj,inj_num)=cc(ij,kj,inj_num)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_DOML)=cc(ij,kj,id_DOML)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_O2)=cc(ij,kj,id_O2)-dt*6.625_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_PO4)=cc(ij,kj,id_PO4)+dt/16.0_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_NH4)=cc(ij,kj,id_NH4)+dt*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                  cc(ij,kj,id_DIC)=cc(ij,kj,id_DIC)+dt*6.625_rk*injection_rate/(dx(ij)*dx(ij)*dz(kj))
!                enddo
!              enddo

!            end if 
!        end if

       !_______Bubble rising_________!
        do i=i_min, i_max
            call calculate_bubble(i, k_max, par_max, model, cc, sink, N_bubbles, &
                dcc, hz, dz, z, t, k_bbl_sed, julianday, dt, freq_float, &
               is_gas, wbio, cc0)
        enddo
        
        !________Check for NaNs (stopping if any found)____________________!
        do ip=1,par_max
            do i=i_min,i_max
                if (any(isnan(cc(i,1:k_max,ip)))) then
                    write(*,*) "Time step within day id = ", id
                    write(*,*) "NaN detected in concentration array, ip = ", ip
                    write(*,*) "Variable name = ", model%state_variables(ip)%name
                    if (show_nan.eq.1) write(*,*) cc(i,1:k_max,ip)
                    if (show_nan_kztCFL.gt.0) then
                        kztCFL(:,ip) = (kzti(i,2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
                        write(*,*) "maxval(kzti(i,:,ip)) = ", maxval(kzti(i,:,ip))
                        write(*,*) "maxval(kztCFL(:,ip)) = ", maxval(kztCFL(:,ip))
                        write(*,*) "fick = ", fick(i,:,ip)
                        if (show_nan_kztCFL.eq.2.and.maxval(kztCFL(:,ip)).gt.0.5_rk) then
                            write(*,*) "(z_L, kz, kztCFL) where kztCFL>0.5 = "
                            do k=1,k_max-1
                                if (kztCFL(k,ip).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kzti(i,k+1,ip), kztCFL(k,ip)
                            end do
                        end if
                    end if
                    if (show_nan_wCFL.gt.0) then
                        wCFL(:,ip) = (abs(wti(i,2:k_max,ip))*dt/freq_sed)/dz(1:k_max-1)
                        write(*,*) "maxval(wti(i,:,ip)) = ", maxval(wti(i,:,ip))
                        write(*,*) "maxval(wCFL(:,ip)) = ", maxval(wCFL(:,ip))
                        write(*,*) "wti = ", wti(i,:,ip)
                        if (show_nan_wCFL.eq.2.and.maxval(wCFL(:,ip)).gt.1.0_rk) then
                            write(*,*) "(z_L, wti, wCFL) where wCFL>1 = "
                            do k=1,k_max-1
                                if (wCFL(k,ip).gt.1.0_rk) write(*,*) z(k)+hz(k)/2, wti(i,k+1,ip), wCFL(k,ip)
                            end do
                        end if
                    end if
                    stop
                end if
            enddo
        end do
               vv=1._rk
	
    ! OUTPUT
	if (mod(int(id*dt),output_step).eq.0) then
           
    fick_per_day = 86400.0_rk * fick
    sink_per_day = 86400.0_rk * sink
	
    if (sediments_units_convert.eq.0) then 
!		write(*,*) "output_step: ", output_step
!		write(*,*) "input_step: ", input_step
!		write(*,*) "id: ", id
        if (ncoutfile_type == 1) then
            call save_netcdf(i_max, k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, i_day+1, id, idt, output_step, input_step)
        else
            call save_netcdf(i_max, k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, julianday, id, idt, output_step, input_step)
        endif
!    else
!        !convert into observarional units in the sediments for dissolved (mass/pore water) and solids (mass/mass) with an exception for biota
!        cc_out(:,:,:)=cc(:,:,:)
!        do ip=1,par_max
!            if (ip.ne.id_Phy.or.ip.ne.id_Het.or.ip.ne.id_Baae.or.ip.ne.id_Baan.or.ip.ne.id_Bhae.or.ip.ne.id_Bhan) then 
!                cc_out(:,k_bbl_sed+1:k_max,:)=cc(:,k_bbl_sed+1:k_max,:)*pF1(:,k_bbl_sed:k_max-1,:)
!            endif
!        enddo 
!        if (ncoutfile_type == 1) then
!            call save_netcdf(i_max, k_max, max(1,julianday), cc_out, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
!            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, i_day+1, id, idt,i_sec_pr)
!        else
!            call save_netcdf(i_max, k_max, max(1,julianday), cc_out, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
!            fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x, julianday, id, idt,i_sec_pr)
!        endif        
    endif

    endif ! END OF OUTPUT
    
    end do ! end of Subloop over timesteps in the course of one day (idt)


    if (julianday == 364) then !Note: saving to ascii
        vv=1._rk
        call saving_state_variables(trim(outfile_name), model_year, julianday, &
                       i_max, k_max, par_max, par_name, z, hz, cc, vv, t, s, kz)
    endif
    
    !Write daily output to screen if required
    do i=i_min,i_max
        if (show_maxmin.eq.1) then
            write(*,*) "maxval(cc(i,:,:),1) = ", maxval(cc(i,:,:),1) !Good to keep an eye on max/min values for each parameter
            write(*,*) "minval(cc(i,:,:),1) = ", minval(cc(i,:,:),1)
        end if
        if (show_kztCFL.gt.0) then
            do ip=1,par_max
                kztCFL(:,ip) = (kzti(i,2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
                
            end do
            if (show_kztCFL.eq.2) then
                write(*,*) "maxval(kzti,2) = ", maxval(kzti,2)
                write(*,*) "maxval(kztCFL,2) = ", maxval(kztCFL,2)
            else
                write(*,*) "maxval(kzti) = ", maxval(kzti)
                write(*,*) "maxval(kztCFL) = ", maxval(kztCFL)
            end if
        end if
        if (show_wCFL.gt.0) then
            do ip=1,par_max
                wCFL(:,ip) = (abs(wti(i,2:k_max,ip))*dt/freq_sed)/dz(1:k_max-1)
            end do
            if (show_wCFL.eq.2) then
                write(*,*) "maxval(wti(i,:,:),2) = ", maxval(wti(i,:,:),2)
                write(*,*) "maxval(wCFL,2) = ", maxval(wCFL,2)
            else
                write(*,*) "maxval(wti) = ", maxval(wti)
                write(*,*) "maxval(wCFL) = ", maxval(wCFL)
            end if
        end if
    enddo
    
     end do   ! end of BIG cycle    

    end subroutine do_brom_transport
!=======================================================================================================================







!=======================================================================================================================
    subroutine clear_brom_transport()

    deallocate(z_w)
    deallocate(dz_w)
    deallocate(hz_w)
    deallocate(t_w)
    deallocate(u_x_w)
    deallocate(s_w)
    deallocate(kz_w)
    deallocate(cc_hmix_w)
    deallocate(z)
    deallocate(dz)
    deallocate(hz)
    deallocate(t)
    deallocate(s)
    deallocate(u_x)
    deallocate(kz)
    deallocate(kl)
    deallocate(cc_hmix)
    deallocate(kz_bio)
    deallocate(kz_mol)
    deallocate(alpha)
    deallocate(phi)
    deallocate(phi1)
    deallocate(phi_inv)
    deallocate(tortuosity)
    deallocate(w_b)
    deallocate(u_b)
    deallocate(wti)
    deallocate(pF1)
    deallocate(pF2)
    deallocate(pWC)
    deallocate(cc)
    deallocate(cc_out)
    deallocate(dcc)
    deallocate(vv)
    deallocate(dVV)
    deallocate(fick)
    deallocate(fick_per_day)
    deallocate(sink)
    deallocate(sink_per_day)
    deallocate(wbio)
    deallocate(surf_flux)
    deallocate(bott_flux)
    deallocate(bc_top)
    deallocate(bc_bottom)
    deallocate(bctype_top)
    deallocate(bctype_bottom)
    deallocate(par_name)
    deallocate(Izt)
    deallocate(pressure)
    deallocate(Eair)
    deallocate(hice)
    deallocate(cc_top)
    deallocate(cc_bottom)
    deallocate(kzti)
    deallocate(kztCFL)
    deallocate(wCFL)
    deallocate(k_wat)
    deallocate(k_sed)
    deallocate(k_sed1)
    deallocate(x)
    deallocate(dx)
    deallocate(hx)
    if (diff_method.gt.0) call clean_tridiagonal()

    end subroutine clear_brom_transport
!=======================================================================================================================


    end module brom_transport
