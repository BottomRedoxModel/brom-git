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


    module io_netcdf

    use netcdf
    use fabm, only: type_model, fabm_get_bulk_diagnostic_data
    use fabm_types, only: attribute_length, rk


    implicit none
    !all is private
    private
    !public functions
    public input_netcdf, init_netcdf, save_netcdf, close_netcdf !input_netcdf, 
    !netCDF file id
    integer               :: nc_id
    integer, allocatable  :: parameter_id(:)
    integer, allocatable  :: parameter_fick_id(:)
    integer, allocatable  :: parameter_sink_id(:)
    integer, allocatable  :: parameter_id_diag(:)
    !parameter_ids
    integer               :: i_id, z_id, z2_id, time_id, Eair_id, hice_id
    integer               :: pH_id, T_id, S_id, Kz_id, Kz_sol_id, Kz_par_id, w_sol_id, w_par_id, u_x_w_id
    integer               :: pCO2_id, Om_Ca_id, Om_Ar_id

    logical               :: first



    contains



!=======================================================================================================================
    subroutine input_netcdf(z_w, dz_w, hz_w, t_w, s_w, kz_w, hmix_rate_w, Eair, use_Eair, hice, use_hice, year, i_water, i_max, &
        days_in_yr, k_wat_bbl, par_name, par_max, bctype_top, bctype_bottom, cc_top, cc_bottom, hmixtype, cc_hmix_w, h_adv, u_x_w)

    !Inputs data from netCDF files

    use io_ascii, only: get_brom_name, get_brom_par
    use ids         !Provides access to variable indices id_O2 etc.

    !Input variables
    integer, intent(in)                         :: use_Eair, use_hice, year, i_water, i_max, days_in_yr, par_max, h_adv
    integer, dimension(:,:), intent(in)         :: bctype_top, bctype_bottom, hmixtype
    character(len=attribute_length), dimension(:), intent(in) :: par_name

    !Output variables
    real(rk), allocatable, dimension(:), intent(out)        :: z_w, dz_w         !Layer midpoint depths and spacing between them
    real(rk), allocatable, dimension(:), intent(out)        :: hz_w              !Layer thicknesses
    real(rk), allocatable, dimension(:,:,:), intent(out)    :: t_w, s_w, kz_w    !Temperature, salinity, and vertical diffusivity
    real(rk), allocatable, dimension(:,:,:), intent(out)    :: hmix_rate_w,u_x_w !Horizontal relaxation rates [day^-1] and horizontal advection
    real(rk), pointer, dimension(:), intent(out)            :: hice              !Ice thickness [m]
    real(rk), allocatable, dimension(:), intent(out)        :: Eair              ! 24-hr average surface downwelling shortwave irradiance in air [W/m2]
    real(rk), dimension(:,:,:), intent(out)                 :: cc_top, cc_bottom !Variable Dirichlet boundary conditions for surface and bottom
    real(rk), allocatable, dimension(:,:,:,:), intent(out)  :: cc_hmix_w         !Horizontal mixing concentrations [uM]

    !Input/output variables
    integer, intent(out)                        :: k_wat_bbl

    !Local variables
    integer                                     :: ncid, i, j, ip, year0, yeari, nyrdays, nyrdaysm1, ni,inext,ilast, iday
    integer                                     :: ndims, idim_z, idim_z2, idim_time, idim_lat, idim_lon, ll_rec(2), ll_sel(2)
    integer                                     :: lat_rec, lon_rec, time_rec, h_rec, h_rec2, istart, iend
    integer                                     :: t_varid, s_varid, Eair_varid, hice_varid, kz_varid, z_varid, z2_varid, time_varid, &
                                                   lat_varid, lon_varid, hmix_rate_varid, u_varid, v_varid
    integer, dimension(nf90_max_var_dims)       :: dimids1, dimids2
    integer, dimension(par_max)                 :: surf_varid, bot_varid, hmix_varid
    integer, allocatable, dimension(:)          :: inds, inds2
    real(rk), allocatable, dimension(:)         :: z_w2, z_w_error, time_temp, Eair_temp, hice_temp, z_temp, z_temp2
    real(rk), allocatable, dimension(:,:)       :: t_temp, s_temp, kz_temp, u_temp, v_temp
    real(rk), allocatable, dimension(:,:,:)     :: Eair_temp2, hice_temp2
    real(rk), allocatable, dimension(:,:,:,:)   :: t_temp2, s_temp2, kz_temp2
    real(rk)                                    :: time0, timei, time00, dind, dind2
    real(rk), dimension(days_in_yr)             :: cc_temp1
    real(rk), allocatable, dimension(:,:)       :: cc_temp2, hmix_rate_temp
    character(len=64) :: ncinfile_name, ncint_name, ncins_name, ncinkz_name, ncinEair_name, ncinhice_name, &
                        ncinz_name, ncinz2_name, ncintime_name, ncinlat_name, ncinlon_name, ncinhmix_rate_name
    character(len=attribute_length), dimension(par_max) :: ncinsurfpar_name, ncinbotpar_name, ncinhmixpar_name
    real(rk)  :: ncinkz_fac, ncinEair_fac, ncinhice_fac, ncinsurfpar_fac(par_max), ncinbotpar_fac(par_max), &
                 ncinhmixpar_fac(par_max), ncinhmix_rate_fac
    integer   :: nc_set_k_wat_bbl, nc_staggered_grid, nc_bottom_to_top, nc_z_increasing_upward, nc_year0, nc_file_source


    !Input and output file names
    ncinfile_name = get_brom_name("ncinfile_name")
    nc_file_source = get_brom_par("nc_file_source", 1.0_rk)
    
    !NetCDF input names of dimensions
    ncinz_name = get_brom_name("ncinz_name")
    ncinz2_name = get_brom_name("ncinz2_name")
    ncintime_name = get_brom_name("ncintime_name")

    !NetCDF input names and scale factors for physical forcing variables
    ncint_name = get_brom_name("ncint_name")
    ncins_name = get_brom_name("ncins_name")
    ncinkz_name = get_brom_name("ncinkz_name")
    ncinkz_fac = get_brom_par("ncinkz_fac",1.0_rk)
    if (use_Eair.eq.1) then
        ncinEair_name = get_brom_name("ncinEair_name")
        ncinEair_fac = get_brom_par("ncinEair_fac",1.0_rk)
    end if
    if (use_hice.eq.1) then
        ncinhice_name = get_brom_name("ncinhice_name")
        ncinhice_fac = get_brom_par("ncinhice_fac",1.0_rk)
    end if
    if (h_adv.ne.1) then
    !NetCDF input names and scale factors for surface forcing variables (this list may need to be appended to include more parameters)
    ncinsurfpar_fac = 1.0_rk
    ncinsurfpar_name(id_NO3) = get_brom_name("ncinNO3s_name")
    ncinsurfpar_fac(id_NO3) = get_brom_par("ncinNO3s_fac",1.0_rk)
    ncinsurfpar_name(id_NH4) = get_brom_name("ncinNH4s_name")
    ncinsurfpar_fac(id_NH4) = get_brom_par("ncinNH4s_fac",1.0_rk)
    ncinsurfpar_name(id_PO4) = get_brom_name("ncinPO4s_name")
    ncinsurfpar_fac(id_PO4) = get_brom_par("ncinPO4s_fac",1.0_rk)
    ncinsurfpar_name(id_Si) = get_brom_name("ncinSis_name")
    ncinsurfpar_fac(id_Si) = get_brom_par("ncinSis_fac",1.0_rk)
    ncinsurfpar_name(id_Alk) = get_brom_name("ncinAlks_name")
    ncinsurfpar_fac(id_Alk) = get_brom_par("ncinAlks_fac",1.0_rk)

    !NetCDF input names and scale factors for bottom forcing variables (this list may need to be appended to include more parameters)
    ncinbotpar_fac = 1.0_rk

    !NetCDF input names and scale factors for horizontal mixing variables (this list may need to be appended to include more parameters)
    ncinhmixpar_fac = 1.0_rk
    ncinhmixpar_name(id_NO3) = get_brom_name("ncinNO3hmix_name")
    ncinhmixpar_fac(id_NO3) = get_brom_par("ncinNO3hmix_fac",1.0_rk)
    ncinhmixpar_name(id_NH4) = get_brom_name("ncinNH4hmix_name")
    ncinhmixpar_fac(id_NH4) = get_brom_par("ncinNH4hmix_fac",1.0_rk)
    ncinhmixpar_name(id_PO4) = get_brom_name("ncinPO4hmix_name")
    ncinhmixpar_fac(id_PO4) = get_brom_par("ncinPO4hmix_fac",1.0_rk)
    ncinhmixpar_name(id_Si) = get_brom_name("ncinSihmix_name")
    ncinhmixpar_fac(id_Si) = get_brom_par("ncinSihmix_fac",1.0_rk)
    ncinhmixpar_name(id_O2) = get_brom_name("ncinO2hmix_name")
    ncinhmixpar_fac(id_O2) = get_brom_par("ncinO2hmix_fac",1.0_rk)

    !NetCDF input for horizontal mixing rate
    ncinhmix_rate_name = get_brom_name("ncinhmix_rate_name")
    ncinhmix_rate_fac = get_brom_par("ncinhmix_rate_fac",1.0_rk)
    endif
    !Other NetCDF parameters
    nc_set_k_wat_bbl = get_brom_par("nc_set_k_wat_bbl")
    nc_staggered_grid = get_brom_par("nc_staggered_grid")
    nc_bottom_to_top = get_brom_par("nc_bottom_to_top")
    nc_z_increasing_upward = get_brom_par("nc_z_increasing_upward")
    nc_year0 = get_brom_par("nc_year0")
    if (nc_staggered_grid.eq.0) then
        write(*,*) "Assuming unstaggered grid, tracers and diffusivity on layer centres"
    else
        write(*,*) "Assuming staggered grid, tracers on layer centres and diffusivity of layer interfaces"
    end if

    !Get netcdf variable ids
    write(*,*) "Opening ", trim(ncinfile_name)
    call check_err(nf90_open(trim(ncinfile_name), NF90_NOWRITE, ncid))
    call check_err(nf90_inq_varid(ncid, trim(ncint_name), t_varid))
    call check_err(nf90_inq_varid(ncid, trim(ncins_name), s_varid))
    call check_err(nf90_inq_varid(ncid, trim(ncinkz_name), kz_varid))
    call check_err(nf90_inq_varid(ncid, trim(ncinz_name), z_varid))
    call check_err(nf90_inq_varid(ncid, trim(ncinz2_name), z2_varid))
    call check_err(nf90_inq_varid(ncid, trim(ncintime_name), time_varid))
    if (h_adv.ne.0)then
    call check_err(nf90_inq_varid(ncid, 'u', u_varid))
    call check_err(nf90_inq_varid(ncid, 'v', v_varid))
    endif
    !use temperature variable to get ids of midpoint depth, time, and possibly lat/lon dimensions
    call check_err(nf90_inquire_variable(ncid, t_varid, dimids = dimids1))
    !if (nc_file_source.ne.3)then
    !    ndims = maxloc(dimids1(1:100),1,mask=(dimids1(1:100).gt.0))
    !else
    !    ndims = 2
    !endif
    select case (nc_file_source) 
      case (1) !ROMS   
        ndims = maxloc(dimids1(1:100),1,mask=(dimids1(1:100).gt.0))       
      case (2) !GETM
        ndims = 4       
      case (3) !FVCOM
        ndims = 2
    end select            
    
    if (ndims.eq.2) write(*,*) "Assuming (temperature,salinity,diffusivity) variables have netCDF dimensions (depth,time)"
    if (ndims.eq.4) write(*,*) "Assuming (temperature,salinity,diffusivity) variables have netCDF dimensions (lon,lat,depth,time) or (lat,lon,depth,time)"
    if (ndims.eq.4) then
        ncinlat_name = get_brom_name("ncinlat_name")
        ncinlon_name = get_brom_name("ncinlon_name")
        call check_err(nf90_inq_varid(ncid, trim(ncinlat_name), lat_varid))
        call check_err(nf90_inq_varid(ncid, trim(ncinlon_name), lon_varid))
    end if
    idim_z = minloc(dimids1,1,mask=dimids1.eq.z_varid)
    select case (nc_file_source) 
        case (1) !ROMS
          idim_time = 2 !minloc(dimids1,1,mask=dimids1.eq.time_varid)  
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_z), len = h_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_time), len = time_rec))
        case (2) !GETM
          idim_time = minloc(dimids1,1,mask=dimids1.eq.time_varid)
          call check_err(nf90_inquire_dimension(ncid, dimids1(3), len = h_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(4), len = time_rec))
        case (3) !FVCOM
          idim_time = 2 !minloc(dimids1,1,mask=dimids1.eq.time_varid)  
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_z), len = h_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_time), len = time_rec))
    end select
    !use diffusivity variable to get length of second input depth variable (possibly same as first)
    call check_err(nf90_inquire_variable(ncid, kz_varid, dimids = dimids2))
    idim_z2 = minloc(dimids2,1,mask=dimids2.eq.z2_varid)
    call check_err(nf90_inquire_dimension(ncid, dimids2(idim_z2), len = h_rec2))
    ll_rec = 1
    ll_sel = 1
    if (ndims.eq.4) then
        idim_lat = minloc(dimids1,1,mask=dimids1.eq.lat_varid)
        idim_lon = minloc(dimids1,1,mask=dimids1.eq.lon_varid)
    select case (nc_file_source) 
        case (1) !ROMS
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_lat), len = lat_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_lon), len = lon_rec))
        case (2) !GETM
          call check_err(nf90_inquire_dimension(ncid, dimids1(1), len = lat_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(2), len = lon_rec))
        case (3) !FVCOM
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_lat), len = lat_rec))
          call check_err(nf90_inquire_dimension(ncid, dimids1(idim_lon), len = lon_rec))
    end select
        ll_rec(idim_lat) = lat_rec
        ll_rec(idim_lon) = lon_rec
        if (lat_rec.gt.1) then
            ll_sel(idim_lat) = get_brom_par("nc_latsel")
            write(*,*) "Selected latitude index = ", ll_sel(idim_lat)
        end if

        if (lon_rec.gt.1) then
            ll_sel(idim_lon) = get_brom_par("nc_lonsel")
            write(*,*) "Selected longitude index = ", ll_sel(idim_lon)
        end if
    end if
    if (use_Eair.eq.1) call check_err(nf90_inq_varid(ncid, trim(ncinEair_name), Eair_varid))
    if (use_hice.eq.1) call check_err(nf90_inq_varid(ncid, trim(ncinhice_name), hice_varid))
    do ip=1,par_max
        if (bctype_top(i_water,ip).eq.3) call check_err(nf90_inq_varid(ncid, trim(ncinsurfpar_name(ip)), surf_varid(ip)))
        if (bctype_bottom(i_water,ip).eq.3) call check_err(nf90_inq_varid(ncid, trim(ncinbotpar_name(ip)), bot_varid(ip)))
        
    select case (nc_file_source) 
        case (1,2) !ROMS, GETM
        if (hmixtype(i_water,ip).eq.1) call check_err(nf90_inq_varid(ncid, trim(ncinhmixpar_name(ip)), hmix_varid(ip)))
        case (3)
        continue
    end select

    end do
    select case (nc_file_source) 
        case (1,2) !ROMS, GETM
        if (maxval(hmixtype(i_water,:)).eq.1) call check_err(nf90_inq_varid(ncid, trim(ncinhmix_rate_name), hmix_rate_varid))
        case (3) !FVCOM
        continue
    end select


    !Set the no. grid points in the water column k_wat_bbl using the netCDF input
    if (nc_set_k_wat_bbl.eq.1) then
        k_wat_bbl = h_rec
        write(*,*) "k_wat_bbl set from netCDF input to ", k_wat_bbl
    else
        write(*,*) "k_wat_bbl set from brom.yaml to ", k_wat_bbl
    end if


    !Allocate temporary variables
    if (ndims.eq.2) then
        allocate(t_temp(h_rec,time_rec))
        allocate(s_temp(h_rec,time_rec))
        allocate(kz_temp(h_rec2,time_rec))
        allocate(u_temp(h_rec2,time_rec))
        allocate(v_temp(h_rec2,time_rec))
        if (use_Eair.eq.1) allocate(Eair_temp(time_rec))
        if (use_hice.eq.1) allocate(hice_temp(time_rec))
    end if
    if (nc_file_source.ne.3) then
      if (ndims.eq.4) then
          allocate(t_temp2(ll_rec(1),ll_rec(2),h_rec,time_rec))
          allocate(s_temp2(ll_rec(1),ll_rec(2),h_rec,time_rec))
          allocate(kz_temp2(ll_rec(1),ll_rec(2),h_rec2,time_rec))
          allocate(kz_temp(h_rec2,time_rec))
          if (use_Eair.eq.1) allocate(Eair_temp2(ll_rec(1),ll_rec(2),time_rec))
          if (use_hice.eq.1) allocate(hice_temp2(ll_rec(1),ll_rec(2),time_rec))
      end if
    endif
    allocate(z_temp(h_rec))
    allocate(time_temp(time_rec))
    if (nc_staggered_grid.eq.1) allocate(z_temp2(h_rec2))
    allocate(cc_temp2(h_rec,days_in_yr))
    allocate(hmix_rate_temp(h_rec,days_in_yr))


    !Allocate permanent variables
    allocate(t_w(i_max,k_wat_bbl,days_in_yr))
    allocate(s_w(i_max,k_wat_bbl,days_in_yr))
    allocate(kz_w(i_max,k_wat_bbl+1,days_in_yr))
    allocate(z_w(k_wat_bbl))
    allocate(z_w2(k_wat_bbl+1))
    allocate(z_w_error(k_wat_bbl))
    allocate(dz_w(k_wat_bbl))
    allocate(hz_w(k_wat_bbl))
    allocate(Eair(days_in_yr))
    allocate(hice(days_in_yr))
    allocate(inds(k_wat_bbl))
    allocate(cc_hmix_w(i_max,par_max,k_wat_bbl,days_in_yr))
    allocate(hmix_rate_w(i_max,k_wat_bbl,days_in_yr))
    allocate(u_x_w(i_max,k_wat_bbl,days_in_yr))

    !Get physical forcings and apply scale factors
    if (ndims.eq.2) then
        call check_err(nf90_get_var(ncid, t_varid, t_temp))
        call check_err(nf90_get_var(ncid, s_varid, s_temp))
        call check_err(nf90_get_var(ncid, kz_varid, kz_temp))
        if (nc_file_source.eq.3) then
            call check_err(nf90_get_var(ncid, u_varid, u_temp))
            call check_err(nf90_get_var(ncid, v_varid, v_temp))
        endif
        kz_temp = ncinkz_fac * kz_temp
        if (use_Eair.eq.1) then
            call check_err(nf90_get_var(ncid, Eair_varid, Eair_temp))
            Eair_temp = ncinEair_fac * Eair_temp
        end if
        if (use_hice.eq.1) then
            call check_err(nf90_get_var(ncid, hice_varid, hice_temp))
            hice_temp = ncinhice_fac * hice_temp
        else
            hice=0.0_rk     
        end if
    end if
    if (nc_file_source.ne.3) then
        if (ndims.eq.4) then
            call check_err(nf90_get_var(ncid, t_varid, t_temp2))
            call check_err(nf90_get_var(ncid, s_varid, s_temp2))
            call check_err(nf90_get_var(ncid, kz_varid, kz_temp))
            kz_temp = ncinkz_fac * kz_temp
            if (use_Eair.eq.1) then
                call check_err(nf90_get_var(ncid, Eair_varid, Eair_temp2))
                Eair_temp2 = ncinEair_fac * Eair_temp2
            end if
            if (use_hice.eq.1) then
                call check_err(nf90_get_var(ncid, hice_varid, hice_temp2))
                hice_temp2 = ncinhice_fac * hice_temp2
            end if
        end if
    endif
    call check_err(nf90_get_var(ncid, z_varid, z_temp))
    call check_err(nf90_get_var(ncid, time_varid, time_temp))
    if (nc_staggered_grid.eq.1) call check_err(nf90_get_var(ncid, z2_varid, z_temp2))

    !Get surface/bottom boundary conditions (dimension [days_in_yr]) and apply scale factors
    cc_top = 0.0_rk
    cc_bottom = 0.0_rk
    do ip=1,par_max
        !Here the surface and bottom boundary conditions are read in where necessary
        !Note that surface and bottom BCs MUST have dimensions (days_in_yr) - they are read in directly with no subsetting
        if (bctype_top(i_water,ip).eq.3) then
            call check_err(nf90_get_var(ncid, surf_varid(ip), cc_temp1))
            cc_top(i_water,ip,1:days_in_yr) = ncinsurfpar_fac(ip) * cc_temp1(1:days_in_yr)
        end if
        if (bctype_bottom(i_water,ip).eq.3) then
            call check_err(nf90_get_var(ncid, bot_varid(ip), cc_temp1))
            cc_bottom(i_water,ip,1:days_in_yr) = ncinbotpar_fac(ip) * cc_temp1(1:days_in_yr)
        end if
    end do

    if (nc_file_source.eq.3) then
        time_temp=time_temp*86400.  ! FVCOM correction ! time is given in days not seconds
    endif


    !!Establish the initial index istart using the netcdf time variable and the chosen year
    !First calculate the earliest year "year0" and corresponding time "time00" [s] at start of this year
    time0 = time_temp(1) !Assume that netcdf time is time in seconds since nc_year0-01-01 00:00:00
    year0 = -1
    yeari = nc_year0
    timei = 0
    do while (year0.eq.-1)
        nyrdays = 365 + merge(1,0,(mod(yeari,4).eq.0))
        if ((time0-timei).ge.0.0_rk.and.(time0-timei).lt.(nyrdays*86400)) then !Found initial year in file
            year0 = yeari
            time00 = timei
        elseif (time0.gt.0.0_rk) then !Advance candidate year0 = yeari forwards by one year
            yeari = yeari + 1
            timei = timei + nyrdays*86400
        elseif (time0.lt.0.0_rk) then !Advance candidate year0 = yeari backwards by one year
            nyrdaysm1 = 365 + merge(1,0,(mod(yeari-1,4).eq.0))
            yeari = yeari - 1
            timei = timei - nyrdaysm1*86400
        end if
    end do
    !Find the starting end finishing time indices (istart, iend) for the selected year "year"
    istart = -1
    iend = -1
    yeari = year0
    timei = time00
    do i=1,time_rec
        nyrdays = 365 + merge(1,0,(mod(yeari,4).eq.0))
        if ((time_temp(i)-timei).ge.(nyrdays*86400)) then !Advance the year and time counters
            yeari = yeari + 1
            timei = timei + nyrdays*86400
            nyrdays = 365 + merge(1,0,(mod(yeari,4).eq.0))
        end if
        if (yeari.eq.year.and.istart.eq.-1) istart = i    !istart = first index where (yeari==year)
        !if (i-istart.eq.days_in_yr) then
        !    iend = i                                      !iend
        !    exit
        !end if
        !if (iend.eq.-1.and.i.eq.time_rec) iend = i        !...or the last index in the file
    end do
    iend=istart+365
    ni = iend-istart+1
    !select case (nc_file_source) 
!    case (1,2) !ROMS, GETM
      if (ni.lt.days_in_yr) then
        write(*,*) "Could not find days_in_yr time inputs starting from 1st day of selected year (ni = ", ni, ", stopping"
        stop
      end if
    !case (3) !FVCOM
    !  if (ni.lt.359) then !(ni.lt.days_in_yr) then
    !    write(*,*) "Could not find days_in_yr time inputs starting from 1st day of selected year (ni = ", ni, ", stopping"
    !    stop
    !  end if
    !end select


    !Check - these results can be validated by importing ocean_time into Matlab and converting to human dates using http://www.epochconverter.com/
    write(*,*) "First year in netCDF input: year0 = ", year0
    write(*,*) "Reference year for netCDF time in seconds: nc_year0 = ", nc_year0
    write(*,*) "NetCDF time in seconds at beginning of first recorded year: time00 = ", time00
    write(*,*) "Target year for BROM model input: year = ", year
    write(*,*) "Calculated starting index: istart = ", istart
    write(*,*) "Calculated finishing index: iend = ", iend
    write(*,*) "Number of input time points: ni = ", ni
    write(*,*) "Initial input time in netCDF units for chosen year: time_temp(istart) = ", time_temp(istart)
    write(*,*) "Final input time in netCDF units for chosen year: time_temp(iend) = ", time_temp(iend)


    !Set the index vectors for subsampling the netCDF input (coarsening the vertical resolution)

    if (nc_set_k_wat_bbl.eq.0) then
        dind = real(h_rec-1,kind=rk)/real(k_wat_bbl,kind=rk)
        do i=1,k_wat_bbl
            inds(i) = nint(1+0.5_rk*dind+(i-1)*dind)
        end do
        inds = min(h_rec, inds)
        if (nc_bottom_to_top.eq.1) inds = inds(k_wat_bbl:1:-1) !Reverse order

        allocate(inds2(k_wat_bbl+1))
        do i=1,k_wat_bbl+1
            inds2(i) = nint(1+(i-1)*dind)
        end do
        inds2 = min(h_rec2, inds2)
        if (nc_bottom_to_top.eq.1) inds2 = inds2(k_wat_bbl+1:1:-1) !Reverse order
    else
        inds(1:h_rec) = (/ (i, i=1, h_rec) /)
        if (nc_bottom_to_top.eq.1) inds = inds(h_rec:1:-1) !Reverse order
        allocate(inds2(h_rec2))
        inds2(1:h_rec2) = (/ (i, i=1, h_rec2) /)
        if (nc_bottom_to_top.eq.1) inds2 = inds2(h_rec2:1:-1) !Reverse order
    end if




    !Set the water column grid geometry
    if (nc_z_increasing_upward.eq.1) then !Flip sign of input if necessary: BROM assumes z increasing DOWNWARD
        z_temp = -1.0_rk*z_temp
        if (nc_staggered_grid.eq.1) z_temp2 = -1.0_rk*z_temp2
    end if

    !Heights of layer midpoints (z_w)
    select case (nc_file_source) 
    case (1,2) !ROMS
      z_w(1:k_wat_bbl) = z_temp(inds)
    case (3) !FVCOM
      z_w(:) = z_temp(:)
    end select

    !Spacing between layer midpoints (dz_w)
    dz_w(1:k_wat_bbl-1) = z_w(2:k_wat_bbl) - z_w(1:k_wat_bbl-1)
    if (nc_file_source.ne.3)    dz_w(k_wat_bbl) = 0.0_rk

    !Layer thicknesses (hz_w)
    if (nc_file_source.ne.3)  then
      if (nc_staggered_grid.eq.0) then
        hz_w(1) = dz_w(1)
        write(*,*) "Assuming the input depths correspond to layer midpoints (min(abs(depth))>0.)"
        write(*,*) "Therefore assuming that surface layer thickness hz_w(1) = spacing between first two grid points dz_w(1)"
        write(*,*) "Given this first layer thickness, subsequent layer thicknesses are inferred"
        do j=2,k_wat_bbl
            hz_w(j) = 2.0_rk*dz_w(j-1) - hz_w(j-1)
        end do
        z_w2(1:k_wat_bbl) = z_w(1:k_wat_bbl) - 0.5_rk*hz_w(1:k_wat_bbl)
        z_w2(k_wat_bbl+1) = z_w(k_wat_bbl) + 0.5_rk*hz_w(k_wat_bbl)
      end if
    else
        hz_w(1) = dz_w(1)
        write(*,*) "Assuming the input depths correspond to layer midpoints (min(abs(depth))>0.)"
        write(*,*) "Therefore assuming that surface layer thickness hz_w(1) = spacing between first two grid points dz_w(1)"
        write(*,*) "Given this first layer thickness, subsequent layer thicknesses are inferred"
        do j=2,k_wat_bbl
            hz_w(j) = 2.0_rk*dz_w(j-1) - hz_w(j-1)
        end do
        z_w2(1:k_wat_bbl) = z_w(1:k_wat_bbl) - 0.5_rk*hz_w(1:k_wat_bbl)
        z_w2(k_wat_bbl+1) = z_w(k_wat_bbl) + 0.5_rk*hz_w(k_wat_bbl)
    endif


    if (nc_file_source.ne.3)  then


    if (nc_staggered_grid.eq.1) then
        write(*,*) "For staggered grid, calculating layer thicknesses using the second netCDF depth variable"
        if (h_rec2.eq.h_rec+1) then       !All interfaces supplied by second depth variable
            z_w2(1:k_wat_bbl+1) = z_temp2(inds2(1:h_rec2))
        else if (h_rec2.eq.h_rec.and.minval(z_temp2).lt.minval(z_temp)) then !Only upper interfaces supplied: infer lowermost interface depth
            z_w2(1:k_wat_bbl) = z_temp2(inds2(1:k_wat_bbl))
            z_w2(k_wat_bbl+1) = maxval(z_temp) + (maxval(z_temp)-maxval(z_temp2))
        else if (h_rec2.eq.h_rec.and.minval(z_temp2).gt.minval(z_temp)) then !Only lower interfaces supplied: infer uppermost interface depth
            z_w2(1) = minval(z_temp) - (minval(z_temp2)-minval(z_temp))
            z_w2(2:k_wat_bbl+1) = z_temp2(inds2(1:k_wat_bbl))
        else if (h_rec2.eq.h_rec-1) then                                     !Only internal interfaces supplied: infer lower and uppermost depths
            z_w2(1) = minval(z_temp) - (minval(z_temp2)-minval(z_temp))
            z_w2(2:k_wat_bbl) = z_temp2(inds2(1:k_wat_bbl-1))
            z_w2(k_wat_bbl+1) = maxval(z_temp) + (maxval(z_temp)-maxval(z_temp2))
        end if
        !Check the depth variables are consistent (midpoints z_w should lie half-way between interfaces z_w2)
        z_w_error = (z_w-z_w2(1:k_wat_bbl))-(z_w2(2:k_wat_bbl+1)-z_w)
        if (maxval(abs(z_w_error)).gt.1.0E-8_rk) then
            write(*,*) "Warning: Input layer midpoint depths are not located exactly halfway between input layer interface depths"
            write(*,*) "maxval(abs(z_w_error)) = ", maxval(abs(z_w_error))
        end if
        !Infer layer thicknesses, midpoint depths and their increments from z_w2
        hz_w = z_w2(2:k_wat_bbl+1) - z_w2(1:k_wat_bbl)
        z_w = z_w2(1:k_wat_bbl) + 0.5_rk*hz_w
        dz_w(1:k_wat_bbl-1) = z_w(2:k_wat_bbl) - z_w(1:k_wat_bbl-1)
        if (maxval(hz_w).eq.0.0_rk) then
            write(*,*) "netCDF grid is apparently NOT staggered: Please check netCDF file or set nc_staggered_grid = 0"
            stop
        end if
    end if

    !Set the water temperature (t_w), salinity (s_w), vertical diffusivity (kz_w),
    !and (if required) the surface irradiance (Eair) and ice thickness (hice)

    if (nc_file_source.ne.3)  istart = istart + 2  !for GETM, ROMS
!    if (nc_file_source.eq.3)  istart = istart -140
    if (ndims.eq.2) then  !Assuming netcdf input dimensions (depth,time) for variables (t,s,kz)
        do i=1,days_in_yr !Loop over days_in_yr
            t_w(i_water,1:k_wat_bbl,i) = t_temp(inds,istart+i-1)
            s_w(i_water,1:k_wat_bbl,i) = s_temp(inds,istart+i-1)
            if (nc_staggered_grid.eq.0) then
                !If not staggered, linearly interpolate to layer interfaces
                kz_w(i_water,1,i)           = 0.0_rk
                kz_w(i_water,2:k_wat_bbl,i) = kz_temp(inds(2:k_wat_bbl),istart+i-1) + 0.5_rk*hz_w(1:k_wat_bbl-1)*&
                    (kz_temp(inds(2:k_wat_bbl),istart+i-1) - kz_temp(inds(1:k_wat_bbl-1),istart+i-1))/dz_w(1:k_wat_bbl-1)
                kz_w(i_water,k_wat_bbl+1,i) = 0.0_rk
            end if
            if (nc_staggered_grid.eq.1) then
                if ((h_rec2.eq.h_rec+1).or.(nc_set_k_wat_bbl.eq.0)) then !All interfaces supplied or subsampled
                    kz_w(i_water,1:k_wat_bbl+1,i) = kz_temp(inds2,istart+i-1)
                else if (h_rec2.eq.h_rec.and.z_w2(1).lt.z_w(1)) then !Only upper interfaces supplied (so duplicate lowermost input)
                    if (i.eq.1) write(*,*) "Only upper interfaces supplied: duplicating lowermost input"
                    kz_w(i_water,1:h_rec,i) = kz_temp(inds2(1:h_rec),istart+i-1)
                    kz_w(i_water,h_rec+1,i) = kz_temp(inds2(h_rec),istart+i-1)
                else if (h_rec2.eq.h_rec.and.z_w2(1).gt.z_w(1)) then !Only lower interfaces supplied (so duplicate uppermost input)
                    if (i.eq.1) write(*,*) "Only lower interfaces supplied: duplicating uppermost input"
                    kz_w(i_water,1,i) = kz_temp(inds2(1),istart+i-1)
                    kz_w(i_water,2:h_rec+1,i) = kz_temp(inds2(1:h_rec),istart+i-1)
                else if (h_rec2.eq.(h_rec-1)) then                   !Only internal interfaces supplied (so duplicate lower and uppermost inputs)
                    if (i.eq.1) write(*,*) "Only internal interfaces supplied: duplicating lowermost and uppermost inputs"
                    kz_w(i_water,1,i) = kz_temp(inds2(1),istart+i-1)
                    kz_w(i_water,2:h_rec,i) = kz_temp(inds2(1:h_rec-1),istart+i-1)
                    kz_w(i_water,h_rec+1,i) = kz_temp(inds2(h_rec-1),istart+i-1)
                end if
            end if
        end do
        Eair = 0.0_rk !Default in case not read from Eair_temp (use_Eair = 0)
        if (use_Eair.eq.1) Eair(1:days_in_yr) = Eair_temp(istart:istart+days_in_yr-1)
        hice = 0.0_rk !Default in case not read from hice_temp (use_hice = 0)
        if (use_hice.eq.1) hice(1:days_in_yr) = hice_temp(istart:istart+days_in_yr-1)
    end if

    if (ndims.eq.4) then  !Assuming netcdf input dimensions (lat/lon,lat/lon,depth,time) for variables (t,s,kz)
        do i=1,days_in_yr !Loop over days_in_yr
            t_w(i_water,1:k_wat_bbl,i) = t_temp2(ll_sel(1),ll_sel(2),inds,istart+i-1)
            s_w(i_water,1:k_wat_bbl,i) = s_temp2(ll_sel(1),ll_sel(2),inds,istart+i-1)
            if (nc_staggered_grid.eq.0) then
                !If not staggered, linearly interpolate to layer interfaces
                kz_w(i_water,1,i)           = 0.0_rk
                kz_w(i_water,2:k_wat_bbl,i) = kz_temp(inds(2:k_wat_bbl),istart+i-1) 
                !kz_w(i_water,2:k_wat_bbl,i) = abs(kz_temp(inds(2:k_wat_bbl),istart+i-1) + 0.5_rk*hz_w(1:k_wat_bbl-1)*&
                !    (kz_temp(inds(2:k_wat_bbl),istart+i-1) - kz_temp(inds(1:k_wat_bbl-1),istart+i-1))/dz_w(1:k_wat_bbl-1))
                kz_w(i_water,k_wat_bbl+1,i) = 0.0_rk
            end if
            if (nc_staggered_grid.eq.1) then
                if ((h_rec2.eq.h_rec+1).or.(nc_set_k_wat_bbl.eq.0)) then !All interfaces supplied or subsampled
                    kz_w(i_water,1:k_wat_bbl+1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2,istart+i-1)
                else if (h_rec2.eq.h_rec.and.z_w2(1).lt.z_w(1)) then !Only upper interfaces supplied (so duplicate lowermost input)
                    if (i.eq.1) write(*,*) "Only upper interfaces supplied: duplicating lowermost input"
                    kz_w(i_water,1:h_rec,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(1:h_rec),istart+i-1)
                    kz_w(i_water,h_rec+1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(h_rec),istart+i-1)
                else if (h_rec2.eq.h_rec.and.z_w2(1).gt.z_w(1)) then !Only lower interfaces supplied (so duplicate uppermost input)
                    if (i.eq.1) write(*,*) "Only lower interfaces supplied: duplicating uppermost input"
                    kz_w(i_water,1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(1),istart+i-1)
                    kz_w(i_water,2:h_rec+1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(1:h_rec),istart+i-1)
                else if (h_rec2.eq.(h_rec-1)) then                   !Only internal interfaces supplied (so duplicate lower and uppermost inputs)
                    if (i.eq.1) write(*,*) "Only internal interfaces supplied: duplicating lowermost and uppermost inputs"
                    kz_w(i_water,1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(1),istart+i-1)
                    kz_w(i_water,2:h_rec,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(1:h_rec-1),istart+i-1)
                    kz_w(i_water,h_rec+1,i) = kz_temp2(ll_sel(1),ll_sel(2),inds2(h_rec-1),istart+i-1)
                end if
            end if
        end do
        Eair = 0.0_rk !Default in case not read from Eair_temp (use_Eair = 0)
        if (use_Eair.eq.1) Eair(1:days_in_yr) = Eair_temp2(ll_sel(1),ll_sel(2),istart:istart+days_in_yr-1)
        hice = 0.0_rk !Default in case not read from hice_temp (use_hice = 0)
        if (use_hice.eq.1) hice(1:days_in_yr) = hice_temp2(ll_sel(1),ll_sel(2),istart:istart+days_in_yr-1)
    end if

    !If kz at top and bottom interfaces are set to zero, set equal to adjacent value (this has been observed in e.g. ROMS input)
    do i=1,days_in_yr
        if (kz_w(i_water,1,i).eq.0.0_rk) then
            kz_w(i_water,1,i) = kz_w(i_water,2,i)
        end if
        if (kz_w(i_water,k_wat_bbl+1,i).eq.0.0_rk) then
            kz_w(i_water,k_wat_bbl+1,i) = kz_w(i_water,k_wat_bbl,i)
        end if
    end do


    !Horizontal mixing inputs
    hmix_rate_w = 0.0_rk
    if (maxval(hmixtype(i_water,:)).eq.1) then
        do ip=1,par_max
            !Here the horizontal mixing variables are read in where necessary
            !Note: horizontal mixing variables MUST have dimensions (h_rec,days_in_yr), with depth indexing in agreement with (t,s) inputs
            if (hmixtype(i_water,ip).eq.1) then
                call check_err(nf90_get_var(ncid, hmix_varid(ip), cc_temp2))
                cc_hmix_w(i_water,ip,1:k_wat_bbl,1:days_in_yr) = ncinhmixpar_fac(ip) * cc_temp2(inds,1:days_in_yr)
            end if
        end do
        !Note: horizontal mixing rate MUST have dimensions (h_rec,days_in_yr), with depth indexing in agreement with (t,s) inputs
        call check_err(nf90_get_var(ncid, hmix_rate_varid, hmix_rate_temp))
        hmix_rate_w(i_water,1:k_wat_bbl,1:days_in_yr) = ncinhmix_rate_fac * hmix_rate_temp(inds,1:days_in_yr)
    end if

   endif


    if (nc_file_source.eq.3)  then  ! FVCOM for Lindesnes
!  correct days, since they start from August 20 2015, i.e.+122)
    do iday=1,243 !Loop over days availabel from FVCOM
            t_w(i_water,1:k_wat_bbl,iday) = t_temp(inds,istart+iday+1+122)
            s_w(i_water,1:k_wat_bbl,iday) = s_temp(inds,istart+iday+1+122)
            kz_w(i_water,1:k_wat_bbl,iday) = kz_temp(inds,istart+iday+1+122)
                        kz_w(i_water,k_wat_bbl+1,iday) = 0.
            u_x_w(i_water,1:k_wat_bbl,iday) = v_temp(inds,istart+iday+1+122) !(u_temp(inds,istart+iday+1)*u_temp(inds,istart+iday+1)+v_temp(inds,istart+iday+1)*v_temp(inds,istart+iday+1))**0.5
    enddo
    do iday=244,365 !Loop over days availabel from FVCOM
            t_w(i_water,1:k_wat_bbl,iday) = t_temp(inds,istart+iday+1-243)
            s_w(i_water,1:k_wat_bbl,iday) = s_temp(inds,istart+iday+1-243)
            kz_w(i_water,1:k_wat_bbl,iday) = kz_temp(inds,istart+iday+1-243)
                        kz_w(i_water,k_wat_bbl+1,iday) = 0.
            u_x_w(i_water,1:k_wat_bbl,iday) = v_temp(inds,istart+iday+1-243) !(u_temp(inds,istart+iday+1)*u_temp(inds,istart+iday+1)+v_temp(inds,istart+iday+1)*v_temp(inds,istart+iday+1))**0.5
    enddo

    !  correct depths array direction
    do j=1,k_wat_bbl
        t_w  (1,j,:) =  t_w(i_water,k_wat_bbl+1-j,:)
        s_w  (1,j,:) =  s_w(i_water,k_wat_bbl+1-j,:)
        kz_w (1,j,:) = kz_w(i_water,k_wat_bbl+1-j,:)
        u_x_w(1,j,:) =u_x_w(i_water,k_wat_bbl+1-j,:)
    enddo


    ! fill all the horizontal columns 
    do i = 1, i_water
        t_w(i,:,:)  =  t_w(1,:,:)
        s_w(i,:,:)  =  s_w(1,:,:)
        kz_w(i,:,:) = max(0.000001,min(0.15,kz_w(1,:,:)))
        u_x_w(i,:,:)= u_x_w(1,:,:)  !convert to m/s from cm/s
    enddo
   else
    u_x_w=0.0_rk     
   endif

    !Free up memory
    if (ndims.eq.2) then
        deallocate(t_temp)
        deallocate(s_temp)
        deallocate(kz_temp)
        deallocate(u_temp)
        deallocate(v_temp)
        if (use_Eair.eq.1) deallocate(Eair_temp)
        if (use_hice.eq.1) deallocate(hice_temp)
    end if
    if (ndims.eq.4) then
        deallocate(t_temp2)
        deallocate(s_temp2)
        deallocate(kz_temp2)
        if (use_Eair.eq.1) deallocate(Eair_temp2)
        if (use_hice.eq.1) deallocate(hice_temp2)
    end if
    deallocate(z_temp)
    deallocate(time_temp)
    deallocate(z_w2)
    deallocate(z_w_error)
    if (nc_staggered_grid.eq.1) deallocate(z_temp2)
    deallocate(cc_temp2)
    deallocate(hmix_rate_temp)

    call check_err(nf90_close(ncid))

    end subroutine input_netcdf
!=======================================================================================================================







!=======================================================================================================================
    subroutine init_netcdf(fn, i_max, k_max, model, use_Eair, use_hice, year)

    !Input variables
    character(len=*), intent(in)     :: fn
    integer, intent(in)              :: k_max, i_max, use_Eair, use_hice, year
    type (type_model),intent(in)     :: model

    !Local variables
    integer                          :: i_dim_id, z_dim_id, z2_dim_id, time_dim_id
    integer                          :: ip, iret, ilast
    integer, parameter               :: time_len = NF90_UNLIMITED
    character(len=4)                 :: yearstr
    integer                          :: dim1d
    integer                          :: dim_ids(3), dim_ids2(3)


    first = .true.
    print *, 'NetCDF version: ', trim(nf90_inq_libvers())
    nc_id = -1
    call check_err(nf90_create(fn, NF90_CLOBBER, nc_id))

    !Define the dimensions
    call check_err(nf90_def_dim(nc_id, "i", i_max, i_dim_id))
    call check_err(nf90_def_dim(nc_id, "z", k_max, z_dim_id))
    call check_err(nf90_def_dim(nc_id, "z2", k_max+1, z2_dim_id))
    call check_err(nf90_def_dim(nc_id, "time", time_len, time_dim_id))

    !Define coordinates
    dim1d = z_dim_id
    call check_err(nf90_def_var(nc_id, "z", NF90_REAL, dim1d, z_id))
    call check_err(nf90_put_att(nc_id, z_id, "positive", "down"))
    call check_err(nf90_put_att(nc_id, z_id, "long_name", "depth at layer midpoints"))
    call check_err(nf90_put_att(nc_id, z_id, "units", "metres"))
    call check_err(nf90_put_att(nc_id, z_id, "axis", "Z"))
    dim1d = z2_dim_id
    call check_err(nf90_def_var(nc_id, "z2", NF90_REAL, dim1d, z2_id))
    call check_err(nf90_put_att(nc_id, z2_id, "positive", "down"))
    call check_err(nf90_put_att(nc_id, z2_id, "long_name", "depth at layer interfaces"))
    call check_err(nf90_put_att(nc_id, z2_id, "units", "metres"))
    call check_err(nf90_put_att(nc_id, z2_id, "axis", "Z"))
    dim1d = i_dim_id
    call check_err(nf90_def_var(nc_id, "i", NF90_REAL, dim1d, i_id))
    call check_err(nf90_put_att(nc_id, i_id, "positive", "right"))
    call check_err(nf90_put_att(nc_id, i_id, "long_name", "distance at horizontal axis"))
    call check_err(nf90_put_att(nc_id, i_id, "units", "metres"))
    call check_err(nf90_put_att(nc_id, i_id, "axis", "X"))
    dim1d = time_dim_id
    call check_err(nf90_def_var(nc_id, "time", NF90_REAL, dim1d, time_id))
    write(yearstr,'(i4)') year
    call check_err(nf90_put_att(nc_id, time_id, "long_name", "time"))
    call check_err(nf90_put_att(nc_id, time_id, "units", "days since "//trim(yearstr)//"-01-01 00:00:00"))
    call check_err(nf90_put_att(nc_id, time_id, "axis", "T"))

    !Define state variables
    dim_ids = (/i_dim_id, z_dim_id, time_dim_id/)
    dim_ids2 = (/i_dim_id, z2_dim_id, time_dim_id/)
    allocate(parameter_id(size(model%state_variables)))
    allocate(parameter_fick_id(size(model%state_variables)))
    allocate(parameter_sink_id(size(model%state_variables)))
    do ip=1,size(model%state_variables)
        ilast = index(model%state_variables(ip)%path,'/',.true.)
        call check_err(nf90_def_var(nc_id, model%state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id(ip)))
        call check_err(nf90_def_var(nc_id, 'fick:'//model%state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_fick_id(ip)))
        call check_err(nf90_def_var(nc_id, 'sink:'//model%state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_sink_id(ip)))
        call check_err(set_attributes(ncid=nc_id, id=parameter_id(ip), units=model%state_variables(ip)%units, &
            long_name=model%state_variables(ip)%long_name, missing_value=model%state_variables(ip)%missing_value))
        call check_err(set_attributes(ncid=nc_id, id=parameter_fick_id(ip), units='mmol/m^2/day', &
            long_name='fick:'//model%state_variables(ip)%long_name,missing_value=model%state_variables(ip)%missing_value))
        call check_err(set_attributes(ncid=nc_id, id=parameter_sink_id(ip), units='mmol/m^2/day', &
            long_name='sink:'//model%state_variables(ip)%long_name,missing_value=model%state_variables(ip)%missing_value))
        call check_err(nf90_put_att(nc_id, parameter_fick_id(ip), "positive", "down"))
        call check_err(nf90_put_att(nc_id, parameter_sink_id(ip), "positive", "down"))
    end do

    !Define diagnostic variables
    allocate(parameter_id_diag(size(model%diagnostic_variables)))
    do ip=1,size(model%diagnostic_variables)
        if (model%diagnostic_variables(ip)%save) then
            ilast = index(model%diagnostic_variables(ip)%path,'/',.true.)
            call check_err(nf90_def_var(nc_id, model%diagnostic_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id_diag(ip)))
            call check_err(set_attributes(ncid=nc_id, id=parameter_id_diag(ip), units=model%diagnostic_variables(ip)%units, &
                long_name=model%diagnostic_variables(ip)%long_name,missing_value=model%diagnostic_variables(ip)%missing_value))
        end if
    end do

    !Define forcing variables used in the run
    call check_err(nf90_def_var(nc_id, "T", NF90_REAL, dim_ids, T_id))
    call check_err(nf90_put_att(nc_id, T_id, "long_name", "temperature"))
    call check_err(nf90_put_att(nc_id, T_id, "units", "degC"))
    call check_err(nf90_def_var(nc_id, "S", NF90_REAL, dim_ids, S_id))
    call check_err(nf90_put_att(nc_id, S_id, "long_name", "salinity"))
    call check_err(nf90_def_var(nc_id, "Kz", NF90_REAL, dim_ids2, Kz_id))
    call check_err(nf90_put_att(nc_id, Kz_id, "long_name", "vertical eddy diffusivity"))
    call check_err(nf90_put_att(nc_id, Kz_id, "units", "m2/s"))
    call check_err(nf90_def_var(nc_id, "Ux", NF90_REAL, dim_ids, u_x_w_id))
    call check_err(nf90_put_att(nc_id, u_x_w_id, "long_name", "horizontal advection"))
    call check_err(nf90_put_att(nc_id, u_x_w_id, "units", "m/s"))
    call check_err(nf90_def_var(nc_id, "Kz_sol", NF90_REAL, dim_ids2, Kz_sol_id))
    call check_err(nf90_put_att(nc_id, Kz_sol_id, "long_name", "total vertical diffusivity of a solute"))
    call check_err(nf90_put_att(nc_id, Kz_sol_id, "units", "m2/s"))
    call check_err(nf90_def_var(nc_id, "Kz_par", NF90_REAL, dim_ids2, Kz_par_id))
    call check_err(nf90_put_att(nc_id, Kz_par_id, "long_name", "total vertical diffusivity of a particulate"))
    call check_err(nf90_put_att(nc_id, Kz_par_id, "units", "m2/s"))
    call check_err(nf90_def_var(nc_id, "w_sol", NF90_REAL, dim_ids2, w_sol_id))
    call check_err(nf90_put_att(nc_id, w_sol_id, "long_name", "total advective velocity of a solute"))
    call check_err(nf90_put_att(nc_id, w_sol_id, "units", "m/s"))
    call check_err(nf90_def_var(nc_id, "w_par", NF90_REAL, dim_ids2, w_par_id))
    call check_err(nf90_put_att(nc_id, w_par_id, "long_name", "total advective velocity of a particulate"))
    call check_err(nf90_put_att(nc_id, w_par_id, "units", "m/s"))
    if (use_Eair.eq.1) then
        call check_err(nf90_def_var(nc_id, "Eair", NF90_REAL, time_dim_id, Eair_id))
        call check_err(nf90_put_att(nc_id, Eair_id, "units", "W/m2"))
    end if
    if (use_hice.eq.1) then
        call check_err(nf90_def_var(nc_id, "hice", NF90_REAL, time_dim_id, hice_id))
        call check_err(nf90_put_att(nc_id, hice_id, "units", "m"))
    end if

    call check_err(nf90_enddef(nc_id))

    end subroutine init_netcdf
!=======================================================================================================================







!=======================================================================================================================
    subroutine save_netcdf(i_max, k_max, julianday, cc, t, s, kz, kzti, wti, model, z, hz, Eair, use_Eair, hice, use_hice, &
        fick_per_day, sink_per_day, ip_sol, ip_par, x, u_x_w, i_day)! i_day here = i_day + 1

    !Input variables
    integer, intent(in)                    :: i_max, k_max, julianday, use_Eair, use_hice, ip_sol, ip_par,i_day
    real(rk), dimension(:,:,:), intent(in) :: cc, t, s, kz, kzti, wti, fick_per_day, sink_per_day, u_x_w
    type (type_model), intent(in)          :: model
    real(rk), dimension(:), intent(in)     :: z, hz, Eair, hice, x

    !Local variables
    integer, dimension(1)                  :: start_z, count_z, start_z2, count_z2, start_time, count_time, start_x, count_x
    integer, dimension(3)                  :: start_cc, count_cc, start_flux, count_flux
    real(rk)                               :: temp_matrix(i_max,k_max), dum(1), z2(k_max+1)
    !Note: The input arguments to nf90_put_var MUST be vectors, even if the length is 1
    !      Removing the dimension(1) or (1) from dum above triggers a spurious error "not finding nf90_put_var"
    integer                                :: ip, i


    !Define nf90_put_var arguments "start" and "count" for z, z2, time, (cc,t,s) and (fick,kz)
    start_z = 1
    count_z = k_max
    start_z2 = 1
    count_z2 = k_max+1
    start_time = i_day !julianday
    count_time = 1
    start_x = 1
    count_x = i_max !number of columns
    start_cc(1:3) = (/1, 1, i_day/) !julianday
    count_cc(1:3) = (/i_max, k_max, 1/)
    start_flux(1:3) = (/1, 1, i_day/) !julianday
    count_flux(1:3) = (/i_max, k_max+1, 1/)

    !At first call only, output depth variable mz = -1*z
    if (first) then
        call check_err(nf90_put_var(nc_id, z_id, z, start_z, count_z))
        z2(1:k_max) = z(1:k_max) - 0.5_rk*hz(1:k_max)
        z2(k_max+1) = z(k_max) + 0.5_rk*hz(k_max)
        call check_err(nf90_put_var(nc_id, z2_id, z2, start_z2, count_z2))
        call check_err(nf90_put_var(nc_id, i_id, x, start_x, count_x))
        !Note: nc_id, z_id and z2_id are available to the entire module and defined in init_netcdf
        first = .false.
    end if

    dum(1) = real(i_day) !julianday
 
    !For all calls output cc, fick, diagnostics and forcings (t,s,kz)
    if (nc_id.ne.-1) then
        call check_err(nf90_put_var(nc_id, time_id, dum, start_time, count_time))
        do ip=1,size(model%state_variables)
            call check_err(nf90_put_var(nc_id, parameter_id(ip), cc(:,:,ip), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, parameter_fick_id(ip), fick_per_day(:,:,ip), start_flux, count_flux))
            call check_err(nf90_put_var(nc_id, parameter_sink_id(ip), sink_per_day(:,:,ip), start_flux, count_flux))
        end do
        do ip=1,size(model%diagnostic_variables)
            if (model%diagnostic_variables(ip)%save) then
                temp_matrix = fabm_get_bulk_diagnostic_data(model,ip)
                if (maxval(abs(temp_matrix)).lt.1.0E37) then
                    !This is to avoid "NetCDF numeric conversion error" for some NetCDF configurations
                    call check_err(nf90_put_var(nc_id, parameter_id_diag(ip), temp_matrix, start_cc, count_cc))
                end if
            end if
        end do
        call check_err(nf90_put_var(nc_id, T_id, t(:,:,julianday), start_cc, count_cc))
        call check_err(nf90_put_var(nc_id, S_id, s(:,:,julianday), start_cc, count_cc))
        call check_err(nf90_put_var(nc_id, Kz_id, kz(:,:,julianday), start_flux, count_flux))
        if (i_max.gt.1) call check_err(nf90_put_var(nc_id, u_x_w_id, u_x_w(:,:,julianday), start_cc, count_cc))
        call check_err(nf90_put_var(nc_id, Kz_sol_id, kzti(:,:,ip_sol), start_flux, count_flux))
        call check_err(nf90_put_var(nc_id, Kz_par_id, kzti(:,:,ip_par), start_flux, count_flux))
        call check_err(nf90_put_var(nc_id, w_sol_id, wti(:,:,ip_sol), start_flux, count_flux))
        call check_err(nf90_put_var(nc_id, w_par_id, wti(:,:,ip_par), start_flux, count_flux))
        if (use_Eair.eq.1) then
            dum(1) = Eair(i_day) !julianday
            call check_err(nf90_put_var(nc_id, Eair_id, dum, start_time, count_time))
        end if
        if (use_hice.eq.1) then
            dum(1) = hice(i_day) !julianday
            call check_err(nf90_put_var(nc_id, hice_id, dum, start_time, count_time))
        end if
        call check_err(nf90_sync(nc_id))
    end if

    end subroutine save_netcdf
!=======================================================================================================================








!=======================================================================================================================
    subroutine close_netcdf()

    implicit none
    if (nc_id.ne.-1) then
        call check_err(nf90_close(nc_id))
        deallocate(parameter_id)
        deallocate(parameter_fick_id)
        deallocate(parameter_sink_id)
        deallocate(parameter_id_diag)
        write (*,'(a)') "finished"
    end if
    nc_id = -1

    end subroutine close_netcdf
!=======================================================================================================================





!=======================================================================================================================
    integer function set_attributes(ncid,id,                         &
                                    units,long_name,                 &
                                    valid_min,valid_max,valid_range, &
                                    scale_factor,add_offset,         &
                                    FillValue,missing_value,         &
                                    C_format,FORTRAN_format)
    !
    ! !DESCRIPTION:
    !  This routine is used to set a number of attributes for
    !  variables. The routine makes heavy use of the {\tt optional} keyword.
    !  The list of recognized keywords is very easy to extend. We have
    !  included a sub-set of the COARDS conventions.
    !
    ! !USES:
    !  IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    integer, intent(in)                     :: ncid,id
    character(len=*), optional              :: units,long_name
    real, optional                          :: valid_min,valid_max
    real, optional                          :: valid_range(2)
    real, optional                          :: scale_factor,add_offset
    double precision, optional              :: FillValue,missing_value
    character(len=*), optional              :: C_format,FORTRAN_format
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding & Hans Burchard
    !
    ! !LOCAL VARIABLES:
    integer                                 :: iret
    real                                    :: vals(2)
    !
    !
    !-----------------------------------------------------------------------
    !
    if (present(units)) then
        iret = nf90_put_att(ncid,id,'units',trim(units))
    end if

    if (present(long_name)) then
        iret = nf90_put_att(ncid,id,'long_name',trim(long_name))
    end if

    if (present(C_format)) then
        iret = nf90_put_att(ncid,id,'C_format',trim(C_format))
    end if

    if (present(FORTRAN_format)) then
        iret = nf90_put_att(ncid,id,'FORTRAN_format',trim(FORTRAN_format))
    end if

    if (present(valid_min)) then
        vals(1) = valid_min
        iret = nf90_put_att(ncid,id,'valid_min',vals(1:1))
    end if

    if (present(valid_max)) then
        vals(1) = valid_max
        iret = nf90_put_att(ncid,id,'valid_max',vals(1:1))
    end if

    if (present(valid_range)) then
        vals(1) = valid_range(1)
        vals(2) = valid_range(2)
        iret = nf90_put_att(ncid,id,'valid_range',vals(1:2))
    end if

    if (present(scale_factor)) then
        vals(1) = scale_factor
        iret = nf90_put_att(ncid,id,'scale_factor',vals(1:1))
    end if

    if (present(add_offset)) then
        vals(1) = add_offset
        iret = nf90_put_att(ncid,id,'add_offset',vals(1:1))
    end if

    if (present(FillValue)) then
        vals(1) = FillValue
        iret = nf90_put_att(ncid,id,'_FillValue',vals(1:1))
    end if

    if (present(missing_value)) then
        vals(1) = missing_value
        iret = nf90_put_att(ncid,id,'missing_value',vals(1:1))
    end if

    set_attributes = 0

    return

    end function set_attributes
!=======================================================================================================================






!=======================================================================================================================
    subroutine check_err(status)

    integer, intent (in) :: status

    if (status .ne. NF90_NOERR) then
        print *, trim(nf90_strerror(status))
        stop
    endif

    end subroutine check_err
!=======================================================================================================================


    end module io_netcdf
