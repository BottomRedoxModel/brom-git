module io_netcdf

!    use netcdf
!    use fabm, only:type_model,fabm_get_bulk_diagnostic_data
!    use fabm_types, only: rk

!    implicit none
!    !all is private
!    private
!    !public functions
!    public input_netcdf !, netcdf_o

!!    type netcdf_o
!!        private
!!        integer                 :: nc_id  ! NetCDF file id
!!        ! parameter_ids
!!        integer, allocatable	:: parameter_id(:)
!!        integer, allocatable    :: parameter_fick_id(:) ! ?
!!		integer, allocatable    :: parameter_sink_id(:) ! ?
!!		integer, allocatable    :: parameter_id_diag(:) ! ?
!!        integer                 :: i_id, depth_id, time_id, Eair_id, hice_id
!!        ! auxiliary
!!        logical                                 :: first
!!	contains
!!        private
!!        procedure, public 		:: init_netcdf
!!        procedure, public 		:: save_netcdf
!!        procedure, public		:: close_netcdf
!!    end type netcdf_o
    
!!    interface netcdf_o
!!        procedure constructor_netcdf_o
!!    end interface

!contains
  
!    ! inputs hydrophysical data (time, depth, t, s, Kz, u, v, ice, light) from netCDF files
!    subroutine input_netcdf(z_w, dz_w, hz_w, t_w, s_w, kz_w, Eair, use_Eair, &
!							hice, use_hice, gargett_a0, gargett_q, use_gargett, &
!							year, i_water, i_max, days_in_yr, k_wat_bbl, u_x_w) ! check the list of arguments AB
		
!		use io_ascii, only: get_brom_name, get_brom_par ! EYA

!		! Input variables
!		integer, intent(in)                         :: use_Eair, use_hice, year, i_water, i_max, days_in_yr, use_Gargett
!		real(rk), intent(in)                        :: gargett_a0, gargett_q
!		! Output variables
!		real(rk), allocatable, dimension(:), intent(out)        :: z_w, dz_w         ! Layer midpoint depths and spacing between them
!		real(rk), allocatable, dimension(:), intent(out)        :: hz_w              ! Layer thicknesses
!		real(rk), allocatable, dimension(:,:,:), intent(out)    :: t_w, s_w, kz_w    ! Temperature, salinity, and vertical diffusivity
!		real(rk), allocatable, dimension(:,:,:), intent(out)    :: u_x_w             ! Horizontal advection [m/s]
!		real(rk), pointer, dimension(:), intent(out)            :: hice              ! Ice thickness [m]
!		real(rk), allocatable, dimension(:), intent(out)        :: Eair              ! 24-hr average surface downwelling shortwave irradiance in air [W/m2]
!		! Input/output variables
!		integer, intent(out)                        :: k_wat_bbl
!		! Local variables
!		integer                 					:: ncid  ! NetCDF file id
!		! variables ids
!		integer                                     :: depth_varid, time_varid, hice_varid, eair_varid, s_varid, t_varid, Kz_varid, u_varid, v_varid
!		! length of z and time coordinates
!		! dimensions ids
!		integer, dimension(nf90_max_var_dims)       :: dimids
!		integer                                     :: h_rec, time_rec
!		real(rk), allocatable, dimension(:)         :: time_temp, z_temp, hice_temp, Eair_temp, dens
!		real(rk), allocatable, dimension(:,:)       :: t_temp, s_temp, kz_temp, u_temp, v_temp ! where is v component?
!		character(len=64) 	:: file_name 
!		integer           	:: nc_year0, nc_set_k_wat_bbl 
!		real(rk)          	:: time0, timei, time00
!		real(rk)          	:: l ! auxiliary parameter

!		! Input and output file names
!		file_name = get_brom_name("ncinfile_name") ! EYA

!		! Other NetCDF parameters
!		nc_set_k_wat_bbl = get_brom_par("nc_set_k_wat_bbl") ! EYA
!		nc_year0 = get_brom_par("nc_year0") ! EYA

!		! open NetCDF file
!        call check_err(nf90_open(file_name, nf90_nowrite, ncid)) ! arguments: path to file, mode, ncid
        
!		! get variables ids
!        call check_err(nf90_inq_varid(ncid, "depth", depth_varid))
!		call check_err(nf90_inq_varid(ncid, "oc_time", time_varid))
!		call check_err(nf90_inq_varid(ncid, "hice", hice_varid))
!		call check_err(nf90_inq_varid(ncid, "eair", eair_varid))
!		call check_err(nf90_inq_varid(ncid, "salt", s_varid))
!		call check_err(nf90_inq_varid(ncid, "temp", t_varid))
!        call check_err(nf90_inq_varid(ncid, "Kz", Kz_varid))
!        call check_err(nf90_inq_varid(ncid, "u", u_varid))
!		call check_err(nf90_inq_varid(ncid, "v", v_varid))
        
!        ! find parameters dimensions and check their sequence 
!        call check_err(nf90_inquire_variable(ncid, t_varid, dimids = dimids))
!        call check_err(nf90_inquire_dimension(ncid, dimids(1), len = h_rec)) ! get length of z coordinates
!        call check_err(nf90_inquire_dimension(ncid, dimids(2), len = time_rec)) ! get length of time coordinate
!        !DEDUG
!		write(*,*) "h_rec = ", h_rec, "time_rec = ", time_rec
	
!	end subroutine input_netcdf
	
!	! ... to be continued
	
!!	function constructor_netcdf_o()
!!        type(netcdf_o), pointer:: constructor_netcdf_o
!!        allocate(constructor_netcdf_o)
!!    end function constructor_netcdf_o
	
!	subroutine check_err(status)

!        integer, intent (in) :: status

!        if (status .ne. NF90_NOERR) then
!            print *, trim(nf90_strerror(status))
!            stop
!        endif
        
!    end subroutine check_err
    
end module io_netcdf
	
