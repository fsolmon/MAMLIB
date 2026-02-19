module mam_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code
!
! adapted from E3SM/CESM/modal_aer_opt.F90 for Geoschem

  use netcdf
  use precision_mod,   only:  r8 => f8
 
  implicit none
  private
  
  ! Public interfaces
  public :: mam_init_opt, mam_aero_sw
  public :: read_species_optics_file
  public :: read_modal_optics_file
  public :: species_optics_properties
  public :: modal_optics_properties
  public :: cleanup_species_properties
  public :: cleanup_modal_properties
  public :: initialize_optics_arrays
  public :: cleanup_optics_arrays
  public ::  mam_optics_diagnostics
  ! Public global variables
  public :: spcopp
  public :: modopp
!  public :: nspecies
!  public :: nmode
  public :: mamoptdiag
  !public :: nswbands, nlwbands conflict with existing in rad_constants 
  
  !==============================================================================
  ! TYPE: species_optics_properties
  ! 
  ! Contains optical properties for a single aerosol species
  ! Read from files like: sulfate_rrtmg_c080918.nc, bcpho_rrtmg_c100508.nc
  !==============================================================================
  type :: species_optics_properties
    
    ! Dimensions
    integer :: nrh              ! Number of RH points
    integer :: nswbands         ! Number of shortwave bands (14 for RRTMG)
    integer :: nlwbands         ! Number of longwave bands (16 for RRTMG)
    
    ! RH grid
    real(r8), allocatable :: rh(:)                    ! (nrh) - RH values [%]
    
    ! Refractive indices (NOT! RH-dependent)
    real(r8), allocatable :: refindex_real_sw(:)    !  nswbands)
    real(r8), allocatable :: refindex_im_sw(:)      !  nswbands)
    real(r8), allocatable :: refindex_real_lw(:)    !  nlwbands)
    real(r8), allocatable :: refindex_im_lw(:)      !  nlwbands)
    
    ! Physical properties
    real(r8) :: specdens        ! Species density [kg/m³]
    real(r8) :: rhcrystal       ! Crystallization RH [%]
    real(r8) :: rhdeliques      ! Deliquescence RH [%]
    
    ! Optical properties (Mie scattering tables)
    real(r8), allocatable :: ext_sw(:,:)              ! (nrh, nswbands) - Extinction
    real(r8), allocatable :: ssa_sw(:,:)              ! (nrh, nswbands) - Single scat. albedo
    real(r8), allocatable :: asm_sw(:,:)              ! (nrh, nswbands) - Asymmetry param
    real(r8), allocatable :: abs_sw(:,:)              ! (nrh, nswbands) - Absorption
    real(r8), allocatable :: ext_lw(:,:)              ! (nrh, nlwbands) - Extinction LW
    real(r8), allocatable :: abs_lw(:,:)              ! (nrh, nlwbands) - Absorption LW
    
  end type species_optics_properties
  
  !==============================================================================
  ! TYPE: modal_optics_properties
  ! 
  ! Contains optical properties for an aerosol mode
  ! Read from files like: mam4_mode1_rrtmg_c141106.nc
  !==============================================================================
  type :: modal_optics_properties
    
    ! Dimensions
    integer :: nswbands         ! Number of shortwave bands (14)
    integer :: nlwbands         ! Number of longwave bands (16)
    integer :: ncoef            ! Number of polynomial coefficients (4)
    integer :: nrefr            ! Number of real refractive index grid points (10)
    integer :: nrefi            ! Number of imag refractive index grid points (10)
    
    ! Mode geometry
    real(r8) :: sigma_logr_aer  ! Geometric standard deviation
    real(r8) :: dgnumlo         ! Lower bound dry diameter [m]
    real(r8) :: dgnumhi         ! Upper bound dry diameter [m]
    
    ! Refractive index grids (for lookup table interpolation)
    ! NetCDF dimension order: (band, refindex)
    ! Fortran declaration: REVERSED for proper memory layout
    real(r8), allocatable :: refrtabsw(:,:)    ! (nrefr, nswbands)
    real(r8), allocatable :: refitabsw(:,:)    ! (nrefi, nswbands)
    real(r8), allocatable :: refrtablw(:,:)    ! (nrefr, nlwbands)
    real(r8), allocatable :: refitablw(:,:)    ! (nrefi, nlwbands)
    
    ! Polynomial coefficient lookup tables
    ! NetCDF dimension order: (band, mode, refindex_im, refindex_real, coef_number)
    ! Fortran declaration: REVERSED for proper memory layout
    real(r8), allocatable :: extpsw(:,:,:,:,:)  ! (ncoef, nrefr, nrefi, 1, nswbands)
    real(r8), allocatable :: abspsw(:,:,:,:,:)  ! (ncoef, nrefr, nrefi, 1, nswbands)
    real(r8), allocatable :: asmpsw(:,:,:,:,:)  ! (ncoef, nrefr, nrefi, 1, nswbands)
    real(r8), allocatable :: extplw(:,:,:,:,:)  ! (ncoef, nrefr, nrefi, 1, nlwbands)
    real(r8), allocatable :: absplw(:,:,:,:,:)  ! (ncoef, nrefr, nrefi, 1, nlwbands)
    
  end type modal_optics_properties
!===============================================================================
! MAM OPTICS DIAGNOSTICS TYPE 
!===============================================================================

type :: mam_optics_diagnostics
   
   ! Extinction for each species across all SW bands (m2/kg air)
   real(r8), allocatable :: vext_sulfate(:,:,:)    ! (pcols, pver, nswbands)
   real(r8), allocatable :: vext_bc(:,:,:)
   real(r8), allocatable :: vext_pom(:,:,:)
   real(r8), allocatable :: vext_soa(:,:,:)
   real(r8), allocatable :: vext_dust(:,:,:)
   real(r8), allocatable :: vext_seasalt(:,:,:)
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
   real(r8), allocatable :: vext_mom(:,:,:)
   real(r8), allocatable :: vext_no3(:,:,:)
   real(r8), allocatable :: vext_nh4(:,:,:)
#elif ( defined MODAL_AERO_4MODE_MOM )
   real(r8), allocatable :: vext_mom(:,:,:)
#endif
   
   ! Single scattering albedo for each species across all SW bands
   real(r8), allocatable :: vssa_sulfate(:,:,:)
   real(r8), allocatable :: vssa_bc(:,:,:)
   real(r8), allocatable :: vssa_pom(:,:,:)
   real(r8), allocatable :: vssa_soa(:,:,:)
   real(r8), allocatable :: vssa_dust(:,:,:)
   real(r8), allocatable :: vssa_seasalt(:,:,:)
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
   real(r8), allocatable :: vssa_mom(:,:,:)
   real(r8), allocatable :: vssa_no3(:,:,:)
   real(r8), allocatable :: vssa_nh4(:,:,:)
#elif ( defined MODAL_AERO_4MODE_MOM )
   real(r8), allocatable :: vssa_mom(:,:,:)
#endif
   
   ! Asymmetry parameter for each species across all SW bands
   real(r8), allocatable :: vasm_sulfate(:,:,:)
   real(r8), allocatable :: vasm_bc(:,:,:)
   real(r8), allocatable :: vasm_pom(:,:,:)
   real(r8), allocatable :: vasm_soa(:,:,:)
   real(r8), allocatable :: vasm_dust(:,:,:)
   real(r8), allocatable :: vasm_seasalt(:,:,:)
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
   real(r8), allocatable :: vasm_mom(:,:,:)
   real(r8), allocatable :: vasm_no3(:,:,:)
   real(r8), allocatable :: vasm_nh4(:,:,:)
#elif ( defined MODAL_AERO_4MODE_MOM )
   real(r8), allocatable :: vasm_mom(:,:,:)
#endif
   
   ! Total mode optical properties across all SW bands
   real(r8), allocatable :: vext_mode(:,:,:)       ! (pcols, pver, nswbands) Total extinction
   real(r8), allocatable :: vssa_mode(:,:,:)        ! Bulk SSA for this mode
   real(r8), allocatable :: vasm_mode(:,:,:)        ! Bulk asymmetry for this mode

   ! Per-mode optical properties across all SW bands
   ! Stored in mamoptdiag(m) : the mode index is implicit in the array subscript.
   ! Dimension order (pcols, pver, nswbands) keeps (i,k) contiguous for the
   ! inner loops and lets isw vary last, matching the loop-nest order.
   real(r8), allocatable :: tauxar(:,:,:)        ! (pcols, pver, nswbands)  opt. depth
   real(r8), allocatable :: ssa(:,:,:)           ! (pcols, pver, nswbands)  single-scat. albedo
   real(r8), allocatable :: g(:,:,:)             ! (pcols, pver, nswbands)  asymmetry param.
   
end type mam_optics_diagnostics

   !==============================================================================
  ! GLOBAL MODULE VARIABLES - Arrays of derived types
  !==============================================================================

  ! Global arrays to store optical properties for all species and modes
  type(species_optics_properties), allocatable :: spcopp(:)  ! (nspecies)
  type(modal_optics_properties), allocatable :: modopp(:)    ! (nmode)
  type(mam_optics_diagnostics), allocatable :: mamoptdiag(:)  ! (nmode)
  !
  integer :: nspecies = 6  ! Number of aerosol species in optic files
  integer :: nmode = 4     ! Number of aerosol modes in optic files 
  integer :: ncoef = 5
  integer :: nswbands,nlwbands
  integer :: iopsulf, iopbc, iopocphi, iopocpho,iopdust,iopsslt  

  complex(r8), allocatable :: crefwsw(:), crefwlw(:)  !
  real(r8) :: xrmin, xrmax

contains
!-------------------------------------------------------------------------------------------------------------`

subroutine mam_init_opt()

!---
! FAB add BLABLA
character(len=50) :: filename
!type(species_optics_properties)  :: spcopt
!type(modal_optics_properties), dimension(4) :: modopt
!---------------
real(r8) :: rmmin, rmmax
integer :: isw

call initialize_optics_arrays(6, 4)

iopsulf = 1
filename = 'sulfate_rrtmg_c080918.nc'
call read_species_optics_file(trim(filename), iopsulf)

iopbc = 2
filename = 'bcpho_rrtmg_c100508.nc'
call read_species_optics_file(trim(filename), iopbc)

iopocphi = 3
filename = 'ocphi_rrtmg_c100508.nc'
call read_species_optics_file(trim(filename), iopocphi)

iopocpho = 4
filename = 'ocpho_rrtmg_c130709.nc'
call read_species_optics_file(trim(filename), iopocpho)

iopdust = 5
filename = 'dust_aeronet_rrtmg_c141106.nc'
call read_species_optics_file(trim(filename), iopdust)

iopsslt = 6
filename = 'ssam_rrtmg_c100508.nc'
call read_species_optics_file(trim(filename), iopsslt)

!------------------------------------------------------

filename = 'mam4_mode1_rrtmg_aeronetdust_c141106.nc'
call read_modal_optics_file(trim(filename), 1)

filename = 'mam4_mode2_rrtmg_aitkendust_c141106.nc'
call read_modal_optics_file(trim(filename), 2)

filename = 'mam4_mode3_rrtmg_aeronetdust_c141106.nc'
call read_modal_optics_file(trim(filename), 3)

filename = 'mam4_mode4_rrtmg_c130628.nc'
call read_modal_optics_file(trim(filename), 4)

!
 rmmin = 0.01e-6_r8
 rmmax = 25.e-6_r8
 xrmin = log(rmmin)
 xrmax = log(rmmax)
! ncoef should be the same for all modes , test it ??

 ncoef = modopp(1)%ncoef
! make sure that nswbands read from the files ( designed for rrtmg) 
! is the same as hard-coded in radprop/ normally yes  
 nswbands = modopp(1)%nswbands
 nlwbands = modopp(1)%nlwbands

 !effective ref index 
 allocate(crefwsw(nswbands))
 allocate(crefwlw(nlwbands))

 ! FAB define pure water ref index, simple but to be re-checked
 ! to be completed for lw
 do isw = 1, nswbands
      crefwsw(isw) = cmplx(1.33_r8, 0.0_r8, kind=r8)
 end do


 ! Allocate and initialize diagnostics 
 call initialize_diagnostics_arrays()


 end subroutine mam_init_opt

subroutine get_species_refind(spectype, refind_sw, refind_lw)
   character(len=*), intent(in) :: spectype
   complex(r8), pointer, intent(out), optional :: refind_sw(:)
   complex(r8), pointer, intent(out), optional :: refind_lw(:)
   
   integer :: iband, ispec_surrogate
   
   ! Determine which species index to use (original or surrogate)
   SELECT CASE (trim(spectype))
      CASE ('sulfate', 'ammonium', 'nitrate')
         ispec_surrogate = iopsulf
      CASE ('black-c')
         ispec_surrogate = iopbc
      CASE ('p-organic', 'm-organic')
         ispec_surrogate = iopocpho
      CASE ('s-organic')
         ispec_surrogate = iopocphi
      CASE ('dust', 'calcium', 'carbonate')
         ispec_surrogate = iopdust
      CASE ('seasalt', 'chloride')
         ispec_surrogate = iopsslt
      CASE DEFAULT
         write(*,*) 'ERROR: Unknown species type: ', trim(spectype)
         stop 'get_species_refind: invalid species'
   END SELECT
  
   ! Fill SW refractive indices if requested
   if (present(refind_sw)) then
      do iband = 1, nswbands
         refind_sw(iband) = cmplx(spcopp(ispec_surrogate)%refindex_real_sw(iband), &
                                  spcopp(ispec_surrogate)%refindex_im_sw(iband), kind=r8)
      end do
   end if
   
   ! Fill LW refractive indices if requested
   if (present(refind_lw)) then
      do iband = 1, nlwbands
         refind_lw(iband) = cmplx(spcopp(ispec_surrogate)%refindex_real_lw(iband), &
                                  spcopp(ispec_surrogate)%refindex_im_lw(iband), kind=r8)
      end do
   end if
   
end subroutine get_species_refind


!-------------------------------------------------------------------------------------


  !==============================================================================
  ! SUBROUTINE: initialize_optics_arrays
  !
  ! PURPOSE: Allocate and initialize global arrays for species and modes
  !
  ! INPUTS:
  !   nspec  - Number of aerosol species
  !   nmod   - Number of aerosol modes
  !
  ! NOTES:
  !   - Call this ONCE before reading any NetCDF files
  !   - This allocates the outer arrays only
  !   - Inner allocatable components are allocated when reading each file
  !==============================================================================
  subroutine initialize_optics_arrays(nspec, nmod)
    
    integer, intent(in) :: nspec  ! Number of species
    integer, intent(in) :: nmod   ! Number of modes
    
    ! Store dimensions as module variables
    nspecies = nspec
    nmode = nmod
    
    ! Allocate the arrays of derived types
    if (allocated(spcopp)) then
      write(*,*) 'WARNING: spcopp already allocated, deallocating first'
      call cleanup_optics_arrays()
    end if
    
    allocate(spcopp(nspecies))
    allocate(modopp(nmode))
    
    write(*,'(A,I3,A)') 'Allocated spcopp for ', nspecies, ' species'
    write(*,'(A,I3,A)') 'Allocated modopp for ', nmode, ' modes'
    
  end subroutine initialize_optics_arrays


  !==============================================================================
  ! SUBROUTINE: cleanup_optics_arrays
  !
  ! PURPOSE: Deallocate global arrays and their components
  !
  !==============================================================================
  subroutine cleanup_optics_arrays()
    
    integer :: i
    
    ! Clean up each species structure
    if (allocated(spcopp)) then
      do i = 1, nspecies
        call cleanup_species_properties(spcopp(i))
      end do
      deallocate(spcopp)
      write(*,*) 'Deallocated spcopp array'
    end if
    
    ! Clean up each mode structure
    if (allocated(modopp)) then
      do i = 1, nmode
        call cleanup_modal_properties(modopp(i))
      end do
      deallocate(modopp)
      write(*,*) 'Deallocated modopp array'
    end if
    
    nspecies = 0
    nmode = 0
    
  end subroutine cleanup_optics_arrays


  !==============================================================================
  ! SUBROUTINE: read_species_optics_file
  !
  ! PURPOSE: Read species-specific aerosol optical properties into global array
  !
  ! INPUTS:
  !   filename - Path to NetCDF file (e.g., 'sulfate_rrtmg_c080918.nc')
  !   ispec    - Index in spcopp array (1 to nspecies)
  !
  ! OUTPUTS:
  !   status   - 0 on success, non-zero on error (optional)
  !
  ! NOTES:
  !   - Must call initialize_optics_arrays() first
  !   - Allocates components of spcopp(ispec) based on NetCDF dimensions
  !   - Follows dimension order exactly as in mam_opt_allocate
  !==============================================================================
  subroutine read_species_optics_file(filename, ispec, status)
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ispec
    integer, intent(out), optional :: status
    
    integer :: ncid, dimid, varid, ierr
    integer :: local_status
    
    local_status = 0
    
    ! Check if arrays are initialized
    if (.not. allocated(spcopp)) then
      write(*,*) 'ERROR: Must call initialize_optics_arrays() first!'
      local_status = -999
      if (present(status)) status = local_status
      return
    end if
    
    ! Check index bounds
    if (ispec < 1 .or. ispec > nspecies) then
      write(*,*) 'ERROR: ispec out of bounds: ', ispec
      local_status = -998
      if (present(status)) status = local_status
      return
    end if
    
    ! Open NetCDF file
    ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR opening file: ', trim(filename)
      write(*,*) 'NetCDF error: ', trim(nf90_strerror(ierr))
      local_status = ierr
      if (present(status)) status = local_status
      return
    end if
    
    write(*,'(A,I3,A,A)') 'Reading species optics file [', ispec, ']: ', trim(filename)
    
    ! Get dimensions from NetCDF
    ierr = nf90_inq_dimid(ncid, 'rh_idx', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension rh_idx not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=spcopp(ispec)%nrh)
    
    ierr = nf90_inq_dimid(ncid, 'sw_band', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension sw_band not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=spcopp(ispec)%nswbands)
    
    ierr = nf90_inq_dimid(ncid, 'lw_band', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension lw_band not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=spcopp(ispec)%nlwbands)
    
    write(*,'(A,3I4)') '  Dimensions: nrh, nswbands, nlwbands = ', &
                       spcopp(ispec)%nrh, spcopp(ispec)%nswbands, spcopp(ispec)%nlwbands
    
    ! ========================================================================
    ! ALLOCATE ARRAYS - Following dimension order from mam_opt_allocate
    ! ========================================================================
    ! NetCDF shows: variable(dim1, dim2) 
    ! Fortran uses: array(dim2, dim1)  <- REVERSED!
    ! ========================================================================
    
    allocate(spcopp(ispec)%rh(spcopp(ispec)%nrh))
    
    ! NetCDF: refindex_real_aer_sw(sw_band) 
    allocate(spcopp(ispec)%refindex_real_sw(spcopp(ispec)%nswbands))
    allocate(spcopp(ispec)%refindex_im_sw(spcopp(ispec)%nswbands))
    
    ! NetCDF: refindex_real_aer_lw(lw_band)
    allocate(spcopp(ispec)%refindex_real_lw( spcopp(ispec)%nlwbands))
    allocate(spcopp(ispec)%refindex_im_lw(spcopp(ispec)%nlwbands))
    
    ! Optical property tables
    allocate(spcopp(ispec)%ext_sw(spcopp(ispec)%nrh, spcopp(ispec)%nswbands))
    allocate(spcopp(ispec)%ssa_sw(spcopp(ispec)%nrh, spcopp(ispec)%nswbands))
    allocate(spcopp(ispec)%asm_sw(spcopp(ispec)%nrh, spcopp(ispec)%nswbands))
    allocate(spcopp(ispec)%abs_sw(spcopp(ispec)%nrh, spcopp(ispec)%nswbands))
    allocate(spcopp(ispec)%ext_lw(spcopp(ispec)%nrh, spcopp(ispec)%nlwbands))
    allocate(spcopp(ispec)%abs_lw(spcopp(ispec)%nrh, spcopp(ispec)%nlwbands))
    
    ! ========================================================================
    ! READ VARIABLES FROM NETCDF
    ! ========================================================================
    
    ! Read RH grid
    ierr = nf90_inq_varid(ncid, 'rh', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%rh)
      write(*,*) '  Read rh'
    else
      write(*,*) '  WARNING: rh variable not found'
    end if
    
    ! Read refractive indices (SW)
    ierr = nf90_inq_varid(ncid, 'refindex_real_aer_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%refindex_real_sw)
      write(*,*) '  Read refindex_real_aer_sw'
    else
      write(*,*) '  WARNING: refindex_real_aer_sw not found'
    end if
    
    ierr = nf90_inq_varid(ncid, 'refindex_im_aer_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%refindex_im_sw)
      write(*,*) '  Read refindex_im_aer_sw'
    else
      write(*,*) '  WARNING: refindex_im_aer_sw not found'
    end if
    
    ! Read refractive indices (LW)
    ierr = nf90_inq_varid(ncid, 'refindex_real_aer_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%refindex_real_lw)
      write(*,*) '  Read refindex_real_aer_lw'
    else
      write(*,*) '  WARNING: refindex_real_aer_lw not found'
    end if
    
    ierr = nf90_inq_varid(ncid, 'refindex_im_aer_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%refindex_im_lw)
      write(*,*) '  Read refindex_im_aer_lw'
    else
      write(*,*) '  WARNING: refindex_im_aer_lw not found'
    end if
    
    ! Read physical properties
    ierr = nf90_inq_varid(ncid, 'density_aer', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%specdens)
      write(*,'(A,F8.2,A)') '  Species density: ', spcopp(ispec)%specdens, ' kg/m³'
    else
      spcopp(ispec)%specdens = 1770.0_r8  ! Default
    end if
    
    ierr = nf90_inq_varid(ncid, 'rhcrystal', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%rhcrystal)
    else
      spcopp(ispec)%rhcrystal = 35.0_r8
    end if
    
    ierr = nf90_inq_varid(ncid, 'rhdeliques', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%rhdeliques)
    else
      spcopp(ispec)%rhdeliques = 80.0_r8
    end if
    
    ! Read optical properties (SW)
    ierr = nf90_inq_varid(ncid, 'ext_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%ext_sw)
      write(*,*) '  Read ext_sw'
    end if
    
    ierr = nf90_inq_varid(ncid, 'ssa_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%ssa_sw)
      write(*,*) '  Read ssa_sw'
    end if
    
    ierr = nf90_inq_varid(ncid, 'asm_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%asm_sw)
      write(*,*) '  Read asm_sw'
    end if
    
    ierr = nf90_inq_varid(ncid, 'abs_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%abs_sw)
      write(*,*) '  Read abs_sw'
    end if
    
    ! Read optical properties (LW)
    ierr = nf90_inq_varid(ncid, 'ext_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%ext_lw)
      write(*,*) '  Read ext_lw'
    end if
    
    ierr = nf90_inq_varid(ncid, 'abs_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, spcopp(ispec)%abs_lw)
      write(*,*) '  Read abs_lw'
    end if
    
    if (local_status == 0) then
      write(*,'(A,I3)') 'Successfully read species optics file into spcopp(', ispec, ')'
    end if
    
999 continue
    ! Close file
    ierr = nf90_close(ncid)
    
    if (present(status)) status = local_status
    
  end subroutine read_species_optics_file


  !==============================================================================
  ! SUBROUTINE: read_modal_optics_file
  !
  ! PURPOSE: Read modal aerosol optical properties into global array
  !
  ! INPUTS:
  !   filename - Path to NetCDF file (e.g., 'mam4_mode1_rrtmg_c141106.nc')
  !   imode    - Index in modopp array (1 to nmode)
  !
  ! OUTPUTS:
  !   status   - 0 on success, non-zero on error (optional)
  !
  ! NOTES:
  !   - Must call initialize_optics_arrays() first
  !   - Allocates components of modopp(imode) based on NetCDF dimensions
  !   - Follows REVERSED dimension order (C -> Fortran conversion)
  !==============================================================================
  subroutine read_modal_optics_file(filename, imode, status)
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: imode
    integer, intent(out), optional :: status
    
    integer :: ncid, dimid, varid, ierr
    integer :: local_status
    
    local_status = 0
    
    ! Check if arrays are initialized
    if (.not. allocated(modopp)) then
      write(*,*) 'ERROR: Must call initialize_optics_arrays() first!'
      local_status = -999
      if (present(status)) status = local_status
      return
    end if
    
    ! Check index bounds
    if (imode < 1 .or. imode > nmode) then
      write(*,*) 'ERROR: imode out of bounds: ', imode
      local_status = -998
      if (present(status)) status = local_status
      return
    end if
    
    ! Open NetCDF file
    ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR opening file: ', trim(filename)
      write(*,*) 'NetCDF error: ', trim(nf90_strerror(ierr))
      local_status = ierr
      if (present(status)) status = local_status
      return
    end if
    
    write(*,'(A,I3,A,A)') 'Reading modal optics file [', imode, ']: ', trim(filename)
    
    ! ========================================================================
    ! GET DIMENSIONS FROM NETCDF
    ! ========================================================================
    
    ierr = nf90_inq_dimid(ncid, 'sw_band', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension sw_band not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=modopp(imode)%nswbands)
    
    ierr = nf90_inq_dimid(ncid, 'lw_band', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension lw_band not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=modopp(imode)%nlwbands)
    
    ierr = nf90_inq_dimid(ncid, 'coef_number', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension coef_number not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=modopp(imode)%ncoef)
    
    ierr = nf90_inq_dimid(ncid, 'refindex_real', dimid)
    if (ierr /= NF90_NOERR) then
      write(*,*) 'ERROR: dimension refindex_real not found'
      local_status = ierr
      goto 999
    end if
    ierr = nf90_inquire_dimension(ncid, dimid, len=modopp(imode)%nrefr)
    
    ierr = nf90_inq_dimid(ncid, 'refindex_im', dimid)
    if (ierr /= NF90_NOERR) then
      ! If not found, assume same as real part
      modopp(imode)%nrefi = modopp(imode)%nrefr
    else
      ierr = nf90_inquire_dimension(ncid, dimid, len=modopp(imode)%nrefi)
    end if

    write(*,'(A,5I4)') '  Dimensions: nswbands, nlwbands, ncoef, nrefr, nrefi = ', &
                       modopp(imode)%nswbands, modopp(imode)%nlwbands, modopp(imode)%ncoef, &
                       modopp(imode)%nrefr, modopp(imode)%nrefi
    
    ! ========================================================================
    ! ALLOCATE ARRAYS - CRITICAL: REVERSE NetCDF dimension order!
    ! ========================================================================
    ! NetCDF (C convention):     variable(band, refindex)
    ! Fortran (column-major):    array(refindex, band)  <- REVERSED!
    ! ========================================================================
    
    ! NetCDF: refindex_real_sw(sw_band, refindex_real) -> Fortran: (nrefr, nswbands)
    allocate(modopp(imode)%refrtabsw(modopp(imode)%nrefr, modopp(imode)%nswbands))
    allocate(modopp(imode)%refitabsw(modopp(imode)%nrefi, modopp(imode)%nswbands))
    allocate(modopp(imode)%refrtablw(modopp(imode)%nrefr, modopp(imode)%nlwbands))
    allocate(modopp(imode)%refitablw(modopp(imode)%nrefi, modopp(imode)%nlwbands))
    
    ! NetCDF: extpsw(sw_band, mode, refindex_im, refindex_real, coef_number)
    ! Fortran: FULLY REVERSED -> (ncoef, nrefr, nrefi, 1, nswbands)
    allocate(modopp(imode)%extpsw(modopp(imode)%ncoef, modopp(imode)%nrefr, &
                                   modopp(imode)%nrefi, 1, modopp(imode)%nswbands))
    allocate(modopp(imode)%abspsw(modopp(imode)%ncoef, modopp(imode)%nrefr, &
                                   modopp(imode)%nrefi, 1, modopp(imode)%nswbands))
    allocate(modopp(imode)%asmpsw(modopp(imode)%ncoef, modopp(imode)%nrefr, &
                                   modopp(imode)%nrefi, 1, modopp(imode)%nswbands))
    allocate(modopp(imode)%extplw(modopp(imode)%ncoef, modopp(imode)%nrefr, &
                                   modopp(imode)%nrefi, 1, modopp(imode)%nlwbands))
    allocate(modopp(imode)%absplw(modopp(imode)%ncoef, modopp(imode)%nrefr, &
                                   modopp(imode)%nrefi, 1, modopp(imode)%nlwbands))
    
    write(*,*) '  Allocated all arrays with REVERSED dimension order'
    
    ! ========================================================================
    ! READ VARIABLES FROM NETCDF
    ! ========================================================================
    
    ! Read mode geometry
    ierr = nf90_inq_varid(ncid, 'sigmag', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%sigma_logr_aer)
      write(*,'(A,F6.2)') '  Geometric std dev: ', modopp(imode)%sigma_logr_aer
    else
      write(*,*) '  WARNING: sigmag not found'
      modopp(imode)%sigma_logr_aer = 1.8_r8  ! Default value
    end if
    
    ierr = nf90_inq_varid(ncid, 'dgnumlo', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%dgnumlo)
      write(*,'(A,ES10.3,A)') '  Diameter lower: ', modopp(imode)%dgnumlo, ' m'
    else
      modopp(imode)%dgnumlo = 0.0_r8
    end if
    
    ierr = nf90_inq_varid(ncid, 'dgnumhi', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%dgnumhi)
      write(*,'(A,ES10.3,A)') '  Diameter upper: ', modopp(imode)%dgnumhi, ' m'
    else
      modopp(imode)%dgnumhi = 0.0_r8
    end if
    
    ! Read refractive index grids
    ierr = nf90_inq_varid(ncid, 'refindex_real_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%refrtabsw)
      write(*,*) '  Read refindex_real_sw grid'
      ! Print first band's grid
      write(*,'(A,10F6.2)') '    Band 1 values: ', modopp(imode)%refrtabsw(:,1)
    else
      write(*,*) '  WARNING: refindex_real_sw not found'
    end if
    
    ierr = nf90_inq_varid(ncid, 'refindex_im_sw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%refitabsw)
      write(*,*) '  Read refindex_im_sw grid'
    else
      write(*,*) '  WARNING: refindex_im_sw not found'
    end if
    
    ierr = nf90_inq_varid(ncid, 'refindex_real_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%refrtablw)
      write(*,*) '  Read refindex_real_lw grid'
    end if
    
    ierr = nf90_inq_varid(ncid, 'refindex_im_lw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%refitablw)
      write(*,*) '  Read refindex_im_lw grid'
    end if
    
    ! Read polynomial coefficient tables (SW)
    ierr = nf90_inq_varid(ncid, 'extpsw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%extpsw)
      write(*,*) '  Read extpsw polynomial table'
    else
      write(*,*) '  ERROR: extpsw not found - this is required!'
      local_status = -1
    end if
    
    ierr = nf90_inq_varid(ncid, 'abspsw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%abspsw)
      write(*,*) '  Read abspsw polynomial table'
    else
      write(*,*) '  ERROR: abspsw not found - this is required!'
      local_status = -1
    end if
    
    ierr = nf90_inq_varid(ncid, 'asmpsw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%asmpsw)
      write(*,*) '  Read asmpsw polynomial table'
    else
      write(*,*) '  ERROR: asmpsw not found - this is required!'
      local_status = -1
    end if
    
    ! Read polynomial coefficient tables (LW)
    ierr = nf90_inq_varid(ncid, 'extplw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%extplw)
      write(*,*) '  Read extplw polynomial table'
    end if
    
    ierr = nf90_inq_varid(ncid, 'absplw', varid)
    if (ierr == NF90_NOERR) then
      ierr = nf90_get_var(ncid, varid, modopp(imode)%absplw)
      write(*,*) '  Read absplw polynomial table'
    end if
    
    if (local_status == 0) then
      write(*,'(A,I3)') 'Successfully read modal optics file into modopp(', imode, ')'
    else
      write(*,*) 'WARNING: Some required variables missing from modal file'
    end if
    
999 continue
    ! Close file
    ierr = nf90_close(ncid)
    
    if (present(status)) status = local_status
    
  end subroutine read_modal_optics_file

  !==============================================================================
  ! SUBROUTINE: cleanup_species_properties
  !
  ! PURPOSE: Deallocate arrays in species_optics_properties structure
  !
  !==============================================================================
  subroutine cleanup_species_properties(props)
    type(species_optics_properties), intent(inout) :: props
    
    if (allocated(props%rh)) deallocate(props%rh)
    if (allocated(props%refindex_real_sw)) deallocate(props%refindex_real_sw)
    if (allocated(props%refindex_im_sw)) deallocate(props%refindex_im_sw)
    if (allocated(props%refindex_real_lw)) deallocate(props%refindex_real_lw)
    if (allocated(props%refindex_im_lw)) deallocate(props%refindex_im_lw)
    if (allocated(props%ext_sw)) deallocate(props%ext_sw)
    if (allocated(props%ssa_sw)) deallocate(props%ssa_sw)
    if (allocated(props%asm_sw)) deallocate(props%asm_sw)
    if (allocated(props%abs_sw)) deallocate(props%abs_sw)
    if (allocated(props%ext_lw)) deallocate(props%ext_lw)
    if (allocated(props%abs_lw)) deallocate(props%abs_lw)
    
  end subroutine cleanup_species_properties

  !==============================================================================
  ! SUBROUTINE: cleanup_modal_properties
  !
  ! PURPOSE: Deallocate arrays in modal_optics_properties structure
  !
  !==============================================================================
  subroutine cleanup_modal_properties(props)
    type(modal_optics_properties), intent(inout) :: props
    
    if (allocated(props%refrtabsw)) deallocate(props%refrtabsw)
    if (allocated(props%refitabsw)) deallocate(props%refitabsw)
    if (allocated(props%refrtablw)) deallocate(props%refrtablw)
    if (allocated(props%refitablw)) deallocate(props%refitablw)
    if (allocated(props%extpsw)) deallocate(props%extpsw)
    if (allocated(props%abspsw)) deallocate(props%abspsw)
    if (allocated(props%asmpsw)) deallocate(props%asmpsw)
    if (allocated(props%extplw)) deallocate(props%extplw)
    if (allocated(props%absplw)) deallocate(props%absplw)
    
  end subroutine cleanup_modal_properties

  
  subroutine initialize_diagnostics_arrays()
    use mam_utils, only: pcols,pver
    integer :: m
    
    if (allocated(mamoptdiag)) then
      write(*,*) 'WARNING: mamoptdiag already allocated'
      call cleanup_diagnostics_arrays()
    end if

    ! Allocate outer array
    allocate(mamoptdiag(nmode))
    ! Allocate inner arrays for each mode
    do m = 1, nmode
      
      ! Sulfate
      allocate(mamoptdiag(m)%vext_sulfate(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_sulfate(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_sulfate(pcols, pver, nswbands))
      
      ! BC
      allocate(mamoptdiag(m)%vext_bc(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_bc(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_bc(pcols, pver, nswbands))
      
      ! POM
      allocate(mamoptdiag(m)%vext_pom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_pom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_pom(pcols, pver, nswbands))
      
      ! SOA
      allocate(mamoptdiag(m)%vext_soa(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_soa(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_soa(pcols, pver, nswbands))
      
      ! Dust
      allocate(mamoptdiag(m)%vext_dust(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_dust(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_dust(pcols, pver, nswbands))
      
      ! Sea salt
      allocate(mamoptdiag(m)%vext_seasalt(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_seasalt(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_seasalt(pcols, pver, nswbands))
      
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
      ! MOM
      allocate(mamoptdiag(m)%vext_mom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_mom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_mom(pcols, pver, nswbands))
      
      ! NO3
      allocate(mamoptdiag(m)%vext_no3(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_no3(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_no3(pcols, pver, nswbands))
      
      ! NH4
      allocate(mamoptdiag(m)%vext_nh4(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_nh4(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_nh4(pcols, pver, nswbands))
#elif ( defined MODAL_AERO_4MODE_MOM )
      ! MOM only
      allocate(mamoptdiag(m)%vext_mom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_mom(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_mom(pcols, pver, nswbands))
#endif
      
      ! Mode totals (all SW bands)
      allocate(mamoptdiag(m)%vext_mode(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vssa_mode(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%vasm_mode(pcols, pver, nswbands))

      ! Per-mode SW optical properties (all bands)
      allocate(mamoptdiag(m)%tauxar(pcols, pver, nswbands))
      allocate(mamoptdiag(m)%ssa   (pcols, pver, nswbands))
      allocate(mamoptdiag(m)%g     (pcols, pver, nswbands))
      ! Initialize to zero
      mamoptdiag(m)%vext_sulfate(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_sulfate(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_sulfate(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_bc(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_bc(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_bc(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_pom(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_pom(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_pom(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_soa(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_soa(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_soa(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_dust(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_dust(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_dust(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_seasalt(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_seasalt(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_seasalt(:,:,:) = 0._r8
      
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
      mamoptdiag(m)%vext_mom(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_mom(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_mom(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_no3(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_no3(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_no3(:,:,:) = 0._r8
      
      mamoptdiag(m)%vext_nh4(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_nh4(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_nh4(:,:,:) = 0._r8
#elif ( defined MODAL_AERO_4MODE_MOM )
      mamoptdiag(m)%vext_mom(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_mom(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_mom(:,:,:) = 0._r8
#endif
      
      mamoptdiag(m)%vext_mode(:,:,:) = 0._r8
      mamoptdiag(m)%vssa_mode(:,:,:) = 0._r8
      mamoptdiag(m)%vasm_mode(:,:,:) = 0._r8

      mamoptdiag(m)%tauxar(:,:,:) = 0._r8
      mamoptdiag(m)%ssa(:,:,:)    = 0._r8
      mamoptdiag(m)%g(:,:,:)      = 0._r8
      
    end do
    
    write(*,*) 'Initialized diagnostics arrays for ', nmode, ' modes'
    
  end subroutine initialize_diagnostics_arrays
  
  
  subroutine cleanup_diagnostics_arrays()
    
    integer :: m
    
    if (.not. allocated(mamoptdiag)) return
    
    do m = 1, nmode
      if (allocated(mamoptdiag(m)%vext_sulfate)) deallocate(mamoptdiag(m)%vext_sulfate)
      if (allocated(mamoptdiag(m)%vssa_sulfate)) deallocate(mamoptdiag(m)%vssa_sulfate)
      if (allocated(mamoptdiag(m)%vasm_sulfate)) deallocate(mamoptdiag(m)%vasm_sulfate)
      
      if (allocated(mamoptdiag(m)%vext_bc)) deallocate(mamoptdiag(m)%vext_bc)
      if (allocated(mamoptdiag(m)%vssa_bc)) deallocate(mamoptdiag(m)%vssa_bc)
      if (allocated(mamoptdiag(m)%vasm_bc)) deallocate(mamoptdiag(m)%vasm_bc)
      
      if (allocated(mamoptdiag(m)%vext_pom)) deallocate(mamoptdiag(m)%vext_pom)
      if (allocated(mamoptdiag(m)%vssa_pom)) deallocate(mamoptdiag(m)%vssa_pom)
      if (allocated(mamoptdiag(m)%vasm_pom)) deallocate(mamoptdiag(m)%vasm_pom)
      
      if (allocated(mamoptdiag(m)%vext_soa)) deallocate(mamoptdiag(m)%vext_soa)
      if (allocated(mamoptdiag(m)%vssa_soa)) deallocate(mamoptdiag(m)%vssa_soa)
      if (allocated(mamoptdiag(m)%vasm_soa)) deallocate(mamoptdiag(m)%vasm_soa)
      
      if (allocated(mamoptdiag(m)%vext_dust)) deallocate(mamoptdiag(m)%vext_dust)
      if (allocated(mamoptdiag(m)%vssa_dust)) deallocate(mamoptdiag(m)%vssa_dust)
      if (allocated(mamoptdiag(m)%vasm_dust)) deallocate(mamoptdiag(m)%vasm_dust)
      
      if (allocated(mamoptdiag(m)%vext_seasalt)) deallocate(mamoptdiag(m)%vext_seasalt)
      if (allocated(mamoptdiag(m)%vssa_seasalt)) deallocate(mamoptdiag(m)%vssa_seasalt)
      if (allocated(mamoptdiag(m)%vasm_seasalt)) deallocate(mamoptdiag(m)%vasm_seasalt)
      
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
      if (allocated(mamoptdiag(m)%vext_mom)) deallocate(mamoptdiag(m)%vext_mom)
      if (allocated(mamoptdiag(m)%vssa_mom)) deallocate(mamoptdiag(m)%vssa_mom)
      if (allocated(mamoptdiag(m)%vasm_mom)) deallocate(mamoptdiag(m)%vasm_mom)
      
      if (allocated(mamoptdiag(m)%vext_no3)) deallocate(mamoptdiag(m)%vext_no3)
      if (allocated(mamoptdiag(m)%vssa_no3)) deallocate(mamoptdiag(m)%vssa_no3)
      if (allocated(mamoptdiag(m)%vasm_no3)) deallocate(mamoptdiag(m)%vasm_no3)
      
      if (allocated(mamoptdiag(m)%vext_nh4)) deallocate(mamoptdiag(m)%vext_nh4)
      if (allocated(mamoptdiag(m)%vssa_nh4)) deallocate(mamoptdiag(m)%vssa_nh4)
      if (allocated(mamoptdiag(m)%vasm_nh4)) deallocate(mamoptdiag(m)%vasm_nh4)
#elif ( defined MODAL_AERO_4MODE_MOM )
      if (allocated(mamoptdiag(m)%vext_mom)) deallocate(mamoptdiag(m)%vext_mom)
      if (allocated(mamoptdiag(m)%vssa_mom)) deallocate(mamoptdiag(m)%vssa_mom)
      if (allocated(mamoptdiag(m)%vasm_mom)) deallocate(mamoptdiag(m)%vasm_mom)
#endif
      
      if (allocated(mamoptdiag(m)%vext_mode)) deallocate(mamoptdiag(m)%vext_mode)
      if (allocated(mamoptdiag(m)%vssa_mode)) deallocate(mamoptdiag(m)%vssa_mode)
      if (allocated(mamoptdiag(m)%vasm_mode)) deallocate(mamoptdiag(m)%vasm_mode)

      if (allocated(mamoptdiag(m)%tauxar)) deallocate(mamoptdiag(m)%tauxar)
      if (allocated(mamoptdiag(m)%ssa))    deallocate(mamoptdiag(m)%ssa)
      if (allocated(mamoptdiag(m)%g))      deallocate(mamoptdiag(m)%g)
    end do
    
    deallocate(mamoptdiag)
    
  end subroutine cleanup_diagnostics_arrays

  
!===============================================================================

subroutine mam_aero_sw(state, tauxar, wa, ga, fa, mamoptdiag) 

!FAB add blalbla



use mam_utils,   only: pcols, pver, top_lev => trop_cloud_top_lev 
use physics_types, only: physics_state
use physics_buffer, only: physics_buffer_desc
use rad_constituents,  only: rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: idx_sw_diag  ! Index for visible band (kept for reference)

! Arguments
type(physics_state), intent(in)  :: state
real(r8), intent(out) :: tauxar(pcols,pver,nswbands)  
real(r8), intent(out) :: wa(pcols,pver,nswbands)      
real(r8), intent(out) :: ga(pcols,pver,nswbands)      
real(r8), intent(out) :: fa(pcols,pver,nswbands)      
type(mam_optics_diagnostics), intent(inout) :: mamoptdiag(:)
! Local variables
type(physics_buffer_desc), pointer :: pbuf(:)

integer :: ncol, lchnk
integer :: nmodes, nspec
integer :: m, k, l, isw, i, nc
integer :: itab(pcols), jtab(pcols)

real(r8), pointer :: dgnumwet(:,:)     
real(r8), pointer :: qaerwat(:,:)      
real(r8), pointer :: specmmr(:,:)      

real(r8) :: sigma_logr_aer             
real(r8) :: radsurf(pcols,pver)        
real(r8) :: logradsurf(pcols,pver)     
real(r8) :: cheb(ncoef,pcols,pver)     

real(r8) :: mass(pcols,pver)           
real(r8) :: air_density(pcols,pver)    

real(r8) :: specdens                   
complex(r8), pointer :: specrefindex(:)  
character(len=32) :: spectype          
real(r8) :: hygro_aer                  

! Volume fractions
real(r8) :: vol(pcols)                 
real(r8) :: dryvol(pcols)              
real(r8) :: watervol(pcols)            
real(r8) :: wetvol(pcols)              

! Complex refractive index
complex(r8) :: crefin(pcols)           
real(r8) :: refr(pcols)                
real(r8) :: refi(pcols)                
real(r8) :: ttab(pcols), utab(pcols)   

! Chebyshev coefficients from bilinear interpolation
real(r8) :: cext(pcols,ncoef)          
real(r8) :: cabs(pcols,ncoef)          
real(r8) :: casm(pcols,ncoef)          

! Parameterized optical properties
real(r8) :: pext(pcols)                
real(r8) :: pabs(pcols)                
real(r8) :: pasm(pcols)                
real(r8) :: palb(pcols)                
real(r8) :: dopaer(pcols)              

! Species-specific scattering and absorption for partitioning
real(r8) :: scatso4(pcols), absso4(pcols), hygroso4(pcols)
real(r8) :: scatbc(pcols), absbc(pcols), hygrobc(pcols)
real(r8) :: scatpom(pcols), abspom(pcols), hygropom(pcols)
real(r8) :: scatsoa(pcols), abssoa(pcols), hygrosoa(pcols)
real(r8) :: scatdust(pcols), absdust(pcols), hygrodust(pcols)
real(r8) :: scatseasalt(pcols), absseasalt(pcols), hygroseasalt(pcols)
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
real(r8) :: scatmom(pcols), absmom(pcols), hygromom(pcols)
real(r8) :: scatno3(pcols), absno3(pcols), hygrono3(pcols)
real(r8) :: scatnh4(pcols), absnh4(pcols), hygronh4(pcols)
#elif ( defined MODAL_AERO_4MODE_MOM )
real(r8) :: scatmom(pcols), absmom(pcols), hygromom(pcols)
#endif

real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
real(r8) :: specrefr, specrefi

! Species optical properties (ext, ssa, asm for each species)
real(r8) :: spec_ext(pcols), spec_ssa(pcols), spec_asm(pcols)

nullify(pbuf)

lchnk = state%lchnk
ncol  = state%ncol

! ========================================================================
! Initialize output variables
! ========================================================================
tauxar(:ncol,:,:) = 0._r8
wa(:ncol,:,:)     = 0._r8
ga(:ncol,:,:)     = 0._r8
fa(:ncol,:,:)     = 0._r8

! Layer 0 is above model top (bottom-up convention)
!tauxar(1:ncol,0,:)  = 0._r8
!wa(1:ncol,0,:)      = 0.925_r8
!ga(1:ncol,0,:)      = 0.850_r8
!fa(1:ncol,0,:)      = 0.7225_r8

! Calculate layer mass and air density
mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

! Allocate species refractive index array
allocate(specrefindex(nswbands))

! ========================================================================
! Get number of modes
! ========================================================================
call rad_cnst_get_info(0, nmodes=nmodes)

! ========================================================================
! Loop over all aerosol modes
! ========================================================================
do m = 1, nmodes

   ! Get wet diameter and water mass for this mode
   dgnumwet => state%dgncur_awet(:,:,m)
   qaerwat  => state%qaerwat(:,:,m)

   ! Get mode geometric standard deviation
   sigma_logr_aer = modopp(m)%sigma_logr_aer

   ! Get number of species in this mode
   call rad_cnst_get_info(0, m, nspec=nspec)

   ! Calculate size parameters
   call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   ! =====================================================================
   ! Loop over shortwave bands
   ! =====================================================================
   do isw = 1, nswbands

      ! Initialize mode-total diagnostics for this band
      mamoptdiag(m)%vext_mode(:ncol,:,isw) = 0._r8
      mamoptdiag(m)%vssa_mode(:ncol,:,isw) = 0._r8
      mamoptdiag(m)%vasm_mode(:ncol,:,isw) = 0._r8

      ! ==================================================================
      ! Loop over vertical levels (BOTTOM-UP: k=1 surface, increasing upward)
      ! ==================================================================
      do k = top_lev, pver ! remember top_lev is set to 1 for mam in gc context( bootom up) .  

         ! Initialize bulk refractive index and volumes
         crefin(:ncol)   = (0._r8, 0._r8)
         dryvol(:ncol)   = 0._r8

         ! Initialize species scattering/absorption buffers (all bands)
         scatso4(:ncol)     = 0._r8
         absso4(:ncol)      = 0._r8
         hygroso4(:ncol)    = 0._r8
         scatbc(:ncol)      = 0._r8
         absbc(:ncol)       = 0._r8
         hygrobc(:ncol)     = 0._r8
         scatpom(:ncol)     = 0._r8
            abspom(:ncol)      = 0._r8
            hygropom(:ncol)    = 0._r8
            scatsoa(:ncol)     = 0._r8
         abssoa(:ncol)      = 0._r8
         hygrosoa(:ncol)    = 0._r8
         scatdust(:ncol)    = 0._r8
         absdust(:ncol)     = 0._r8
         hygrodust(:ncol)   = 0._r8
         scatseasalt(:ncol) = 0._r8
         absseasalt(:ncol)  = 0._r8
         hygroseasalt(:ncol)= 0._r8
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
         scatmom(:ncol)     = 0._r8
         absmom(:ncol)      = 0._r8
         hygromom(:ncol)    = 0._r8
         scatno3(:ncol)     = 0._r8
         absno3(:ncol)      = 0._r8
         hygrono3(:ncol)    = 0._r8
         scatnh4(:ncol)     = 0._r8
         absnh4(:ncol)      = 0._r8
         hygronh4(:ncol)    = 0._r8
#elif ( defined MODAL_AERO_4MODE_MOM )
         scatmom(:ncol)     = 0._r8
         absmom(:ncol)      = 0._r8
         hygromom(:ncol)    = 0._r8
#endif

         ! ===============================================================
         ! Loop over species - accumulate volume-weighted refractive index
         ! ===============================================================
         do l = 1, nspec
            
            ! Get species mass mixing ratio
            call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, specmmr)
            
            ! Get species properties
            call rad_cnst_get_aer_props(0, m, l, density_aer=specdens, &
                                        hygro_aer=hygro_aer,          &
                                        spectype=spectype)

            ! Get species complex refractive index
            call get_species_refind(spectype, refind_sw=specrefindex)

            ! Calculate volume fraction and accumulate bulk refractive index
            do i = 1, ncol
               vol(i)    = specmmr(i,k)/specdens
               dryvol(i) = dryvol(i) + vol(i)
               crefin(i) = crefin(i) + vol(i)*specrefindex(isw)
            end do

            ! ============================================================
            ! Accumulate species-specific scattering/absorption (all bands)
            ! ============================================================
            specrefr = real(specrefindex(isw))
            specrefi = aimag(specrefindex(isw))

            ! Accumulate species-specific scattering and absorption
            ! These will be used to partition the total optical properties
            
#if ( defined MOSAIC_SPECIES )
            if ((trim(spectype) == 'dust') .or. (trim(spectype) == 'carbonate') .or. &
                (trim(spectype) == 'calcium')) then
               do i = 1, ncol
                  scatdust(i)  = scatdust(i) + vol(i)*specrefr
                  absdust(i)   = absdust(i) - vol(i)*specrefi
                  hygrodust(i) = hygrodust(i) + vol(i)*hygro_aer
               end do
            end if
#else
            if (trim(spectype) == 'dust') then
               do i = 1, ncol
                  scatdust(i)  = vol(i)*specrefr
                  absdust(i)   = -vol(i)*specrefi
                  hygrodust(i) = vol(i)*hygro_aer
               end do
            end if
#endif
            
            if (trim(spectype) == 'sulfate') then
               do i = 1, ncol
                  scatso4(i)  = vol(i)*specrefr
                  absso4(i)   = -vol(i)*specrefi
                  hygroso4(i) = vol(i)*hygro_aer
               end do
            end if

#if ( defined MOSAIC_SPECIES )
               if (trim(spectype) == 'nitrate') then
                  do i = 1, ncol
                     scatno3(i)  = vol(i)*specrefr
                     absno3(i)   = -vol(i)*specrefi
                     hygrono3(i) = vol(i)*hygro_aer
                  end do
               end if
               
               if (trim(spectype) == 'ammonium') then
                  do i = 1, ncol
                     scatnh4(i)  = vol(i)*specrefr
                     absnh4(i)   = -vol(i)*specrefi
                     hygronh4(i) = vol(i)*hygro_aer
                  end do
               end if
#endif
            
            if (trim(spectype) == 'black-c') then
               do i = 1, ncol
                  scatbc(i)  = vol(i)*specrefr
                  absbc(i)   = -vol(i)*specrefi
                  hygrobc(i) = vol(i)*hygro_aer
               end do
            end if
            
            if (trim(spectype) == 'p-organic') then
               do i = 1, ncol
                  scatpom(i)  = vol(i)*specrefr
                  abspom(i)   = -vol(i)*specrefi
                  hygropom(i) = vol(i)*hygro_aer
               end do
            end if
            
            if (trim(spectype) == 's-organic') then
               do i = 1, ncol
                  scatsoa(i)  = vol(i)*specrefr
                  abssoa(i)   = -vol(i)*specrefi
                  hygrosoa(i) = vol(i)*hygro_aer
               end do
            end if

#if ( defined MOSAIC_SPECIES )
            if ((trim(spectype) == 'seasalt') .or. (trim(spectype) == 'chloride')) then
               do i = 1, ncol
                  scatseasalt(i)  = scatseasalt(i) + vol(i)*specrefr
                  absseasalt(i)   = absseasalt(i) - vol(i)*specrefi
                  hygroseasalt(i) = hygroseasalt(i) + vol(i)*hygro_aer
               end do
            end if
#else
            if (trim(spectype) == 'seasalt') then
               do i = 1, ncol
                  scatseasalt(i)  = vol(i)*specrefr
                  absseasalt(i)   = -vol(i)*specrefi
                  hygroseasalt(i) = vol(i)*hygro_aer
               end do
            end if
#endif

#if ( defined MODAL_AERO_4MODE_MOM )
            if (trim(spectype) == 'm-organic') then
               do i = 1, ncol
                  scatmom(i)  = vol(i)*specrefr
                  absmom(i)   = -vol(i)*specrefi
                  hygromom(i) = vol(i)*hygro_aer
               end do
            end if
#endif

         end do ! species loop

         ! ===============================================================
         ! Add water contribution to bulk refractive index
         ! ===============================================================
         do i = 1, ncol
            watervol(i) = qaerwat(i,k)/rhoh2o
            wetvol(i)   = watervol(i) + dryvol(i)
            
            if (watervol(i) < 0._r8) then
               watervol(i) = 0._r8
               wetvol(i)   = dryvol(i)
            end if

            ! Volume mixing with water
            crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
            crefin(i) = crefin(i)/max(wetvol(i), 1.e-60_r8)
            
            refr(i) = real(crefin(i))
            refi(i) = abs(aimag(crefin(i)))
         end do

         ! ===============================================================
         ! Interpolate Chebyshev coefficients based on refractive index
         ! ===============================================================
         itab(:ncol) = 0
         
         call binterp(modopp(m)%extpsw(:,:,:,1,isw), ncol, ncoef, &
                      modopp(m)%nrefr, modopp(m)%nrefi, &
                      refr, refi, modopp(m)%refrtabsw(:,isw), &
                      modopp(m)%refitabsw(:,isw), &
                      itab, jtab, ttab, utab, cext)
         
         call binterp(modopp(m)%abspsw(:,:,:,1,isw), ncol, ncoef, &
                      modopp(m)%nrefr, modopp(m)%nrefi, &
                      refr, refi, modopp(m)%refrtabsw(:,isw), &
                      modopp(m)%refitabsw(:,isw), &
                      itab, jtab, ttab, utab, cabs)
         
         call binterp(modopp(m)%asmpsw(:,:,:,1,isw), ncol, ncoef, &
                      modopp(m)%nrefr, modopp(m)%nrefi, &
                      refr, refi, modopp(m)%refrtabsw(:,isw), &
                      modopp(m)%refitabsw(:,isw), &
                      itab, jtab, ttab, utab, casm)

         ! ===============================================================
         ! Calculate parameterized optical properties
         ! ===============================================================
         do i = 1, ncol

            if (logradsurf(i,k) <= xrmax) then
               pext(i) = 0.5_r8*cext(i,1)
               do nc = 2, ncoef
                  pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc)
               end do
               pext(i) = exp(pext(i))
            else
               pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o)
            end if

            pext(i) = pext(i)*wetvol(i)*rhoh2o

            pabs(i) = 0.5_r8*cabs(i,1)
            do nc = 2, ncoef
               pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc)
            end do
            pabs(i) = pabs(i)*wetvol(i)*rhoh2o
            pabs(i) = max(0._r8, pabs(i))
            pabs(i) = min(pext(i), pabs(i))

            pasm(i) = 0.5_r8*casm(i,1)
            do nc = 2, ncoef
               pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc)
            end do

            palb(i) = 1._r8 - pabs(i)/max(pext(i), 1.e-40_r8)
            dopaer(i) = pext(i)*mass(i,k)

         end do

         ! ===============================================================
         ! Partition optical properties by species across all SW bands
         ! ===============================================================
            
            do i = 1, ncol
               
               if (wetvol(i) > 1.e-40_r8) then
                  
                  ! Add water contribution
                  scath2o = watervol(i)*real(crefwsw(isw))
                  absh2o  = -watervol(i)*aimag(crefwsw(isw))
                  
                  ! Calculate totals
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
                  sumscat = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                            scatdust(i) + scatseasalt(i) + scath2o + &
                            scatmom(i) + scatno3(i) + scatnh4(i)
                  sumabs  = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                            absdust(i) + absseasalt(i) + absh2o + &
                            absmom(i) + absno3(i) + absnh4(i)
                  sumhygro = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                             hygrodust(i) + hygroseasalt(i) + &
                             hygromom(i) + hygrono3(i) + hygronh4(i)
#elif ( defined MODAL_AERO_4MODE_MOM )
                  sumscat = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                            scatdust(i) + scatseasalt(i) + scath2o + scatmom(i)
                  sumabs  = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                            absdust(i) + absseasalt(i) + absh2o + absmom(i)
                  sumhygro = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                             hygrodust(i) + hygroseasalt(i) + hygromom(i)
#else
                  sumscat = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                            scatdust(i) + scatseasalt(i) + scath2o
                  sumabs  = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                            absdust(i) + absseasalt(i) + absh2o
                  sumhygro = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                             hygrodust(i) + hygroseasalt(i)
#endif
                  
                  ! Partition optical depth by species
                  ! Each species gets its fraction of scattering/absorption plus
                  ! its fraction of water uptake
                  
                  ! Sulfate
                  spec_ext(i) = ((scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat*palb(i) + &
                                 (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  spec_asm(i) = pasm(i)  ! Assume same asymmetry for all species
                  mamoptdiag(m)%vext_sulfate(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_sulfate(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_sulfate(i,k,isw) = spec_asm(i)
                  
                  ! BC
                  spec_ext(i) = ((scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat*palb(i) + &
                                 (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_bc(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_bc(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_bc(i,k,isw) = pasm(i)
                  
                  ! POM
                  spec_ext(i) = ((scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat*palb(i) + &
                                 (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_pom(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_pom(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_pom(i,k,isw) = pasm(i)
                  
                  ! SOA
                  spec_ext(i) = ((scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat*palb(i) + &
                                 (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_soa(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_soa(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_soa(i,k,isw) = pasm(i)
                  
                  ! Dust
                  spec_ext(i) = ((scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat*palb(i) + &
                                 (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_dust(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_dust(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_dust(i,k,isw) = pasm(i)
                  
                  ! Sea salt
                  spec_ext(i) = ((scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat*palb(i) + &
                                 (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_seasalt(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_seasalt(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_seasalt(i,k,isw) = pasm(i)

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
                  ! MOM
                  spec_ext(i) = ((scatmom(i) + scath2o*hygromom(i)/sumhygro)/sumscat*palb(i) + &
                                 (absmom(i) + absh2o*hygromom(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatmom(i) + scath2o*hygromom(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_mom(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_mom(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_mom(i,k,isw) = pasm(i)
                  
                  ! NO3
                  spec_ext(i) = ((scatno3(i) + scath2o*hygrono3(i)/sumhygro)/sumscat*palb(i) + &
                                 (absno3(i) + absh2o*hygrono3(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatno3(i) + scath2o*hygrono3(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_no3(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_no3(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_no3(i,k,isw) = pasm(i)
                  
                  ! NH4
                  spec_ext(i) = ((scatnh4(i) + scath2o*hygronh4(i)/sumhygro)/sumscat*palb(i) + &
                                 (absnh4(i) + absh2o*hygronh4(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatnh4(i) + scath2o*hygronh4(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_nh4(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_nh4(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_nh4(i,k,isw) = pasm(i)
#elif ( defined MODAL_AERO_4MODE_MOM )
                  ! MOM
                  spec_ext(i) = ((scatmom(i) + scath2o*hygromom(i)/sumhygro)/sumscat*palb(i) + &
                                 (absmom(i) + absh2o*hygromom(i)/sumhygro)/sumabs*(1._r8-palb(i))) * &
                                dopaer(i)
                  spec_ssa(i) = ((scatmom(i) + scath2o*hygromom(i)/sumhygro)/sumscat) * palb(i) / &
                                max(spec_ext(i)/dopaer(i), 1.e-40_r8)
                  mamoptdiag(m)%vext_mom(i,k,isw) = spec_ext(i) / mass(i,k)
                  mamoptdiag(m)%vssa_mom(i,k,isw) = spec_ssa(i)
                  mamoptdiag(m)%vasm_mom(i,k,isw) = pasm(i)
#endif
                  
                  ! Mode totals
                  mamoptdiag(m)%vext_mode(i,k,isw) = dopaer(i) / mass(i,k)
                  mamoptdiag(m)%vssa_mode(i,k,isw) = palb(i)
                  mamoptdiag(m)%vasm_mode(i,k,isw) = pasm(i)
                  
               end if ! wetvol > 0
               
            end do ! ncol

         ! ===============================================================
         ! Store per-mode optical properties (no accumulation here)
         ! ===============================================================
         do i = 1, ncol
            mamoptdiag(m)%tauxar(i,k,isw) = dopaer(i)
            mamoptdiag(m)%ssa(i,k,isw)    = palb(i)
            mamoptdiag(m)%g(i,k,isw)      = pasm(i)
         end do

      end do ! vertical level loop

   end do ! shortwave band loop

end do ! mode loop

! ========================================================================
! Sum per-mode contributions  ->  total column optical properties
!   tauxar  = sum_m  tau_m
!   wa      = sum_m  tau_m * ssa_m            (= sum_m  wa_m)
!   ga      = sum_m  tau_m * ssa_m * g_m      (= sum_m  ga_m)
!   fa      = sum_m  tau_m * ssa_m * g_m^2    (= sum_m  fa_m)
! ========================================================================
do isw = 1, nswbands
   do k =  top_lev, pver
      do i = 1, ncol
         do m = 1, nmodes
            tauxar(i,k,isw) = tauxar(i,k,isw) + mamoptdiag(m)%tauxar(i,k,isw)
            wa(i,k,isw)     = wa(i,k,isw)     + mamoptdiag(m)%tauxar(i,k,isw) &
                                               * mamoptdiag(m)%ssa(i,k,isw)
            ga(i,k,isw)     = ga(i,k,isw)     + mamoptdiag(m)%tauxar(i,k,isw) &
                                               * mamoptdiag(m)%ssa(i,k,isw) &
                                               * mamoptdiag(m)%g(i,k,isw)
            fa(i,k,isw)     = fa(i,k,isw)     + mamoptdiag(m)%tauxar(i,k,isw) &
                                               * mamoptdiag(m)%ssa(i,k,isw) &
                                               * mamoptdiag(m)%g(i,k,isw) &
                                               * mamoptdiag(m)%g(i,k,isw)
         end do
      end do
   end do
end do

! Cleanup
deallocate(specrefindex)

end subroutine mam_aero_sw

  
  !==========================================================================  
subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   use mam_utils, only : pver, top_lev => trop_cloud_top_lev ! 1 in mam/gc context !

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: explnsigma
   real(r8) :: xrad(ncol) ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)
! top_lev = 1 in gc context 
   do k = top_lev, pver
      do i = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(i,k) = 0.5_r8*dgnumwet(i,k)*explnsigma
         logradsurf(i,k) = log(radsurf(i,k))
         ! normalize size parameter
         xrad(i) = max(logradsurf(i,k),xrmin)
         xrad(i) = min(xrad(i),xrmax)
         xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,i,k) = 1._r8
         cheb(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_size_parameters



!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)
      
!     bilinear interpolation of table
!
      use mam_utils, only: pcols, iulog
 
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*) 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp


end module mam_opt
