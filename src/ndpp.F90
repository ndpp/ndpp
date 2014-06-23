module ndpp_class

  use ace,         only: read_ace_table
  use ace_header
  use array_merge
  use chi
  use constants
  use dict_header
  use error,       only: fatal_error, warning
  use global
  use output,      only: write_message, header, time_stamp, print_ascii_array
  use scatt
  use string,      only: lower_case, starts_with, ends_with, to_str
  use thin,        only: thin_grid
  use timer_header

#ifdef MPI
  use mpi
#endif

#ifdef OPENMP
  use omp_lib
#endif

#ifdef HDF5
  use hdf5_interface
#endif

  implicit none
  private
  public  :: nuclearDataPreProc

  type :: nuclearDataPreProc
    ! Status of whether or not the preprocessor has been initialized
    logical              :: is_init = .false.
    ! Path to cross-section xml file to create library for
    character(len=255)   :: path_cross_sections
    ! Number of listings in cross_sections.xml
    integer              :: n_listings    = 0
    ! Start and stop indices of files for me to parse
    integer              :: list_stt      = 0
    integer              :: list_stp      = 0
    ! Incoming energy points for elastic scattering
    real(8), allocatable :: Ein_el(:)
    ! Incoming energy points for inelastic scattering
    real(8), allocatable :: Ein_inel(:)
    ! Energy group structure
    real(8), allocatable :: energy_bins(:)
    ! Number of energy groups
    integer              :: energy_groups = 0
    ! Output library name
    character(len=255)   :: lib_name
    ! Flag to denote library write format
    integer              :: lib_format    = ASCII
    ! Scattering data output type (currently only Legendre or Histogram)
    integer              :: scatt_type    = SCATT_TYPE_LEGENDRE
    ! Scattering data output size (Number of Legendre Orders or Number of Bins)
    integer              :: scatt_order   = SCATT_ORDER_DEFAULT
    ! Whether or not to include neutron multiplication in the distributions
    ! (i.e., whether or not nu-scatter or simply scatter is desired)
    logical              :: nuscatter     = NUSCATTER_DEFAULT
    ! Number of angular bins to use during f_{n,MT} conversion
    integer              :: mu_bins       = MU_BINS_DEFAULT
    ! Flag to integrate chi or not
    logical              :: integrate_chi = INTEGRATE_CHI_DEFAULT
    ! Flag to use free-gas treatment or not
    logical              :: use_freegas   = USE_FREEGAS_DEFAULT
    ! The energy cutoff (divided by kT) to turn off the free gas treatment
    ! (unless overridden by per-nuclide cutoffs)
    real(8)              :: freegas_cutoff = FREEGAS_THRESHOLD_DEFAULT
    ! Minimum group to group transfer probability to bother printing
    real(8)              :: print_tol = PRINT_TOL_DEFAULT
    ! Tolerance on the union energy grid thinning
    real(8)              :: thin_tol
    ! Total Time
    type(Timer)          :: time_total
    ! Initialization Time
    type(Timer)          :: time_initialize
    ! Time to Read XS Data
    type(Timer)          :: time_read_xs
    ! Time to Perform Pre-Processing Calculations
    type(Timer)          :: time_preproc
    ! Time to Perform Scattering Pre-Processing Calculations
    type(Timer)          :: time_scatt_preproc
    ! Time to Perform Chi Pre-Processing Calculations
    type(Timer)          :: time_chi_preproc
    ! Time to Print the Output Data Libraries
    type(Timer)          :: time_print

    ! List of type-bound procedures
    contains
      ! Initialization
      procedure, pass :: init  => init_ndpp
      ! Clear
      procedure, pass :: clear => clear_ndpp
      ! Pre-Process Data for all the nuclides
      ! this will do scattering and chi accordingly
      procedure, pass :: preprocess => preprocess_ndpp
      ! Print the runtime information
      procedure, pass :: print_runtime => print_runtime_ndpp
  end type nuclearDataPreProc

  contains

!===============================================================================
!===============================================================================
! nuclearDataPreProc CLASS OPERATIONS
!===============================================================================
!===============================================================================

!===============================================================================
! INIT_NDPP initializes the data preprocessor object
!===============================================================================

    subroutine init_ndpp(self)

      use xml_data_ndpp_t

      class(nuclearDataPreProc), intent(inout)   :: self ! data preprocessor to initialize
      !character(MAX_FILE_LEN), intent(in) :: path_input ! Path to input file

      logical :: file_exists
      character(MAX_FILE_LEN) :: env_variable
      character(MAX_LINE_LEN) :: filename
      integer :: g

      ! Display output message
      message = "Reading Nuclear Data Pre-Processor XML file..."
      call write_message(5)

      ! Start total timer and initialization timer
      call timer_start(self % time_total)
      call timer_start(self % time_initialize)

      ! Check if ndpp.xml exists
      filename = trim(path_input) // "ndpp.xml"

      inquire(FILE=filename, EXIST=file_exists)
      if (.not. file_exists) then
        message = "Data Pre-Processing XML file '" // trim(filename) // &
                  "' does not exist!"
        call fatal_error()
      end if

      ! Initialize XML scalar variables
      cross_sections_ = ''
      integrate_chi_  = ''
      library_name_   = ''
      output_format_  = ''
      scatt_type_     = 'legendre'
      scatt_order_    = SCATT_ORDER_DEFAULT
      nuscatter_ = ''
      thinning_tol_   = THIN_TOL_DEFAULT
      mu_bins_        = MU_BINS_DEFAULT
      freegas_cutoff_ = FREEGAS_THRESHOLD_DEFAULT
      threads_        = THREADS_DEFAULT

      ! Parse ndpp.xml file
      call read_xml_file_ndpp_t(filename)

      ! Find cross_sections.xml file -- the first place to look is the
      ! settings.xml file. If no file is found there, then we check the
      ! CROSS_SECTIONS environment variable

      if (len_trim(cross_sections_) == 0) then
        ! No cross_sections.xml file specified in settings.xml, check environment
        ! variable
        call get_environment_variable("CROSS_SECTIONS", env_variable)
        if (len_trim(env_variable) == 0) then
          message = "No cross_sections.xml file was specified in " // &
                    "ndpp.xml or in the CROSS_SECTIONS " // &
                    "environment variable."
          call fatal_error()
        else
          self % path_cross_sections = trim(env_variable)
        end if
      else
        self % path_cross_sections = trim(cross_sections_)
      end if

      ! Work with OpenMP threading information
#ifdef OPENMP
      if (threads_ == THREADS_DEFAULT) then
        omp_threads = omp_get_max_threads()
      else
        omp_threads = threads_
      end if
#else
      omp_threads = 1
#endif

      ! Read cross_sections.xml
      call read_cross_sections_xml(self)

      ! Get integrate_chi flag if provided, if not, default provided by class
      call lower_case(integrate_chi_)
      if (integrate_chi_ == '') then
        self % integrate_chi = INTEGRATE_CHI_DEFAULT
      elseif (integrate_chi_ == 'false') then
        self % integrate_chi = .false.
      elseif (integrate_chi_ /= 'true') then
        message = "Value for <integrate_chi> provided, but does not match " // &
                  "TRUE or FALSE. Using default of TRUE."
        call warning()
      end if

      ! Get energy groups and bins
      if (associated(energy_bins_)) then
        ! Lets check to see that the bins are all in increasing order, and positive
        do g = 1, size(energy_bins_) - 1
          if (energy_bins_(g) < ZERO) then
            message = "Invalid energy group structure specified in ndpp.xml; " // &
                      "Groups boundaries be positive."
            call fatal_error()
          end if
          if (energy_bins_(g) >= energy_bins_(g + 1)) then
            message = "Invalid energy group structure specified in ndpp.xml; " // &
                      "Group boundaries must be in increasing order."
            call fatal_error()
          end if
        end do
        ! Ensure that self % energy_bins(1) is 0
        if (energy_bins_(1) /= ZERO) then
          message = "Invalid Lower Energy Boundary: Bottom of Lowest Group " // &
                    " Must be Zero!"
          call fatal_error()
        end if
        self % energy_groups = size(energy_bins_) - 1
        allocate(self % energy_bins(size(energy_bins_)))
        self % energy_bins = energy_bins_
      else
        message = "No energy group structure was specified in ndpp.xml."
        call fatal_error()
      end if

      ! Get the output type, if none is provided, the default is set by the class.
      call lower_case(output_format_)
      if (len_trim(output_format_) > 0) then ! one is provided, make sure it is correct
        if (output_format_ == 'ascii') then
          self % lib_format = ASCII
        elseif (output_format_ == 'binary') then
          self % lib_format = BINARY
        elseif (output_format_ == 'hdf5') then
#ifdef HDF5
          self % lib_format = H5
#else
          message = "Value of HDF5 provided for <output_format>. " // &
                    "NDPP must be compiled with HDF5 enabled."
          call fatal_error()
#endif
        elseif (output_format_ == 'none') then
          self % lib_format = NO_OUT
        elseif (output_format_ == 'human') then
          self % lib_format = HUMAN
          message = "Value of HUMAN provided for <output_format>, " // &
                    "Beware that this output type is incompatible with Monte " // &
                    "Carlo codes."
          call warning()
        else ! incorrect value, print warning, but use default.
          message = "Value for <output_format> provided, but does not match " // &
                    "ASCII, BINARY, HDF5, HUMAN, or NONE. Using default of ASCII."
          call warning()
          self % lib_format = ASCII
        end if
      end if

      ! Get lib_name, if none provided, and not using HDF5,
      ! use the number of groups as ".g###"
      ! The distinction is because HDF5 data all goes in to one file.
      if (len_trim(library_name_) == 0) then
        if (self % lib_format == H5) then
          library_name_ = 'g' // trim(adjustl(to_str(self % energy_groups))) &
                          // '.h5'
        else
          library_name_ = '.g' // trim(adjustl(to_str(self % energy_groups)))
        end if
      end if
      self % lib_name = library_name_

      ! Get scattering type information, if not provided, default provided by class
      call lower_case(scatt_type_)
      if (scatt_type_ == '') then
        self % scatt_type = SCATT_TYPE_DEFAULT
      elseif (scatt_type_ == 'legendre') then
        self % scatt_type = SCATT_TYPE_LEGENDRE
      elseif (scatt_type_ == 'tabular') then
        self % scatt_type = SCATT_TYPE_TABULAR
      end if

      ! Get scattering order information.
      if (scatt_order_ > 0) then
        if ((self % scatt_type == SCATT_TYPE_LEGENDRE) .and. &
          (scatt_order_ > MAX_LEGENDRE_ORDER)) then
          message = "Invalid negative or zero scatt_order value specified in " // &
                    "ndpp.xml."
          call fatal_error()
        else
          self % scatt_order = scatt_order_
        end if
      else
        message = "Invalid negative or zero scatt_order value specified in " // &
                  "ndpp.xml."
        call fatal_error()
      end if

      ! Get nuscatter information
      call lower_case(nuscatter_)
      if (nuscatter_ == '') then
        self % nuscatter = NUSCATTER_DEFAULT
      elseif (nuscatter_ == 'false') then
        self % nuscatter = .false.
      elseif (nuscatter_ == 'true') then
        self % nuscatter = .true.
      else
        message = "Value for <nuscatter> provided, but does not match " // &
                  "TRUE or FALSE. Using default of FALSE."
        call warning()
      end if

      ! Get mu_bins information
      if (mu_bins_ > 1) then
        self % mu_bins = mu_bins_
      else
        message = "Invalid mu_bins value specified in " // &
                  "ndpp.xml. Mu_bins must be two or greater."
        call fatal_error()
      end if

      ! Now get the free-gas threshol
      ! If the user entered a negative, that means they want freegas for the
      ! entire energy range. They can also set the cutoff to a
      ! very large number. The result of this if-block will be the largest
      ! number possible if freegas_cutoff < 0, and the provided value if not.
      if (freegas_cutoff_ == INFINITE_FREEGAS_CUTOFF) then
        self % freegas_cutoff = INFINITY
      else if (freegas_cutoff_ >= ZERO) then
        self % freegas_cutoff = freegas_cutoff_
      else
        message = "Invalid negative value of <freegas_cutoff> specified in " // &
                "ndpp.xml. Specify -1 if no cutoff is desired; all other " // &
                "values are invalid."
        call fatal_error()
      end if

      ! Get grid thinning information
      if (thinning_tol_ < ZERO) then
        message = "Invalid thinning tolerance provided, setting to default" // &
                   " of no thinning."
        call warning()
      end if
      ! Convert from percent to fraction and store
      self % thin_tol = 0.01_8 * thinning_tol_

      ! Get printing tolerance information
      if (print_tol_ > ZERO) then
        self % print_tol = print_tol_
      else
        message = "Invalid printing tolerance provided, setting to default."
        call warning()
      end if

      call partition_work(self % n_listings, self % list_stt, self % list_stp)

      self % is_init = .true.

      ! Stop initialization timer
      call timer_stop(self % time_initialize)

    end subroutine init_ndpp

!===============================================================================
! CLEAR_DPP clears the data preprocessor object
!===============================================================================

    subroutine clear_ndpp(self)
      class(nuclearDataPreProc), intent(inout) :: self ! data preprocessor to clear

      ! Revert to the default/uninitialized values
      self % path_cross_sections = ""
      self % n_listings     = 0
      self % list_stt       = 0
      self % list_stp       = 0
      if (allocated(self % Ein_el)) deallocate(self % Ein_el)
      if (allocated(self % Ein_inel)) deallocate(self % Ein_inel)
      if (allocated(self % energy_bins)) deallocate(self % energy_bins)
      self % energy_groups  = 0
      self % lib_name       = ''
      self % lib_format     = ASCII
      self % scatt_type     = SCATT_TYPE_LEGENDRE
      self % scatt_order    = SCATT_ORDER_DEFAULT
      self % nuscatter      = NUSCATTER_DEFAULT
      self % mu_bins        = MU_BINS_DEFAULT
      self % integrate_chi  = INTEGRATE_CHI_DEFAULT
      self % use_freegas    = USE_FREEGAS_DEFAULT
      self % freegas_cutoff = FREEGAS_THRESHOLD_DEFAULT
      self % print_tol      = PRINT_TOL_DEFAULT
      self % thin_tol       = ZERO

      ! Reset the timers
      call timer_reset(self % time_total)
      call timer_reset(self % time_initialize)
      call timer_reset(self % time_read_xs)
      call timer_reset(self % time_preproc)
      call timer_reset(self % time_scatt_preproc)
      call timer_reset(self % time_chi_preproc)

      ! Complete the clear
      self % is_init = .false.

    end subroutine clear_ndpp

!===============================================================================
! PREPROCESS_DPP Performs the pre-processing of each nuclide requested and
! writes the results to the output library.
!===============================================================================

    subroutine preprocess_ndpp(self)
      class(nuclearDataPreProc), intent(inout) :: self ! data preprocessor to preprocess

      type(Nuclide), pointer    :: nuc => null()  ! Nuclide cross-sections
      type(SAlphaBeta), pointer :: sab => null()  ! S(a,b) tables
      integer                   :: i_listing      ! index of xs_listing
      character(MAX_FILE_LEN)   :: nuc_lib_name   ! nuclidic library's filename
      integer                   :: g              ! Energy group index
      real(8)                   :: thin_compr     ! Thinning compression fraction
      real(8)                   :: thin_err       ! Thinning compression error
      character(MAX_LINE_LEN)   :: xmllib_line    ! XML library tag for nuclides
      character(MAX_LINE_LEN)   :: msg_prepend    ! Text to print before current working lib
      integer, allocatable      :: group_index_el(:) ! Group locations in self % Ein_el
      integer, allocatable      :: group_index_inel(:) ! Group locations in self % Ein_inel
      ! Scattering specific data
      real(8), allocatable :: el_mat(:,:,:)     ! elastic matrix moments,
                                                ! order x g_out x E_in
      real(8), allocatable :: inel_mat(:,:,:)   ! inelastic matrix moments,
                                                ! order x g_out x E_in
      real(8), allocatable :: nuinel_mat(:,:,:) ! nu-inelastic matrix moments,
                                                ! order x g_out x E_in
      real(8), allocatable :: normalization(:)  ! Inelastic normalization data
      ! Chi specific data
      real(8), allocatable   :: Ein_chi(:)   ! List of energy points for chi
      real(8), allocatable   :: chi_t(:,:)   ! grp x E_in chi tot values
      real(8), allocatable   :: chi_p(:,:)   ! grp x E_in chi prompt values
      real(8), allocatable   :: chi_d(:,:,:) ! grp x E_in x precursor chi delayed

      ! S(a,b) specific data
      integer :: islash  ! location in name of slash

      ! Before we begin, write the metadata for the ndpp_lib.xml file.
      ! The ndpp_lib.xml file is NDPPs version of an
      ! xsdata/xsdir/cross_sections.xml file, which describes the nuclear datas
      ! location in a set of files.
      ! We will write the header and metadata here and then, after each nuclide
      ! is complete, print that nuclide's entry.
      ! Finally, we will close the xml file.
      ! We also will write the HDF5 metafile (the library) which will contain
      ! the data for all ACE libraries to be processed. This is different than
      ! the ASCII or binary files; in those cases, each nuclide has its own file.

      if (master) then
        call timer_start(self % time_print)
        call print_ndpp_lib_xml_header(self % n_listings, self % energy_bins, &
          self % lib_format, self % scatt_type, self % scatt_order, &
          self % mu_bins, self % nuscatter, self % integrate_chi, &
          self % print_tol, self % thin_tol)
        call timer_stop(self % time_print)

      ! Display output message
        message = "Beginning Pre-Processing..."
        call write_message(5)

#ifdef OPENMP
        message = "Using " // trim(to_str(omp_threads)) // " OpenMP Threads"
        call write_message(5)
#endif

#ifdef MPI
        message = "Using " // trim(to_str(n_procs)) // " MPI Processes"
        call write_message(5)
#endif

        ! Start PreProcessor Timer
        call timer_start(self % time_preproc)
      end if

#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

      do i_listing = self % list_stt, self % list_stp
        if (mpi_enabled) then
          msg_prepend = "Processor " // trim(adjustl(to_str(rank))) // ":"
        else
          msg_prepend = ""
        end if
        message = trim(adjustl(msg_prepend)) // " Performing Pre-Processing For " // &
          trim(adjustl(xs_listings(i_listing) % name))
        call write_message()
        if (xs_listings(i_listing) % type == ACE_NEUTRON) then
          ! ===================================================================
          ! PERFORM CONTINUOUS ENERGY LIBRARY CALCULATIONS
          allocate(nuclides(1))

          ! Read the ACE library
          call timer_start(self % time_read_xs)
          call read_ace_table(1, i_listing)
          nuc => nuclides(1)
          call timer_stop(self % time_read_xs)

          nuc_lib_name = trim(adjustl(nuc % name)) // &
            trim(adjustl(self % lib_name))

          ! Set the free gas cutoff point (in MeV now, not MeV/kT)
          ! for this nuclide
          if (xs_listings(i_listing) % freegas_cutoff /= &
            GLOBAL_FREEGAS_CUTOFF) then
            if (xs_listings(i_listing) % freegas_cutoff == &
              INFINITE_FREEGAS_CUTOFF) then
              nuc % freegas_cutoff = INFINITY
            else
              nuc % freegas_cutoff = xs_listings(i_listing) % freegas_cutoff * &
                nuc % kT
            end if
          else
            if (self % freegas_cutoff == INFINITY) then
              nuc % freegas_cutoff = INFINITY
            else
              nuc % freegas_cutoff = self % freegas_cutoff * nuc % kT
            end if
          end if
          ! Also set xs_listings()%freegas_cutoff to the new nuclide value
          ! so we can print our library by passing in just xs_listing()
          xs_listings(i_listing) % freegas_cutoff = nuc % freegas_cutoff

          ! Setup output for nuclear data library
          call init_library(self, nuc_lib_name, nuc % name, nuc % kT, &
                            fiss=nuc % fissionable, sab=.false.)

          ! display message
          if (.not. mpi_enabled) then
            message = "....Performing Scattering Integration"
            call write_message(6)
          end if

          ! Integrate Scattering Distributions
          ! Calc_scatt will also update self % Ein_* to include data points necessary
          ! to improve inelastic level scattering interpolation. Cannot do this
          ! a priori since calc_scatt reads in the reactions.
          call timer_start(self % time_scatt_preproc)
          call calc_scatt(nuc, self % energy_bins, self % scatt_type, &
            self % scatt_order, self % mu_bins, self % nuscatter, self % Ein_el, &
            self % Ein_inel, el_mat, inel_mat, nuinel_mat, normalization)

          ! Thin the grid, unless thin_tol is zero
          if (self % thin_tol > ZERO) then
            !For elastic
            call thin_grid(self % Ein_el, el_mat, self % energy_bins, &
                           self % thin_tol, thin_compr, thin_err)
            if (.not. mpi_enabled) then
              ! Report results of thinning
              message = "....Completed Elastic Thinning, Reduced Storage By " // &
                        trim(to_str(100.0_8 * thin_compr)) // "%"
              call write_message(6)
              message = "....Maximum Elastic Thinning Error Was " // &
                        trim(to_str(100.0_8 * thin_err)) // "%"
              call write_message(6)
            end if
            ! And for inelastic
            call thin_grid(self % Ein_inel, inel_mat, self % energy_bins, &
                           self % thin_tol, thin_compr, thin_err, nuinel_mat, &
                           normalization)
            if (.not. mpi_enabled) then
              ! Report results of thinning
              message = "....Completed Elastic Thinning, Reduced Storage By " // &
                        trim(to_str(100.0_8 * thin_compr)) // "%"
              call write_message(6)
              message = "....Maximum Elastic Thinning Error Was " // &
                        trim(to_str(100.0_8 * thin_err)) // "%"
              call write_message(6)
            end if
          end if

          ! Get the energy group boundary indices in self % Ein_* for printing
          allocate(group_index_el(size(self % energy_bins)))
          do g = 1, size(self % energy_bins - 1)
            if (self % energy_bins(g) < self % Ein_el(1)) then
              group_index_el(g) = 1
            else if (self % energy_bins(g) >= self % Ein_el(size(self % Ein_el))) then
              group_index_el(g) = size(self % Ein_el)
            else
              group_index_el(g) = binary_search(self % Ein_el, size(self % Ein_el), &
                                             self % energy_bins(g))
            end if
          end do
          ! Set final group (w/ rounding error this could not include top E_bin pt,
          ! so do it manually (and avoid a search to boot)
          group_index_el(size(self % energy_bins)) = size(self % Ein_el)

          allocate(group_index_inel(size(self % energy_bins)))
          do g = 1, size(self % energy_bins - 1)
            if (self % energy_bins(g) < self % Ein_inel(1)) then
              group_index_inel(g) = 1
            else if (self % energy_bins(g) >= self % Ein_inel(size(self % Ein_inel))) then
              group_index_inel(g) = size(self % Ein_inel)
            else
              group_index_inel(g) = binary_search(self % Ein_inel, size(self % Ein_inel), &
                                             self % energy_bins(g))
            end if
          end do
          ! Set final group (w/ rounding error this could not include top E_bin pt,
          ! so do it manually (and avoid a search to boot)
          group_index_inel(size(self % energy_bins)) = size(self % Ein_inel)

          ! Print the results to file
          call timer_start(self % time_print)
          if (self % nuscatter) then
            call print_scatt(self % lib_format, group_index_el, group_index_inel, &
                             self % Ein_el, self % Ein_inel, self % print_tol, &
                             self % thin_tol, el_mat, inel_mat, normalization, &
                             nuinel_mat)
          else
            call print_scatt(self % lib_format, group_index_el, group_index_inel, &
                             self % Ein_el, self % Ein_inel, self % print_tol, &
                             self % thin_tol, el_mat, inel_mat, normalization)
          end if
          call timer_stop(self % time_print)

          if (allocated(el_mat)) then
            deallocate(el_mat)
          end if
          if (allocated(inel_mat)) then
            deallocate(inel_mat)
          end if
          if (allocated(nuinel_mat)) then
            deallocate(nuinel_mat)
          end if
          call timer_stop(self % time_scatt_preproc)

          ! Integrate Chi
          if (self % integrate_chi) then
            if (.not. mpi_enabled) then
              ! display message
              message = "....Performing Fission Neutron Energy Integration"
              call write_message(6)
            end if

            call timer_start(self % time_chi_preproc)
            if (nuc % fissionable) then
              call calc_chi(nuc, self % energy_bins, Ein_chi, chi_t, chi_p, chi_d)

              ! Print the results to file
              call timer_start(self % time_print)
              call print_chi(nuc % name, self % lib_format, Ein_chi, chi_t, &
                chi_p, chi_d)
              call timer_stop(self % time_print)

              ! Deallocate data created for chi
              deallocate(Ein_chi, chi_t, chi_p, chi_d)
            end if
          end if

          call timer_stop(self % time_chi_preproc)

          nullify(nuc)
          call nuclides(1) % clear()
          deallocate(nuclides)

        else if (xs_listings(i_listing) % type == ACE_THERMAL) then
          ! ===================================================================
          ! PERFORM THERMAL SCATTERING LIBRARY CALCULATIONS
          ! Get the data and then point to it with sab
          allocate(sab_tables(1))
          call timer_start(self % time_read_xs)
          call read_ace_table(1, i_listing)
          sab => sab_tables(1)
          call timer_stop(self % time_read_xs)

          nuc_lib_name = trim(adjustl(sab % name)) // &
            trim(adjustl(self % lib_name))

          xs_listings(i_listing) % alias = xs_listings(i_listing) % name

          ! Convert old S(a,b) filename (from MCNP RSICC distributed data)
          ! to new format (i.e., u/o2.10t to u-o2.10t)
          islash = index(nuc_lib_name, '/')
          if (islash > 1) then
            nuc_lib_name(islash:islash) = '-'
          end if

          ! Setup output for nuclear data library
          call init_library(self, nuc_lib_name, sab % name, sab % kT, sab=.true.)

          if (.not. mpi_enabled) then
            ! display message
            message = "....Performing Scattering Integration"
            call write_message(6)
          end if


          ! Integrate Scattering Distributions
          call timer_start(self % time_scatt_preproc)

          ! Create energy grid to use (inelastic_e_in, elastic_e_in,
          ! energy_bins)
          call sab_egrid(sab, self % energy_bins, self % Ein_el)
          ! Finally add in one point above energy_bins to give MC code something to
          ! interpolate to if Ein==E_bins(size(E_bins))
          call add_one_more_point(self % Ein_el)

          call calc_scattsab(sab, self % energy_bins, self % scatt_type, &
                             self % scatt_order, el_mat, self % mu_bins, &
                             self % Ein_el)

          ! Thin the grid, unless thin_tol is zero
          if (self % thin_tol > ZERO) then
            call thin_grid(self % Ein_el, el_mat, self % energy_bins, &
                           self % thin_tol, thin_compr, thin_err)
            if (.not. mpi_enabled) then
              ! Report results of thinning
              message = "....Completed Thinning, Reduced Storage By " // &
                        trim(to_str(100.0_8 * thin_compr)) // "%"
              call write_message(6)
              message = "....Maximum Thinning Error Was " // &
                        trim(to_str(100.0_8 * thin_err)) // "%"
              call write_message(6)
            end if
          end if

          ! Get the energy group boundary indices in self % Ein for printing
          allocate(group_index_el(size(self % energy_bins)))
          do g = 1, size(self % energy_bins)
            if (self % energy_bins(g) < self % Ein_el(1)) then
              group_index_el(g) = 1
            else if (self % energy_bins(g) >= self % Ein_el(size(self % Ein_el))) then
              group_index_el(g) = size(self % Ein_el)
            else
              group_index_el(g) = binary_search(self % Ein_el, size(self % Ein_el), &
                                             self % energy_bins(g))
            end if
          end do
          ! Set final group (w/ rounding error this could not include top E_bin pt,
          ! so do it manually (and avoid a search to boot)
          group_index_el(size(self % energy_bins)) = size(self % Ein_el)

          ! Print the results to file
          call timer_start(self % time_print)
          call print_scatt(self % lib_format, group_index_el, group_index_inel, &
                           self % Ein_el, self % Ein_inel, self % print_tol, &
                           self % thin_tol, el_mat, inel_mat, normalization)
          call timer_stop(self % time_print)

          if (allocated(el_mat)) then
            deallocate(el_mat)
          end if
          call timer_stop(self % time_scatt_preproc)

          nullify(sab)
          deallocate(sab_tables)
        else
          message = "Invalid Entry in cross_sections listings.  " // &
            "NDPP does not support dosimetry Tables! Entry will be ignored."
          call warning()
        end if

        ! Write this nuclide to the ndpp_lib.xml file
        if (master) then
          call timer_start(self % time_print)

          call print_ndpp_lib_xml_nuclide(nuc_lib_name, &
            self % lib_format, xs_listings(i_listing))
#ifdef MPI
        else
          xmllib_line = write_ndpp_lib_xml_nuclide(nuc_lib_name, &
            self % lib_format, xs_listings(i_listing))
          call MPI_SEND(xmllib_line, MAX_LINE_LEN, MPI_CHARACTER, 0, &
                        i_listing, MPI_COMM_WORLD, mpi_err)
#endif
        end if
        if (master) then
          call timer_stop(self % time_print)
        end if

        ! Close the file or HDF5 group
        call finalize_library(self % lib_format)

        deallocate(self % Ein_el)
        deallocate(self % Ein_inel)
        deallocate(group_index_el)
        deallocate(group_index_inel)
      end do

#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

      ! Now we have to receive all those xmllib_line vals and print them
      if (master) then
        do i_listing = self % list_stp + 1, self % n_listings
          call MPI_RECV(xmllib_line, MAX_LINE_LEN, MPI_CHARACTER, &
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
                        MPI_STATUS_IGNORE, mpi_err)
          write(UNIT_NDPP, '(A)') xmllib_line
        end do
      end if
#endif

      if (master) then
        call timer_stop(self % time_preproc)

        ! Close the ndpp_lib.xml file
        call timer_start(self % time_print)
        call print_ndpp_lib_xml_closer(self % lib_format)
#ifdef HDF5
        ! Close the hdf5 file
        if (self % lib_format == H5) then
          call hdf5_file_close(hdf5_output_file)
        end if
#endif

        call timer_stop(self % time_print)
      end if

    end subroutine preprocess_ndpp

!===============================================================================
! PRINT_RUNTIME_NDPP displays the total time elapsed for the entire run, for
! initialization, and for computation
!===============================================================================

    subroutine print_runtime_ndpp(self)
      use, intrinsic :: ISO_FORTRAN_ENV

      class(nuclearDataPreProc), intent(inout) :: self ! data preprocessor to preprocess
      integer :: ou = OUTPUT_UNIT

      call timer_stop(self % time_total)

      if (master) then
        ! display header block
        call header("Timing Statistics")

        ! display time elapsed for various sections
        if (.not. mpi_enabled) then
          write(ou,100) "Total time for initialization", self % time_initialize % elapsed
          write(ou,100) "Total time for data pre-processing", self % time_preproc % elapsed
          write(ou,100) "  Reading cross sections", self % time_read_xs % elapsed
          write(ou,100) "  Time for scattering integration", &
            self % time_scatt_preproc % elapsed
          write(ou,100) "  Time for chi integration", &
            self % time_chi_preproc % elapsed
        end if
        write(ou,100) "Total time elapsed", self % time_total % elapsed

        ! format for write statements
        100 format (1X,A,T36,"= ",ES11.4," seconds")
      end if

    end subroutine print_runtime_ndpp

!===============================================================================
!===============================================================================
! SUBROUTINES/FUNCTIONS TO SUPPORT nuclearDataPreProc CLASS
!===============================================================================
!===============================================================================

!===============================================================================
! PARTITION_WORK sets the boundaries of the work to be performed by each distrib.
! memory process.
!===============================================================================

  subroutine partition_work(n_listings, stt, stp)
    integer, intent(in)    :: n_listings ! Total # of work items
    integer, intent(inout) :: stt  ! Starting index for this processor to work
    integer, intent(inout) :: stp  ! Final index for this processor to work on

    integer :: work_per ! Work per process

    work_per = n_listings / n_procs

    stt = 1 + rank * work_per
    stp = (rank + 1) * work_per

    if (rank == n_procs - 1) then
      stp = n_listings
    end if

  end subroutine partition_work


!===============================================================================
! PRINT_NDPP_LIB_XML_HEADER prints the metadata for this run of ndpp to
! ndpp_lib.xml in the path of the output libraries.
!===============================================================================

    subroutine print_ndpp_lib_xml_header(n_listings, energy_bins, lib_format, &
      scatt_type, scatt_order, mu_bins, nuscatter, integrate_chi, print_tol, &
      thin_tol)

      integer, intent(in)              :: n_listings     ! Number of entries
      real(8), allocatable, intent(in) :: energy_bins(:) ! Energy group structure
      integer, intent(in)              :: lib_format     ! Library type
      integer, intent(in)              :: scatt_type     ! Representation of scattering data
      integer, intent(in)              :: scatt_order    ! Order of scattering data
      integer, intent(in)              :: mu_bins        ! Number of angular bins to use
      logical, intent(in)              :: nuscatter      ! Flag on if nuscatter data is included
      logical, intent(in)              :: integrate_chi  ! Flag on if chi data is included
      real(8), intent(in)              :: print_tol      ! Minimum g'->g transfer to bother printing
      real(8), intent(in)              :: thin_tol       ! Tolerance on the union energy grid thinning

      character(2) :: indent = '  '

      ! This file is unnecessary if no output is chosen, therefore, exit if that
      ! is the case
      if (lib_format == NO_OUT) return

      ! Open file for writing
      open(FILE="ndpp_lib.xml", UNIT=UNIT_NDPP, STATUS='replace', ACTION='write')
      ! Write xml version header
      write(UNIT_NDPP, '(A)') '<?xml version="1.0"?>'
      ! Create ndpp_lib xml object (open tag)
      write(UNIT_NDPP, '(A)') '<ndpp_lib>'
      ! Write metadata
      write(UNIT_NDPP, '(A)') indent // '<directory> ' // trim(path_input) // &
        '  </directory>'
      if (lib_format == ASCII) then
        write(UNIT_NDPP, '(A)') indent // '<filetype> ascii </filetype>'
      else if (lib_format == BINARY) then
        write(UNIT_NDPP, '(A)') indent // '<filetype> binary </filetype>'
      else if (lib_format == H5) then
        write(UNIT_NDPP, '(A)') indent // '<filetype> hdf5 </filetype>'
      else if (lib_format == HUMAN) then
        write(UNIT_NDPP, '(A)') indent // '<filetype> human </filetype>'
      end if
      write(UNIT_NDPP, '(A)') indent // '<entries> ' // trim(to_str(n_listings))// &
        '  </entries>'
      !!! Skipping over record length for now, not sure I need it.
      if (nuscatter) then
        write(UNIT_NDPP, '(A)') indent // '<nuscatter> true </nuscatter>'
      else
        write(UNIT_NDPP, '(A)') indent // '<nuscatter> false </nuscatter>'
      end if
      if (integrate_chi) then
        write(UNIT_NDPP, '(A)') indent // '<chi_present> true </chi_present>'
      else
        write(UNIT_NDPP, '(A)') indent // '<chi_present> false </chi_present>'
      end if
      write(UNIT_NDPP, '(A)') indent // '<scatt_type> ' // &
        trim(to_str(scatt_type)) // ' </scatt_type>'
      write(UNIT_NDPP, '(A)') indent // '<scatt_order> ' // &
        trim(to_str(scatt_order)) // ' </scatt_order>'
      write(UNIT_NDPP, '(A)') indent // '<print_tol> ' // &
        trim(to_str(print_tol)) // ' </print_tol>'
      write(UNIT_NDPP, '(A)') indent // '<thin_tol> ' // &
        trim(to_str(thin_tol)) // ' </thin_tol>'
      write(UNIT_NDPP, '(A)') indent // '<mu_bins> ' // &
        trim(to_str(mu_bins)) // ' </mu_bins>'
      write(UNIT_NDPP, '(A)') indent // '<energy_bins>'
      call print_ascii_array(energy_bins, UNIT_NDPP)
      write(UNIT_NDPP, '(A)') indent // '</energy_bins>'
      write(UNIT_NDPP, '(A)') ! blank line

    end subroutine print_ndpp_lib_xml_header

!===============================================================================
! PRINT_NDPP_LIB_XML_NUCLIDE prints the entry for each nuclide to the UNIT_NUC
! file.
!===============================================================================

    subroutine print_ndpp_lib_xml_nuclide(filename, lib_format, nuc)
      character(*), intent(in)    :: filename   ! Output filename
      integer, intent(in)         :: lib_format ! Library type
      type(XsListing), intent(in) :: nuc        ! The nuclide to print

      character(2) :: indent = '  '

      ! This file is unnecessary if no output is chosen, therefore, exit if that
      ! is the case
      if (lib_format == NO_OUT) return

      if (nuc % metastable) then ! include metastable attribute
        write(UNIT_NDPP, '(A)') indent // '<ndpp_table alias="' // &
          trim(nuc % alias) // '" awr="' // trim(to_str(nuc % awr)) // &
          '" location="1" name="' // trim(nuc % name) // '" path="' // &
          trim(filename) // '" temperature="' // trim(to_str(nuc % kT)) // &
          '" zaid="' // trim(to_str(nuc % zaid)) // '" metastable= "1" ' // &
          'freegas_cutoff="' // trim(to_str(nuc % freegas_cutoff)) // '"/>'
      else
        write(UNIT_NDPP, '(A)') indent // '<ndpp_table alias="' // &
          trim(nuc % alias) // '" awr="' // trim(to_str(nuc % awr)) // &
          '" location="1" name="' // trim(nuc % name) // '" path="' // &
          trim(filename) // '" temperature="' // trim(to_str(nuc % kT)) // &
          '" zaid="' // trim(to_str(nuc % zaid)) // '" ' // &
          'freegas_cutoff="' // trim(to_str(nuc % freegas_cutoff)) // '"/>'
      end if

    end subroutine print_ndpp_lib_xml_nuclide

!===============================================================================
! WRITE_NDPP_LIB_XML_NUCLIDE prints the entry for each nuclide but to a string
! as opposed to the UNIT_NUC file like PRINT_NDPP_LIB_XML_NUCLIDE.
!===============================================================================

    function write_ndpp_lib_xml_nuclide(filename, lib_format, nuc) result(line)
      character(*), intent(in)    :: filename   ! Output filename
      integer, intent(in)         :: lib_format ! Library type
      type(XsListing), intent(in) :: nuc        ! The nuclide to print

      character(MAX_LINE_LEN) :: line
      character(2) :: indent = '  '

      ! This file is unnecessary if no output is chosen, therefore, exit if that
      ! is the case
      if (lib_format == NO_OUT) return

      if (nuc % metastable) then ! include metastable attribute
        write(line, '(A)') indent // '<ndpp_table alias="' // &
          trim(nuc % alias) // '" awr="' // trim(to_str(nuc % awr)) // &
          '" location="1" name="' // trim(nuc % name) // '" path="' // &
          trim(filename) // '" temperature="' // trim(to_str(nuc % kT)) // &
          '" zaid="' // trim(to_str(nuc % zaid)) // '" metastable= "1" ' // &
          'freegas_cutoff="' // trim(to_str(nuc % freegas_cutoff)) // '"/>'
      else
        write(line, '(A)') indent // '<ndpp_table alias="' // &
          trim(nuc % alias) // '" awr="' // trim(to_str(nuc % awr)) // &
          '" location="1" name="' // trim(nuc % name) // '" path="' // &
          trim(filename) // '" temperature="' // trim(to_str(nuc % kT)) // &
          '" zaid="' // trim(to_str(nuc % zaid)) // '" ' // &
          'freegas_cutoff="' // trim(to_str(nuc % freegas_cutoff)) // '"/>'
      end if

    end function write_ndpp_lib_xml_nuclide

!===============================================================================
! PRINT_NDPP_LIB_XML_CLOSER prints the final tag of ndpp_lib.xml and closes the
! file.
!===============================================================================

    subroutine print_ndpp_lib_xml_closer(lib_format)
      integer, intent(in) :: lib_format     ! Library type
      ! This file is unnecessary if no output is chosen, therefore, exit if that
      ! is the case
      if (lib_format == NO_OUT) return

      write(UNIT_NDPP, '(A)') '</ndpp_lib>'
      close(UNIT_NDPP)

    end subroutine print_ndpp_lib_xml_closer

!===============================================================================
! READ_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the ACE cross sections that may be used. This is
! a near-duplicate of the same routine in input_xml.F90, but is modified to used
! the nuclearDataPreProc structure where it can.
!===============================================================================

    subroutine read_cross_sections_xml(this_ndpp)

      use xml_data_cross_sections_t

      class(nuclearDataPreProc), intent(inout) :: this_ndpp ! data preprocessor to use

      integer :: i           ! loop index
      integer :: filetype    ! default file type
      integer :: recl        ! default record length
      integer :: entries     ! default number of entries
      logical :: file_exists ! does cross_sections.xml exist?
      character(MAX_WORD_LEN)  :: directory ! directory with cross sections
      type(XsListing), pointer :: listing => null()

      ! Check if cross_sections.xml exists
      inquire(FILE=this_ndpp % path_cross_sections, EXIST=file_exists)
      if (.not. file_exists) then
         ! Could not find cross_sections.xml file
         message = "Cross sections XML file '" // trim(this_ndpp % path_cross_sections) // &
              "' does not exist!"
         call fatal_error()
      end if

      message = "Reading cross sections XML file..."
      call write_message(5)

      ! Initialize variables that may go unused
      directory_ = ""
      filetype_ = ""
      record_length_ = 0
      entries_ = 0

      ! Parse cross_sections.xml file
      call read_xml_file_cross_sections_t(this_ndpp % path_cross_sections)

      if (len_trim(directory_) > 0) then
         ! Copy directory information if present
         directory = trim(directory_)
      else
         ! If no directory is listed in cross_sections.xml, by default select the
         ! directory in which the cross_sections.xml file resides
         i = index(this_ndpp % path_cross_sections, "/", BACK=.true.)
         directory = this_ndpp % path_cross_sections(1:i)
      end if

      ! determine whether binary/ascii
      if (filetype_ == 'ascii') then
         filetype = ASCII
      elseif (filetype_ == 'binary') then
         filetype = BINARY
      elseif (len_trim(filetype_) == 0) then
         filetype = ASCII
      else
         message = "Unknown filetype in cross_sections.xml: " // trim(filetype_)
         call fatal_error()
      end if

      ! copy default record length and entries for binary files
      recl = record_length_
      entries = entries_

      ! Allocate xs_listings array
      if (.not. associated(ace_tables_)) then
         message = "No ACE table listings present in cross_sections.xml file!"
         call fatal_error()
      else
         this_ndpp % n_listings = size(ace_tables_)
         allocate(xs_listings(this_ndpp % n_listings))
      end if

      do i = 1, this_ndpp % n_listings
         listing => xs_listings(i)

         ! copy a number of attributes
         listing % name       = trim(ace_tables_(i) % name)
         listing % alias      = trim(ace_tables_(i) % alias)
         listing % zaid       = ace_tables_(i) % zaid
         listing % awr        = ace_tables_(i) % awr
         listing % kT         = ace_tables_(i) % temperature
         listing % location   = ace_tables_(i) % location
         if ((ace_tables_(i) % freegas_cutoff < ZERO) .AND. &
            ((ace_tables_(i) % freegas_cutoff /= INFINITE_FREEGAS_CUTOFF) .AND. &
            (ace_tables_(i) % freegas_cutoff /= GLOBAL_FREEGAS_CUTOFF))) then
           message = "Invalid value of freegas_cutoff element in cross_sections.xml file!"
           call fatal_error()
         else
           listing % freegas_cutoff   = ace_tables_(i) % freegas_cutoff
         end if

         ! determine type of cross section
         if (ends_with(listing % name, 'c')) then
            listing % type = ACE_NEUTRON
         elseif (ends_with(listing % name, 't')) then
            listing % type = ACE_THERMAL
         end if

         ! set filetype, record length, and number of entries
         listing % filetype = filetype
         listing % recl     = recl
         listing % entries  = entries

         ! determine metastable state
         if (ace_tables_(i) % metastable == 0) then
            listing % metastable = .false.
         else
            listing % metastable = .true.
         end if

         ! determine path of cross section table
         if (starts_with(ace_tables_(i) % path, '/')) then
            listing % path = ace_tables_(i) % path
         else
            if (ends_with(directory,'/')) then
               listing % path = trim(directory) // trim(ace_tables_(i) % path)
            else
               listing % path = trim(directory) // '/' // trim(ace_tables_(i) % path)
            end if
         end if

         ! create dictionary entry for both name and alias
         call xs_listing_dict % add_key(listing % name, i)
         call xs_listing_dict % add_key(listing % alias, i)
      end do

    end subroutine read_cross_sections_xml

!===============================================================================
! INIT_LIBRARY writes the header for the NDPP output library; it accepts the
! options and nuclidic information and writes the header lines as appropriate.
!===============================================================================

    subroutine init_library(this_ndpp, filename, name, kT, fiss, sab)
      class(nuclearDataPreProc), intent(in) :: this_ndpp ! NDPP data
      character(*), intent(in)              :: filename  ! output filename
      character(*), intent(in)              :: name      ! Library name
      real(8), intent(in)                   :: kT        ! Library Temperature
      logical, optional, intent(in)         :: fiss      ! Is it fissionable?
      logical, optional, intent(in)         :: sab       ! Is it an S(a,b) table?

      character(MAX_LINE_LEN) :: line
      integer                 :: chi_present_int, nuscatter_int
      logical                 :: fissionable
#ifdef HDF5
      character(MAX_FILE_LEN) :: h_filename
      integer                 :: period_loc
#endif

      ! Deal with fissionable value
      if (.not. present(fiss)) then
        fissionable = .false.
      else
        fissionable = fiss
      end if

      ! First convert the logical value of nuscatter & Chi Present to an integer.
      if (this_ndpp % nuscatter .and. (.not. sab)) then
        nuscatter_int = 1
      else
        nuscatter_int = 0
      end if
      if(this_ndpp % integrate_chi .AND. fissionable) then
        chi_present_int = 1
      else
        chi_present_int = 0
      end if

      if ((this_ndpp % lib_format == ASCII) .or. &
        (this_ndpp % lib_format == HUMAN)) then

        ! Open file for writing
        open(FILE=filename, UNIT=UNIT_NUC, STATUS='replace', ACTION='write')

        ! Write header information:
        ! Nuclide Name, Temperature, Run Date
        line = ''
        write(line,'(A20,1PE20.12,I20,A20)') name, kT, this_ndpp % energy_groups
        write(UNIT_NUC,'(A)') trim(line)
        ! Energy Bin Structure
        call print_ascii_array(this_ndpp % energy_bins, UNIT_NUC)
        ! Scattering Type (Legendre/Hist), Order of this Type, Nu-Scatter,
        ! Chi Present
        line= ''
        write(line,'(I20,I20,I20,I20)') this_ndpp % scatt_type, &
          this_ndpp % scatt_order, nuscatter_int, chi_present_int
        write(UNIT_NUC,'(A)') trim(line)
        ! Do the same, and on a new line, for this_ndpp % mu_bins & thin_tol
        line= ''
        write(line,'(I20,1PE20.12)') this_ndpp % mu_bins, this_ndpp % thin_tol
        write(UNIT_NUC,'(A)') trim(line)

      else if (this_ndpp % lib_format == BINARY) then

        ! Open file for writing
        open(FILE=filename, UNIT=UNIT_NUC, STATUS='replace', ACTION='write', &
          ACCESS = 'stream')

        ! Write header information:
        ! Nuclide Name, Temperature, Run Date
        write(UNIT_NUC) name
        write(UNIT_NUC) kT
        write(UNIT_NUC) this_ndpp % energy_groups

        ! Energy Bin Structure
        write(UNIT_NUC) this_ndpp % energy_bins

        ! Scattering Type (Legendre/Hist), Order of this Type, Nu-Scatter,
        ! Chi Present
        write(UNIT_NUC) this_ndpp % scatt_type
        write(UNIT_NUC) this_ndpp % scatt_order
        write(UNIT_NUC) nuscatter_int
        write(UNIT_NUC) chi_present_int

        ! Write mu_bins, thin_tol
        write(UNIT_NUC) this_ndpp % mu_bins
        write(UNIT_NUC) this_ndpp % thin_tol
#ifdef HDF5
      else if (this_ndpp % lib_format == H5) then
        h_filename = trim(adjustl(name))
        period_loc = scan(h_filename, '.')
        h_filename(period_loc : period_loc) = '_'
        h_filename = "/" // trim(adjustl(h_filename))

        call hdf5_file_create(h_filename, h5_file)

        call hdf5_open_group(h5_file, 'metadata', temp_group)

        ! Write name, kt, energy_groups, energy_bins,
        ! scatt_type, scatt_order, nuscatter, integrate_chi, thin_tol, mu_bins
        ! Now we can print (these should be attributes, but for we need
        ! to first incorporate more routines for this in hdf5_interface)
        call hdf5_write_string(temp_group, 'name', trim(adjustl(name)), &
          len(trim(adjustl(nuc % name))))
        call hdf5_write_double(temp_group, 'kT', kT)
        call hdf5_write_integer(temp_group, 'energy_groups', &
                                this_ndpp % energy_groups)
        call hdf5_write_double_1Darray(temp_group, 'energy_bins', &
          this_ndpp % energy_bins, this_ndpp % energy_groups + 1)
        call hdf5_write_integer(temp_group, 'scatt_type', &
                                this_ndpp % scatt_type)
        call hdf5_write_integer(temp_group, 'scatt_order', &
                                this_ndpp % scatt_order)
        call hdf5_write_integer(temp_group, 'nuscatter', nuscatter_int)
        call hdf5_write_integer(temp_group, 'integrate_chi', chi_present_int)
        call hdf5_write_double(temp_group, 'thin_tol', this_ndpp % thin_tol)
        call hdf5_write_integer(temp_group, 'mu_bins', this_ndpp % mu_bins)

        call hdf5_close_group(temp_group)
#endif
      end if

    end subroutine init_library

!===============================================================================
! FINALIZE_LIBRARY closes the opened library accordingly.
!===============================================================================
    subroutine finalize_library(lib_format)
      integer, intent(in) :: lib_format

      if (lib_format /= H5) then
        close(UNIT_NUC)
#ifdef HDF5
      else
        call hdf5_file_close(h5_file)
#endif
      end if
    end subroutine finalize_library

  end module ndpp_class

