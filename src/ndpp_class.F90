module ndpp_class
  !!! Right now I dont require that the upper bound of the energy bins matches the top of the data.
  !!! e.g., E_max could be 2.0 MeV. But, since the data is still tabulated above that, the E_in
  !!! data in the ACE files will still be collected, but the code is only concerned with the outgoing
  !!! energy groups. I think this is OK, and it should be up to the user to have energy_in in
  !!! OpenMC match energy_out, but... do I want to do anything about that here???
  use ace,              only: read_ace_table
  use ace_header
  use constants
  use dict_header
  use error,            only: fatal_error, warning
  use global,           only: message, path_input, master, xs_listings, &
                              xs_listing_dict, nuclides
  use ndpp_chi
  use output,           only: write_message, header, time_stamp, print_ascii_array
  use ndpp_scatt
  use string,           only: lower_case, starts_with, ends_with, to_str
  use timer_header

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
    ! Energy group structure
    real(8), allocatable :: energy_bins(:)
    ! Number of energy groups
    integer              :: energy_groups = 0
    ! Output library name
    character(len=255)   :: lib_name
    ! Flag to denote if the output library is ASCII or Binary (!!!HDF5 in future?)
    integer              :: lib_format    = ASCII
    ! Scattering data output type (currently only Legendre or Histogram)
    integer              :: scatt_type    = SCATT_TYPE_LEGENDRE
    ! Scattering data output size (Number of Legendre Orders or Number of Bins)
    integer              :: scatt_order   = 5
    ! Number of angular bins to use during f_{n,MT} conversion
    integer              :: mu_bins = 2001
    ! Flag to integrate chi or not
    logical              :: integrate_chi = .true.
    ! Minimum group to group transfer to bother printing
    real(8)              :: print_tol = 1.0E-8_8
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
      integrate_chi_  = 'true'
      library_name_   = ''
      output_format_  = ''
      scatt_type_     = 'legendre'
      scatt_order_    = 5
      thinning_tol_   = 0.2
      mu_bins_        = 2001
      
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
      
      ! Read cross_sections.xml
      call read_cross_sections_xml(self)
      
      ! Get integrate_chi flag if provided
      call lower_case(integrate_chi_)
      if (integrate_chi_ == 'false') then
        self % integrate_chi = .false.
      elseif (integrate_chi_ /= 'true') then
        message = "Value for <integrate_chi> provided, but does not match " // &
                  "TRUE or FALSE. Using default of TRUE."
        call warning()
      end if
      
      ! Get energy groups and bins
      if (associated(energy_bins_)) then
        self % energy_groups = size(energy_bins_) - 1
        allocate(self % energy_bins(size(energy_bins_)))
        self % energy_bins = energy_bins_
      else
        message = "No energy group structure was specified in ndpp.xml."
        call fatal_error()
      end if
      
      ! Get lib_name, if none provided use the number of groups as ".g###"
      if (len_trim(library_name_) == 0) then
        library_name_ = '.g' // trim(adjustl(to_str(self % energy_groups)))
      end if
      self % lib_name = library_name_
      
      ! Get the output type, if none is provided, the default is set by the class.
      call lower_case(output_format_)
      if (len_trim(output_format_) > 0) then ! one is provided, make sure it is correct
        if (output_format_ == 'ascii') then
          self % lib_format = ASCII
        elseif (output_format_ == 'binary') then
          self % lib_format = BINARY
        elseif (output_format_ == 'hdf5') then
          self % lib_format = HDF5
        else ! incorrect value, print warning, but use default.
          message = "Value for <output_format> provided, but does not match " // &
                    "ASCII, BINARY, or HDF5. Using default of ASCII."
          call warning()
          self % lib_format = ASCII
        end if
      end if
      
      ! Get scattering type information.
      call lower_case(scatt_type_)
      if (scatt_type_ == 'tabular') then
        self % scatt_type = SCATT_TYPE_TABULAR
      elseif (scatt_type_ /= 'legendre') then
        message = "Invalid scattering type provided, setting to default of LEGENDRE."
        call warning()
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
      
      ! Get mu_bins information
      if (mu_bins_ > 1) then
        self % mu_bins = mu_bins_
      else
        message = "Invalid mu_bins value specified in " // &
                  "ndpp.xml. Mu_bins must be two or greater."
        call fatal_error()
      end if
      
      ! Get grid thinning information
      if (thinning_tol_ > ZERO) then
        ! Convert from percent to fraction and store
        self % thin_tol = 0.01_8 * thinning_tol_ 
      else 
        message = "Invalid thinning tolerance provided, setting to default of 0.2%."
        call warning()
      end if
      
      ! Get printing tolerance information
      if (print_tol_ > ZERO) then
        self % print_tol = print_tol_ 
      else 
        message = "Invalid printing tolerance provided, setting to default."
        call warning()
      end if
      
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
!       if (associated(self % xs_listings)) deallocate(self % xs_listings)
!       nullify(self % xs_listings)
!       if (associated(self % xs_listing_dict)) deallocate(self % xs_listing_dict)
!       nullify(self % xs_listing_dict)
      self % n_listings    = 0
      if (allocated(self % energy_bins)) deallocate(self % energy_bins)
      self % energy_groups = 0
      self % lib_name  = ''
      self % lib_format    = ASCII
      self % scatt_type    = SCATT_TYPE_LEGENDRE
      self % scatt_order   = 5
      self % integrate_chi = .true.
      self % mu_bins       = 2001
      
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

      type(Nuclide), pointer :: nuc => null() ! Nuclide cross-sections
      integer                :: i_listing     ! index of xs_listing
      ! Scattering specific data
      real(8), allocatable   :: scatt_mat(:,:,:) !scattering matrix moments, 
                                                 ! order x g_out x E_in                                              
      ! Chi specific data
      real(8), allocatable   :: chi_t(:,:)  ! grp x E_in chi tot values on a union grid
      real(8), allocatable   :: e_t_grid(:) ! List of energy points for chi tot
      real(8), allocatable   :: chi_p(:,:)  ! grp x E_in chi prompt values on a union grid
      real(8), allocatable   :: e_p_grid(:) ! List of energy points for chi prompt
      real(8), allocatable   :: chi_d(:,:)  ! grp x E_in chi delayed values on a union grid
      real(8), allocatable   :: e_d_grid(:) ! List of energy points for chi delayed
      
      ! For each xs library requested, this routine will:
      ! 1) Read the ACE file
      ! 2) Calculate the scattering information as follows:
      !   2.a) For each scattering MT of the xs library:
      !     2.a.i) For each angular distribution Ein:
      !       2.a.i.1) Avg nuclide temperature adjustment???
      !       2.a.i.2) Convert from CM-LAB if necessary
      !       2.a.i.3) For each outgoing energy energy *group*
      !         2.a.i.3.a) Do the conversion, integration, etc. (the hard part...)
      !         2.a.i.3.b) Multiply the end-result by the microscopic x/s for this MT
      !     2.a.ii) Sum each (MT,E_in,g_out) set of results to a total for each incoming E
      !   2.b) Normalize?
      ! 3) If requested, calculate the Chi data as follows:
      !   3.a) For each Fission MT
      !     3.a.i) For each incoming energy:
      !       3.a.i.1) Integrate the chi(MT,Eout) over each energy group, store.
      !       3.a.i.1) Multiply end result by microscopic x/s for this MT
      !     3.a.ii) Sum each (MT,E_in,g_out) results to a total for each incoming E
      !   3.b) Normalize?
      ! 4) Print the data according to the requested filename and format type
      
      ! Display output message
      message = "Beginning Pre-Processing..."
      call write_message(5)
      
      ! Start PreProcessor Timer
      call timer_start(self % time_preproc)
      
      
      do i_listing = 1, self % n_listings
        if (xs_listings(i_listing) % type == ACE_NEUTRON) then
          allocate(nuclides(1))
          
          ! Read the ACE library
          call timer_start(self % time_read_xs)
          call read_ace_table(1, i_listing)
          nuc => nuclides(1)
          call timer_stop(self % time_read_xs)
          
          ! Setup output library
          call init_library(self, nuc)
          
          ! Integrate Scattering Distributions
          call timer_start(self % time_scatt_preproc)
          call calc_scatt(nuc, self % energy_bins, self % scatt_type, &
            self % scatt_order, scatt_mat, self % mu_bins, self % thin_tol)
            
          ! Print the results to file
          call timer_start(self % time_print)
          call print_scatt(self % lib_format, scatt_mat, nuc % energy, &
            self % print_tol)
          call timer_stop(self % time_print)
          
          if (allocated(scatt_mat)) then
            deallocate(scatt_mat)
          end if
          call timer_stop(self % time_scatt_preproc)
          
          ! Integrate Chi
          if (self % integrate_chi) then
            call timer_start(self % time_chi_preproc)
            if (nuc % fissionable) then
              call calc_chis(nuc, self % energy_bins, self % energy_groups, chi_t, &
                chi_p, chi_d, e_t_grid, e_p_grid, e_d_grid, self % thin_tol)
              
              ! Print the results to file
              call timer_start(self % time_print)
              call print_chi(self % lib_format, chi_t, chi_p, chi_d, e_t_grid, &
                e_p_grid, e_d_grid)
              call timer_stop(self % time_print)
            end if
          end if
          
          if (allocated(chi_t)) deallocate(chi_t)
          if (allocated(e_t_grid)) deallocate(e_t_grid)
          if (allocated(chi_p)) deallocate(chi_p)
          if (allocated(e_p_grid)) deallocate(e_p_grid)
          if (allocated(chi_d)) deallocate(chi_d)
          if (allocated(e_d_grid)) deallocate(e_d_grid)
          
          call timer_stop(self % time_chi_preproc)
        else
          message = "Invalid Entry in cross_sections listings.  " // &
            "NDPP does not support S(A,B) tables and Dosimetry Tables!"
          call warning()
        end if
        ! Close the file/HDF5 object
        !!!call finalize_output()
        ! Deallocate the nuclide and its member data.
        if (allocated(nuclides)) then
          call nuclides(1) % clear()
          deallocate(nuclides)
        end if
        close(UNIT_NUC)
      end do
      
      call timer_stop(self % time_preproc)
      
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
      
      ! display header block
      call header("Timing Statistics")

      ! display time elapsed for various sections
      write(ou,100) "Total time for initialization", self % time_initialize % elapsed
      write(ou,100) "  Reading cross sections", self % time_read_xs % elapsed
      write(ou,100) "Total time for data pre-processing", self % time_preproc % elapsed
      write(ou,100) "  Time performing scattering calculations", &
        self % time_scatt_preproc % elapsed
      write(ou,100) "  Time performing chi integration calculations", &
        self % time_chi_preproc % elapsed
      write(ou,100) "Total time elapsed", self % time_total % elapsed
      
      ! format for write statements
  100 format (1X,A,T36,"= ",ES11.4," seconds")    

    end subroutine print_runtime_ndpp
  
!===============================================================================
!===============================================================================
! SUBROUTINES/FUNCTIONS TO SUPPORT nuclearDataPreProc CLASS
!===============================================================================
!===============================================================================
  
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
    
    subroutine init_library(this_ndpp, nuc)
      type(nuclearDataPreProc), intent(in) :: this_ndpp ! NDPP data
      type(Nuclide), pointer, intent(in)   :: nuc       ! Nuclide data
      
      character(MAX_LINE_LEN) :: line
      character(MAX_FILE_LEN) :: filename
      integer                 :: chi_present_int
      
      if (this_ndpp % lib_format == ASCII) then
        ! Create filename for output library
        filename = trim(adjustl(nuc % name)) // trim(adjustl(this_ndpp % lib_name))
        
        ! Open file for writing
        open(FILE=filename, UNIT=UNIT_NUC, STATUS='replace', ACTION='write')
        
        ! Write header information:
        ! Nuclide Name, Temperature, Run Date
        line = ''
        write(line,'(A20,1PE20.12,A20,A20)') nuc % name, nuc % kT, ' ', time_stamp()
        write(UNIT_NUC,'(A)') trim(line)
        ! Energy Bin Structure
        call print_ascii_array(this_ndpp % energy_bins, UNIT_NUC)
        ! Scattering Type (Legendre/Hist), Order of this Type, Chi Present, Thinning Tolerance
        ! First convert the logical value of Chi Present to an integer. It seemas as if a type-cast
        ! is not in the standard, so the next if-then  block will explicitly do the cast.
        if(this_ndpp % integrate_chi .AND. nuc % fissionable) then 
          chi_present_int = 1
        else
          chi_present_int = 0
        end if
        ! Now print the results
        line= ''
        write(line,'(I20,I20,I20,1PE20.12)') this_ndpp % scatt_type, &
          this_ndpp % scatt_order, chi_present_int, this_ndpp % thin_tol
        write(UNIT_NUC,'(A)') trim(line)
        ! Do the same, and on a new line, for this_ndpp % mu_bins
        line= ''
        write(line,'(I20)') this_ndpp % mu_bins
        write(UNIT_NUC,'(A)') trim(line)
      else if (this_ndpp % lib_format == BINARY) then
        ! Create filename for output library
        filename = trim(adjustl(nuc % name)) // trim(adjustl(this_ndpp % lib_name))
        
        ! Open file for writing
        open(FILE=filename, UNIT=UNIT_NUC, STATUS='replace', ACTION='write', &
          ACCESS = 'stream')
        
        ! Write header information:
        ! Nuclide Name, Temperature, Run Date
        write(UNIT_NUC) nuc % name
        write(UNIT_NUC) nuc % kT
        write(UNIT_NUC) time_stamp()
        
        ! Energy Bin Structure
        write(UNIT_NUC) this_ndpp % energy_bins
        
        ! Scattering Type (Legendre/Hist), Order of this Type, Chi Present, 
        ! Thinning Tolerance
        ! First convert the logical value of Chi Present to an integer. It seemas as if a type-cast
        ! is not in the standard, so the next if-then  block will explicitly do the cast.
        if(this_ndpp % integrate_chi .AND. nuc % fissionable) then 
          chi_present_int = 1
        else
          chi_present_int = 0
        end if
        ! Now print the results
        write(UNIT_NUC) this_ndpp % scatt_type
        write(UNIT_NUC) this_ndpp % scatt_order
        write(UNIT_NUC) chi_present_int
        write(UNIT_NUC) this_ndpp % thin_tol
        
        ! Write mu_bins
        write(UNIT_NUC) this_ndpp % mu_bins
          
      else if (this_ndpp % lib_format == HDF5) then
        !!! TBI
      end if
      
    end subroutine init_library
  
  
  end module ndpp_class
  
