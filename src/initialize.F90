module initialize

  use constants
  use error,            only: fatal_error
  use global,           only: message, path_input, master
  use string,           only: starts_with, ends_with
  use output,           only: title, header, print_usage, print_version

  implicit none

contains

!===============================================================================
! INITIALIZE_RUN takes care of all initialization tasks, i.e. reading
! from command line, reading xml input files, etc.
!===============================================================================

  subroutine init_run()


    ! Read command line arguments
    call read_command_line()

    if (master) then
      ! Display title and initialization header
       call title()
       call header("INITIALIZATION", level=1)
    end if

  end subroutine init_run

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line()

    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    integer :: last_flag ! index of last flag
    character(MAX_FILE_LEN) :: pwd      ! present working directory
    character(MAX_WORD_LEN), allocatable :: argv(:) ! command line arguments

    ! Get working directory
    call GET_ENVIRONMENT_VARIABLE("PWD", pwd)

    ! Check number of command line arguments and allocate argv
    argc = COMMAND_ARGUMENT_COUNT()

    ! Allocate and retrieve command arguments
    allocate(argv(argc))
    do i = 1, argc
      call GET_COMMAND_ARGUMENT(i, argv(i))
    end do

    ! Process command arguments
    last_flag = 0
    i = 1
    do while (i <= argc)
      ! Check for flags
      if (starts_with(argv(i), "-")) then
        select case (argv(i))
        case ('-?', '-help', '--help')
          call print_usage()
          stop
        case ('-v', '-version', '--version')
          call print_version()
          stop
          message = "Unknown command line option: " // argv(i)
          call fatal_error()
        end select

        last_flag = i
      end if

      ! Increment counter
      i = i + 1
    end do

    ! Determine directory where XML input files are
    if (argc > 0 .and. last_flag < argc) then
      path_input = argv(last_flag + 1)
      ! Need to add working directory if the given path is a relative path
      if (.not. starts_with(path_input, "/")) then
        path_input = trim(pwd) // "/" // trim(path_input)
      end if
    else
      path_input = pwd
    end if

    ! Add slash at end of directory if it isn't there
    if (.not. ends_with(path_input, "/")) then
      path_input = trim(path_input) // "/"
    end if

    ! Free memory from argv
    deallocate(argv)

    ! TODO: Check that directory exists

  end subroutine read_command_line

end module initialize
