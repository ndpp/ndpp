module output

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header,      only: Nuclide, Reaction, UrrData
  use constants
  use dict_header
  use endf,            only: reaction_name
  use global,          only: master, message, verbosity
  use string,          only: upper_case, to_str

  implicit none

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title()

    write(UNIT=OUTPUT_UNIT, FMT='(/11(A/))') &
         '                 888b    888 8888888b.  8888888b.  8888888b.  ', &
         '                 8888b   888 888  "Y88b 888   Y88b 888   Y88b ', &
         '                 88888b  888 888    888 888    888 888    888 ', &
         '                 888Y88b 888 888    888 888   d88P 888   d88P ', &
         '                 888 Y88b888 888    888 8888888P"  8888888P"  ', &
         '                 888  Y88888 888    888 888        888        ', &
         '                 888   Y8888 888  .d88P 888        888        ', &
         '                 888    Y888 8888888P"  888        888        ', &
         '____________________________________________________________________________'

    ! Write version information
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '     Developed At:  University of Michigan and'
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '                    Massachusetts Institute of Technology'
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Version:",7X,I1,".",I1,".",I1)') &
         VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Git SHA1:",6X,A)') GIT_SHA1
#endif

    ! Write the date and time
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Date/Time:",5X,A)') &
         time_stamp()

#ifdef MPI
    ! Write number of processors
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"MPI Processes:",1X,A)') &
         trim(to_str(n_procs))
#endif

  end subroutine title

!===============================================================================
! TIME_STAMP returns the current date and time in a formatted string
!===============================================================================

  function time_stamp() result(current_time)

    character(19) :: current_time ! ccyy-mm-dd hh:mm:ss
    character(8)  :: date_        ! ccyymmdd
    character(10) :: time_        ! hhmmss.sss

    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end function time_stamp

!===============================================================================
! HEADER displays a header block according to a specified level. If no level is
! specified, it is assumed to be a minor header block (H3).
!===============================================================================

  subroutine header(msg, unit, level)

    character(*), intent(in) :: msg ! header message
    integer, optional :: unit       ! unit to write to
    integer, optional :: level      ! specified header level

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    integer :: header_level ! actual header level
    character(MAX_LINE_LEN) :: line

    ! set default level
    if (present(level)) then
      header_level = level
    else
      header_level = 3
    end if

    ! set default unit
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! determine how many times to repeat '=' character
    n = (63 - len_trim(msg))/2
    m = n
    if (mod(len_trim(msg),2) == 0) m = m + 1

    ! convert line to upper case
    line = msg
    call upper_case(line)

    ! print header based on level
    select case (header_level)
    case (1)
      write(UNIT=unit_, FMT='(/3(1X,A/))') repeat('=', 75), & 
           repeat('=', n) // '>     ' // trim(line) // '     <' // &
           repeat('=', m), repeat('=', 75)
    case (2)
      write(UNIT=unit_, FMT='(/2(1X,A/))') trim(line), repeat('-', 75)
    case (3)
      write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
           trim(line) // '     <' // repeat('=', m)
    end select

  end subroutine header

!===============================================================================
! PRINT_VERSION shows the current version as well as copright and license
! information
!===============================================================================

  subroutine print_version()

    if (master) then
      write(UNIT=OUTPUT_UNIT, FMT='(1X,A,1X,I1,".",I1,".",I1)') &
           "NDPP version", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
      write(UNIT=OUTPUT_UNIT, FMT=*) "Copyright (c) 2011-2012 &
           &Massachusetts Institute of Technology"
      write(UNIT=OUTPUT_UNIT, FMT=*) "MIT/X license at &
           &<http://mit-crpg.github.com/openmc/license.html>"
    end if

  end subroutine print_version

!===============================================================================
! PRINT_USAGE displays information about command line usage of NDPP
!===============================================================================

  subroutine print_usage()

    if (master) then
      write(OUTPUT_UNIT,*) 'Usage: ndpp [options] [directory]'
      write(OUTPUT_UNIT,*)
      write(OUTPUT_UNIT,*) 'Options:'
      write(OUTPUT_UNIT,*) '  -v, --version   Show version information'
      write(OUTPUT_UNIT,*) '  -?, --help      Show this message'
    end if

  end subroutine print_usage

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the 
! standard output stream.
!===============================================================================

  subroutine write_message(level)

    integer, optional :: level ! verbosity level

    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message

    ! Set length of line
    line_wrap = 80

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    if (.not. present(level) .or. level <= verbosity) then
      ! Determine length of message
      length = len_trim(message)

      i_start = 0
      do
        if (length - i_start < line_wrap - 1) then
          ! Remainder of message will fit on line
          write(ou, fmt='(1X,A)') message(i_start+1:length)
          exit

        else
          ! Determine last space in current line
          i_end = i_start + index(message(i_start+1:i_start+line_wrap), &
               ' ', BACK=.true.)

          ! Write up to last space
          write(ou, fmt='(1X,A)') message(i_start+1:i_end-1)

          ! Advance starting position
          i_start = i_end
          if (i_start > length) exit
        end if
      end do
    end if

  end subroutine write_message

  subroutine print_ascii_array(array, ou)
    real(8), intent(in) :: array(:)
    integer, intent(in) :: ou
    
    integer                 :: i_array
    character(MAX_LINE_LEN) :: line
    
    i_array = 1
    do
      line = ''
      if ((size(array) - i_array) >= 3) then
        write(line, '(1PE20.12,1PE20.12,1PE20.12,1PE20.12)') &
          array(i_array), array(i_array + 1), array(i_array + 2), &
          array(i_array + 3)
        i_array = i_array + 4
        write(ou, '(A)') trim(line)
      else
        select case(size(array) - i_array)
        case (0)
          write(line, '(1PE20.12)') array(i_array)
          write(ou, '(A)') trim(line)
        case (1)
          write(line, '(1PE20.12,1PE20.12)') array(i_array), array(i_array + 1)
          write(ou, '(A)') trim(line)
        case (2)
          write(line, '(1PE20.12,1PE20.12,1PE20.12)') array(i_array), &
            array(i_array + 1), array(i_array + 2)
          write(ou, '(A)') trim(line)
        end select
        exit
      end if
    end do
  end subroutine print_ascii_array
  
end module output
