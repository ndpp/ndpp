program ndpp

  use ndpp_class,   only: nuclearDataPreProc
  use global,       only: free_memory
  use initialize,   only: init_run
  use output

#ifdef MPI
  use mpi
#endif

  implicit none

  type(nuclearDataPreProc) :: my_ndpp

  call init_run()

  call my_ndpp % init()

  call my_ndpp % preprocess()

  call my_ndpp % print_runtime()

  call my_ndpp % clear()

  call free_memory()

#ifdef MPI
  call MPI_FINALIZE(mpi_err)
#endif

end program ndpp
