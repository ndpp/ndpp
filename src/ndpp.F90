program ndpp

  use ndpp_class,   only: nuclearDataPreProc
  use global,       only: free_memory
  use ndpp_init,    only: init_run
  use output
  
  implicit none
  
  type(nuclearDataPreProc) :: my_ndpp
  
  call init_run()
  
  call my_ndpp % init()
  
  call my_ndpp % preprocess()
  
  call my_ndpp % print_runtime()
  
  call my_ndpp % clear()
  
  call free_memory()
  
end program ndpp
