! note of unit numbers
! 7: input file
! 8: log file
! 9: data file
! 10: error
! 11: warning file: warning.txt

program main

  use mt19937
  use Ising_MCMC
  use Jackknife

  implicit none

  ! model setting
  integer :: L ! system size
  integer :: N_T ! number of temperatures
  double precision, allocatable, dimension(:) :: T ! temperature
  character(128) init_cond ! initial configuration. "random" or "ferromag"
  character(128) algorithm ! update algorithm (You can choose only "Metropolis" currently)
  integer thermalization_step ! number of MC steps for thermalization(discarded)
  integer MC_step ! number of MC steps
  double precision beta ! inverse temperature

  ! MCMC variables
  ! double precision E ! energy
  ! double precision M ! total magnetization
  ! double precision M2 ! squared total magnetization

  ! data analysis parameters
  character(128) data_analysis_method
  integer :: blocksize

  ! array to store data
  double precision,allocatable,dimension(:) :: E_arr
  double precision,allocatable,dimension(:) :: M_arr
  double precision,allocatable,dimension(:) :: M2_arr

  ! physical observables
  double precision :: E_ave, E_err ! average total energy and its error
  double precision :: M_ave, M_err ! average total magnetization and its error
  double precision :: M2_ave, M2_err ! average total squared magnetization and its error
  double precision :: Cv, Cv_err ! heat capacity and its error
  double precision :: chi, chi_err ! magnetic susceptibility and its error

  ! other variables
  integer seed_MT ! seed of the Mersenne twister

  integer i, j, k, ret_val
  

  ! open log file
  open(8, file = "result/Ising_MCMC.log", status = "replace")
  
  ! read input file
  call read_input(L, seed_MT, init_cond, algorithm, thermalization_step, MC_step, data_analysis_method, &
  blocksize, N_T, T, write_flg = .true.) ! write parameters to the log file if write_flg = .true.
  
  ! set seed of the Mersenne twister
  call sgrnd(seed_MT)

  write(8, *) "calculation: "
  ! allocate E_arr, M_arr
  allocate(E_arr(MC_step))
  allocate(M_arr(MC_step))
  allocate(M2_arr(MC_step))

  ! run MCMC simulation
  open(9, file = "result/T_dependence.txt", status = "replace")
  write(9, *) "#1: T, 2: E_ave/(L*L), 3: E_err/(L*L), 4: M_ave/(L*L), &
  5, M_err/(L*L), 6: M2_ave/(L^4), 7: M2_err/(L^4), 8: Cv/(L*L), 9: Cv_err/(L*L), &
  10: chi/(L*L), 11: chi_err/(L*L)"
  write(9, *) "# all the printed physical quantities are &
  normalized as the quantity per site."

  do i = 1, N_T
    print *, "T = ", T(i)

    ! run MCMC simulation
    ret_val = run_MCMC(L, T(i), init_cond, algorithm, thermalization_step, MC_step, E_arr, M_arr, M2_arr)

    ! data analysis
    call jackknife_ave(MC_step, E_arr, blocksize, E_ave, E_err) ! calculate average energy and its error
    call jackknife_ave(MC_step, M_arr, blocksize, M_ave, M_err) ! calculate average magnetization and its error
    call jackknife_ave(MC_step, M2_arr, blocksize, M2_ave, M2_err) ! calculate average squared magnetization and its error

    call jackknife_stddv(MC_step, E_arr, blocksize, Cv, Cv_err) ! calculate heat capacity and its error
    Cv = Cv/(T(i)**2)
    Cv_err = Cv_err/(T(i)**2)

    call jackknife_stddv(MC_step, M_arr, blocksize, chi, chi_err) ! calculate magnetic susceptibility and its error
    chi = chi/T(i)
    chi_err = chi_err/T(i)

    ! write result
    write(9, *) T(i), E_ave/(L*L), E_err/(L*L), M_ave/(L*L), M_err/(L*L), &
     M2_ave/(L*L*L*L), M2_err/(L*L*L*L), Cv/(L*L), Cv_err/(L*L), chi/(L*L), chi_err/(L*L)
    write(8, *) T(i), E_ave/(L*L), E_err/(L*L), M_ave/(L*L), M_err/(L*L), &
     M2_ave/(L*L*L*L), M2_err/(L*L*L*L), Cv/(L*L), Cv_err/(L*L), chi/(L*L), chi_err/(L*L)
  end do 

  close(9)

  ! close log file
  close(8)
  
  ! deallocate arrays
  deallocate(E_arr)
  deallocate(M_arr)
  deallocate(M2_arr)


contains
  !=========================================================
  ! read simulation setting from the input file "T_dependence.inp"
  ! write the 
  !=========================================================

  subroutine read_input(L, seed_MT, init_cond, algorithm, thermalization_step, MC_step, data_analysis, &
    blocksize, N_T, T, write_flg)

    integer,intent(out) :: L
    integer, intent(out) :: seed_MT
    character(128), intent(out) :: init_cond
    character(128), intent(out) :: algorithm
    integer, intent(out) :: thermalization_step
    integer, intent(out) :: MC_step

    character(128), intent(out) :: data_analysis
    integer, intent(out) :: blocksize
    integer, intent(out) :: N_T
    double precision, allocatable, dimension(:), intent(out) :: T

    logical, intent(in) :: write_flg
    character(128) tmp_string

    ! read input file
    open(7, file = "T_dependence.inp", status = "old")

    read(7, *) L
    !read(7, *) T
    read(7, *) seed_MT
    read(7, *) init_cond, tmp_string
    read(7, *) algorithm

    read(7, *) thermalization_step
    read(7, *) MC_step

    read(7, *) data_analysis
    read(7, *) blocksize
    read(7, *) N_T
    allocate(T(N_T))
    read(7, *) T

    close(7)

    ! write log file
    if(write_flg) then
      write(8, *) "simulation setting: "
      write(8, *) "L: ", L
      write(8, *) "initial condition: ", init_cond
      write(8, *) "algorithm: ", algorithm
      write(8, *) "thermalization_step: ", thermalization_step
      write(8, *) "MC_Step: ", MC_step
      write(8, *) "data analysis method: ", data_analysis_method
      write(8, *) "Jackknife blocksizes: ", blocksize
      write(8, *) "number of temperatures: ", N_T
      write(8, *) "temperatures: ", T
      write(8, *) "seed of the Mersenne twister: ", seed_MT
    end if

  end subroutine read_input

end program main