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
  double precision :: T ! temperature
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
  integer :: N_blocksizes
  integer,allocatable, dimension(:) :: blocksizes

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
  call read_input(L, T, seed_MT, init_cond, algorithm, thermalization_step, MC_step, data_analysis_method, &
  N_blocksizes, blocksizes, write_flg = .true.) ! write parameters to the log file if write_flg = .true.
  
  ! set seed of the Mersenne twister
  call sgrnd(seed_MT)

  write(8, *) "calculation: "
  ! allocate E_arr, M_arr
  allocate(E_arr(MC_step))
  allocate(M_arr(MC_step))
  allocate(M2_arr(MC_step))

  ! run MCMC simulation
  ret_val = run_MCMC(L, T, init_cond, algorithm, thermalization_step, MC_step, E_arr, M_arr, M2_arr)

  ! write MCMC result
  open(9, file = "result/MCMC_sequence.txt", status = "replace")
  write(9, *) "#1: MC_step, 2: Energy E, 3: total magnetization M, 4: squared total magnetization M2"
  write(9, *) "# Note that the quantities are not normalized as the quantities per site."
  do i = 1, MC_step
    write(9,*) i, E_arr(i), M_arr(i), M2_arr(i)
  end do
  close(9)

  ! data analysis
  open(9, file = "result/jackknife_result.txt")
  write(9, *) "#1: block size, 2: E_ave/(L^2), 3: E_err/(L^2), 4: M_ave/(L^2), 5, M_err/(L^2), &
  6: M2_ave/(L^4), 7: M2_err/(L^4), 8: Cv/(L^2), 9: Cv_err/(L^2), 10: chi/(L^2), 11: chi_err/(L^2)"
  write(9, *) "the printed quantities are normalized as the quantity per site."
  do i = 1, N_blocksizes
    call jackknife_ave(MC_step, E_arr, blocksizes(i), E_ave, E_err) ! calculate average energy and its error
    call jackknife_ave(MC_step, M_arr, blocksizes(i), M_ave, M_err) ! calculate average magnetization and its error
    call jackknife_ave(MC_step, M2_arr, blocksizes(i), M2_ave, M2_err) ! calculate average squared magnetization and its error

    call jackknife_stddv(MC_step, E_arr, blocksizes(i), Cv, Cv_err) ! calculate heat capacity and its error
    Cv = Cv/(T**2)
    Cv_err = Cv_err/(T**2)

    call jackknife_stddv(MC_step, M_arr, blocksizes(i), chi, chi_err) ! calculate magnetic susceptibility and its error
    chi = chi/T
    chi_err = chi_err/T

    write(9, *) blocksizes(i), E_ave/(L*L), E_err/(L*L), M_ave/(L*L), &
    M_err/(L*L), M2_ave/(L**4), M2_err/(L**4), Cv/(L*L), Cv_err/(L*L), chi/(L*L), chi_err/(L*L)
  end do
  close(9)

  ! write final
  write(8, *) "#1: L, 2: T, 3: E_ave, 4: E_err, 5: M_ave, 6, M_err, 7: M2_ave, 8: M2_err, 9: Cv, 10: Cv_err, 11: chi, 12: chi_err"
  write(8, *) L, T, E_ave, E_err, M_ave, M_err, M2_ave, M2_err, Cv, Cv_err, chi, chi_err

  ! print result
  write(8, *) "just for check: "
  write(8, *) "average energy: ", sum(E_arr)/MC_step
  write(8, *) "average magnetization: ", sum(M_arr)/MC_step
  write(8, *) "average squared magnetization: ", sum(M2_arr)/MC_step

  ! close log file
  close(8)
  
  ! deallocate arrays
  deallocate(E_arr)
  deallocate(M_arr)
  deallocate(M2_arr)


contains
  !=========================================================
  ! read simulation setting from the input file "single_T.inp"
  ! write the 
  !=========================================================

  subroutine read_input(L, T, seed_MT, init_cond, algorithm, thermalization_step, MC_step, data_analysis, &
    N_blocksizes, blocksizes, write_flg)

    integer,intent(out) :: L
    double precision, intent(out) :: T
    integer, intent(out) :: seed_MT
    character(128), intent(out) :: init_cond
    character(128), intent(out) :: algorithm
    integer, intent(out) :: thermalization_step
    integer, intent(out) :: MC_step

    character(128), intent(out) :: data_analysis
    integer, intent(out) :: N_blocksizes
    integer, allocatable, dimension(:), intent(out) :: blocksizes

    logical, intent(in) :: write_flg
    character(128) tmp_string

    ! read input file
    open(7, file = "single_T.inp", status = "old")

    read(7, *) L
    read(7, *) T
    read(7, *) seed_MT
    read(7, *) init_cond, tmp_string
    read(7, *) algorithm

    read(7, *) thermalization_step
    read(7, *) MC_step

    read(7, *) data_analysis
    read(7, *) N_blocksizes
    allocate(blocksizes(N_blocksizes))
    read(7, *) blocksizes

    close(7)

    ! write log file
    if(write_flg) then
      write(8, *) "simulation setting: "
      write(8, *) "L: ", L
      write(8, *) "T: ", T
      write(8, *) "initial condition: ", init_cond
      write(8, *) "algorithm: ", algorithm
      write(8, *) "thermalization_step: ", thermalization_step
      write(8, *) "MC_Step: ", MC_step
      write(8, *) "data analysis method: ", data_analysis_method
      write(8, *) "number of Jackknife blocksizes: ", N_blocksizes
      write(8, *) "Jackknife blocksizes: ", blocksizes
      write(8, *) "seed of the Mersenne twister: ", seed_MT
    end if

  end subroutine read_input

end program main