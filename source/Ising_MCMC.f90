! note of unit numbers
! 7: input file
! 8: log file
! 9: data file
! 10: error
! 11: warning file: warning.txt

module Ising_MCMC

  use mt19937


  implicit none
  

  contains

  !=========================================================
  ! run MCMC simulation for fixed simulation setting.
  ! store the simulation result in E_arr, M_arr, M2_arr
  ! 
  ! L : perform simulation on square lattice of L*L sites
  ! T : temperature
  ! algorithm : Metropolis
  ! thermalization_step : number of discarded MCMC steps for thermalization
  ! MC_step : number of MCMC steps
  !
  ! E_arr : array to store the total energy at each MC step
  ! M_arr : array to store the total magnetization at each MC step
  ! M2_arr : array to store the squared total magnetization at each MC step
  !=========================================================

  integer function run_MCMC(L, T, init_cond, algorithm, thermalization_step, MC_step, E_arr, M_arr, M2_arr)

    ! input variables
    integer, intent(in) :: L
    double precision, intent(in) :: T
    character(128), intent(in) :: init_cond
    character(128), intent(in) :: algorithm
    integer, intent(in) :: thermalization_step
    integer, intent(in) :: MC_step

    ! output variables
    double precision, intent(out) :: E_arr(:)
    double precision, intent(out) :: M_arr(:)
    double precision, intent(out) :: M2_arr(:)

    ! internal variables
    double precision beta ! inverse temperature
    integer, allocatable, dimension(:) :: S ! spin configuration
    integer,allocatable,dimension(:,:) :: neighbor_list ! neighbor list
    double precision E, M

    ! other variables
    integer i

    ! initial setting 
    beta = 1.0/T
    allocate(S(L*L))
    allocate(neighbor_list(L*L,4))

    call set_neighbor_list(L, neighbor_list, .false.)

    call init_Ising(L, init_cond, neighbor_list, E, M, S, .false.)

    do i = 1, thermalization_step
      call update_state_Metropolis(L, neighbor_list, beta, S, E, M)
    end do
  
    ! MCMC
    do i = 1, MC_step
      call update_state_Metropolis(L, neighbor_list, beta, S, E, M)
      E_arr(i) = E
      M_arr(i) = M
      M2_arr(i) = M*M
    end do

    deallocate(S)
    deallocate(neighbor_list)

    run_MCMC = 0
    return

  end function run_MCMC

  !=========================================================
  ! set neighbor_list
  !=========================================================

  subroutine set_neighbor_list(L, neighbor_list, flg_print)

    integer, intent(in) :: L
    integer, intent(out) :: neighbor_list(:,:)
    logical, intent(in) :: flg_print

    ! internal variables
    integer i

    ! set neighbor_list
    do i = 1, L*L
      ! set left
      if(mod(i,L) == 1) then
        neighbor_list(i,1) = i+L-1
      else
        neighbor_list(i,1) = i-1
      end if
      ! set right
      if(mod(i,L) == 0) then
        neighbor_list(i,2) = i-L+1
      else
        neighbor_list(i,2) = i+1
      end if
      ! set top
      if(i <= L) then
        neighbor_list(i,3) = i + L*(L-1)
      else
        neighbor_list(i,3) = i-L
      end if
      ! set bottom
      if(i > L*(L-1)) then
        neighbor_list(i,4) = i - L*(L-1)
      else
        neighbor_list(i,4) = i+L
      end if 
    end do

    if(flg_print) then
      print *, "neighbor_list: "
      do i = 1, L*L
        print *, i, ":", neighbor_list(i,:)
      end do
    end if

  end subroutine set_neighbor_list

  !=========================================================
  ! initialize Ising model
  !=========================================================

  subroutine init_Ising(L, init_cond, neighbor_list, E, M, S, flg_print)

    ! input arguments
    integer,intent(in) :: L
    character(128), intent(in) :: init_cond
    integer,intent(in) :: neighbor_list(:,:)
    logical,intent(in) :: flg_print
    ! output arguments
    double precision, intent(out) :: E
    double precision, intent(out) :: M
    integer,intent(out) :: S(:)
    

    ! other variables
    integer i,j
    double precision rndnum

    ! set random spin configuration to S
    !if(init_cond(:) == "ferromag") then
    if(index(init_cond, "ferromag") == 1) then
      ! set ferromagnetic initial condition
      do i = 1, L*L
        S(i) = 1
      end do
    else 
      ! warning if init_cond is not "random"
      !if(init_cond /= "random")then
      if(index(init_cond, "random") /= 1) then
        print *, "Warning in subroutine Init_Ising: input initial condition is invalid." // &
        "The initial condition is set randomly by default."

        open(11, file = "result/warning.txt", position = "append")
        write(11, *) "Warning in subroutine Init_Ising: input initial condition is invalid." // &
        "The initial condition is set randomly by default."
        close(11)
      end if

      ! set random initial condition
      do i = 1, L*L
        rndnum = grnd()
        if(rndnum < 0.5) then
          S(i) = -1
        else 
          S(i) = 1
        end if
      end do
    end if

    ! calculate energy and magnetization
    E = 0.0
    M = 0.0
    do i = 1, L*L
      M = M + S(i)
      do j = 1, 4
        E = E - S(i) * S(neighbor_list(i,j)) * 0.5
      end do
    end do

    ! print result
    if(flg_print) then

      print *, "The result of subroutine init_Ising."

      print *, "initial spin configuration S: "
      do i = 1, L
      print *, S(i*L-L+1:i*L)
      end do

      print *, "initial energy E: ", E
      print *, "initial total magnetization M: ", M

    end if
  end subroutine init_Ising

  !=========================================================
  ! update spin configuration S by Metropolis update
  !=========================================================

  subroutine update_state_Metropolis(L, neighbor_list, beta, S, E, M)

    ! arguments
    integer,intent(in) :: L
    integer,intent(in) :: neighbor_list(:,:)
    double precision,intent(in) :: beta
    integer,intent(inout) :: S(:)
    double precision, intent(inout) :: E
    double precision, intent(inout) :: M

    ! other variables
    integer i, j
    double precision rndnum
    double precision dE

    ! select site randomly
    !call random_number(rndnum)
    rndnum = grnd()
    i = min(int(rndnum*L*L) + 1, L*L)

    ! local update
    dE = 0.0
    do j = 1,4
      dE = dE + 2.0*S(i)*S(neighbor_list(i,j))
    end do

    !call random_number(rndnum)
    rndnum = grnd()
    if(rndnum < exp(- beta*dE)) then
      S(i) = -S(i)
      E = E + dE
      M = M + 2.0*S(i)
    end if

  end subroutine update_state_Metropolis

end module Ising_MCMC