! note of unit numbers
! 7: input file
! 8: log file
! 9: data file
! 10: error
! 11: warning file: warning.txt

module Jackknife
  implicit none

  contains
    !=========================================================
  ! calculate average and its error of array of length N
  ! using Jackknife analysis with blocksize
  !=========================================================

  subroutine jackknife_ave(N, array, blocksize, ave, err)

    ! arguments
    integer, intent(in) :: N
    double precision, dimension(:), intent(in) :: array
    integer, intent(in) :: blocksize
    double precision, intent(out) :: ave
    double precision, intent(out) :: err

    ! interior variables
    integer M
    double precision, allocatable, dimension(:) :: array_m
    double precision err2
    integer i

    ! warning
    if(mod(N, blocksize) /= 0) then
      open(11, file = "result/warning.txt", position = "append")
      write(11, *) "Warning in subroutine jackknife_ave : N cannot be devided by blocksize."
      close(11)
    end if

    ! jackknife analysis
    M = N/blocksize
    allocate(array_m(M))
    do i = 1, M
      array_m(i) = sum(array((i-1)*blocksize+1: i*blocksize))/blocksize
    end do

    ave = sum(array_m)/M
    
    err2 = 0.0
    do i = 1, M
      err2 = err2 + 1.0/(M*(M-1)) * (array_m(i) - ave)**2
    end do
    err = sqrt(err2)

    deallocate(array_m)

  end subroutine jackknife_ave

  !=========================================================
  ! calculate standard deviation and its error of array of length N
  ! using Jackknife analysis with blocksize
  !=========================================================

  subroutine jackknife_stddv(N, array, blocksize, stddv, err)

    ! arguments
    integer, intent(in) :: N
    double precision, dimension(:), intent(in) :: array
    integer, intent(in) :: blocksize
    double precision, intent(out) :: stddv
    double precision, intent(out) :: err

    ! interior variables
    integer M
    double precision, allocatable, dimension(:) :: array_m, array2_m
    double precision array_ave, array2_ave
    double precision err2
    integer i

    ! warning
    if(mod(N, blocksize) /= 0) then
      open(11, file = "result/warning.txt", position = "append")
      write(11, *) "Warning in subroutine jackknife_stddv : N cannot be devided by blocksize."
      close(11)
    end if

    ! jackknife analysis
    M = N/blocksize
    allocate(array_m(M))
    allocate(array2_m(M))
    do i = 1, M
      array_m(i) = sum(array((i-1)*blocksize+1: i*blocksize))/blocksize
      array2_m(i) = sum( (array((i-1)*blocksize+1: i*blocksize))**2 )/blocksize
    end do

    array_ave = sum(array_m)/M
    array2_ave = sum(array2_m)/M
    stddv = array2_ave - array_ave**2

    err2 = 0.0
    do i = 1, M
      err2 = err2 + 1.0/(M*(M-1.0)) * (array2_ave - array2_m(i) - (2.0*M-1.0)/(M-1.0)*array_ave**2 &
      + (2.0*M)/(M-1.0)*array_ave*array_m(i) - 1.0/(M-1.0)*array_m(i)**2)**2
    end do

    err = sqrt(err2)

  end subroutine jackknife_stddv



end module Jackknife
