module math_setup  
  implicit none
  integer, parameter :: int_kind = 4
  integer, parameter :: real_kind = kind(1.d0)
  integer, parameter :: STDOUT = 10
  real (real_kind), parameter :: PI_ = 4d0*atan(1d0)

contains

  integer function choose(j,k)
    implicit none
    integer :: j,k,i2
    integer (kind=8) :: tj,tk,tchoose
    tchoose = 1
    tj = j
    tk = k
    do i2 = tj-tk+1,tj
       tchoose = tchoose * i2
    end do
    do i2 = 2, tk
       tchoose = tchoose / i2
    end do
    choose = tchoose
  end function choose

end module math_setup
