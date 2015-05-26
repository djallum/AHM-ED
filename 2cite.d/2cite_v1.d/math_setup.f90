module math_setup  
  implicit none

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
