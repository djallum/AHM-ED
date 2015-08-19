program math

  ! This program contains a subroutine that calculates the choose function of two variables. This is the typical
  ! choose function from statistics. 

  ! %---------------------------------------------%
  ! |  simple choose function from statistics:    |
  ! |                                             |
  ! |                         i!                  |
  ! |       choose(i,j) =  --------               |
  ! |                      j!(i-j)!               |
  ! %---------------------------------------------%
  
	implicit none
	    
    integer :: j,k,i2
    integer (kind=8) :: tj,tk,tchoose

   	read(*,*), tj,tk
    tchoose = 1
    do i2 = tj-tk+1,tj
       tchoose = tchoose * i2
    end do
    do i2 = 2, tk
       tchoose = tchoose / i2
    end do

    write(*,*) tchoose


end program math