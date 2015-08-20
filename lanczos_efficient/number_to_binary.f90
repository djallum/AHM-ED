program test

  implicit none

  integer :: a
  character (16) :: ac,convert_to_binary


  a = 1051
  ac = convert_to_binary(a)
  print *, a,ac

end program test

character(16) function convert_to_binary(a)
  
  implicit none
  
  integer :: a
  character (16) :: ac
  integer :: i

  ac=''
  do i=0,15
     if (btest(a,i)) then
        ac = '1'//trim(ac)
     else
        ac = '0'//trim(ac)
     end if
  end do
     
  convert_to_binary = ac
  
end function convert_to_binary
  
