program test

use random_gen

integer :: i
real, dimension(8) :: hi

call gen_seed()
call random_number(hi)

do i=1,8,1
   hi(i) = hi(i) - 0.5
end do

write(*,*), hi

end program test
