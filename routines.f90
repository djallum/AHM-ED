module routines

implicit none 

real :: H00
real, dimension(3,3) :: H10, H01 

contains 

subroutine random_gen_seed()

integer :: seed_size, clock, i
integer, allocatable, dimension(:) :: seed

call random_seed(size=seed_size)
allocate(seed(seed_size))
call system_clock(count = clock)
seed=clock + 37*(/(i-1,i=1,seed_size)/)

call random_seed(put=seed)

deallocate(seed)

end subroutine random_gen_seed

subroutine site_potentials(delta,E)

real, intent(in) :: delta
real, intent(out) :: E(2)
real :: random
integer :: i ! counter

do i=1,3
   call random_number(random)
   E(i) = random - 0.5          !centering the random numbers about 0
end do

end subroutine site_potentials

subroutine make_hamiltonians(E,t)

real, intent(in) :: E(3)
real, intent(in) :: t
integer :: i ! counter

H00 = 0

H10 = t

do i=1,3
   H10(i,i) = E(i)
end do

H01 = H10 

end subroutine make_hamiltonians

end module routines
