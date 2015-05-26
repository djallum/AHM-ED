module routines

implicit none 

real :: H00, H30, H03
real, dimension(3,3) :: H10, H01, H20, H02 
real, dimension(9,9) :: H11

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

subroutine make_hamiltonians(E,t,U)

real, intent(in) :: E(3)
real, intent(in) :: t
real, intent(in) :: U
integer :: i,j ! counter

H00 = 0

H10 = t

do i=1,3
   H10(i,i) = E(i)
end do

H01 = H10 

H11(1,2) = t;  H11(1,3) = -t; H11(1,4) = t;  H11(1,5) = -t
H11(2,4) = t;  H11(2,5) = t;  H11(2,8) = -t
H11(3,4) = -t; H11(3,6) = t;  H11(3,7) = -t
H11(4,7) = t;  H11(4,8) = -t
H11(5,7) = t;  H11(5,8) = t
H11(6,8) = t;  H11(6,9) = t 
H11(7,9) = t
H11(8,9) = -t

do i=1,9
   do j=1,9
      if(i>j) then
         H11(i,j) = H11(j,i)
      end if
   end do
end do

H11(1,1) = 2*E(1) + U
H11(2,2) = E(1) + E(2); H11(3,3) = E(1) + E(2)
H11(4,4) = 2*E(2) + U
H11(5,5) = E(1) + E(3); H11(6,6) = E(1) + E(3)
H11(7,7) = E(2) + E(3); H11(8,8) = E(2) + E(3)
H11(9,9) = 2*E(3) + U

end subroutine make_hamiltonians

end module routines
