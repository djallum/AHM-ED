module parameters

!type :: parameters_list

!real :: t      ! hopping potential
!real :: U      ! on site repulsion
!real :: mu     ! chemical potential
!real :: delta  ! widtth of disorder potential
!real :: npair  ! number of pairs
!end type parameters_list

!type(parameters_list) :: parameters

contains

subroutine read_parameters()

integer :: status
!parameters_list :: parameters

open(UNIT=10, FILE="input.dat", STATUS="OLD", ACTION="read", IOSTAT=status) 
if (status /= 0) then
   write(*,*) 'Error opening input file'
   stop
end if

end subroutine read_parameters

!subroutine evaluate_parameters()

!end subroutine evaluate_parameters

subroutine site_potentials(npairs,nsites,E,delta)              ! randomly creates site potentials and stores them in the array E

real :: random1, random2, temp
real, intent(out) :: E(npairs,nsites)
real, intent(in) :: delta

do i=1,npairs
   call random_number(random1)
   call random_number(random2)
   random1 = random1 - 0.5
   random2 = random2 - 0.5
   E(i,1) = delta*max(random1,random2)
   E(i,2) = delta*min(random1,random2)
end do
   
end subroutine site_potentials

end module parameters

