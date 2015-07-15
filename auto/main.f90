program main

	use routines

	implicit none

	integer :: i
	call num_sites()

	!do i=1,4**nsites
	!	write(*,*) fock_states(1,i), fock_states(2,i)
	!end do

	call transformations()

	write(*,*) PES_up(1,:)
write(*,*) phase_PES_up(1,:)

end program
