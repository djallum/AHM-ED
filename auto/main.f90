program main

	use routines

	implicit none

	real(dp), parameter :: t = 1
	real(dp) :: E(nsites)
	real(dp), parameter :: U=2
  	real(dp), parameter :: mu = U/2
  	real(dp), parameter :: delta=12

	integer :: i,j

	call num_sites()

	!do i=1,4**nsites
	!	write(*,*) fock_states(1,i), fock_states(2,i)
	!end do

	call random_gen_seed()
	call transformations()
	call make_neighbours()
	call make_hamiltonian2(t)

	!do i=0,nsites
	!	do j=0,nsites
	!		write(*,*) i,j,msize(i,j),mblock(i,j)	 
	!	end do
	!end do

	!write(*,*) SUM(msize)
	!write(*,*) H00
	!do i=1,2
	!	write(*,*) H10(i,:)
	!end do

	call site_potentials(delta,E)
	E(1)=1; E(2)=-2
	call solve_hamiltonian2(E,U,mu)
	
end program
