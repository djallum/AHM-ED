module routines

	implicit none

  	!--------------Other Parameters and Variables-------------------------------
  	integer, parameter :: nsites = 2             ! must be set to two
	integer, parameter :: nstates = 4**nsites    ! total number of states
	real, dimension(nstates) :: grand_potential           ! grand potentials (eigenenergies - mu*number electrons)
	real :: grand_potential_ground                        ! the lowest grand potential energy
	real, dimension(nstates,nstates) :: eigenvectors      ! the eigenvectors found from LAPACK
	real, dimension(nstates) :: v_ground                  ! the ground state eigenvector
	integer, dimension(nsites,nstates) :: PES_dn, PES_up                ! lookup tables for PES of ground states vector transformation
	integer, dimension(nsites,nstates) :: IPES_dn, IPES_up              ! lookup tables for IPES of ground states vector transformation
	integer, dimension(nsites,nstates) :: phase_PES_dn, phase_PES_up    ! to get the anticommutation signs right
	integer, dimension(nsites,nstates) :: phase_IPES_dn, phase_IPES_up  ! to get the anticommutation signs right

contains

!************************************************************************************

subroutine random_gen_seed()

	integer :: seed_size         ! the size of the random generators seed
	integer :: clock             ! the system clock (used to seed random generator) 
	integer :: i                 ! counter for loop
	integer, allocatable, dimension(:) :: seed    ! will store our seed

	call random_seed(size = seed_size)       ! find the size of the random generators seed
	allocate(seed(seed_size))                ! make our seed that same size
	call system_clock(count = clock)         ! find the system time
	!seed=clock + 37*(/(i-1,i=1,seed_size)/) ! create our seed using the clock
	seed=3
	call random_seed(put=seed)               ! seed the random generator with our seed

	deallocate(seed)

end subroutine random_gen_seed

!************************************************************************************

subroutine site_potentials(delta,E)

	real, dimension(nsites), intent(out) :: E  ! site potentials
	real, intent(in) :: delta                  ! width of disorder
	real :: random1, random2                   ! random numbers

	call random_number(random1)
	call random_number(random2)
	random1 = random1 - 0.5            ! centering the random numbers about 0
	random2 = random2 - 0.5
	E(1) = delta*max(random1,random2)  ! multiplying the random numbers by the width of disorder
	E(2) = delta*min(random1,random2)

end subroutine site_potentials

!************************************************************************************

subroutine hamiltonian(t,U,mu,E)

	implicit none

	real, intent(in) :: t, U, mu
	real, dimension(nsites), intent(in) :: E             ! site potentials
	real :: H00, H10(2,2), H20, H11(4,4), H21(2,2), H22  ! hamiltonian submatrices
	real :: W00, W10(2), W20, W11(4), W21(2), W22        ! eigenvalues of hamiltonian sub matrices
	integer :: i                                         ! counter for loops

	!------for LAPACK------------
  	integer :: INFO = 0
  	integer :: LWORK
 	real, allocatable, dimension(:) :: WORK

 	!-------------zero variables------------------------------------

 	H00 = 0; H10 = 0; H20 = 0; H11 = 0; H21 = 0; H22 = 0
 	eigenvectors = 0

 	!-------------------------H00------------------------------------------

 	W00 = H00
 	grand_potential(1) = W00
 	eigenvectors(1,1) = 1

 	!-------------------H10 and H01 (same)---------------------------------

 	! enter precalculated hamiltonian
 	H10(1,1) = E(1); H10(1,2) = t
 	H10(2,1) = t;    H10(2,2) = E(2)

 	! solve the eigenvalues and eigenvectors
  	LWORK = 10
  	allocate(WORK(LWORK))

  	call ssyev('v','u',2,H10,2,W10,WORK,LWORK,INFO)
  	if (INFO /= 0) then
     	write(*,*) 'Problem with Lapack for H10 matrix. Error code:', INFO
     	stop
  	end if 

  	do i=1,2
     	grand_potential(i+1) = W10(i) - mu*1  ! grand potentials of H10
     	grand_potential(i+3) = W10(i) - mu*1  ! grand potentials of H01
  	end do

  	do i=1,2
     	eigenvectors(i+1,2:3) = H10(1:2,i)    ! eigenvectors of H10
     	eigenvectors(i+3,4:5) = H10(1:2,i)    ! eigenvectors of H01
  	end do

  	deallocate(WORK)

  	!-------------------H20 and H02 (same)---------------------------------

  	H20 = E(1) + E(2)

  	W20 = H20

  	grand_potential(6) = W20 - 2*mu   ! grand potentials of H20
  	grand_potential(7) = W20 - 2*mu   ! grand potentials of H02

  	eigenvectors(6,6) = 1        ! eigenvectors of H20
  	eigenvectors(7,7) = 1        ! eigenvectors of H02

  	!-------------------------H11------------------------------------------

  	! enter precalculated hamiltonian
 	H11(1,1) = E(1) + E(2); H11(1,3) = t;  H11(1,4) = t
 	H11(2,2) = E(1) + E(2); H11(2,3) = -t; H11(2,4) = -t
 	H11(3,1) = t;           H11(3,2) = -t; H11(3,3) = 2*E(1) + U
 	H11(4,1) = t;           H11(4,2) = -t; H11(4,4) = 2*E(2) + U

 	! solve the eigenvalues and eigenvectors
  	LWORK = 15
  	allocate(WORK(LWORK))

  	call ssyev('v','u',4,H11,4,W11,WORK,LWORK,INFO)
  	if (INFO /= 0) then
     	write(*,*) 'Problem with Lapack for H11 matrix. Error code:', INFO
     	stop
  	end if 

  	do i=1,4
     	grand_potential(i+7) = W11(i) - mu*2  ! grand potentials of H11
  	end do

  	do i=1,4
     	eigenvectors(i+7,8:11) = H11(1:4,i)    ! eigenvectors of H11
  	end do

  	deallocate(WORK)

  	!-------------------H21 and H12 (same)---------------------------------

  	! enter precalculated hamiltonian
 	H21(1,1) = 2*E(1) + E(2) + U; H21(1,2) = -t
 	H21(2,1) = -t;                H21(2,2) = E(1) + 2*E(2) + U

 	! solve the eigenvalues and eigenvectors
  	LWORK = 10
  	allocate(WORK(LWORK))

  	call ssyev('v','u',2,H21,2,W21,WORK,LWORK,INFO)
  	if (INFO /= 0) then
     	write(*,*) 'Problem with Lapack for H21 matrix. Error code:', INFO
     	stop
  	end if 

  	do i=1,2
     	grand_potential(i+11) = W21(i) - mu*3  ! grand potentials of H21
     	grand_potential(i+13) = W21(i) - mu*3  ! grand potentials of H12
  	end do

  	do i=1,2
     	eigenvectors(i+11,12:13) = H21(1:2,i)    ! eigenvectors of H21
     	eigenvectors(i+13,14:15) = H21(1:2,i)    ! eigenvectors of H12
  	end do

  	deallocate(WORK)

  	!-------------------------H22------------------------------------------

  	H22 = 2*(E(1) + E(2) + U)

  	W22 = H22

  	grand_potential(16) = W22 - 4*mu

  	eigenvectors(16,16) = 1

end subroutine hamiltonian

!************************************************************************************

subroutine transformations()

	PES_up = 0; PES_dn = 0
	IPES_up = 0; IPES_dn = 0
	phase_PES_up = 0; phase_PES_dn = 0
	phase_IPES_up = 0; phase_IPES_dn = 0

	PES_up(1,1) = 2; phase_PES_up(1,1) = 1
	PES_up(1,3) = 6; phase_PES_up(1,3) = -1 
	PES_up(1,4) = 10; phase_PES_up(1,4) = -1
	PES_up(1,5) = 8;  phase_PES_up(1,5) = -1 
	PES_up(1,7) = 14;  phase_PES_up(1,7) = 1
	PES_up(1,9) = 12; phase_PES_up(1,9) = 1
	PES_up(1,11) = 13; phase_PES_up(1,11) = 1
	PES_up(1,15) = 16; phase_PES_up(1,15) = -1

	PES_dn(1,1) = 4; phase_PES_dn(1,1) = 1
	PES_dn(1,2) = 10; phase_PES_dn(1,2) = 1
	PES_dn(1,3) = 9; phase_PES_dn(1,3) = -1 
	PES_dn(1,5) = 7; phase_PES_dn(1,4) = -1
	PES_dn(1,6) = 12;  phase_PES_dn(1,6) = -1
	PES_dn(1,8) = 14; phase_PES_dn(1,8) = -1
	PES_dn(1,11) = 15; phase_PES_dn(1,11) = 1
	PES_dn(1,13) = 16; phase_PES_dn(1,13) = 1

	IPES_up(1,2) = 1; phase_IPES_up(1,2) = 1
	IPES_up(1,6) = 3; phase_IPES_up(1,6) = -1 
	IPES_up(1,8) = 5; phase_IPES_up(1,8) = -1
	IPES_up(1,10) = 4; phase_IPES_up(1,10) = -1 
	IPES_up(1,12) = 9;  phase_IPES_up(1,12) = 1
	IPES_up(1,13) = 11; phase_IPES_up(1,13) = 1
	IPES_up(1,14) = 7; phase_IPES_up(1,14) = 1
	IPES_up(1,16) = 15; phase_IPES_up(1,16) = -1           !found by hand

	IPES_dn(1,4) = 1; phase_IPES_dn(1,4) = 1
	IPES_dn(1,7) = 5; phase_IPES_dn(1,7) = -1
	IPES_dn(1,9) = 3; phase_IPES_dn(1,9) = -1 
	IPES_dn(1,10) = 2; phase_IPES_dn(1,10) = 1
	IPES_dn(1,12) = 6; phase_IPES_dn(1,12) = -1
	IPES_dn(1,14) = 8; phase_IPES_dn(1,14) = -1
	IPES_dn(1,15) = 11; phase_IPES_dn(1,15) = 1
	IPES_dn(1,16) = 13; phase_IPES_dn(1,16) = 1

	PES_up(2,1) = 3; phase_PES_up(2,1) = 1
	PES_up(2,2) = 6; phase_PES_up(2,2) = 1 
	PES_up(2,4) = 9; phase_PES_up(2,4) = 1
	PES_up(2,5) = 11; phase_PES_up(2,5) = -1 
	PES_up(2,7) = 15; phase_PES_up(2,7) = -1
	PES_up(2,8) = 13; phase_PES_up(2,8) = -1
	PES_up(2,10) = 12; phase_PES_up(2,10) = 1
	PES_up(2,14) = 16; phase_PES_up(2,14) = -1

	PES_dn(2,1) = 5; phase_PES_dn(2,1) = 1
	PES_dn(2,2) = 8; phase_PES_dn(2,2) = 1
	PES_dn(2,3) = 11; phase_PES_dn(2,3) = 1 
	PES_dn(2,4) = 7; phase_PES_dn(2,4) = 1
	PES_dn(2,6) = 13; phase_PES_dn(2,6) = 1
	PES_dn(2,9) = 15; phase_PES_dn(2,9) = 1
	PES_dn(2,10) = 14; phase_PES_dn(2,10) = 1
	PES_dn(2,12) = 16; phase_PES_dn(2,12) = 1

	IPES_up(2,3) = 1; phase_IPES_up(2,3) = 1
	IPES_up(2,6) = 2; phase_IPES_up(2,6) = 1 
	IPES_up(2,9) = 4; phase_IPES_up(2,9) = 1
	IPES_up(2,11) = 5; phase_IPES_up(2,11) = -1 
	IPES_up(2,12) = 10; phase_IPES_up(2,12) = 1
	IPES_up(2,13) = 8; phase_IPES_up(2,13) = -1
	IPES_up(2,15) = 7; phase_IPES_up(2,15) = -1
	IPES_up(2,16) = 14; phase_IPES_up(2,16) = -1           !found by hand

	IPES_dn(2,5) = 1; phase_IPES_dn(2,5) = 1
	IPES_dn(2,7) = 4; phase_IPES_dn(2,7) = 1
	IPES_dn(2,8) = 2; phase_IPES_dn(2,8) = 1 
	IPES_dn(2,11) = 3; phase_IPES_dn(2,11) = 1
	IPES_dn(2,13) = 6; phase_IPES_dn(2,13) = 1
	IPES_dn(2,14) = 10; phase_IPES_dn(2,14) = 1
	IPES_dn(2,15) = 9; phase_IPES_dn(2,15) = 1
	IPES_dn(2,16) = 12; phase_IPES_dn(2,16) = 1

end subroutine transformations

end module routines
