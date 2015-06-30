program main

! for one iteration of the two site problem. All the code is contained without modules for the test. (150519)
! fock state basis in order is (+ is up, - is down) :
! 00,0+,+0,0-,-0,++,--,-+,+-,02,20,+2,2+,-2,2-,22

implicit none

integer, parameter :: nsites = 2
real, parameter :: t=-0.5 ! hopping
real, parameter :: U=5 ! on-site interaction
real, parameter :: mu=0 ! chemical potential (half-filling)
real, parameter :: delta =5
real :: E(nsites)
real :: H0=0,H1(2,2)=0,H2(2,2)=0,H3=0,H4=0,H5(4,4)=0,H6(2,2)=0,H7(2,2)=0,H8=0  !hamiltonian sub matrices
real :: W0=0,W1(2)=0.,W2(2)=0,W3=0.,W4=0,W5(4)=0,W6(2)=0,W7(2)=0,W8=0  !matrices for eigenvectors
real :: omega(16) ! grand esemble energies (eigen_energies - mu*number_electrons)
real :: omega_ground=0 ! the minimum grand ensemble energy
real :: eigenvector(16,16)=0 ! the 16 eigenvectors
real :: v_ground(16)=0  !ground state eigenvector
integer :: location(1)  !will store the location in the omega array of the lowest energy
integer :: PES_up1(16,16),PES_up2(16,16), PES_down1(16,16), PES_down2 ! matrices that operate on eigenstate for photo emmisions (up2 mean remove up spin from cite 2 etc.)
integer :: IPES_up1,IPES_up2,IPES_down1, IPES_down2 ! matrices for inverse photo emmisions (up2 mean add up spin on cite 2 etc.)

!-----(for lapack)----------------
integer :: INFO = 0
integer :: LWORK = 16
integer :: WORK(16)
!-----(for random_seed)---------------
integer :: seed_size, clock, i
integer, allocatable, dimension(:) :: seed
!-----(for site_potentials)-----------
real :: random1, random2

!-----(random_seed)--------------------

call random_seed(size=seed_size)
allocate(seed(seed_size))
call system_clock(count = clock)
seed=clock + 37*(/(i-1,i=1,seed_size)/)

call random_seed(put=seed)

deallocate(seed)

!-----(site_potentials)--------------------------

call random_number(random1)
call random_number(random2)
random1 = random1 - 0.5          !centering the random numbers about 0
random2 = random2 - 0.5
E(1) = delta*max(random1,random2)
E(2) = delta*min(random1,random2)

!-------(making matrices and solving eigenstates and eigenvaleus)-----------------------
H0 = 0
W0 = H0
omega(1) = W0
eigenvector(1,1) = 1

!---one electrons states--------------
H1(1,1) = E(1)
H1(2,2) = E(2)
H1(2,1) = t
H1(1,2) = t

H2 = H1

call ssyev('v','u',2,H1,2,W1,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
   stop
end if

omega(2) = W1(1) - mu*1
omega(3) = W1(2) - mu*1 
eigenvector(2,2:3) = H1(1:2,1)
eigenvector(3,2:3) = H1(1:2,2)

call ssyev('v','u',2,H2,2,W2,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H2 matrix. Error code:', INFO
   stop
end if

omega(4) = W2(1) - mu*1
omega(5) = W2(2) - mu*1
eigenvector(4,4:5) = H2(1:2,1)
eigenvector(5,4:5) = H2(1:2,2)

!---two electrons states--------------

H3 = E(1) + E(2)
H4 = E(1) + E(2)
W3 = H3
W4 = H4

H5(1,1) = E(1) + E(2)
H5(2,2) = E(1) + E(2)
H5(3,3) = 2*E(1) + U
H5(4,4) = 2*E(2) + U
H5(1,3) = t
H5(1,4) = t
H5(2,3) = -t
H5(2,4) = -t
H5(3,1) = t
H5(3,2) = -t
H5(4,1) = t
H5(4,2) = -t


call ssyev('v','u',4,H5,4,W5,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H5 matrix. Error code:', INFO
   stop
end if

omega(6) = W3 - mu*2
omega(7) = W4 - mu*2
omega(8) = W5(1) - mu*2
omega(9) = W5(2) - mu*2
omega(10) = W5(3) - mu*2
omega(11) = W5(4) - mu*2
eigenvector(6,6) = 1
eigenvector(7,7) = 1
eigenvector(8,8:11) = H5(1:4,1)
eigenvector(9,8:11) = H5(1:4,2)
eigenvector(10,8:11) = H5(1:4,3)
eigenvector(11,8:11) = H5(1:4,4)

!-----three electron states---------------

H6(1,1) = 2*E(1) + E(2) + U
H6(2,2) = E(1) + 2*E(2) + U
H6(2,1) = -t
H6(1,2) = -t

H7 = H6

call ssyev('v','u',2,H6,2,W6,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H6 matrix. Error code:', INFO
   stop
end if

omega(12) = W6(1) - mu*3 
omega(13) = W6(2) - mu*3 
eigenvector(12,12:13) = H6(1:2,1)
eigenvector(13,12:13) = H6(1:2,2)

call ssyev('v','u',2,H7,2,W7,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H7 matrix. Error code:', INFO
   stop
end if

omega(14) = W7(1) - mu*3 
omega(15) = W7(2) - mu*3
eigenvector(14,14:15) = H7(1:2,1)
eigenvector(15,14:15) = H7(1:2,2)

!------four electron states-------------------------

H8 = 2*E(1) + 2*E(2) + 2*U
W8 = H8

omega(16) = W8 - mu*4
eigenvector(16,16) = 1

!-----find ground state energy and eigenvector--------

omega_ground = minval(omega)   ! find the lowest grand ensemble energy

!find the corresponding eigenvector

location = minloc(omega)  !find the location of the lowest energy  
v_ground = eigenvector(location(1),:) !set v ground to the eigenvector corresponding to the lowest energy



end program main
