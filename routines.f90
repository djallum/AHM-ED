module routines

implicit none 

real, dimension(3) :: W10=0, W20=0, W31=0, W32=0
real , dimension(9) :: W11=0, W21=0, W12=0, W22=0
real :: W00=0, W30=0, W03=0                                          ! matrices for eigenvalues
real :: omega(64)=0                                              ! grand potentials (eigen_energies - mu*number_electrons)
real :: omega_ground=0                                           ! the lowest grand ensemble energy
real :: eigenvectors(64,64)=0                                    ! the 16 eigenvectors
real :: v_ground(64)=0                                           ! the ground state eigenvector
integer, dimension(3,64) :: PES_down=0, PES_up=0, IPES_down=0, IPES_up=0  !matrices for PES and IPES 
integer, dimension(3,64) :: phase_PES_down=0, phase_PES_up=0, phase_IPES_down=0, phase_IPES_up=0  !do get anticommutation sign right

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
   E(i) = delta*(random - 0.5)          !centering the random numbers about 0
end do

end subroutine site_potentials

subroutine hamiltonian(E,t,U,mu)

real, intent(in) :: E(3)
real, intent(in) :: mu 
real, intent(in) :: t
real, intent(in) :: U
real :: H00, H30, H03, H33
real, dimension(3,3) :: H10, H01, H20, H31, H32
real, dimension(9,9) :: H11, H21, H12, H22
integer :: i,j ! counter

!------for lapack------------
integer :: INFO = 0
integer :: LWORK
integer, allocatable, dimension(:) :: WORK

!--------zero variables for each loop---------------

W10=0; W20=0; W31=0; W32=0; W11=0; W21=0; W12=0; W22=0
H10=0; H20=0; H31=0; H32=0; H11=0; H21=0; H12=0; H22=0

!--------zero electron state-----------------------------
H00 = 0

W00 = H00 ! eigenvalues

omega(1) = W00 - mu*0     ! grand potential

eigenvectors(1,1) = 1     ! eigenvector matrix

!-----------H10 and H01 (only H10 since H10 is same)------------------------------

H10 = t !the off-diagonal terms (the ones on diagonal are replaced in next step)

! do the on diagoanl terms (replace the ts that were there)
do i=1,3
   H10(i,i) = E(i)
end do

! solve the eigenvalues and eigenectors
LWORK = 10
allocate(WORK(10))

call ssyev('v','u',3,H10,3,W10,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H10 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,3
   omega(i+1) = W10(i) - mu*1  ! grandpotentials of H10
   omega(i+4) = W10(i) - mu*1  ! grandpotentials of H01
end do

do i=1,3
   eigenvectors(i+1,2:4) = H10(1:3,i)    ! eigenvectors of H10
   eigenvectors(i+4,5:7) = H10(1:3,i)    ! eigenvectors of H01
end do 

!--------H11---------------------------------------

! enter the off diagonal terms for the upper half of the matrix
H11(1,2) = t;  H11(1,3) = -t; H11(1,4) = t;  H11(1,5) = -t
H11(2,4) = t;  H11(2,5) = t;  H11(2,8) = -t
H11(3,4) = -t; H11(3,6) = t;  H11(3,7) = -t
H11(4,7) = t;  H11(4,8) = -t
H11(5,7) = t;  H11(5,8) = t
H11(6,8) = t;  H11(6,9) = t 
H11(7,9) = t
H11(8,9) = -t

!make it symetric
do i=1,9
   do j=1,9
      if(i>j) then
         H11(i,j) = H11(j,i)
      end if
   end do
end do

! enter the on diagonal terms
H11(1,1) = 2*E(1) + U
H11(2,2) = E(1) + E(2); H11(3,3) = E(1) + E(2)
H11(4,4) = 2*E(2) + U
H11(5,5) = E(1) + E(3); H11(6,6) = E(1) + E(3)
H11(7,7) = E(2) + E(3); H11(8,8) = E(2) + E(3)
H11(9,9) = 2*E(3) + U

LWORK = 30
allocate(WORK(30))

call ssyev('v','u',9,H11,9,W11,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H11 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,9
   omega(i+7) = W11(i) - mu*2  ! grandpotentials of H11
end do

do i=1,9
   eigenvectors(i+7,8:16) = H11(1:9,i)  ! eigenvectors of H11
end do

!------H20 and H02 (only H20 since H02 is the same)------------------------

H20 = 0

H20(1,1) = E(1) + E(2); H20(2,2) = E(1) + E(3); H20(3,3) = E(2) + E(3) ! on diagonal terms

H20(1,2) = t; H20(1,3) = -t; H20(2,3) = t   ! off diagonal upper matrix  

! make it symetric
do i=1,3
   do j=1,3
      if(i>j) then
         H20(i,j) = H20(j,i)
      end if
   end do
end do

! lapack portion to solve the eigenvalues and eigenvectors

LWORK = 10
allocate(WORK(10))

call ssyev('v','u',3,H20,3,W20,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,3
   omega(i+16) = W20(i) - mu*2    ! grandpotentials of H20
   omega(i+19) = W20(i) - mu*2    ! grandpotentials of H02
end do

do i=1,3
   eigenvectors(i+16,17:19) = H20(1:3,i)    ! eigenvectors of H20
   eigenvectors(i+19,20:22) = H20(1:3,i)    ! eigenvectors of H02
end do 

!-------H21 and H12 (H12 = -H21  (for off diagonal entries))------------------

H21 = 0; H12 = 0

! off diagonal entries of upper matrix
H21(1,4) = -t; H21(1,5) = t;  H21(1,7) = t;  H21(1,9) = -t
H21(2,5) = -t; H21(2,6) = t;  H21(2,7) = -t; H21(2,8) = t
H21(3,4) = t;  H21(3,6) = -t; H21(3,8) = -t; H21(3,9) = t
H21(4,5) = -t; H21(4,6) = t
H21(5,8) = t
H21(6,7) = -t
H21(7,9) = t
H21(8,9) = -t

!make it symetric
do i=1,9
   do j=1,9
      if(i>j) then
         H21(i,j) = H21(j,i)
      end if
   end do
end do

H12 = -H21

! on diagonal entries
H21(1,1) = E(1) + E(2) + E(3); H12(1,1) = E(1) + E(2) + E(3)
H21(2,2) = E(1) + E(2) + E(3); H12(2,2) = E(1) + E(2) + E(3)
H21(3,3) = E(1) + E(2) + E(3); H12(3,3) = E(1) + E(2) + E(3)
H21(4,4) = 2*E(1) + E(2) + U;  H12(4,4) = 2*E(1) + E(2) + U
H21(5,5) = E(1) + 2*E(2) + U;  H12(5,5) = E(1) + 2*E(2) + U
H21(6,6) = 2*E(1) + E(3) + U;  H12(6,6) = 2*E(1) + E(3) + U
H21(7,7) = E(1) + 2*E(3) + U;  H12(7,7) = E(1) + 2*E(3) + U
H21(8,8) = 2*E(2) + E(3) + U;  H12(8,8) = 2*E(2) + E(3) + U
H21(9,9) = E(2) + 2*E(3) + U;  H12(9,9) = E(2) + 2*E(3) + U

!lapack section

LWORK = 30
allocate(WORK(30))

call ssyev('v','u',9,H21,9,W21,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H21 matrix. Error code:', INFO
   stop
end if 

call ssyev('v','u',9,H12,9,W12,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H12 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,9
   omega(i+22) = W21(i) - mu*3  ! grandpotentials of H21
   omega(i+31) = W12(i) - mu*3  ! grandpotentials of H12
end do

do i=1,9
   eigenvectors(i+22,23:31) = H21(1:9,i)  ! eigenvectors of H21
   eigenvectors(i+31,32:40) = H12(1:9,i)  ! eigenvectors of H12
end do

!------H30 and H03 (same)--------------------------------------------------

H30 = E(1) + E(2) + E(3)

omega(41) = H30 - 3*mu
omega(42) = omega(41)

eigenvectors(41,41) = 1
eigenvectors(42,42) = 1

!------H31 and H13 (they are the same)-------------------------------------------------

H31(1,1) = 2*E(1) + E(2) + E(3); H31(2,2) = E(1) + 2*E(2) + E(3); H31(3,3) = E(1) + E(2) + 2*E(3) ! on diagonal terms

H31(1,2) = -t; H31(1,3) = t; H31(2,3) = -t   ! off diagonal upper matrix  

! make it symetric
do i=1,3
   do j=1,3
      if(i>j) then
         H31(i,j) = H31(j,i)
      end if
   end do
end do

! lapack portion to solve the eigenvalues and eigenvectors

LWORK = 10
allocate(WORK(10))

call ssyev('v','u',3,H31,3,W31,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H31 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,3
   omega(i+42) = W31(i) - mu*4    ! grandpotentials of H31
   omega(i+45) = W31(i) - mu*4    ! grandpotentials of H13
end do

do i=1,3
   eigenvectors(i+42,43:45) = H31(1:3,i)    ! eigenvectors of H31
   eigenvectors(i+45,46:48) = H31(1:3,i)    ! eigenvectors of H13
end do 

!------H22 (9x9 matrix)-----------------------------------------

! off diagonal entries of upper triangular part of matrix
H22(1,3) = -t; H22(1,6) = t;  H22(1,7) = t;  H22(1,8) = t
H22(2,4) = -t; H22(2,5) = t;  H22(2,7) = -t; H22(2,8) = -t
H22(3,5) = -t; H22(3,7) = t;  H22(3,9) = t
H22(4,6) = -t; H22(4,7) = -t; H22(4,9) = -t
H22(5,8) = t;  H22(5,9) = t;
H22(6,8) = -t; H22(6,9) = -t

!make it symetric
do i=1,9
   do j=1,9
      if(i>j) then
         H22(i,j) = H22(j,i)
      end if
   end do
end do

! on diagonal entries
H22(1,1) = 2*E(1) + E(2) + E(3) + U
H22(2,2) = 2*E(1) + E(2) + E(3) + U
H22(3,3) = E(1) + 2*E(2) + E(3) + U
H22(4,4) = E(1) + 2*E(2) + E(3) + U
H22(5,5) = E(1) + E(2) + 2*E(3) + U 
H22(6,6) = E(1) + E(2) + 2*E(3) + U  
H22(7,7) = 2*E(1) + 2*E(2) + 2*U  
H22(8,8) = 2*E(1) + 2*E(3) + 2*U 
H22(9,9) = 2*E(2) + 2*E(3) + 2*U

!lapack section

LWORK = 30
allocate(WORK(30))

call ssyev('v','u',9,H22,9,W22,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H22 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,9
   omega(i+48) = W22(i) - mu*4  ! grandpotentials of H11
end do

do i=1,9
   eigenvectors(i+48,49:57) = H22(1:9,i)  ! eigenvectors of H11
end do

!------H32 and H23 (same)--------------------------

H32 = -t ! all off diagonal entries are -t the on diagonal are replaced in next step

H32(1,1) = 2*E(1) + 2*E(2) + E(3) + 2*U ! on diagonal entries
H32(2,2) = 2*E(1) + E(2) + 2*E(3) + 2*U 
H32(3,3) = E(1) + 2*E(2) + 2*E(3) + 2*U 

! lapack portion to solve the eigenvalues and eigenvectors

LWORK = 10
allocate(WORK(10))

call ssyev('v','u',3,H32,3,W32,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H32 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,3
   omega(i+57) = W32(i) - mu*5    ! grandpotentials of H20
   omega(i+60) = W32(i) - mu*5    ! grandpotentials of H02
end do

do i=1,3
   eigenvectors(i+57,58:60) = H32(1:3,i)    ! eigenvectors of H20
   eigenvectors(i+60,61:63) = H32(1:3,i)    ! eigenvectors of H02
end do 

!---------H33--------------------------------

H33 = 2*E(1) + 2*E(2) + 2*E(3) + 3*U 

omega(64) = H33 - mu*6

eigenvectors(64,64) = 1

end subroutine hamiltonian

subroutine transformations()

integer :: i,j ! counters

PES_up = 0; PES_down = 0
IPES_up = 0; IPES_down = 0
phase_PES_up = 0; phase_PES_down = 0
phase_IPES_up = 0; phase_IPES_down = 0

PES_up(1,1) = 2;   phase_PES_up(1,1) = 1  ! PES_up(i,j) = k means a PES from cite i of state k makes state j
PES_up(1,3) = 17;  phase_PES_up(1,3) = -1 
PES_up(1,4) = 18;  phase_PES_up(1,4) = -1
PES_up(1,5) = 14;  phase_PES_up(1,5) = -1 
PES_up(1,6) = 8;   phase_PES_up(1,6) = -1
PES_up(1,7) = 9;   phase_PES_up(1,7) = -1
PES_up(1,10) = 23; phase_PES_up(1,10) = 1
PES_up(1,11) = 26; phase_PES_up(1,11) = 1
PES_up(1,12) = 28; phase_PES_up(1,12) = 1
PES_up(1,13) = 24; phase_PES_up(1,13) = 1 
PES_up(1,15) = 27; phase_PES_up(1,15) = 1
PES_up(1,16) = 29; phase_PES_up(1,16) = 1 
PES_up(1,19) = 41; phase_PES_up(1,19) = 1
PES_up(1,20) = 35; phase_PES_up(1,20) = 1
PES_up(1,21) = 37; phase_PES_up(1,21) = 1
PES_up(1,22) = 34; phase_PES_up(1,22) = 1
PES_up(1,25) = 43; phase_PES_up(1,25) = -1 
PES_up(1,30) = 44; phase_PES_up(1,30) = -1
PES_up(1,31) = 45; phase_PES_up(1,31) = -1
PES_up(1,32) = 50; phase_PES_up(1,32) = -1
PES_up(1,33) = 49; phase_PES_up(1,33) = -1
PES_up(1,36) = 55; phase_PES_up(1,36) = -1
PES_up(1,38) = 56; phase_PES_up(1,38) = -1
PES_up(1,39) = 51; phase_PES_up(1,39) = -1
PES_up(1,40) = 53; phase_PES_up(1,40) = -1 
PES_up(1,42) = 46; phase_PES_up(1,42) = -1
PES_up(1,47) = 61; phase_PES_up(1,47) = 1 
PES_up(1,48) = 62; phase_PES_up(1,48) = 1
PES_up(1,52) = 58; phase_PES_up(1,52) = 1
PES_up(1,54) = 59; phase_PES_up(1,54) = 1 
PES_up(1,57) = 60; phase_PES_up(1,57) = 1
PES_up(1,63) = 64; phase_PES_up(1,63) = -1 

PES_down(1,1) = 5;   phase_PES_down(1,1) = 1
PES_down(1,2) = 14;  phase_PES_down(1,2) = 1
PES_down(1,3) = 11;  phase_PES_down(1,3) = -1 
PES_down(1,4) = 12;  phase_PES_down(1,4) = -1
PES_down(1,6) = 20;  phase_PES_down(1,6) = -1
PES_down(1,7) = 21;  phase_PES_down(1,7) = -1
PES_down(1,8) = 35;  phase_PES_down(1,8) = -1
PES_down(1,9) = 37;  phase_PES_down(1,9) = -1
PES_down(1,10) = 33; phase_PES_down(1,10) = 1
PES_down(1,13) = 32; phase_PES_down(1,13) = 1
PES_down(1,15) = 36; phase_PES_down(1,15) = 1 
PES_down(1,16) = 38; phase_PES_down(1,16) = 1
PES_down(1,17) = 26; phase_PES_down(1,17) = -1
PES_down(1,18) = 28; phase_PES_down(1,18) = -1
PES_down(1,19) = 25; phase_PES_down(1,19) = 1
PES_down(1,22) = 42; phase_PES_down(1,22) = 1
PES_down(1,23) = 49; phase_PES_down(1,23) = 1
PES_down(1,24) = 50; phase_PES_down(1,24) = 1
PES_down(1,27) = 55; phase_PES_down(1,27) = 1 
PES_down(1,29) = 56; phase_PES_down(1,29) = 1
PES_down(1,30) = 52; phase_PES_down(1,30) = -1
PES_down(1,31) = 54; phase_PES_down(1,31) = -1
PES_down(1,34) = 46; phase_PES_down(1,34) = 1
PES_down(1,39) = 47; phase_PES_down(1,39) = -1
PES_down(1,40) = 48; phase_PES_down(1,40) = -1
PES_down(1,41) = 43; phase_PES_down(1,41) = 1
PES_down(1,44) = 58; phase_PES_down(1,44) = -1 
PES_down(1,45) = 59; phase_PES_down(1,45) = -1
PES_down(1,51) = 61; phase_PES_down(1,51) = -1
PES_down(1,53) = 62; phase_PES_down(1,53) = -1
PES_down(1,57) = 63; phase_PES_down(1,57) = 1
PES_down(1,60) = 64; phase_PES_down(1,60) = 1

PES_up(2,1) = 3;   phase_PES_up(2,1) = 1  ! PES_up(i,j) = k means a PES from cite i of state k makes state j
PES_up(2,2) = 17;  phase_PES_up(2,2) = 1 
PES_up(2,4) = 19;  phase_PES_up(2,4) = -1
PES_up(2,5) = 11;  phase_PES_up(2,5) = 1
PES_up(2,6) = 15;  phase_PES_up(2,6) = -1
PES_up(2,7) = 10;  phase_PES_up(2,7) = -1
PES_up(2,8) = 27;  phase_PES_up(2,8) = -1
PES_up(2,9) = 23;  phase_PES_up(2,9) = -1
PES_up(2,12) = 25; phase_PES_up(2,12) = -1
PES_up(2,13) = 30; phase_PES_up(2,13) = 1 
PES_up(2,14) = 26; phase_PES_up(2,14) = 1
PES_up(2,16) = 31; phase_PES_up(2,16) = 1 
PES_up(2,18) = 41; phase_PES_up(2,18) = -1
PES_up(2,20) = 36; phase_PES_up(2,20) = -1
PES_up(2,21) = 33; phase_PES_up(2,21) = -1
PES_up(2,22) = 39; phase_PES_up(2,22) = 1
PES_up(2,24) = 44; phase_PES_up(2,24) = 1 
PES_up(2,28) = 43; phase_PES_up(2,28) = -1
PES_up(2,29) = 45; phase_PES_up(2,29) = 1
PES_up(2,32) = 52; phase_PES_up(2,32) = 1
PES_up(2,34) = 51; phase_PES_up(2,34) = 1
PES_up(2,35) = 55; phase_PES_up(2,35) = -1
PES_up(2,37) = 49; phase_PES_up(2,37) = -1
PES_up(2,38) = 54; phase_PES_up(2,38) = 1
PES_up(2,40) = 57; phase_PES_up(2,40) = -1 
PES_up(2,42) = 47; phase_PES_up(2,42) = 1
PES_up(2,46) = 61; phase_PES_up(2,46) = 1 
PES_up(2,48) = 63; phase_PES_up(2,48) = -1
PES_up(2,50) = 58; phase_PES_up(2,50) = 1
PES_up(2,53) = 60; phase_PES_up(2,53) = -1 
PES_up(2,56) = 59; phase_PES_up(2,56) = 1
PES_up(2,62) = 64; phase_PES_up(2,62) = -1 

PES_down(2,1) = 6;   phase_PES_down(2,1) = 1
PES_down(2,2) = 8;   phase_PES_down(2,2) = 1
PES_down(2,3) = 15;  phase_PES_down(2,3) = 1 
PES_down(2,4) = 13;  phase_PES_down(2,4) = -1
PES_down(2,5) = 20;  phase_PES_down(2,5) = 1
PES_down(2,7) = 22;  phase_PES_down(2,7) = -1
PES_down(2,9) = 34;  phase_PES_down(2,9) = -1
PES_down(2,10) = 39; phase_PES_down(2,10) = -1
PES_down(2,11) = 36; phase_PES_down(2,11) = 1
PES_down(2,12) = 32; phase_PES_down(2,12) = -1
PES_down(2,14) = 35; phase_PES_down(2,14) = 1
PES_down(2,16) = 40; phase_PES_down(2,16) = 1
PES_down(2,17) = 27; phase_PES_down(2,17) = 1
PES_down(2,18) = 24; phase_PES_down(2,18) = -1
PES_down(2,19) = 30; phase_PES_down(2,19) = -1
PES_down(2,21) = 42; phase_PES_down(2,21) = -1
PES_down(2,23) = 51; phase_PES_down(2,23) = -1
PES_down(2,25) = 52; phase_PES_down(2,25) = -1
PES_down(2,26) = 55; phase_PES_down(2,26) = 1
PES_down(2,28) = 50; phase_PES_down(2,28) = -1
PES_down(2,29) = 53; phase_PES_down(2,29) = 1
PES_down(2,31) = 57; phase_PES_down(2,31) = 1
PES_down(2,33) = 47; phase_PES_down(2,33) = -1
PES_down(2,37) = 46; phase_PES_down(2,37) = -1
PES_down(2,38) = 48; phase_PES_down(2,38) = 1
PES_down(2,41) = 44; phase_PES_down(2,41) = -1
PES_down(2,43) = 58; phase_PES_down(2,43) = -1 
PES_down(2,45) = 60; phase_PES_down(2,45) = 1
PES_down(2,49) = 61; phase_PES_down(2,49) = -1
PES_down(2,54) = 63; phase_PES_down(2,54) = 1
PES_down(2,56) = 62; phase_PES_down(2,56) = 1
PES_down(2,59) = 64; phase_PES_down(2,59) = 1

PES_up(3,1) = 4;   phase_PES_up(3,1) = 1  ! PES_up(i,j) = k means a PES from cite i of state k makes state j
PES_up(3,2) = 18;  phase_PES_up(3,2) = 1 
PES_up(3,3) = 19;  phase_PES_up(3,3) = 1
PES_up(3,5) = 12;  phase_PES_up(3,5) = 1 
PES_up(3,6) = 13;  phase_PES_up(3,6) = 1
PES_up(3,7) = 16;  phase_PES_up(3,7) = -1
PES_up(3,8) = 24;  phase_PES_up(3,8) = 1
PES_up(3,9) = 29;  phase_PES_up(3,9) = -1
PES_up(3,10) = 31; phase_PES_up(3,10) = -1
PES_up(3,11) = 25; phase_PES_up(3,11) = 1 
PES_up(3,14) = 28; phase_PES_up(3,14) = 1
PES_up(3,15) = 30; phase_PES_up(3,15) = 1 
PES_up(3,17) = 41; phase_PES_up(3,17) = 1
PES_up(3,20) = 32; phase_PES_up(3,20) = 1
PES_up(3,21) = 38; phase_PES_up(3,21) = -1
PES_up(3,22) = 40; phase_PES_up(3,22) = -1
PES_up(3,23) = 45; phase_PES_up(3,23) = -1 
PES_up(3,26) = 43; phase_PES_up(3,26) = 1
PES_up(3,27) = 44; phase_PES_up(3,27) = 1 
PES_up(3,33) = 54; phase_PES_up(3,33) = -1
PES_up(3,34) = 53; phase_PES_up(3,34) = -1
PES_up(3,35) = 50; phase_PES_up(3,35) = 1
PES_up(3,36) = 49; phase_PES_up(3,36) = 1
PES_up(3,37) = 52; phase_PES_up(3,37) = -1
PES_up(3,39) = 57; phase_PES_up(3,39) = -1 
PES_up(3,42) = 48; phase_PES_up(3,42) = -1
PES_up(3,46) = 62; phase_PES_up(3,46) = -1 
PES_up(3,47) = 63; phase_PES_up(3,47) = -1
PES_up(3,49) = 59; phase_PES_up(3,49) = -1
PES_up(3,51) = 60; phase_PES_up(3,51) = -1 
PES_up(3,55) = 58; phase_PES_up(3,55) = 1
PES_up(3,61) = 64; phase_PES_up(3,61) = -1 

PES_down(3,1) = 7;   phase_PES_down(3,1) = 1
PES_down(3,2) = 9;   phase_PES_down(3,2) = 1
PES_down(3,3) = 10;  phase_PES_down(3,3) = 1 
PES_down(3,4) = 16;  phase_PES_down(3,4) = 1
PES_down(3,5) = 21;  phase_PES_down(3,5) = 1
PES_down(3,6) = 22;  phase_PES_down(3,6) = 1
PES_down(3,8) = 34;  phase_PES_down(3,8) = 1
PES_down(3,11) = 33; phase_PES_down(3,11) = 1
PES_down(3,12) = 38; phase_PES_down(3,12) = 1
PES_down(3,13) = 40; phase_PES_down(3,13) = 1
PES_down(3,14) = 37; phase_PES_down(3,14) = 1 
PES_down(3,15) = 39; phase_PES_down(3,15) = 1
PES_down(3,17) = 23; phase_PES_down(3,17) = 1
PES_down(3,18) = 29; phase_PES_down(3,18) = 1
PES_down(3,19) = 31; phase_PES_down(3,19) = 1
PES_down(3,20) = 42; phase_PES_down(3,20) = 1
PES_down(3,24) = 53; phase_PES_down(3,24) = 1
PES_down(3,25) = 54; phase_PES_down(3,25) = 1
PES_down(3,26) = 49; phase_PES_down(3,26) = 1 
PES_down(3,27) = 51; phase_PES_down(3,27) = 1
PES_down(3,28) = 56; phase_PES_down(3,28) = 1
PES_down(3,30) = 57; phase_PES_down(3,30) = 1
PES_down(3,32) = 48; phase_PES_down(3,32) = 1
PES_down(3,35) = 46; phase_PES_down(3,35) = 1
PES_down(3,36) = 47; phase_PES_down(3,36) = 1
PES_down(3,41) = 45; phase_PES_down(3,41) = 1
PES_down(3,43) = 59; phase_PES_down(3,43) = 1 
PES_down(3,44) = 60; phase_PES_down(3,44) = 1
PES_down(3,50) = 62; phase_PES_down(3,50) = 1
PES_down(3,52) = 63; phase_PES_down(3,52) = 1
PES_down(3,55) = 61; phase_PES_down(3,55) = 1
PES_down(3,58) = 64; phase_PES_down(3,58) = 1


cites: do i=1,3  ! calculating the IPES matrices 
   do j=1,64
      if (PES_down(i,j) /= 0) then
         phase_IPES_down(i,PES_down(i,j)) = phase_PES_down(i,j)
         IPES_down(i,PES_down(i,j)) = j
      end if
      if (PES_up(i,j) /= 0) then
         IPES_up(i,PES_up(i,j)) = j
         phase_IPES_up(i,PES_up(i,j)) = phase_PES_up(i,j)
      end if
   end do
end do cites

end subroutine transformations

end module routines
