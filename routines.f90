module routines

implicit none 

real, dimension(3) :: W10, W20, W31, W32
real , dimension(9) :: W11, W21, W12, W22
real :: W00, W30, W03                                          ! matrices for eigenvalues
real :: omega(64)                                              ! grand potentials (eigen_energies - mu*number_electrons)
real :: omega_ground                                           ! the lowest grand ensemble energy
real :: eigenvectors(64,64)                                    ! the 16 eigenvectors
real :: v_ground(64)                                           ! the ground state eigenvector

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

H00 = 0

W00 = H00 ! eigenvalues

omega(1) = W00 - mu*0   ! grand potential

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

H11 = 0

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
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
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
   eigenvectors(i+16,17:19) = H10(1:3,i)    ! eigenvectors of H20
   eigenvectors(i+19,20:22) = H10(1:3,i)    ! eigenvectors of H02
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
H21(1,1) = E(1) + E(2) + E(3); H21(1,1) = E(1) + E(2) + E(3)
H21(2,2) = E(1) + E(2) + E(3); H21(2,2) = E(1) + E(2) + E(3)
H21(3,3) = E(1) + E(2) + E(3); H21(3,3) = E(1) + E(2) + E(3)
H21(4,4) = 2*E(1) + E(2) + U;  H21(4,4) = 2*E(1) + E(2) + U
H21(5,5) = E(1) + 2*E(2) + U;  H21(5,5) = E(1) + 2*E(2) + U
H21(6,6) = 2*E(1) + E(3) + U;  H21(6,6) = 2*E(1) + E(3) + U
H21(7,7) = E(1) + 2*E(3) + U;  H21(7,7) = E(1) + 2*E(3) + U
H21(8,8) = 2*E(2) + E(3) + U;  H21(8,8) = 2*E(2) + E(3) + U
H21(9,9) = E(2) + 2*E(3) + U;  H21(9,9) = E(2) + 2*E(3) + U

!lapack section

LWORK = 30
allocate(WORK(30))

call ssyev('v','u',9,H21,9,W21,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
   stop
end if 

call ssyev('v','u',9,H12,9,W12,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
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

H31 = 0

H31(1,1) = 2*E(1) + E(2) + E(3); H31(2,2) = E(1) + 2*E(2) + E(3); H31(3,3) = E(1) + E(2) + 2*E(3) ! on diagonal terms

H31(1,2) = -t; H20(1,3) = t; H20(2,3) = -t   ! off diagonal upper matrix  

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
   omega(i+42) = W31(i) - mu*4    ! grandpotentials of H20
   omega(i+45) = W31(i) - mu*4    ! grandpotentials of H02
end do

do i=1,3
   eigenvectors(i+42,43:45) = H31(1:3,i)    ! eigenvectors of H20
   eigenvectors(i+45,46:48) = H31(1:3,i)    ! eigenvectors of H02
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
   write(*,*) 'Problem with Lapack for H31 matrix. Error code:', INFO
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

omega(64) = H33 - mu*5

eigenvectors(64,64) = 1

end subroutine hamiltonian

end module routines
