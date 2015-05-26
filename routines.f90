module routines

implicit none 

real, dimension(3) :: W10, W20
real , dimension(9) :: W11
real :: W00                                                    ! matrices for eigenvalues
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
real :: H00, H30, H03
real, dimension(3,3) :: H10, H01, H20, H02 
real, dimension(9,9) :: H11
integer :: i,j ! counter

!------for lapack------------
integer :: INFO = 0
integer :: LWORK
integer, allocatable, dimension(:) :: WORK

H00 = 0

W00 = H00 ! eigenvalues

omega(1) = W00 - mu*0   ! grand potential

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
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
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

do i=1,9
   write(*,*), H11(i,:)
end do

LWORK = 30
allocate(WORK(30))

call ssyev('v','u',9,H11,9,W11,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
   stop
end if 

do i=1,9
   omega(i+7) = W11(i) - mu*2  ! grandpotentials of H11
end do

do i=1,9
   eigenvectors(i+7,8:16) = H11(1:9,i)  ! eigenvectors of H11
end do

deallocate(WORK)

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
   eigenvectors(i+16,16:18) = H10(1:3,i)    ! eigenvectors of H20
   eigenvectors(i+19,19:21) = H10(1:3,i)    ! eigenvectors of H02
end do 

!-------make H21 and H12------------------




end subroutine hamiltonian

end module routines
