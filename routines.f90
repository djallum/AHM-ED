module routines

implicit none

real :: Grand_potential(256)=0                    ! grand potentials (eigenenergies - mu*number electrons)
real :: Grand_potential_ground=0                 ! the lowest grand ensemble energy
real :: eigenvectors(256,256)=0                    ! the 64 eigenvectors
real :: v_ground(256)=0                           ! the ground state eigenvector
integer, dimension(4,256) :: PES_down=0, PES_up=0, IPES_down=0, IPES_up=0  !matrices for PES and IPES 
integer, dimension(4,256) :: phase_PES_down=0, phase_PES_up=0, phase_IPES_down=0, phase_IPES_up=0  !do get anticommutation sign right

contains

subroutine random_gen_seed()

implicit none

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

implicit none

real, intent(in) :: delta
real, intent(out) :: E(4)
real :: random
integer :: i ! counter

do i=1,4
   call random_number(random)           ! gives a random number between 0 and 1
   E(i) = delta*(random - 0.5)          ! centers the random numbers about 0
end do

end subroutine site_potentials

subroutine hamiltonian(E,t,U,mu)

implicit none

real, intent(in) :: E(4)
real, intent(in) :: mu 
real, intent(in) :: t
real, intent(in) :: U
real :: H00, W00
real, dimension(4,4) :: H10, W10, H30, H03
real, dimension(6,6) :: H20, W20
real, dimension(16,16) :: H11, W11
real, dimension(24,24) :: H21, W21, H12, W12
integer :: i,j ! counter

!------for lapack------------
integer :: INFO = 0
integer :: LWORK
integer, allocatable, dimension(:) :: WORK

!----------zero variables for each loop---------------------------------------

H00=0; H10=0; H20=0; H11=0; H21=0; H12=0
W00=0; W10=0; H20=0; W11=0; W21=0; W12=0

!-------------zero electron states---------------------------------------------

W00=H00

eigenvectors(1,1) = 1

Grand_potential(1) = W00 - mu*0

!-------one electron states (H10 & H01 are same)--------------------------

H10 = t        ! all off main diagonal terms are t (main diagonal terms are replaced in next step)

do i=1,4
   H10(i,i) = E(i)      ! main diagonal terms
end do

! solve the eigenvalues and eigenvectors
LWORK = 12
allocate(WORK(12))

call ssyev('v','u',4,H10,4,W10,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H10 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,4
   Grand_potential(i+1) = W10(i) - mu*1  ! grand potentials of H10
   Grand_potential(i+5) = W10(i) - mu*1  ! grand potentials of H01
end do

do i=1,4
   eigenvectors(i+1,2:5) = H10(1:4,i)    ! eigenvectors of H10
   eigenvectors(i+5,6:9) = H10(1:4,i)    ! eigenvectors of H01
end do 

!----------------------H20 and H02 (same)------------------------------------

! off main diagonal terms
H20(1,2) = t;  H20(1,3) = t;  H20(1,4) = -t; H20(1,5) = -t
H20(2,3) = t;  H20(2,4) = t;  H20(1,6) = -t
H20(3,5) = t;  H20(3,6) = t
H20(4,5) = t;  H20(4,6) = -t
H20(5,6) = t

H20(1,1) = E(1) + E(2); H20(2,2) = E(1) + E(3)
H20(3,3) = E(1) + E(4); H20(4,4) = E(2) + E(3)
H20(5,5) = E(2) + E(4); H20(6,6) = E(3) + E(4) 

! make it symetric
do i=1,6
   do j=1,6
      if(i>j) then
         H20(i,j) = H20(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 18
allocate(WORK(18))

call ssyev('v','u',6,H20,6,W20,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H20 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,6
   Grand_potential(i+9) = W20(i) - mu*2  ! grand potentials of H20
   Grand_potential(i+15) = W20(i) - mu*2  ! grand potentials of H02
end do

do i=1,6
   eigenvectors(i+9,10:15) = H20(1:6,i)    ! eigenvectors of H20
   eigenvectors(i+15,16:21) = H20(1:6,i)    ! eigenvectors of H02
end do 

!---------------------H11----------------------------------------

H11(1,3) = t;    H11(1,5) = t;    H11(1,8) = -t;  H11(1,10) = -t; H11(1,13) = t;  H11(1,14) = t
H11(2,4) = t;    H11(2,6) = t;    H11(2,7) = -t;  H11(2,9) = -t;  H11(2,13) = -t; H11(2,14) = -t
H11(3,5) = t;    H11(3,7) = t;    H11(3,12) = -t; H11(3,13) = t;  H11(3,15) = t
H11(4,6) = t;    H11(4,8) = t;    H11(4,11) = -t; H11(4,13) = -t; H11(4,15) = -t
H11(5,9) = t;    H11(5,11) = t;   H11(5,13) = t;  H11(5,16) = t
H11(6,10) = t;   H11(6,12) = t;   H11(6,13) = -t; H11(6,16) = -t
H11(7,9) = t;    H11(7,12) = -t;  H11(7,14) = t;  H11(7,15) = t
H11(8,10) = t;   H11(8,11) = -t;  H11(8,14) = -t; H11(8,15) = -t
H11(9,11) = t;   H11(9,14) = t;   H11(9,16) = t
H11(10,12) = t;  H11(10,14) = -t; H11(10,16) = -t
H11(11,15) = t;  H11(11,16) = t
H11(12,15) = -t: H11(12,16) = -t

H11(1,1) = E(1) + E(2); H11(2,2) = E(1) + E(2)
H11(3,3) = E(1) + E(3); H11(4,4) = E(1) + E(3)
H11(5,5) = E(1) + E(4); H11(6,6) = E(1) + E(4)
H11(7,7) = E(2) + E(3); H11(8,8) = E(2) + E(3)
H11(9,9) = E(2) + E(4); H11(10,10) = E(2) + E(4)
H11(11,11) = E(3) + E(4); H11(12,12) = E(3) + E(4)
H11(13,13) = 2*E(1) + U; H11(14,14) = 2*E(2) + U
H11(15,15) = 2*E(3) + U; H11(16,16) = 2*E(4) + U

! make it symetric
do i=1,16
   do j=1,16
      if(i>j) then
         H11(i,j) = H11(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 48
allocate(WORK(48))

call ssyev('v','u',16,H11,16,W11,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H11 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,16
   Grand_potential(i+21) = W11(i) - mu*2  ! grand potentials of H11
end do

do i=1,16
   eigenvectors(i+21,22:37) = H11(1:16,i)    ! eigenvectors of H11
end do 

!-------------H30 & H03 (same)-----------------------------------

H30(1,2) = t; H30(1,3) = -t; H30(1,4) = t
H30(2,3) = t; H30(2,4) = -t;
H30(3,4) = t

H30(1,1) = E(1) + E(2) + E(3); H30(2,2) = E(1) + E(2) + E(4)
H30(3,3) = E(1) + E(3) + E(4); H30(4,4) = E(2) + E(3) + E(4)

! make it symetric
do i=1,4
   do j=1,4
      if(i>j) then
         H30(i,j) = H30(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 12
allocate(WORK(12))

call ssyev('v','u',4,H30,4,W30,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H30 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,4
   Grand_potential(i+37) = W30(i) - mu*3  ! grand potentials of H30
   Grand_potential(i+41) = W30(i) - mu*3  ! grand potentials of H03
end do

do i=1,4
   eigenvectors(i+37,38:41) = H30(1:4,i)    ! eigenvectors of H30
   eigenvectors(i+41,42:45) = H30(1:4,i)    ! eigenvectors of H03
end do 

!---------------------H21 & H12--------------------------------------------

H21(1,4) = t;    H21(1,8) = -t;   H21(1,11) = t;   H21(1,13) = -t; H21(1,14) = t;  H21(1,16) = t;  H21(1,20) = -t
H21(2,5) = t;    H21(2,7) = -t;   H21(2,12) = t;   H21(2,14) = -t; H21(2,15) = t;  H21(2,16) = -t; H21(2,19) = t
H21(3,6) = t;    H21(3,9) = -t;   H21(3,10) = t;   H21(3,13) = t;  H21(3,15) = -t; H21(3,19) = -t; H21(3,20) = t
H21(4,7) = t;    H21(4,10) = -t;  H21(4,13) = -t;  H21(4,14) = t;  H21(4,18) = t   H21(4,22) = -t
H21(5,8) = t;    H21(5,12) = -t;  H21(5,14) = -t;  H21(5,17) = t;  H21(5,18) =-t;  H21(5,21) = t
H21(6,9) = t;    H21(6,11) = -t;  H21(6,13) = t;   H21(6,17) = -t; H21(6,21) = -t; H21(6,22) = t
H21(7,10) = t;   H21(7,15) = -t;  H21(7,16) = t;   H21(7,18) = t;  H21(7,24) = -t
H21(8,11) = t;   H21(8,16) = -t;  H21(8,17) = t;   H21(8,18) = -t; H21(8,23) = t
H21(9,12) = t;   H21(9,15) = t;   H21(9,17) = -t;  H21(9,23) = -t; H21(9,24) = t
H21(10,19) = -t; H21(10,20) = t;  H21(10,22) = t;  H21(10,24) = -t
H21(11,20) = -t; H21(11,21) = t;  H21(11,22) = -t; H21(11,23) = t
H21(12,19) = t;  H21(12,21) = -t; H21(12,23) = -t; H21(12,24) = t
H21(13,14) = -t; H21(13,15) = t;  H21(13,17) = t
H21(14,19) = t;  H21(14,21) = t
H21(15,16) = -t; H21(15,17) = t
H21(16,20) = t;  H21(16,23) = t
H21(17,18) = -t
H21(18,22) = t;  H21(18,24) = t
H21(19,20) = -t; H21(19,21) = t
H21(20,23) = t
H21(21,22) = -t
H21(22,24) = t
H21(23,24) = -t

H21(1,1) = E(1) + E(2) + E(3);  H21(2,2) = E(1) + E(2) + E(3);  H21(3,3) = E(1) + E(2) + E(3) 
H21(4,4) = E(1) + E(2) + E(4);  H21(5,5) = E(1) + E(2) + E(4);  H21(6,6) = E(1) + E(2) + E(4)
H21(7,7) = E(1) + E(3) + E(4);  H21(8,8) = E(1) + E(3) + E(4);  H21(9,9) = E(1) + E(3) + E(4)
H21(10,10) = E(2) + E(3) + E(4);  H21(11,11) = E(2) + E(3) + E(4);  H21(12,12) = E(2) + E(3) + E(4)
H21(13,13) = 2*E(1) + E(2) + U; H21(14,14) = 2*E(2) + E(1) + U
H21(15,15) = 2*E(1) + E(3) + U; H21(16,16) = 2*E(3) + E(1) + U
H21(17,17) = 2*E(1) + E(4) + U; H21(18,18) = 2*E(4) + E(1) + U
H21(19,19) = 2*E(2) + E(3) + U; H21(20,20) = 2*E(3) + E(2) + U
H21(21,21) = 2*E(2) + E(4) + U; H21(22,22) = 2*E(4) + E(2) + U
H21(23,23) = 2*E(3) + E(4) + U; H21(24,24) = 2*E(4) + E(3) + U

H12 = H21 ! not exactly but will make modifications in the next step


! make it symetric and do changes to upper right and lower left quadrants of H12 matrix
do i=1,4
   do j=1,4
      if(i>12 .and. j<12) then
         H12(i,j) = -H21(i,j)
      end if
      if(i<12 .and. j>12) then
         H12(i,j) = -H21(i,j)
      end if
      if(i>j) then
         H21(i,j) = H21(j,i)
         H12(i,j) = H12(j,i)
      end if
   end do
end do

end subroutine hamiltonian
end module routines
