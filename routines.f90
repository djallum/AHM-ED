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
real :: H00, W00, H40, W40, H44
real :: W10(4), W30(4), W41(4), W43(4), W20(6), W42(6), W11(16), W31(16), W13(16)
real :: W33(16), W21(24), W12(24), W32(24), W23(24), W22(36)
real, dimension(4,4) :: H10=0, H30=0, H41=0, H43=0
real, dimension(6,6) :: H20=0, H42=0
real, dimension(16,16) :: H11=0, H31=0, H13=0, H33=0
real, dimension(24,24) :: H21=0, H12=0, H32=0, H23=0
real, dimension(36,36) :: H22=0
integer :: i,j ! counter

!------for lapack------------
integer :: INFO = 0
integer :: LWORK
integer, allocatable, dimension(:) :: WORK

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
H11(12,15) = -t; H11(12,16) = -t

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

do i=1,16
   write(*,'(16F7.2)') (H11(i,j), j=1,16)
end do
! solve the eigenvalues and eigenvectors
LWORK = 50
allocate(WORK(50))

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
H21(4,7) = t;    H21(4,10) = -t;  H21(4,13) = -t;  H21(4,14) = t;  H21(4,18) = t;  H21(4,22) = -t
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
do i=1,24
   do j=1,24
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

! solve the eigenvalues and eigenvectors
LWORK = 75
allocate(WORK(75))

call ssyev('v','u',24,H21,24,W21,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H21 matrix. Error code:', INFO
   stop
end if 

call ssyev('v','u',24,H12,24,W12,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H12 matrix. Error code:', INFO
   stop
end if

deallocate(WORK)

do i=1,24
   Grand_potential(i+45) = W21(i) - mu*3  ! grand potentials of H21
   Grand_potential(i+69) = W12(i) - mu*3  ! grand potentials of H12
end do

do i=1,24
   eigenvectors(i+45,46:69) = H21(1:24,i)    ! eigenvectors of H21
   eigenvectors(i+69,70:93) = H12(1:24,i)    ! eigenvectors of H12
end do 

!----------------H40 & H04 (same)------------------------------------

H40 = E(1) + E(2) + E(3) + E(4)

W40 = H40

Grand_potential(94) = W40 - 4*mu   ! grand potential of H40
Grand_potential(95) = W40 - 4*mu   ! grand potential of H04

eigenvectors(94,94) = 1  ! eigenvector of H40
eigenvectors(95,95) = 1  ! eigenvector of H04

!---------------H31 & H13------------------------------------------

H31(1,5) = t; H31(1,6) = -t; H31(1,7) = t; H31(1,10) = t; H31(1,13) = -t; H31(1,16) = t
H31(2,7) = -t; H31(2,8) = -t; H31(2,9) = t; H31(2,10) = -t; H31(2,12) = t; H31(2,15) = -t
H31(3,6) = t; H31(3,9) = -t; H31(3,11) = t; H31(3,12) = -t; H31(3,13) = t; H31(3,14) = t
H31(4,5) = -t; H31(4,8) = t; H31(4,11) = -t; H31(4,14) = -t; H31(4,15) = t; H31(4,16) = -t
H31(5,6) = -t; H31(5,7) = t; H31(5,8) = t; H31(5,11) = -t
H31(6,7) = -t; H31(6,9) = t; H31(6,14) = -t
H31(7,12) = t; H31(7,15) = -t
H31(8,9) = -t; H31(8,10) = t; H31(8,11) = t
H31(9,10) = -t; H31(9,14) = t
H31(10,13) = t; H31(10,16) = -t
H31(11,12) = -t; H31(11,13) = t
H31(12,13) = -t; H31(12,15) = t
H31(13,16) = t
H31(14,15) = -t; H31(14,16) = t
H31(15,16) = -t

H31(1,1) = E(1) + E(2) + E(3) + E(4); H31(2,2) = H31(1,1); H31(3,3) = H31(1,1); H31(4,4) = H31(1,1)
H31(5,5) = 2*E(1) + E(2) + E(3) + U; H31(6,6) = E(1) + 2*E(2) + E(3) + U; H31(7,7) = E(1) + E(2) + 2*E(3) + U
H31(8,8) = 2*E(1) + E(2) + E(4) + U; H31(9,9) = E(1) + 2*E(2) + E(4) + U; H31(10,10) = E(1) + E(2) + 2*E(4) + U
H31(11,11) = 2*E(1) + E(3) + E(4) + U; H31(12,12) = E(1) + 2*E(3) + E(4) + U; H31(13,13) = E(1) + E(3) + 2*E(4) + U
H31(14,14) = 2*E(2) + E(3) + E(4) + U; H31(15,15) = E(2) + 2*E(3) + E(4) + U; H31(16,16) = E(2) + E(3) + 2*E(4) + U

H13 = H31

! make it symetric and do changes to first four rows of H13
do i=1,16
   do j=1,16
      if(i>5 .and. i/=j) then
         H13(i,j) = -H31(i,j)
      end if
      if(i>j) then
         H31(i,j) = H31(j,i)
         H13(i,j) = H13(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 48
allocate(WORK(48))

call ssyev('v','u',16,H31,16,W31,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H31 matrix. Error code:', INFO
   stop
end if 

call ssyev('v','u',16,H13,16,W13,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H13 matrix. Error code:', INFO
   stop
end if

deallocate(WORK)

do i=1,16
   Grand_potential(i+95) = W31(i) - mu*4  ! grand potentials of H31
   Grand_potential(i+111) = W13(i) - mu*4  ! grand potentials of H13
end do

do i=1,16
   eigenvectors(i+95,96:111) = H31(1:16,i)    ! eigenvectors of H31
   eigenvectors(i+111,112:127) = H13(1:16,i)    ! eigenvectors of H13
end do 

!-----------------------H22---------------------------------------------

H22(1,7) = t; H22(1,9) = -t; H22(1,13) = -t; H22(1,15) = t; H22(1,21) = t; H22(1,23) = -t; H22(1,27) = -t; H22(1,29) = t
H22(2,8) = t; H22(2,11) = t; H22(2,15) = -t; H22(2,17) = t; H22(2,19) = t; H22(2,21) = -t; H22(2,25) = t; H22(2,30) = t
H22(3,9) = t; H22(3,11) = -t; H22(3,14) = -t; H22(3,17) = -t; H22(3,20) = t; H22(3,23) = t; H22(3,26) = t; H22(3,28) = -t
H22(4,10) = -t; H22(4,12) = t; H22(4,13) = t; H22(4,18) = t; H22(4,19) = -t; H22(4,24) = -t; H22(4,25) = -t; H22(4,27) = t
H22(5,7) = -t; H22(5,12) = -t; H22(5,16) = t; H22(5,18) = -t; H22(5,20) = -t; H22(5,22) = t; H22(5,26) = -t; H22(5,29) = -t
H22(6,8) = -t; H22(6,10) = t; H22(6,14) = t;H22(6,16) = -t; H22(6,22) = -t; H22(6,24) = t; H22(6,28) = t; H22(6,30) = -t
H22(7,9) = -t; H22(7,12) = t; H22(7,13) = t; H22(7,20) = -t; H22(7,31) = t; H22(7,32) = t
H22(8,10) = -t; H22(8,11) = t; H22(8,14) = t; H22(8,19) = -t; H22(8,31) = -t; H22(8,32) = -t
H22(9,11) = -t; H22(9,15) = t; H22(9,26) = -t; H22(9,31) = t; H22(9,34) = t
H22(10,12) = -t; H22(10,16) = t; H22(10,25) = -t; H22(10,31) = -t; H22(10,34) = -t
H22(11,21) = t; H22(11,28) = -t; H22(11,32) = t; H22(11,34) = t
H22(12,22) = t; H22(12,27) = -t; H22(12,32) = -t; H22(12,34) = -t
H22(13,15) = -t; H22(13,18) = t; H22(13,19) = t; H22(13,31) = t; H22(13,33) = t
H22(14,16) = -t; H22(14,17) = t; H22(14,20) = t; H22(14,31) = -t; H22(14,33) = -t
H22(15,17) = -t; H22(15,25) = t; H22(15,31) = t; H22(15,35) = t
H22(16,18) = -t; H22(16,26) = t; H22(16,31) = -t; H22(16,35) = -t
H22(17,23) = t; H22(17,30) = -t; H22(17,33) = t; H22(17,35) = t
H22(18,24) = t; H22(18,29) = -t; H22(18,33) = -t; H22(18,35) = -t
H22(19,21) = -t; H22(19,24) = t; H22(19,32) = t; H22(19,33) = t
H22(20,22) = -t; H22(20,23) = t; H22(20,32) = -t; H22(20,33) = -t
H22(21,23) = -t; H22(21,27) = t; H22(21,32) = t; H22(21,36) = t
H22(22,24) = -t; H22(22,28) = t; H22(22,32) = -t; H22(22,36) = -t
H22(23,29) = t; H22(23,33) = t; H22(23,36) = t
H22(24,30) = t; H22(24,33) = -t; H22(24,36) = -t
H22(25,27) = -t; H22(25,30) = t; H22(25,34) = t; H22(25,35) = t
H22(26,28) = -t; H22(26,29) = t; H22(26,34) = -t; H22(26,35) = -t
H22(27,29) = -t; H22(27,34) = t; H22(27,36) = t
H22(28,30) = -t; H22(28,34) = -t; H22(28,36) = -t
H22(29,35) = t; H22(29,36) = t
H22(30,35) = -t; H22(30,36) = -t

H22(1,1) = E(1) + E(2) + E(3) + E(4); H22(3,3) = H22(1,1); H22(5,5) = H22(1,1)
H22(7,7) = 2*E(1) + E(2) + E(3) + U; H22(9,9) = E(1) + 2*E(2) + E(3) + U; H22(11,11) = E(1) + E(2) + 2*E(3) + U
H22(13,13) = 2*E(1) + E(2) + E(4) + U; H22(15,15) = E(1) + 2*E(2) + E(4) + U; H22(17,17) = E(1) + E(2) + 2*E(4) + U
H22(19,19) = 2*E(1) + E(3) + E(4) + U; H22(21,21) = E(1) + 2*E(3) + E(4) + U; H22(23,23) = E(1) + E(3) + 2*E(4) + U
H22(25,25) = 2*E(2) + E(3) + E(4) + U; H22(27,27) = E(2) + 2*E(3) + E(4) + U; H22(29,29) = E(2) + E(3) + 2*E(4) + U
H22(31,31) = 2*E(1) + 2*E(2) + 2*U; H22(32,32) = 2*E(1) + 2*E(3) + 2*U; H22(33,33) = 2*E(1) + 2*E(4) + 2*U
H22(34,34) = 2*E(2) + 2*E(3) + 2*U; H22(35,35) = 2*E(2) + 2*E(4) + 2*U; H22(36,36) = 2*E(3) + 2*E(4) + 2*U

do i=1,29,2
   H22(i+1,i+1) = H22(i,i)
end do

! make it symetric
do i=1,36
   do j=1,36
      if(i>j) then
         H22(i,j) = H22(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 108
allocate(WORK(108))

call ssyev('v','u',36,H22,36,W22,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H22 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,36
   Grand_potential(i+127) = W22(i) - mu*4  ! grand potentials of H22
end do

do i=1,36
   eigenvectors(i+127,128:163) = H22(1:36,i)    ! eigenvectors of H22
end do 

!-------------H41 & H14 (same)-------------------------------------

H41(1,2) = -t; H41(1,3) = t; H41(1,4) = -t
H41(2,3) = t; H41(2,4) = t
H41(3,4) = -t

H41(1,1) = 2*E(1) + E(2) + E(3) + E(4) + U; H41(2,2) = E(1) + 2*E(2) + E(3) + E(4) + U
H41(3,3) = E(1) + E(2) + 2*E(3) + E(4) + U; H41(4,4) = E(1) + E(2) + E(3) + 2*E(4) + U

! make it symetric
do i=1,4
   do j=1,4
      if(i>j) then
         H41(i,j) = H41(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 12
allocate(WORK(12))

call ssyev('v','u',4,H41,4,W41,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H41 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,4
   Grand_potential(i+163) = W41(i) - mu*5  ! grand potentials of H41
   Grand_potential(i+167) = W41(i) - mu*5  ! grand potentials of H14
end do

do i=1,4
   eigenvectors(i+163,164:167) = H41(1:4,i)    ! eigenvectors of H41
   eigenvectors(i+167,168:171) = H41(1:4,i)    ! eigenvectors of H14
end do 

!--------------H32 & H23----------------------------------

H32(1,4) = -t; H32(1,7) = t; H32(1,12) = -t; H32(1,13) = -t; H32(1,14) = t; H32(1,17) = t; H32(1,20) = -t
H32(2,5) = -t; H32(2,9) = t; H32(2,10) = -t; H32(2,14) = -t; H32(2,16) = t; H32(2,17) = -t; H32(2,19) = t
H32(3,6) = -t; H32(3,8) = t; H32(3,11) = -t; H32(3,13) = t; H32(3,16) = -t; H32(3,19) = -t; H32(3,20) = t
H32(4,7) = -t; H32(4,11) = t; H32(4,13) = -t; H32(4,15) = t; H32(4,18) = t; H32(4,23) = -t
H32(5,8) = -t; H32(5,10) = t; H32(5,15) = -t; H32(5,16) = t; H32(5,18) = -t; H32(5,22) = t
H32(6,9) = -t; H32(6,12) = t; H32(6,13) = t; H32(6,16) = -t; H32(6,22) = -t; H32(6,23) = t
H32(7,10) = -t; H32(7,14) = -t; H32(7,15) = t; H32(7,21) = t; H32(7,24) = -t
H32(8,11) = -t; H32(8,15) = -t; H32(8,19) = t; H32(8,21) = -t; H32(8,22) = t
H32(9,12) = -t; H32(9,14) = t; H32(9,19) = -t; H32(9,22) = -t; H32(9,24) = t
H32(10,17) = -t; H32(10,18) = t; H32(10,21) = t; H32(10,24) = -t
H32(11,18) = -t; H32(11,20) = t; H32(11,21) = -t; H32(11,23) = t
H32(12,17) = t; H32(12,20) = -t; H32(12,23) = -t; H32(12,24) = t
H32(13,14) = -t; H32(13,15) = -t; H32(13,16) = t
H32(14,15) = -t; H32(14,19) = t
H32(15,22) = t
H32(16,17) = -t; H32(17,18) = -t
H32(17,18) = -t; H32(17,20) = t
H32(18,23) = t
H32(19,20) = -t; H32(19,21) = -t
H32(20,21) = -t
H32(21,24) = t
H32(22,23) = -t; H32(22,24) = -t
H32(23,24) = -t

H32(1,1) = 2*E(1) + E(2) + E(3) + E(4) + U; H32(4,4) = E(1) + 2*E(2) + E(3) + E(4) + U
H32(7,7) = E(1) + E(2) + 2*E(3) + E(4) + U; H32(10,10) = E(1) + E(2) + E(3) + 2*E(4) + U
H32(13,13) = 2*E(1) + 2*E(2) + E(3) + 2*U; H32(14,14) = 2*E(1) + E(2) + 2*E(3) + 2*U; H32(15,15) = E(1) + 2*E(2) + 2*E(3) + 2*U
H32(16,16) = 2*E(1) + 2*E(2) + E(4) + 2*U; H32(17,17) = 2*E(1) + E(2) + 2*E(4) + 2*U; H32(18,18) = E(1) + 2*E(2) + 2*E(4) + 2*U
H32(19,19) = 2*E(1) + 2*E(3) + E(4) + 2*U; H32(20,20) = 2*E(1) + E(3) + 2*E(4) + 2*U; H32(21,21) = E(1) + 2*E(3) + 2*E(4) + 2*U
H32(22,22) = 2*E(2) + 2*E(3) + E(4) + 2*U; H32(23,23) = 2*E(2) + E(3) + 2*E(4) + 2*U; H32(24,24) = E(2) + 2*E(3) + 2*E(4) + 2*U

do i=1,10,3
   H32(i+1,i+1) = H32(i,i)
   H32(i+2,i+2) = H32(i,i)
end do

! make it symetric and do changes to upper right and lower left quadrants of H12 matrix
do i=1,24
   do j=1,24
      if(i>12 .and. j<12) then
         H23(i,j) = -H32(i,j)
      end if
      if(i<12 .and. j>12) then
         H23(i,j) = -H32(i,j)
      end if
      if(i>j) then
         H32(i,j) = H32(j,i)
         H23(i,j) = H23(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 75
allocate(WORK(75))

call ssyev('v','u',24,H32,24,W32,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H32 matrix. Error code:', INFO
   stop
end if 

call ssyev('v','u',24,H23,24,W23,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H23 matrix. Error code:', INFO
   stop
end if

deallocate(WORK)

do i=1,24
   Grand_potential(i+171) = W32(i) - mu*5  ! grand potentials of H32
   Grand_potential(i+195) = W23(i) - mu*5  ! grand potentials of H23
end do

do i=1,24
   eigenvectors(i+171,172:195) = H32(1:24,i)    ! eigenvectors of H32
   eigenvectors(i+195,196:219) = H23(1:24,i)    ! eigenvectors of H23
end do 

!----------------------H42 & H24 (same)------------------------------------------

H42(1,2) = -t; H42(1,3) = t; H32(1,4) = -t; H32(1,5) = t
H42(2,3) = -t; H32(2,4) = -t; H32(2,6) = t
H42(3,5) = -t; H42(3,6) = t
H42(4,5) = -t; H42(4,6) = -t
H42(5,6) = -t

H42(1,1) = 2*E(1) + 2*E(2) + E(3) + E(4) + 2*U; H42(2,2) = 2*E(1) + E(2) + 2*E(3) + E(4) + 2*U
H42(3,3) = 2*E(1) + E(2) + E(3) + 2*E(4) + 2*U; H42(4,4) = E(1) + 2*E(2) + 2*E(3) + E(4) + 2*U
H42(5,5) = E(1) + 2*E(2) + E(3) + 2*E(4) + 2*U; H42(6,6) = E(1) + E(2) + 2*E(3) + 2*E(4) + 2*U

! make it symetric
do i=1,6
   do j=1,6
      if(i>j) then
         H42(i,j) = H42(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 18
allocate(WORK(18))

call ssyev('v','u',6,H42,6,W42,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H42 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,6
   Grand_potential(i+219) = W42(i) - mu*6  ! grand potentials of H42
   Grand_potential(i+225) = W42(i) - mu*6  ! grand potentials of H24
end do

do i=1,6
   eigenvectors(i+219,220:225) = H42(1:6,i)    ! eigenvectors of H42
   eigenvectors(i+225,226:231) = H42(1:6,i)    ! eigenvectors of H24
end do 

!---------------------------H33-----------------------------------------

H33(1,3) = -t; H33(1,6) = t; H33(1,7) = -t; H33(1,10) = t; H33(1,13) = t; H33(1,14) = t
H33(2,4) = -t; H33(2,5) = t; H33(2,8) = -t; H33(2,9) = t; H33(2,13) = -t; H33(2,14) = -t
H33(3,5) = -t; H33(3,7) = -t; H33(3,12) = t; H33(3,13) = t; H33(3,15) = t
H33(4,6) = -t; H33(4,8) = -t; H33(4,11) = t; H33(4,13) = -t; H33(4,15) = -t
H33(5,9) = -t; H33(5,12) = t; H33(5,14) = t; H33(5,15) = t
H33(6,10) = -t; H33(6,11) = t; H33(6,14) = -t; H33(6,15) = -t
H33(7,9) = -t; H33(7,11) = -t; H33(7,13) = t; H33(7,16) = t
H33(8,10) = -t; H33(8,12) = -t; H33(8,13) = -t; H33(8,16) = -t
H33(9,11) = -t; H33(9,14) = t; H33(9,16) = t
H33(10,12) = -t; H33(10,14) = -t; H33(10,16) = -t
H33(11,15) = t; H33(11,16) = t
H33(12,15) = -t; H33(12,16) = -t

H33(1,1) = 2*E(1) + 2*E(2) + E(3) + E(4) + 2*U; H33(2,2) = 2*E(1) + 2*E(2) + E(3) + E(4) + 2*U
H33(3,3) = 2*E(1) + E(2) + 2*E(3) + E(4) + 2*U; H33(4,4) = 2*E(1) + E(2) + 2*E(3) + E(4) + 2*U
H33(5,5) = 2*E(1) + E(2) + E(3) + 2*E(4) + 2*U; H33(6,6) = 2*E(1) + E(2) + E(3) + 2*E(4) + 2*U
H33(7,7) = E(1) + 2*E(2) + 2*E(3) + E(4) + 2*U; H33(8,8) = E(1) + 2*E(2) + 2*E(3) + E(4) + 2*U
H33(9,9) = E(1) + 2*E(2) + E(3) + 2*E(4) + 2*U; H33(10,10) = E(1) + 2*E(2) + E(3) + 2*E(4) + 2*U
H33(11,11) = E(1) + E(2) + 2*E(3) + 2*E(4) + 2*U; H33(12,12) = E(1) + E(2) + 2*E(3) + 2*E(4) + 2*U
H33(13,13) =  2*E(1) + 2*E(2) + 2*E(3) + 3*U; H33(14,14) =  2*E(1) + 2*E(2) + 2*E(4) + 3*U
H33(15,15) =  2*E(1) + 2*E(3) + 2*E(4) + 3*U; H33(16,16) =  2*E(2) + 2*E(3) + 2*E(4) + 3*U

! make it symetric
do i=1,16
   do j=1,16
      if(i>j) then
         H33(i,j) = H33(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 48
allocate(WORK(48))

call ssyev('v','u',16,H33,16,W33,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H33 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,16
   Grand_potential(i+231) = W33(i) - mu*6  ! grand potentials of H33
end do

do i=1,16
   eigenvectors(i+231,232:247) = H33(1:16,i)    ! eigenvectors of H33
end do 

!--------------------H43 & H34 (same)--------------------------

H43(1,2) = -t; H43(1,3) = t; H43(1,4) = -t
H43(2,3) = -t; H43(2,4) = t
H43(3,4) = -t

H43(1,1) = 2*E(1) + 2*E(2) + 2*E(3) + E(4) + 3*U; H43(2,2) = 2*E(1) + 2*E(2) + E(3) + 2*E(4) + 3*U
H43(1,1) = 2*E(1) + E(2) + 2*E(3) + 2*E(4) + 3*U; H43(2,2) = E(1) + 2*E(2) + 2*E(3) + 2*E(4) + 3*U

! make it symetric
do i=1,4
   do j=1,4
      if(i>j) then
         H43(i,j) = H43(j,i)
      end if
   end do
end do

! solve the eigenvalues and eigenvectors
LWORK = 15
allocate(WORK(15))

call ssyev('v','u',4,H43,4,W43,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H43 matrix. Error code:', INFO
   stop
end if 

deallocate(WORK)

do i=1,4
   Grand_potential(i+247) = W43(i) - mu*7  ! grand potentials of H43
   Grand_potential(i+251) = W43(i) - mu*7  ! grand potentials of H34
end do

do i=1,4
   eigenvectors(i+247,248:251) = H43(1:4,i)    ! eigenvectors of H43
   eigenvectors(i+251,252:255) = H43(1:4,i)    ! eigenvectors of H34
end do 

!------------H44--------------------------

H44 =  2*E(2) + 2*E(2) + 2*E(3) + 2*E(4) + 4*U

Grand_potential(256) = H44 - mu*8     ! grand potentials of H44

eigenvectors(256,256) = 1              ! eigenvectors of H44

end subroutine hamiltonian

end module routines
