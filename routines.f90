module routines

! file before debugging has begun

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

H44 =  2*E(1) + 2*E(2) + 2*E(3) + 2*E(4) + 4*U

Grand_potential(256) = H44 - mu*8     ! grand potentials of H44

eigenvectors(256,256) = 1              ! eigenvectors of H44

end subroutine hamiltonian

subroutine transformations()

implicit none

integer :: i,j ! counters

PES_up = 0; PES_down = 0
IPES_up = 0; IPES_down = 0
phase_PES_up = 0; phase_PES_down = 0
phase_IPES_up = 0; phase_IPES_down = 0

PES_up(1,1) = 2; phase_PES_up(1,1) = 1
PES_up(1,3) = 10; phase_PES_up(1,3) = -1
PES_up(1,4) = 11; phase_PES_up(1,4) = -1
PES_up(1,5) = 12; phase_PES_up(1,5) = -1
PES_up(1,6) = 34; phase_PES_up(1,6) = -1
PES_up(1,7) = 22; phase_PES_up(1,7) = -1
PES_up(1,8) = 24; phase_PES_up(1,8) = -1
PES_up(1,9) = 26; phase_PES_up(1,9) = -1
PES_up(1,13) = 38; phase_PES_up(1,13) = 1
PES_up(1,14) = 39; phase_PES_up(1,14) = 1
PES_up(1,15) = 40; phase_PES_up(1,15) = 1
PES_up(1,16) = 82; phase_PES_up(1,16) = 1
PES_up(1,17) = 84; phase_PES_up(1,17) = 1
PES_up(1,18) = 86; phase_PES_up(1,18) = 1
PES_up(1,19) = 72; phase_PES_up(1,19) = 1
PES_up(1,20) = 75; phase_PES_up(1,20) = 1
PES_up(1,21) = 78; phase_PES_up(1,21) = 1
PES_up(1,23) = 58; phase_PES_up(1,23) = 1
PES_up(1,25) = 60; phase_PES_up(1,25) = 1
PES_up(1,27) = 62; phase_PES_up(1,27) = 1
PES_up(1,28) = 46; phase_PES_up(1,28) = 1
PES_up(1,29) = 47; phase_PES_up(1,29) = 1
PES_up(1,30) = 49; phase_PES_up(1,30) = 1
PES_up(1,31) = 50; phase_PES_up(1,31) = 1
PES_up(1,32) = 52; phase_PES_up(1,32) = 1
PES_up(1,33) = 53; phase_PES_up(1,33) = 1
PES_up(1,35) = 59; phase_PES_up(1,35) = 1
PES_up(1,36) = 61; phase_PES_up(1,36) = 1
PES_up(1,37) = 63; phase_PES_up(1,37) = 1
PES_up(1,41) = 94; phase_PES_up(1,41) = -1
PES_up(1,42) = 116; phase_PES_up(1,42) = -1
PES_up(1,43) = 119; phase_PES_up(1,43) = -1
PES_up(1,44) = 122; phase_PES_up(1,44) = -1
PES_up(1,45) = 115; phase_PES_up(1,45) = -1
PES_up(1,48) = 100; phase_PES_up(1,48) = -1
PES_up(1,51) = 103; phase_PES_up(1,51) = -1
PES_up(1,54) = 106; phase_PES_up(1,54) = -1
PES_up(1,55) = 96; phase_PES_up(1,55) = -1
PES_up(1,56) = 97; phase_PES_up(1,56) = -1
PES_up(1,57) = 98; phase_PES_up(1,57) = -1
PES_up(1,64) = 101; phase_PES_up(1,64) = -1
PES_up(1,65) = 102; phase_PES_up(1,65) = -1
PES_up(1,66) = 104; phase_PES_up(1,66) = -1
PES_up(1,67) = 105; phase_PES_up(1,67) = -1
PES_up(1,68) = 107; phase_PES_up(1,68) = -1
PES_up(1,69) = 108; phase_PES_up(1,69) = -1
PES_up(1,70) = 135; phase_PES_up(1,70) = -1
PES_up(1,71) = 134; phase_PES_up(1,71) = -1
PES_up(1,73) = 141; phase_PES_up(1,73) = -1
PES_up(1,74) = 140; phase_PES_up(1,74) = -1
PES_up(1,76) = 147; phase_PES_up(1,76) = -1
PES_up(1,77) = 146; phase_PES_up(1,77) = -1
PES_up(1,79) = 130; phase_PES_up(1,79) = -1
PES_up(1,80) = 129; phase_PES_up(1,80) = -1
PES_up(1,81) = 128; phase_PES_up(1,81) = -1
PES_up(1,83) = 158; phase_PES_up(1,83) = -1
PES_up(1,85) = 159; phase_PES_up(1,85) = -1
PES_up(1,87) = 160; phase_PES_up(1,87) = -1
PES_up(1,88) = 136; phase_PES_up(1,88) = -1
PES_up(1,89) = 138; phase_PES_up(1,89) = -1
PES_up(1,90) = 142; phase_PES_up(1,90) = -1
PES_up(1,91) = 144; phase_PES_up(1,91) = -1
PES_up(1,92) = 148; phase_PES_up(1,92) = -1
PES_up(1,93) = 150; phase_PES_up(1,93) = -1
PES_up(1,95) = 168; phase_PES_up(1,95) = 1
PES_up(1,99) = 164; phase_PES_up(1,99) = 1
PES_up(1,109) = 165; phase_PES_up(1,109) = 1
PES_up(1,110) = 166; phase_PES_up(1,110) = 1
PES_up(1,111) = 167; phase_PES_up(1,111) = 1
PES_up(1,112) = 196; phase_PES_up(1,112) = 1
PES_up(1,113) = 197; phase_PES_up(1,113) = 1
PES_up(1,114) = 198; phase_PES_up(1,114) = 1
PES_up(1,117) = 208; phase_PES_up(1,117) = 1
PES_up(1,118) = 209; phase_PES_up(1,118) = 1
PES_up(1,120) = 211; phase_PES_up(1,120) = 1
PES_up(1,121) = 212; phase_PES_up(1,121) = 1
PES_up(1,123) = 214; phase_PES_up(1,123) = 1
PES_up(1,124) = 215; phase_PES_up(1,124) = 1
PES_up(1,125) = 201; phase_PES_up(1,125) = 1
PES_up(1,126) = 204; phase_PES_up(1,126) = 1
PES_up(1,127) = 207; phase_PES_up(1,127) = 1
PES_up(1,131) = 172; phase_PES_up(1,131) = 1
PES_up(1,132) = 173; phase_PES_up(1,132) = 1
PES_up(1,133) = 174; phase_PES_up(1,133) = 1
PES_up(1,137) = 184; phase_PES_up(1,137) = 1
PES_up(1,139) = 185; phase_PES_up(1,139) = 1
PES_up(1,143) = 187; phase_PES_up(1,143) = 1
PES_up(1,145) = 188; phase_PES_up(1,145) = 1
PES_up(1,149) = 190; phase_PES_up(1,149) = 1
PES_up(1,151) = 191; phase_PES_up(1,151) = 1
PES_up(1,152) = 175; phase_PES_up(1,152) = 1
PES_up(1,153) = 176; phase_PES_up(1,153) = 1
PES_up(1,154) = 178; phase_PES_up(1,154) = 1
PES_up(1,155) = 179; phase_PES_up(1,155) = 1
PES_up(1,156) = 181; phase_PES_up(1,156) = 1
PES_up(1,157) = 182; phase_PES_up(1,157) = 1
PES_up(1,161) = 186; phase_PES_up(1,161) = 1
PES_up(1,162) = 189; phase_PES_up(1,162) = 1
PES_up(1,163) = 192; phase_PES_up(1,163) = 1
PES_up(1,169) = 226; phase_PES_up(1,169) = -1
PES_up(1,170) = 227; phase_PES_up(1,170) = -1
PES_up(1,171) = 228; phase_PES_up(1,171) = -1
PES_up(1,177) = 220; phase_PES_up(1,177) = -1
PES_up(1,180) = 221; phase_PES_up(1,180) = -1
PES_up(1,183) = 222; phase_PES_up(1,183) = -1
PES_up(1,193) = 223; phase_PES_up(1,193) = -1
PES_up(1,194) = 224; phase_PES_up(1,194) = -1
PES_up(1,195) = 225; phase_PES_up(1,195) = -1
PES_up(1,199) = 233; phase_PES_up(1,199) = -1
PES_up(1,200) = 232; phase_PES_up(1,200) = -1
PES_up(1,202) = 235; phase_PES_up(1,202) = -1
PES_up(1,203) = 234; phase_PES_up(1,203) = -1
PES_up(1,205) = 237; phase_PES_up(1,205) = -1
PES_up(1,206) = 236; phase_PES_up(1,206) = -1
PES_up(1,210) = 244; phase_PES_up(1,210) = -1
PES_up(1,213) = 245; phase_PES_up(1,213) = -1
PES_up(1,216) = 246; phase_PES_up(1,216) = -1
PES_up(1,217) = 238; phase_PES_up(1,217) = -1
PES_up(1,218) = 240; phase_PES_up(1,218) = -1
PES_up(1,219) = 242; phase_PES_up(1,219) = -1
PES_up(1,229) = 252; phase_PES_up(1,229) = 1
PES_up(1,230) = 253; phase_PES_up(1,230) = 1
PES_up(1,231) = 254; phase_PES_up(1,231) = 1
PES_up(1,239) = 248; phase_PES_up(1,239) = 1
PES_up(1,241) = 249; phase_PES_up(1,241) = 1
PES_up(1,243) = 250; phase_PES_up(1,243) = 1
PES_up(1,247) = 251; phase_PES_up(1,247) = 1
PES_up(1,255) = 256; phase_PES_up(1,255) = -1

PES_down(1,1) = 6; phase_PES_down(1,1) = 1
PES_down(1,2) = 34; phase_PES_down(1,2) = 1
PES_down(1,3) = 23; phase_PES_down(1,3) = -1
PES_down(1,4) = 25; phase_PES_down(1,4) = -1
PES_down(1,5) = 27; phase_PES_down(1,5) = -1
PES_down(1,7) = 16; phase_PES_down(1,7) = -1
PES_down(1,8) = 17; phase_PES_down(1,8) = -1
PES_down(1,9) = 18; phase_PES_down(1,9) = -1
PES_down(1,10) = 58; phase_PES_down(1,10) = -1
PES_down(1,11) = 60; phase_PES_down(1,11) = -1
PES_down(1,12) = 62; phase_PES_down(1,12) = -1
PES_down(1,13) = 48; phase_PES_down(1,13) = 1
PES_down(1,14) = 51; phase_PES_down(1,14) = 1
PES_down(1,15) = 54; phase_PES_down(1,15) = 1
PES_down(1,19) = 42; phase_PES_down(1,19) = 1
PES_down(1,20) = 43; phase_PES_down(1,20) = 1
PES_down(1,21) = 44; phase_PES_down(1,21) = 1
PES_down(1,22) = 82; phase_PES_down(1,22) = -1
PES_down(1,24) = 84; phase_PES_down(1,24) = -1
PES_down(1,26) = 86; phase_PES_down(1,26) = -1
PES_down(1,28) = 71; phase_PES_down(1,28) = 1
PES_down(1,29) = 70; phase_PES_down(1,29) = 1
PES_down(1,30) = 74; phase_PES_down(1,30) = 1
PES_down(1,31) = 73; phase_PES_down(1,31) = 1
PES_down(1,32) = 77; phase_PES_down(1,32) = 1
PES_down(1,33) = 76; phase_PES_down(1,33) = 1
PES_down(1,35) = 83; phase_PES_down(1,35) = 1
PES_down(1,36) = 85; phase_PES_down(1,36) = 1
PES_down(1,37) = 87; phase_PES_down(1,37) = 1
PES_down(1,38) = 100; phase_PES_down(1,38) = 1
PES_down(1,39) = 103; phase_PES_down(1,39) = 1
PES_down(1,40) = 106; phase_PES_down(1,40) = 1
PES_down(1,41) = 99; phase_PES_down(1,41) = -1
PES_down(1,45) = 95; phase_PES_down(1,45) = -1
PES_down(1,46) = 134; phase_PES_down(1,46) = 1
PES_down(1,47) = 135; phase_PES_down(1,47) = 1
PES_down(1,49) = 140; phase_PES_down(1,49) = 1
PES_down(1,50) = 141; phase_PES_down(1,50) = 1
PES_down(1,52) = 146; phase_PES_down(1,52) = 1
PES_down(1,53) = 147; phase_PES_down(1,53) = 1
PES_down(1,55) = 131; phase_PES_down(1,55) = -1
PES_down(1,56) = 132; phase_PES_down(1,56) = -1
PES_down(1,57) = 133; phase_PES_down(1,57) = -1
PES_down(1,59) = 158; phase_PES_down(1,59) = 1
PES_down(1,61) = 159; phase_PES_down(1,61) = 1
PES_down(1,63) = 160; phase_PES_down(1,63) = 1
PES_down(1,64) = 137; phase_PES_down(1,64) = -1
PES_down(1,65) = 139; phase_PES_down(1,65) = -1
PES_down(1,66) = 143; phase_PES_down(1,66) = -1
PES_down(1,67) = 145; phase_PES_down(1,67) = -1
PES_down(1,68) = 149; phase_PES_down(1,68) = -1
PES_down(1,69) = 151; phase_PES_down(1,69) = -1
PES_down(1,72) = 116; phase_PES_down(1,72) = 1
PES_down(1,75) = 119; phase_PES_down(1,75) = 1
PES_down(1,78) = 122; phase_PES_down(1,78) = 1
PES_down(1,79) = 112; phase_PES_down(1,79) = -1
PES_down(1,80) = 113; phase_PES_down(1,80) = -1
PES_down(1,81) = 114; phase_PES_down(1,81) = -1
PES_down(1,88) = 117; phase_PES_down(1,88) = -1
PES_down(1,89) = 118; phase_PES_down(1,89) = -1
PES_down(1,90) = 120; phase_PES_down(1,90) = -1
PES_down(1,91) = 121; phase_PES_down(1,91) = -1
PES_down(1,92) = 123; phase_PES_down(1,92) = -1
PES_down(1,93) = 124; phase_PES_down(1,93) = -1
PES_down(1,94) = 164; phase_PES_down(1,94) = -1
PES_down(1,96) = 172; phase_PES_down(1,96) = -1
PES_down(1,97) = 173; phase_PES_down(1,97) = -1
PES_down(1,98) = 174; phase_PES_down(1,98) = -1
PES_down(1,101) = 184; phase_PES_down(1,101) = -1
PES_down(1,102) = 185; phase_PES_down(1,102) = -1
PES_down(1,104) = 187; phase_PES_down(1,104) = -1
PES_down(1,105) = 188; phase_PES_down(1,105) = -1
PES_down(1,107) = 190; phase_PES_down(1,107) = -1
PES_down(1,108) = 191; phase_PES_down(1,108) = -1
PES_down(1,109) = 177; phase_PES_down(1,109) = 1
PES_down(1,110) = 180; phase_PES_down(1,110) = 1
PES_down(1,111) = 183; phase_PES_down(1,111) = 1
PES_down(1,115) = 168; phase_PES_down(1,115) = -1
PES_down(1,125) = 169; phase_PES_down(1,125) = 1
PES_down(1,126) = 170; phase_PES_down(1,126) = 1
PES_down(1,127) = 171; phase_PES_down(1,127) = 1
PES_down(1,128) = 198; phase_PES_down(1,128) = -1
PES_down(1,129) = 197; phase_PES_down(1,129) = -1
PES_down(1,130) = 196; phase_PES_down(1,130) = -1
PES_down(1,136) = 208; phase_PES_down(1,136) = -1
PES_down(1,138) = 209; phase_PES_down(1,138) = -1
PES_down(1,142) = 211; phase_PES_down(1,142) = -1
PES_down(1,144) = 212; phase_PES_down(1,144) = -1
PES_down(1,148) = 214; phase_PES_down(1,148) = -1
PES_down(1,150) = 215; phase_PES_down(1,150) = -1
PES_down(1,152) = 200; phase_PES_down(1,152) = 1
PES_down(1,153) = 199; phase_PES_down(1,153) = 1
PES_down(1,154) = 203; phase_PES_down(1,154) = 1
PES_down(1,155) = 202; phase_PES_down(1,155) = 1
PES_down(1,156) = 206; phase_PES_down(1,156) = 1
PES_down(1,157) = 205; phase_PES_down(1,157) = 1
PES_down(1,161) = 210; phase_PES_down(1,161) = 1
PES_down(1,162) = 213; phase_PES_down(1,162) = 1
PES_down(1,163) = 216; phase_PES_down(1,163) = 1
PES_down(1,165) = 220; phase_PES_down(1,165) = 1
PES_down(1,166) = 221; phase_PES_down(1,166) = 1
PES_down(1,167) = 222; phase_PES_down(1,167) = 1
PES_down(1,175) = 232; phase_PES_down(1,175) = 1
PES_down(1,176) = 233; phase_PES_down(1,176) = 1
PES_down(1,178) = 234; phase_PES_down(1,178) = 1
PES_down(1,179) = 235; phase_PES_down(1,179) = 1
PES_down(1,181) = 236; phase_PES_down(1,181) = 1
PES_down(1,182) = 237; phase_PES_down(1,182) = 1
PES_down(1,186) = 244; phase_PES_down(1,186) = 1
PES_down(1,189) = 245; phase_PES_down(1,189) = 1
PES_down(1,192) = 246; phase_PES_down(1,192) = 1
PES_down(1,193) = 239; phase_PES_down(1,193) = -1
PES_down(1,194) = 241; phase_PES_down(1,194) = -1
PES_down(1,195) = 243; phase_PES_down(1,195) = -1
PES_down(1,201) = 226; phase_PES_down(1,201) = 1
PES_down(1,204) = 227; phase_PES_down(1,204) = 1
PES_down(1,207) = 228; phase_PES_down(1,207) = 1
PES_down(1,217) = 229; phase_PES_down(1,217) = -1
PES_down(1,218) = 230; phase_PES_down(1,218) = -1
PES_down(1,219) = 231; phase_PES_down(1,219) = -1
PES_down(1,223) = 248; phase_PES_down(1,223) = -1
PES_down(1,224) = 249; phase_PES_down(1,224) = -1
PES_down(1,225) = 250; phase_PES_down(1,225) = -1
PES_down(1,238) = 252; phase_PES_down(1,238) = -1
PES_down(1,240) = 253; phase_PES_down(1,240) = -1
PES_down(1,242) = 254; phase_PES_down(1,242) = -1
PES_down(1,247) = 255; phase_PES_down(1,247) = 1
PES_down(1,251) = 256; phase_PES_down(1,251) = 1

PES_up(2,1) = 3; phase_PES_up(2,1) = 1
PES_up(2,2) = 10; phase_PES_up(2,2) = 1
PES_up(2,4) = 13; phase_PES_up(2,4) = -1
PES_up(2,5) = 14; phase_PES_up(2,5) = -1
PES_up(2,6) = 23; phase_PES_up(2,6) = 1
PES_up(2,7) = 35; phase_PES_up(2,7) = -1
PES_up(2,8) = 28; phase_PES_up(2,8) = -1
PES_up(2,9) = 30; phase_PES_up(2,9) = -1
PES_up(2,11) = 38; phase_PES_up(2,11) = -1
PES_up(2,12) = 39; phase_PES_up(2,12) = -1
PES_up(2,15) = 41; phase_PES_up(2,15) = 1
PES_up(2,16) = 83; phase_PES_up(2,16) = -1
PES_up(2,17) = 71; phase_PES_up(2,17) = -1
PES_up(2,18) = 74; phase_PES_up(2,18) = -1
PES_up(2,19) = 88; phase_PES_up(2,19) = 1
PES_up(2,20) = 90; phase_PES_up(2,20) = 1
PES_up(2,21) = 81; phase_PES_up(2,21) = 1
PES_up(2,22) = 59; phase_PES_up(2,22) = -1
PES_up(2,24) = 46; phase_PES_up(2,24) = -1
PES_up(2,25) = 48; phase_PES_up(2,25) = -1
PES_up(2,26) = 49; phase_PES_up(2,26) = -1
PES_up(2,27) = 51; phase_PES_up(2,27) = -1
PES_up(2,29) = 64; phase_PES_up(2,29) = 1
PES_up(2,31) = 66; phase_PES_up(2,31) = 1
PES_up(2,32) = 55; phase_PES_up(2,32) = 1
PES_up(2,33) = 56; phase_PES_up(2,33) = 1
PES_up(2,34) = 58; phase_PES_up(2,34) = 1
PES_up(2,36) = 65; phase_PES_up(2,36) = 1
PES_up(2,37) = 67; phase_PES_up(2,37) = 1
PES_up(2,40) = 94; phase_PES_up(2,40) = 1
PES_up(2,42) = 117; phase_PES_up(2,42) = 1
PES_up(2,43) = 120; phase_PES_up(2,43) = 1
PES_up(2,44) = 114; phase_PES_up(2,44) = 1
PES_up(2,45) = 125; phase_PES_up(2,45) = -1
PES_up(2,47) = 101; phase_PES_up(2,47) = 1
PES_up(2,50) = 104; phase_PES_up(2,50) = 1
PES_up(2,52) = 96; phase_PES_up(2,52) = 1
PES_up(2,53) = 97; phase_PES_up(2,53) = 1
PES_up(2,54) = 99; phase_PES_up(2,54) = 1
PES_up(2,57) = 106; phase_PES_up(2,57) = -1
PES_up(2,60) = 100; phase_PES_up(2,60) = -1
PES_up(2,61) = 102; phase_PES_up(2,61) = 1
PES_up(2,62) = 103; phase_PES_up(2,62) = -1
PES_up(2,63) = 105; phase_PES_up(2,63) = 1
PES_up(2,68) = 110; phase_PES_up(2,68) = -1
PES_up(2,69) = 111; phase_PES_up(2,69) = -1
PES_up(2,70) = 137; phase_PES_up(2,70) = 1
PES_up(2,72) = 136; phase_PES_up(2,72) = 1
PES_up(2,73) = 143; phase_PES_up(2,73) = 1
PES_up(2,75) = 142; phase_PES_up(2,75) = 1
PES_up(2,76) = 132; phase_PES_up(2,76) = 1
PES_up(2,77) = 131; phase_PES_up(2,77) = 1
PES_up(2,78) = 128; phase_PES_up(2,78) = 1
PES_up(2,79) = 153; phase_PES_up(2,79) = -1
PES_up(2,80) = 152; phase_PES_up(2,80) = -1
PES_up(2,82) = 158; phase_PES_up(2,82) = -1
PES_up(2,84) = 134; phase_PES_up(2,84) = -1
PES_up(2,85) = 139; phase_PES_up(2,85) = 1
PES_up(2,86) = 140; phase_PES_up(2,86) = -1
PES_up(2,87) = 145; phase_PES_up(2,87) = 1
PES_up(2,89) = 161; phase_PES_up(2,89) = -1
PES_up(2,91) = 162; phase_PES_up(2,91) = -1
PES_up(2,92) = 154; phase_PES_up(2,92) = -1
PES_up(2,93) = 156; phase_PES_up(2,93) = -1
PES_up(2,95) = 169; phase_PES_up(2,95) = -1
PES_up(2,98) = 165; phase_PES_up(2,98) = -1
PES_up(2,106) = 164; phase_PES_up(2,106) = 1
PES_up(2,107) = 166; phase_PES_up(2,107) = -1
PES_up(2,108) = 167; phase_PES_up(2,108) = -1
PES_up(2,112) = 199; phase_PES_up(2,112) = -1
PES_up(2,113) = 200; phase_PES_up(2,113) = -1
PES_up(2,115) = 201; phase_PES_up(2,115) = -1
PES_up(2,116) = 208; phase_PES_up(2,116) = 1
PES_up(2,118) = 210; phase_PES_up(2,118) = -1
PES_up(2,119) = 211; phase_PES_up(2,119) = 1
PES_up(2,121) = 213; phase_PES_up(2,121) = -1
PES_up(2,122) = 198; phase_PES_up(2,122) = 1
PES_up(2,123) = 203; phase_PES_up(2,123) = -1
PES_up(2,124) = 206; phase_PES_up(2,124) = -1
PES_up(2,126) = 217; phase_PES_up(2,126) = 1
PES_up(2,127) = 218; phase_PES_up(2,127) = 1
PES_up(2,129) = 175; phase_PES_up(2,129) = -1
PES_up(2,130) = 176; phase_PES_up(2,130) = -1
PES_up(2,133) = 177; phase_PES_up(2,133) = -1
PES_up(2,135) = 184; phase_PES_up(2,135) = 1
PES_up(2,138) = 186; phase_PES_up(2,138) = -1
PES_up(2,141) = 187; phase_PES_up(2,141) = 1
PES_up(2,144) = 189; phase_PES_up(2,144) = -1
PES_up(2,146) = 172; phase_PES_up(2,146) = 1
PES_up(2,147) = 173; phase_PES_up(2,147) = 1
PES_up(2,148) = 178; phase_PES_up(2,148) = -1
PES_up(2,149) = 180; phase_PES_up(2,149) = -1
PES_up(2,150) = 181; phase_PES_up(2,150) = -1
PES_up(2,151) = 183; phase_PES_up(2,151) = -1
PES_up(2,155) = 193; phase_PES_up(2,155) = 1
PES_up(2,157) = 194; phase_PES_up(2,157) = 1
PES_up(2,159) = 185; phase_PES_up(2,159) = 1
PES_up(2,160) = 188; phase_PES_up(2,160) = 1
PES_up(2,163) = 195; phase_PES_up(2,163) = 1
PES_up(2,168) = 226; phase_PES_up(2,168) = -1
PES_up(2,170) = 229; phase_PES_up(2,170) = 1
PES_up(2,171) = 230; phase_PES_up(2,171) = 1
PES_up(2,174) = 220; phase_PES_up(2,174) = -1
PES_up(2,179) = 223; phase_PES_up(2,179) = 1
PES_up(2,182) = 224; phase_PES_up(2,182) = 1
PES_up(2,190) = 221; phase_PES_up(2,190) = -1
PES_up(2,191) = 222; phase_PES_up(2,191) = -1
PES_up(2,192) = 225; phase_PES_up(2,192) = 1
PES_up(2,196) = 233; phase_PES_up(2,196) = -1
PES_up(2,197) = 232; phase_PES_up(2,197) = -1
PES_up(2,202) = 239; phase_PES_up(2,202) = 1
PES_up(2,204) = 238; phase_PES_up(2,204) = 1
PES_up(2,205) = 241; phase_PES_up(2,205) = 1
PES_up(2,207) = 240; phase_PES_up(2,207) = 1
PES_up(2,209) = 244; phase_PES_up(2,209) = 1
PES_up(2,212) = 245; phase_PES_up(2,212) = -1
PES_up(2,214) = 234; phase_PES_up(2,214) = -1
PES_up(2,215) = 236; phase_PES_up(2,215) = -1
PES_up(2,216) = 243; phase_PES_up(2,216) = 1
PES_up(2,219) = 247; phase_PES_up(2,219) = -1
PES_up(2,227) = 252; phase_PES_up(2,227) = 1
PES_up(2,228) = 253; phase_PES_up(2,228) = 1
PES_up(2,231) = 255; phase_PES_up(2,231) = -1
PES_up(2,235) = 248; phase_PES_up(2,235) = 1
PES_up(2,237) = 249; phase_PES_up(2,237) = 1
PES_up(2,242) = 251; phase_PES_up(2,242) = -1
PES_up(2,246) = 250; phase_PES_up(2,246) = 1
PES_up(2,254) = 256; phase_PES_up(2,254) = -1

PES_down(2,1) = 7; phase_PES_down(2,1) = 1
PES_down(2,2) = 22; phase_PES_down(2,2) = 1
PES_down(2,3) = 35; phase_PES_down(2,3) = 1
PES_down(2,4) = 29; phase_PES_down(2,4) = -1
PES_down(2,5) = 31; phase_PES_down(2,5) = -1
PES_down(2,6) = 16; phase_PES_down(2,6) = 1
PES_down(2,8) = 19; phase_PES_down(2,8) = -1
PES_down(2,9) = 20; phase_PES_down(2,9) = -1
PES_down(2,10) = 59; phase_PES_down(2,10) = 1
PES_down(2,11) = 47; phase_PES_down(2,11) = -1
PES_down(2,12) = 50; phase_PES_down(2,12) = -1
PES_down(2,13) = 64; phase_PES_down(2,13) = -1
PES_down(2,14) = 66; phase_PES_down(2,14) = -1
PES_down(2,15) = 57; phase_PES_down(2,15) = 1
PES_down(2,17) = 42; phase_PES_down(2,17) = -1
PES_down(2,18) = 43; phase_PES_down(2,18) = -1
PES_down(2,21) = 45; phase_PES_down(2,21) = 1
PES_down(2,23) = 83; phase_PES_down(2,23) = 1
PES_down(2,24) = 72; phase_PES_down(2,24) = -1
PES_down(2,25) = 70; phase_PES_down(2,25) = -1
PES_down(2,26) = 75; phase_PES_down(2,26) = -1
PES_down(2,27) = 73; phase_PES_down(2,27) = -1
PES_down(2,28) = 88; phase_PES_down(2,28) = -1
PES_down(2,30) = 90; phase_PES_down(2,30) = -1
PES_down(2,32) = 80; phase_PES_down(2,32) = 1
PES_down(2,33) = 79; phase_PES_down(2,33) = 1
PES_down(2,34) = 82; phase_PES_down(2,34) = 1
PES_down(2,36) = 89; phase_PES_down(2,36) = 1
PES_down(2,37) = 91; phase_PES_down(2,37) = 1
PES_down(2,38) = 101; phase_PES_down(2,38) = -1
PES_down(2,39) = 104; phase_PES_down(2,39) = -1
PES_down(2,40) = 98; phase_PES_down(2,40) = 1
PES_down(2,41) = 109; phase_PES_down(2,41) = 1
PES_down(2,44) = 95; phase_PES_down(2,44) = 1
PES_down(2,46) = 136; phase_PES_down(2,46) = -1
PES_down(2,48) = 137; phase_PES_down(2,48) = -1
PES_down(2,49) = 142; phase_PES_down(2,49) = -1
PES_down(2,51) = 143; phase_PES_down(2,51) = -1
PES_down(2,52) = 120; phase_PES_down(2,52) = 1
PES_down(2,53) = 130; phase_PES_down(2,53) = 1
PES_down(2,54) = 133; phase_PES_down(2,54) = 1
PES_down(2,55) = 152; phase_PES_down(2,55) = 1
PES_down(2,56) = 153; phase_PES_down(2,56) = 1
PES_down(2,58) = 158; phase_PES_down(2,58) = 1
PES_down(2,60) = 135; phase_PES_down(2,60) = -1
PES_down(2,61) = 138; phase_PES_down(2,61) = 1
PES_down(2,62) = 141; phase_PES_down(2,62) = -1
PES_down(2,63) = 144; phase_PES_down(2,63) = 1
PES_down(2,65) = 161; phase_PES_down(2,65) = 1
PES_down(2,67) = 162; phase_PES_down(2,67) = 1
PES_down(2,68) = 155; phase_PES_down(2,68) = -1
PES_down(2,69) = 157; phase_PES_down(2,69) = -1
PES_down(2,71) = 117; phase_PES_down(2,71) = -1
PES_down(2,74) = 120; phase_PES_down(2,74) = -1
PES_down(2,76) = 112; phase_PES_down(2,76) = 1
PES_down(2,77) = 113; phase_PES_down(2,77) = 1
PES_down(2,78) = 115; phase_PES_down(2,78) = 1
PES_down(2,81) = 125; phase_PES_down(2,81) = 1
PES_down(2,84) = 116; phase_PES_down(2,84) = -1
PES_down(2,85) = 118; phase_PES_down(2,85) = 1
PES_down(2,86) = 119; phase_PES_down(2,86) = -1
PES_down(2,87) = 121; phase_PES_down(2,87) = 1
PES_down(2,92) = 126; phase_PES_down(2,92) = -1
PES_down(2,93) = 127; phase_PES_down(2,93) = -1
PES_down(2,94) = 165; phase_PES_down(2,94) = 1
PES_down(2,96) = 175; phase_PES_down(2,96) = 1
PES_down(2,97) = 176; phase_PES_down(2,97) = 1
PES_down(2,99) = 177; phase_PES_down(2,99) = 1
PES_down(2,100) = 184; phase_PES_down(2,100) = -1
PES_down(2,102) = 186; phase_PES_down(2,102) = 1
PES_down(2,103) = 187; phase_PES_down(2,103) = -1
PES_down(2,105) = 189; phase_PES_down(2,105) = 1
PES_down(2,106) = 174; phase_PES_down(2,106) = 1
PES_down(2,107) = 179; phase_PES_down(2,107) = -1
PES_down(2,108) = 182; phase_PES_down(2,108) = -1
PES_down(2,110) = 193; phase_PES_down(2,110) = -1
PES_down(2,111) = 194; phase_PES_down(2,111) = -1
PES_down(2,114) = 169; phase_PES_down(2,114) = 1
PES_down(2,122) = 168; phase_PES_down(2,122) = 1
PES_down(2,123) = 170; phase_PES_down(2,123) = -1
PES_down(2,124) = 171; phase_PES_down(2,124) = -1
PES_down(2,128) = 201; phase_PES_down(2,128) = 1
PES_down(2,131) = 200; phase_PES_down(2,131) = 1
PES_down(2,132) = 199; phase_PES_down(2,132) = 1
PES_down(2,134) = 208; phase_PES_down(2,134) = -1
PES_down(2,139) = 210; phase_PES_down(2,139) = 1
PES_down(2,140) = 211; phase_PES_down(2,140) = -1
PES_down(2,145) = 213; phase_PES_down(2,145) = 1
PES_down(2,146) = 197; phase_PES_down(2,146) = 1
PES_down(2,147) = 196; phase_PES_down(2,147) = 1
PES_down(2,148) = 204; phase_PES_down(2,148) = -1
PES_down(2,149) = 202; phase_PES_down(2,149) = -1
PES_down(2,150) = 207; phase_PES_down(2,150) = -1
PES_down(2,151) = 200; phase_PES_down(2,151) = -1
PES_down(2,154) = 217; phase_PES_down(2,154) = -1
PES_down(2,156) = 218; phase_PES_down(2,156) = -1
PES_down(2,159) = 209; phase_PES_down(2,159) = 1
PES_down(2,160) = 212; phase_PES_down(2,160) = 1
PES_down(2,163) = 219; phase_PES_down(2,163) = 1
PES_down(2,164) = 220; phase_PES_down(2,164) = 1
PES_down(2,166) = 223; phase_PES_down(2,166) = -1
PES_down(2,167) = 224; phase_PES_down(2,167) = -1
PES_down(2,172) = 232; phase_PES_down(2,172) = 1
PES_down(2,173) = 233; phase_PES_down(2,173) = 1
PES_down(2,178) = 238; phase_PES_down(2,178) = -1
PES_down(2,180) = 239; phase_PES_down(2,180) = -1
PES_down(2,181) = 240; phase_PES_down(2,181) = -1
PES_down(2,183) = 241; phase_PES_down(2,183) = -1
PES_down(2,185) = 244; phase_PES_down(2,185) = 1
PES_down(2,188) = 245; phase_PES_down(2,188) = 1
PES_down(2,190) = 235; phase_PES_down(2,190) = -1
PES_down(2,191) = 237; phase_PES_down(2,191) = -1
PES_down(2,192) = 242; phase_PES_down(2,192) = 1
PES_down(2,195) = 247; phase_PES_down(2,195) = 1
PES_down(2,198) = 226; phase_PES_down(2,198) = 1
PES_down(2,203) = 229; phase_PES_down(2,203) = -1
PES_down(2,206) = 230; phase_PES_down(2,206) = -1
PES_down(2,214) = 227; phase_PES_down(2,214) = -1
PES_down(2,215) = 228; phase_PES_down(2,215) = -1
PES_down(2,216) = 231; phase_PES_down(2,216) = 1
PES_down(2,221) = 248; phase_PES_down(2,221) = -1
PES_down(2,222) = 249; phase_PES_down(2,222) = -1
PES_down(2,225) = 251; phase_PES_down(2,225) = 1
PES_down(2,234) = 252; phase_PES_down(2,234) = -1
PES_down(2,236) = 253; phase_PES_down(2,236) = -1
PES_down(2,243) = 255; phase_PES_down(2,243) = 1
PES_down(2,246) = 254; phase_PES_down(2,246) = 1
PES_down(2,250) = 256; phase_PES_down(2,250) = 1

PES_up(3,1) = 4; phase_PES_up(3,1) = 1
PES_up(3,2) = 11; phase_PES_up(3,2) = 1
PES_up(3,3) = 13; phase_PES_up(3,3) = 1
PES_up(3,5) = 15; phase_PES_up(3,5) = -1
PES_up(3,6) = 25; phase_PES_up(3,6) = 1
PES_up(3,7) = 29; phase_PES_up(3,7) = 1
PES_up(3,8) = 36; phase_PES_up(3,8) = -1
PES_up(3,9) = 32; phase_PES_up(3,9) = -1
PES_up(3,10) = 38; phase_PES_up(3,10) = 1
PES_up(3,12) = 40; phase_PES_up(3,12) = -1
PES_up(3,14) = 41; phase_PES_up(3,14) = -1
PES_up(3,16) = 70; phase_PES_up(3,16) = 1
PES_up(3,17) = 85; phase_PES_up(3,17) = -1
PES_up(3,18) = 77; phase_PES_up(3,18) = -1
PES_up(3,19) = 89; phase_PES_up(3,19) = -1
PES_up(3,20) = 80; phase_PES_up(3,20) = -1
PES_up(3,21) = 92; phase_PES_up(3,21) = 1
PES_up(3,22) = 47; phase_PES_up(3,22) = 1
PES_up(3,23) = 48; phase_PES_up(3,23) = 1
PES_up(3,24) = 61; phase_PES_up(3,24) = -1
PES_up(3,26) = 52; phase_PES_up(3,26) = -1
PES_up(3,27) = 54; phase_PES_up(3,27) = -1
PES_up(3,28) = 65; phase_PES_up(3,28) = -1
PES_up(3,30) = 55; phase_PES_up(3,30) = -1
PES_up(3,31) = 57; phase_PES_up(3,31) = -1
PES_up(3,33) = 68; phase_PES_up(3,33) = 1
PES_up(3,34) = 60; phase_PES_up(3,34) = 1
PES_up(3,35) = 64; phase_PES_up(3,35) = 1
PES_up(3,37) = 69; phase_PES_up(3,37) = 1
PES_up(3,39) = 94; phase_PES_up(3,39) = -1
PES_up(3,42) = 118; phase_PES_up(3,42) = -1
PES_up(3,43) = 113; phase_PES_up(3,43) = -1
PES_up(3,44) = 123; phase_PES_up(3,44) = 1
PES_up(3,45) = 126; phase_PES_up(3,45) = 1
PES_up(3,46) = 102; phase_PES_up(3,46) = -1
PES_up(3,49) = 96; phase_PES_up(3,49) = -1
PES_up(3,50) = 98; phase_PES_up(3,50) = -1
PES_up(3,51) = 99; phase_PES_up(3,51) = -1
PES_up(3,53) = 107; phase_PES_up(3,53) = 1
PES_up(3,56) = 110; phase_PES_up(3,56) = 1
PES_up(3,58) = 100; phase_PES_up(3,58) = 1
PES_up(3,59) = 101; phase_PES_up(3,59) = 1
PES_up(3,62) = 106; phase_PES_up(3,62) = -1
PES_up(3,63) = 108; phase_PES_up(3,63) = 1
PES_up(3,66) = 109; phase_PES_up(3,66) = -1
PES_up(3,67) = 111; phase_PES_up(3,67) = 1
PES_up(3,71) = 139; phase_PES_up(3,71) = -1
PES_up(3,72) = 138; phase_PES_up(3,72) = -1
PES_up(3,73) = 133; phase_PES_up(3,73) = -1
PES_up(3,74) = 131; phase_PES_up(3,74) = -1
PES_up(3,75) = 129; phase_PES_up(3,75) = -1
PES_up(3,76) = 149; phase_PES_up(3,76) = 1
PES_up(3,78) = 148; phase_PES_up(3,78) = 1
PES_up(3,79) = 155; phase_PES_up(3,79) = 1
PES_up(3,81) = 154; phase_PES_up(3,81) = 1
PES_up(3,82) = 135; phase_PES_up(3,82) = 1
PES_up(3,83) = 137; phase_PES_up(3,83) = 1
PES_up(3,84) = 159; phase_PES_up(3,84) = -1
PES_up(3,86) = 146; phase_PES_up(3,86) = -1
PES_up(3,87) = 151; phase_PES_up(3,87) = 1
PES_up(3,88) = 161; phase_PES_up(3,88) = -1
PES_up(3,90) = 152; phase_PES_up(3,90) = -1
PES_up(3,91) = 157; phase_PES_up(3,91) = 1
PES_up(3,93) = 163; phase_PES_up(3,93) = -1
PES_up(3,95) = 170; phase_PES_up(3,95) = 1
PES_up(3,97) = 166; phase_PES_up(3,97) = 1
PES_up(3,103) = 164; phase_PES_up(3,103) = -1
PES_up(3,104) = 165; phase_PES_up(3,104) = -1
PES_up(3,105) = 167; phase_PES_up(3,105) = 1
PES_up(3,112) = 202; phase_PES_up(3,112) = 1
PES_up(3,114) = 203; phase_PES_up(3,114) = 1
PES_up(3,115) = 204; phase_PES_up(3,115) = 1
PES_up(3,116) = 209; phase_PES_up(3,116) = -1
PES_up(3,117) = 210; phase_PES_up(3,117) = -1
PES_up(3,119) = 197; phase_PES_up(3,119) = -1
PES_up(3,120) = 200; phase_PES_up(3,120) = -1
PES_up(3,121) = 205; phase_PES_up(3,121) = 1 
PES_up(3,122) = 214; phase_PES_up(3,122) = 1
PES_up(3,124) = 216; phase_PES_up(3,124) = -1
PES_up(3,125) = 217; phase_PES_up(3,125) = 1
PES_up(3,127) = 219; phase_PES_up(3,127) = -1
PES_up(3,128) = 178; phase_PES_up(3,128) = 1
PES_up(3,130) = 179; phase_PES_up(3,130) = 1
PES_up(3,132) = 180; phase_PES_up(3,132) = 1
PES_up(3,134) = 185; phase_PES_up(3,134) = -1
PES_up(3,136) = 186; phase_PES_up(3,136) = -1
PES_up(3,140) = 172; phase_PES_up(3,140) = -1
PES_up(3,141) = 174; phase_PES_up(3,141) = -1
PES_up(3,142) = 175; phase_PES_up(3,142) = -1
PES_up(3,143) = 177; phase_PES_up(3,143) = -1
PES_up(3,144) = 182; phase_PES_up(3,144) = 1
PES_up(3,145) = 183; phase_PES_up(3,145) = 1
PES_up(3,147) = 190; phase_PES_up(3,147) = 1
PES_up(3,150) = 192; phase_PES_up(3,150) = -1
PES_up(3,153) = 193; phase_PES_up(3,153) = 1
PES_up(3,156) = 195; phase_PES_up(3,156) = -1
PES_up(3,158) = 184; phase_PES_up(3,158) = 1
PES_up(3,160) = 191; phase_PES_up(3,160) = 1
PES_up(3,162) = 194; phase_PES_up(3,162) = 1
PES_up(3,168) = 227; phase_PES_up(3,168) = 1
PES_up(3,169) = 229; phase_PES_up(3,169) = 1
PES_up(3,171) = 231; phase_PES_up(3,171) = -1
PES_up(3,173) = 221; phase_PES_up(3,173) = 1
PES_up(3,176) = 223; phase_PES_up(3,176) = 1
PES_up(3,181) = 225; phase_PES_up(3,181) = -1
PES_up(3,187) = 220; phase_PES_up(3,187) = -1
PES_up(3,188) = 222; phase_PES_up(3,188) = 1
PES_up(3,189) = 224; phase_PES_up(3,189) = 1
PES_up(3,196) = 235; phase_PES_up(3,196) = 1
PES_up(3,198) = 234; phase_PES_up(3,198) = 1
PES_up(3,199) = 239; phase_PES_up(3,199) = 1
PES_up(3,201) = 238; phase_PES_up(3,201) = 1
PES_up(3,206) = 243; phase_PES_up(3,206) = -1
PES_up(3,207) = 242; phase_PES_up(3,207) = -1
PES_up(3,208) = 244; phase_PES_up(3,208) = -1
PES_up(3,211) = 232; phase_PES_up(3,211) = -1
PES_up(3,212) = 237; phase_PES_up(3,212) = 1
PES_up(3,213) = 241; phase_PES_up(3,213) = 1
PES_up(3,215) = 246; phase_PES_up(3,215) = -1
PES_up(3,218) = 247; phase_PES_up(3,218) = -1
PES_up(3,226) = 252; phase_PES_up(3,226) = 1
PES_up(3,228) = 254; phase_PES_up(3,228) = -1
PES_up(3,230) = 255; phase_PES_up(3,230) = -1
PES_up(3,233) = 248; phase_PES_up(3,233) = 1
PES_up(3,236) = 250; phase_PES_up(3,236) = -1
PES_up(3,240) = 251; phase_PES_up(3,240) = -1
PES_up(3,245) = 249; phase_PES_up(3,245) = 1
PES_up(3,253) = 256; phase_PES_up(3,253) = -1

PES_down(3,1) = 8; phase_PES_down(3,1) = 1
PES_down(3,2) = 24; phase_PES_down(3,2) = 1
PES_down(3,3) = 28; phase_PES_down(3,3) = 1
PES_down(3,4) = 36; phase_PES_down(3,4) = 1
PES_down(3,5) = 33; phase_PES_down(3,5) = -1
PES_down(3,6) = 17; phase_PES_down(3,6) = 1
PES_down(3,7) = 19; phase_PES_down(3,7) = 1
PES_down(3,9) = 21; phase_PES_down(3,9) = -1
PES_down(3,10) = 46; phase_PES_down(3,10) = 1
PES_down(3,11) = 61; phase_PES_down(3,11) = 1
PES_down(3,12) = 53; phase_PES_down(3,12) = -1
PES_down(3,13) = 65; phase_PES_down(3,13) = 1
PES_down(3,14) = 56; phase_PES_down(3,14) = -1
PES_down(3,15) = 68; phase_PES_down(3,15) = -1
PES_down(3,16) = 42; phase_PES_down(3,16) = 1
PES_down(3,18) = 44; phase_PES_down(3,18) = -1
PES_down(3,20) = 45; phase_PES_down(3,20) = -1
PES_down(3,22) = 72; phase_PES_down(3,22) = -1
PES_down(3,23) = 71; phase_PES_down(3,23) = 1
PES_down(3,25) = 85; phase_PES_down(3,25) = 1
PES_down(3,26) = 78; phase_PES_down(3,26) = -1
PES_down(3,27) = 76; phase_PES_down(3,27) = -1
PES_down(3,29) = 89; phase_PES_down(3,29) = 1
PES_down(3,30) = 81; phase_PES_down(3,30) = -1
PES_down(3,31) = 79; phase_PES_down(3,31) = -1
PES_down(3,32) = 92; phase_PES_down(3,32) = -1
PES_down(3,34) = 84; phase_PES_down(3,34) = 1
PES_down(3,35) = 88; phase_PES_down(3,35) = 1
PES_down(3,37) = 93; phase_PES_down(3,37) = 1
PES_down(3,38) = 102; phase_PES_down(3,38) = 1
PES_down(3,39) = 97; phase_PES_down(3,39) = -1
PES_down(3,40) = 107; phase_PES_down(3,40) = -1
PES_down(3,41) = 110; phase_PES_down(3,41) = -1
PES_down(3,43) = 95; phase_PES_down(3,43) = -1
PES_down(3,47) = 138; phase_PES_down(3,47) = 1
PES_down(3,48) = 139; phase_PES_down(3,48) = 1
PES_down(3,49) = 128; phase_PES_down(3,49) = -1
PES_down(3,50) = 130; phase_PES_down(3,50) = -1
PES_down(3,51) = 132; phase_PES_down(3,51) = -1
PES_down(3,52) = 148; phase_PES_down(3,52) = -1
PES_down(3,54) = 149; phase_PES_down(3,54) = -1
PES_down(3,55) = 154; phase_PES_down(3,55) = -1
PES_down(3,57) = 155; phase_PES_down(3,57) = -1
PES_down(3,58) = 134; phase_PES_down(3,58) = 1
PES_down(3,59) = 136; phase_PES_down(3,59) = 1
PES_down(3,60) = 159; phase_PES_down(3,60) = 1
PES_down(3,62) = 147; phase_PES_down(3,62) = -1
PES_down(3,63) = 150; phase_PES_down(3,63) = 1
PES_down(3,64) = 161; phase_PES_down(3,64) = 1
PES_down(3,66) = 153; phase_PES_down(3,66) = -1
PES_down(3,67) = 156; phase_PES_down(3,67) = 1
PES_down(3,69) = 163; phase_PES_down(3,69) = 1
PES_down(3,70) = 118; phase_PES_down(3,70) = 1
PES_down(3,73) = 112; phase_PES_down(3,73) = -1
PES_down(3,74) = 114; phase_PES_down(3,74) = -1
PES_down(3,75) = 115; phase_PES_down(3,75) = -1
PES_down(3,77) = 123; phase_PES_down(3,77) = -1
PES_down(3,80) = 126; phase_PES_down(3,80) = -1
PES_down(3,82) = 116; phase_PES_down(3,82) = 1
PES_down(3,83) = 117; phase_PES_down(3,83) = 1
PES_down(3,86) = 122; phase_PES_down(3,86) = -1
PES_down(3,87) = 124; phase_PES_down(3,87) = 1
PES_down(3,90) = 125; phase_PES_down(3,90) = -1
PES_down(3,91) = 127; phase_PES_down(3,91) = 1
PES_down(3,94) = 166; phase_PES_down(3,94) = -1
PES_down(3,96) = 178; phase_PES_down(3,96) = -1
PES_down(3,98) = 179; phase_PES_down(3,98) = -1
PES_down(3,99) = 180; phase_PES_down(3,99) = -1
PES_down(3,100) = 185; phase_PES_down(3,100) = 1
PES_down(3,101) = 186; phase_PES_down(3,101) = 1
PES_down(3,103) = 173; phase_PES_down(3,103) = -1
PES_down(3,104) = 176; phase_PES_down(3,104) = -1
PES_down(3,105) = 181; phase_PES_down(3,105) = 1
PES_down(3,106) = 190; phase_PES_down(3,106) = -1
PES_down(3,108) = 192; phase_PES_down(3,108) = 1
PES_down(3,109) = 193; phase_PES_down(3,109) = -1
PES_down(3,111) = 195; phase_PES_down(3,111) = 1
PES_down(3,113) = 170; phase_PES_down(3,113) = -1
PES_down(3,119) = 168; phase_PES_down(3,119) = -1
PES_down(3,120) = 169; phase_PES_down(3,120) = -1
PES_down(3,121) = 171; phase_PES_down(3,121) = 1
PES_down(3,129) = 204; phase_PES_down(3,129) = -1
PES_down(3,131) = 203; phase_PES_down(3,131) = -1
PES_down(3,133) = 202; phase_PES_down(3,133) = -1
PES_down(3,135) = 209; phase_PES_down(3,135) = 1
PES_down(3,137) = 210; phase_PES_down(3,137) = 1
PES_down(3,140) = 198; phase_PES_down(3,140) = -1
PES_down(3,141) = 196; phase_PES_down(3,141) = -1
PES_down(3,142) = 201; phase_PES_down(3,142) = -1
PES_down(3,143) = 199; phase_PES_down(3,143) = -1
PES_down(3,144) = 207; phase_PES_down(3,144) = 1
PES_down(3,145) = 206; phase_PES_down(3,145) = 1
PES_down(3,146) = 214; phase_PES_down(3,146) = -1
PES_down(3,151) = 216; phase_PES_down(3,151) = 1
PES_down(3,152) = 217; phase_PES_down(3,152) = -1
PES_down(3,157) = 219; phase_PES_down(3,157) = 1
PES_down(3,158) = 208; phase_PES_down(3,158) = 1
PES_down(3,160) = 215; phase_PES_down(3,160) = 1
PES_down(3,162) = 218; phase_PES_down(3,162) = 1
PES_down(3,164) = 221; phase_PES_down(3,164) = -1
PES_down(3,165) = 223; phase_PES_down(3,165) = -1
PES_down(3,167) = 224; phase_PES_down(3,167) = 1
PES_down(3,172) = 234; phase_PES_down(3,172) = -1
PES_down(3,174) = 235; phase_PES_down(3,174) = -1
PES_down(3,175) = 238; phase_PES_down(3,175) = -1
PES_down(3,177) = 239; phase_PES_down(3,177) = -1
PES_down(3,182) = 242; phase_PES_down(3,182) = 1
PES_down(3,183) = 243; phase_PES_down(3,183) = 1
PES_down(3,184) = 244; phase_PES_down(3,184) = 1
PES_down(3,187) = 233; phase_PES_down(3,187) = -1
PES_down(3,188) = 236; phase_PES_down(3,188) = 1
PES_down(3,189) = 240; phase_PES_down(3,189) = 1
PES_down(3,191) = 246; phase_PES_down(3,191) = 1
PES_down(3,194) = 247; phase_PES_down(3,194) = 1
PES_down(3,197) = 227; phase_PES_down(3,197) = -1
PES_down(3,200) = 229; phase_PES_down(3,200) = -1
PES_down(3,205) = 231; phase_PES_down(3,205) = 1
PES_down(3,211) = 226; phase_PES_down(3,211) = -1
PES_down(3,212) = 228; phase_PES_down(3,212) = 1
PES_down(3,213) = 230; phase_PES_down(3,213) = 1
PES_down(3,220) = 248; phase_PES_down(3,220) = -1
PES_down(3,222) = 250; phase_PES_down(3,222) = 1
PES_down(3,224) = 251; phase_PES_down(3,224) = 1
PES_down(3,232) = 252; phase_PES_down(3,232) = -1
PES_down(3,237) = 254; phase_PES_down(3,237) = 1
PES_down(3,241) = 255; phase_PES_down(3,241) = 1
PES_down(3,245) = 253; phase_PES_down(3,245) = 1
PES_down(3,249) = 256; phase_PES_down(3,249) = 1

PES_up(4,1) = 5; phase_PES_up(4,1) = 1
PES_up(4,2) = 12; phase_PES_up(4,2) = 1
PES_up(4,3) = 14; phase_PES_up(4,3) = 1
PES_up(4,4) = 15; phase_PES_up(4,4) = 1
PES_up(4,6) = 27; phase_PES_up(4,6) = 1
PES_up(4,7) = 31; phase_PES_up(4,7) = 1
PES_up(4,8) = 33; phase_PES_up(4,8) = 1
PES_up(4,9) = 37; phase_PES_up(4,9) = -1
PES_up(4,10) = 39; phase_PES_up(4,10) = 1
PES_up(4,11) = 40; phase_PES_up(4,11) = 1
PES_up(4,13) = 41; phase_PES_up(4,13) = 1
PES_up(4,16) = 73; phase_PES_up(4,16) = 1
PES_up(4,17) = 76; phase_PES_up(4,17) = 1
PES_up(4,18) = 87; phase_PES_up(4,18) = -1
PES_up(4,19) = 79; phase_PES_up(4,19) = 1
PES_up(4,20) = 91; phase_PES_up(4,20) = -1
PES_up(4,21) = 93; phase_PES_up(4,21) = -1
PES_up(4,22) = 50; phase_PES_up(4,22) = 1
PES_up(4,23) = 51; phase_PES_up(4,23) = 1
PES_up(4,24) = 53; phase_PES_up(4,24) = 1
PES_up(4,25) = 54; phase_PES_up(4,25) = 1
PES_up(4,26) = 63; phase_PES_up(4,26) = -1
PES_up(4,28) = 56; phase_PES_up(4,28) = 1
PES_up(4,29) = 57; phase_PES_up(4,29) = 1
PES_up(4,30) = 67; phase_PES_up(4,30) = -1
PES_up(4,32) = 69; phase_PES_up(4,32) = -1
PES_up(4,34) = 62; phase_PES_up(4,34) = 1
PES_up(4,35) = 66; phase_PES_up(4,35) = 1
PES_up(4,36) = 68; phase_PES_up(4,36) = 1
PES_up(4,38) = 94; phase_PES_up(4,38) = 1
PES_up(4,42) = 112; phase_PES_up(4,42) = 1
PES_up(4,43) = 121; phase_PES_up(4,43) = -1
PES_up(4,44) = 124; phase_PES_up(4,44) = -1
PES_up(4,45) = 127; phase_PES_up(4,45) = -1
PES_up(4,46) = 97; phase_PES_up(4,46) = 1
PES_up(4,47) = 98; phase_PES_up(4,47) = 1
PES_up(4,48) = 99; phase_PES_up(4,48) = 1
PES_up(4,49) = 105; phase_PES_up(4,49) = -1
PES_up(4,52) = 108; phase_PES_up(4,52) = -1
PES_up(4,55) = 111; phase_PES_up(4,55) = -1
PES_up(4,58) = 103; phase_PES_up(4,58) = 1
PES_up(4,59) = 104; phase_PES_up(4,59) = 1
PES_up(4,60) = 106; phase_PES_up(4,60) = 1
PES_up(4,61) = 107; phase_PES_up(4,61) = 1
PES_up(4,64) = 109; phase_PES_up(4,64) = 1
PES_up(4,65) = 110; phase_PES_up(4,65) = 1
PES_up(4,70) = 133; phase_PES_up(4,70) = 1
PES_up(4,71) = 132; phase_PES_up(4,71) = 1
PES_up(4,72) = 130; phase_PES_up(4,72) = 1
PES_up(4,74) = 145; phase_PES_up(4,74) = -1
PES_up(4,75) = 144; phase_PES_up(4,75) = -1
PES_up(4,77) = 151; phase_PES_up(4,77) = -1
PES_up(4,78) = 150; phase_PES_up(4,78) = -1
PES_up(4,80) = 157; phase_PES_up(4,80) = -1
PES_up(4,81) = 156; phase_PES_up(4,81) = -1
PES_up(4,82) = 141; phase_PES_up(4,82) = 1
PES_up(4,83) = 143; phase_PES_up(4,83) = 1
PES_up(4,84) = 147; phase_PES_up(4,84) = 1
PES_up(4,85) = 149; phase_PES_up(4,85) = 1
PES_up(4,86) = 160; phase_PES_up(4,86) = -1
PES_up(4,88) = 153; phase_PES_up(4,88) = 1
PES_up(4,89) = 155; phase_PES_up(4,89) = 1
PES_up(4,90) = 162; phase_PES_up(4,90) = -1
PES_up(4,92) = 163; phase_PES_up(4,92) = -1
PES_up(4,95) = 171; phase_PES_up(4,95) = -1
PES_up(4,96) = 167; phase_PES_up(4,96) = -1
PES_up(4,100) = 164; phase_PES_up(4,100) = 1
PES_up(4,101) = 165; phase_PES_up(4,101) = 1
PES_up(4,102) = 166; phase_PES_up(4,102) = 1
PES_up(4,113) = 205; phase_PES_up(4,113) = -1
PES_up(4,114) = 206; phase_PES_up(4,114) = -1
PES_up(4,115) = 207; phase_PES_up(4,115) = -1
PES_up(4,116) = 196; phase_PES_up(4,116) = 1
PES_up(4,117) = 199; phase_PES_up(4,117) = 1
PES_up(4,118) = 202; phase_PES_up(4,118) = 1
PES_up(4,119) = 212; phase_PES_up(4,119) = -1
PES_up(4,120) = 213; phase_PES_up(4,120) = -1
PES_up(4,122) = 215; phase_PES_up(4,122) = -1
PES_up(4,123) = 216; phase_PES_up(4,123) = -1
PES_up(4,125) = 218; phase_PES_up(4,125) = -1
PES_up(4,126) = 219; phase_PES_up(4,126) = -1
PES_up(4,128) = 181; phase_PES_up(4,128) = -1
PES_up(4,129) = 182; phase_PES_up(4,129) = -1
PES_up(4,131) = 183; phase_PES_up(4,131) = -1
PES_up(4,134) = 173; phase_PES_up(4,134) = 1
PES_up(4,135) = 174; phase_PES_up(4,135) = 1
PES_up(4,136) = 176; phase_PES_up(4,136) = 1
PES_up(4,137) = 177; phase_PES_up(4,137) = 1
PES_up(4,138) = 179; phase_PES_up(4,138) = 1
PES_up(4,139) = 180; phase_PES_up(4,139) = 1
PES_up(4,140) = 188; phase_PES_up(4,140) = -1
PES_up(4,142) = 189; phase_PES_up(4,142) = -1
PES_up(4,146) = 191; phase_PES_up(4,146) = -1
PES_up(4,148) = 192; phase_PES_up(4,148) = -1
PES_up(4,152) = 194; phase_PES_up(4,152) = -1
PES_up(4,154) = 195; phase_PES_up(4,154) = -1
PES_up(4,158) = 187; phase_PES_up(4,158) = 1
PES_up(4,159) = 190; phase_PES_up(4,159) = 1
PES_up(4,161) = 193; phase_PES_up(4,161) = 1
PES_up(4,168) = 228; phase_PES_up(4,168) = -1
PES_up(4,169) = 230; phase_PES_up(4,169) = -1
PES_up(4,170) = 231; phase_PES_up(4,170) = -1
PES_up(4,172) = 222; phase_PES_up(4,172) = -1
PES_up(4,175) = 224; phase_PES_up(4,175) = -1
PES_up(4,178) = 225; phase_PES_up(4,178) = -1
PES_up(4,184) = 220; phase_PES_up(4,184) = 1
PES_up(4,185) = 221; phase_PES_up(4,185) = 1
PES_up(4,186) = 223; phase_PES_up(4,186) = 1
PES_up(4,197) = 237; phase_PES_up(4,197) = -1
PES_up(4,198) = 236; phase_PES_up(4,198) = -1
PES_up(4,200) = 241; phase_PES_up(4,200) = -1
PES_up(4,201) = 240; phase_PES_up(4,201) = -1
PES_up(4,203) = 243; phase_PES_up(4,203) = -1
PES_up(4,204) = 242; phase_PES_up(4,204) = -1
PES_up(4,208) = 233; phase_PES_up(4,208) = 1
PES_up(4,209) = 235; phase_PES_up(4,209) = 1
PES_up(4,210) = 239; phase_PES_up(4,210) = 1
PES_up(4,211) = 245; phase_PES_up(4,211) = -1
PES_up(4,214) = 246; phase_PES_up(4,214) = -1
PES_up(4,217) = 247; phase_PES_up(4,217) = -1
PES_up(4,226) = 253; phase_PES_up(4,226) = -1
PES_up(4,227) = 254; phase_PES_up(4,227) = -1
PES_up(4,229) = 255; phase_PES_up(4,229) = -1
PES_up(4,232) = 249; phase_PES_up(4,232) = -1
PES_up(4,234) = 250; phase_PES_up(4,234) = -1
PES_up(4,238) = 251; phase_PES_up(4,238) = -1
PES_up(4,244) = 248; phase_PES_up(4,244) = 1
PES_up(4,252) = 256; phase_PES_up(4,252) = -1

PES_down(4,1) = 9; phase_PES_down(4,1) = 1
PES_down(4,2) = 26; phase_PES_down(4,2) = 1
PES_down(4,3) = 30; phase_PES_down(4,3) = 1
PES_down(4,4) = 32; phase_PES_down(4,4) = 1
PES_down(4,5) = 37; phase_PES_down(4,5) = 1
PES_down(4,6) = 18; phase_PES_down(4,6) = 1
PES_down(4,7) = 20; phase_PES_down(4,7) = 1
PES_down(4,8) = 21; phase_PES_down(4,8) = 1
PES_down(4,10) = 49; phase_PES_down(4,10) = 1
PES_down(4,11) = 52; phase_PES_down(4,11) = 1
PES_down(4,12) = 63; phase_PES_down(4,12) = 1
PES_down(4,13) = 55; phase_PES_down(4,13) = 1
PES_down(4,14) = 67; phase_PES_down(4,14) = 1
PES_down(4,15) = 69; phase_PES_down(4,15) = 1
PES_down(4,16) = 43; phase_PES_down(4,16) = 1
PES_down(4,17) = 44; phase_PES_down(4,17) = 1
PES_down(4,19) = 45; phase_PES_down(4,19) = 1
PES_down(4,22) = 75; phase_PES_down(4,22) = 1
PES_down(4,23) = 74; phase_PES_down(4,23) = 1
PES_down(4,24) = 78; phase_PES_down(4,24) = 1
PES_down(4,25) = 77; phase_PES_down(4,25) = 1
PES_down(4,27) = 87; phase_PES_down(4,27) = 1
PES_down(4,28) = 81; phase_PES_down(4,28) = 1
PES_down(4,29) = 80; phase_PES_down(4,29) = 1
PES_down(4,31) = 91; phase_PES_down(4,31) = 1
PES_down(4,33) = 93; phase_PES_down(4,33) = 1
PES_down(4,34) = 86; phase_PES_down(4,34) = 1
PES_down(4,35) = 90; phase_PES_down(4,35) = 1
PES_down(4,36) = 92; phase_PES_down(4,36) = 1
PES_down(4,38) = 96; phase_PES_down(4,38) = 1
PES_down(4,39) = 105; phase_PES_down(4,39) = 1
PES_down(4,40) = 108; phase_PES_down(4,40) = 1
PES_down(4,41) = 111; phase_PES_down(4,41) = 1
PES_down(4,42) = 95; phase_PES_down(4,42) = 1
PES_down(4,46) = 128; phase_PES_down(4,46) = 1
PES_down(4,47) = 129; phase_PES_down(4,47) = 1
PES_down(4,48) = 131; phase_PES_down(4,48) = 1
PES_down(4,50) = 144; phase_PES_down(4,50) = 1
PES_down(4,51) = 145; phase_PES_down(4,51) = 1
PES_down(4,53) = 150; phase_PES_down(4,53) = 1
PES_down(4,54) = 151; phase_PES_down(4,54) = 1
PES_down(4,56) = 156; phase_PES_down(4,56) = 1
PES_down(4,57) = 157; phase_PES_down(4,57) = 1
PES_down(4,58) = 140; phase_PES_down(4,58) = 1
PES_down(4,59) = 142; phase_PES_down(4,59) = 1
PES_down(4,60) = 146; phase_PES_down(4,60) = 1
PES_down(4,61) = 148; phase_PES_down(4,61) = 1
PES_down(4,62) = 160; phase_PES_down(4,62) = 1
PES_down(4,64) = 152; phase_PES_down(4,64) = 1
PES_down(4,65) = 154; phase_PES_down(4,65) = 1
PES_down(4,66) = 162; phase_PES_down(4,66) = 1
PES_down(4,68) = 163; phase_PES_down(4,68) = 1
PES_down(4,70) = 113; phase_PES_down(4,70) = 1
PES_down(4,71) = 114; phase_PES_down(4,71) = 1
PES_down(4,72) = 115; phase_PES_down(4,72) = 1
PES_down(4,73) = 121; phase_PES_down(4,73) = 1
PES_down(4,76) = 124; phase_PES_down(4,76) = 1
PES_down(4,79) = 127; phase_PES_down(4,79) = 1
PES_down(4,82) = 119; phase_PES_down(4,82) = 1
PES_down(4,83) = 120; phase_PES_down(4,83) = 1
PES_down(4,84) = 122; phase_PES_down(4,84) = 1
PES_down(4,85) = 123; phase_PES_down(4,85) = 1
PES_down(4,88) = 125; phase_PES_down(4,88) = 1
PES_down(4,89) = 126; phase_PES_down(4,89) = 1
PES_down(4,94) = 167; phase_PES_down(4,94) = 1
PES_down(4,97) = 181; phase_PES_down(4,97) = 1
PES_down(4,98) = 182; phase_PES_down(4,98) = 1
PES_down(4,99) = 183; phase_PES_down(4,99) = 1
PES_down(4,100) = 172; phase_PES_down(4,100) = 1
PES_down(4,101) = 175; phase_PES_down(4,101) = 1
PES_down(4,102) = 178; phase_PES_down(4,102) = 1
PES_down(4,103) = 188; phase_PES_down(4,103) = 1
PES_down(4,104) = 189; phase_PES_down(4,104) = 1
PES_down(4,106) = 191; phase_PES_down(4,106) = 1
PES_down(4,107) = 192; phase_PES_down(4,107) = 1
PES_down(4,109) = 194; phase_PES_down(4,109) = 1
PES_down(4,110) = 195; phase_PES_down(4,110) = 1
PES_down(4,112) = 171; phase_PES_down(4,112) = 1
PES_down(4,116) = 168; phase_PES_down(4,116) = 1
PES_down(4,117) = 169; phase_PES_down(4,117) = 1
PES_down(4,118) = 170; phase_PES_down(4,118) = 1
PES_down(4,130) = 207; phase_PES_down(4,130) = 1
PES_down(4,132) = 206; phase_PES_down(4,132) = 1
PES_down(4,133) = 205; phase_PES_down(4,133) = 1
PES_down(4,134) = 198; phase_PES_down(4,134) = 1
PES_down(4,135) = 197; phase_PES_down(4,135) = 1
PES_down(4,136) = 201; phase_PES_down(4,136) = 1
PES_down(4,137) = 200; phase_PES_down(4,137) = 1
PES_down(4,138) = 204; phase_PES_down(4,138) = 1
PES_down(4,139) = 203; phase_PES_down(4,139) = 1
PES_down(4,141) = 212; phase_PES_down(4,141) = 1
PES_down(4,143) = 213; phase_PES_down(4,143) = 1
PES_down(4,147) = 215; phase_PES_down(4,147) = 1
PES_down(4,149) = 216; phase_PES_down(4,149) = 1
PES_down(4,153) = 218; phase_PES_down(4,153) = 1
PES_down(4,155) = 219; phase_PES_down(4,155) = 1
PES_down(4,158) = 211; phase_PES_down(4,158) = 1
PES_down(4,159) = 214; phase_PES_down(4,159) = 1
PES_down(4,161) = 217; phase_PES_down(4,161) = 1
PES_down(4,164) = 222; phase_PES_down(4,164) = 1
PES_down(4,165) = 224; phase_PES_down(4,165) = 1
PES_down(4,166) = 225; phase_PES_down(4,166) = 1
PES_down(4,173) = 236; phase_PES_down(4,173) = 1
PES_down(4,174) = 237; phase_PES_down(4,174) = 1
PES_down(4,176) = 240; phase_PES_down(4,176) = 1
PES_down(4,177) = 241; phase_PES_down(4,177) = 1
PES_down(4,179) = 242; phase_PES_down(4,179) = 1
PES_down(4,180) = 243; phase_PES_down(4,180) = 1
PES_down(4,184) = 232; phase_PES_down(4,184) = 1
PES_down(4,185) = 234; phase_PES_down(4,185) = 1
PES_down(4,186) = 238; phase_PES_down(4,186) = 1
PES_down(4,187) = 245; phase_PES_down(4,187) = 1
PES_down(4,190) = 246; phase_PES_down(4,190) = 1
PES_down(4,193) = 247; phase_PES_down(4,193) = 1
PES_down(4,196) = 228; phase_PES_down(4,196) = 1
PES_down(4,199) = 230; phase_PES_down(4,199) = 1
PES_down(4,202) = 231; phase_PES_down(4,202) = 1
PES_down(4,208) = 226; phase_PES_down(4,208) = 1
PES_down(4,209) = 227; phase_PES_down(4,209) = 1
PES_down(4,210) = 229; phase_PES_down(4,210) = 1
PES_down(4,220) = 249; phase_PES_down(4,220) = 1
PES_down(4,221) = 250; phase_PES_down(4,221) = 1
PES_down(4,223) = 251; phase_PES_down(4,223) = 1
PES_down(4,233) = 253; phase_PES_down(4,233) = 1
PES_down(4,235) = 254; phase_PES_down(4,235) = 1
PES_down(4,239) = 255; phase_PES_down(4,239) = 1
PES_down(4,244) = 252; phase_PES_down(4,244) = 1
PES_down(4,248) = 256; phase_PES_down(4,248) = 1

sites: do j=1,4  ! calculating the IPES matrices 
   do i=1,256
      if (PES_down(j,i) /= 0) then
         phase_IPES_down(j,PES_down(j,i)) = phase_PES_down(j,i)
         IPES_down(j,PES_down(j,i)) = i
      end if
      if (PES_up(j,i) /= 0) then
         IPES_up(j,PES_up(j,i)) = i
         phase_IPES_up(j,PES_up(j,i)) = phase_PES_up(j,i)
      end if
   end do
end do sites

end subroutine transformations

end module routines
