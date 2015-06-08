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
integer :: i,j ! counter

!------for lapack------------
integer :: INFO = 0
integer :: LWORK
integer, allocatable, dimension(:) :: WORK

!----------zero variables for each loop---------------------------------------

H00=0; H10=0; H20=0; H11=0;
W00=0; W10=0; H20=0; W11=0;

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

end module routines
