module hamiltonian

contains

!subroutine fock_state(fock)

!integer, allocatable, dimension(:,:), intent (in out) :: fock   ! array holding all the fock states
!integer :: nstates ! total number of states
!integer, intent(in) :: nsites ! total number of sites
!integer :: state !counter that goes through all the integers representing the states
!integer :: position  !puts the states in order by increasing integer value within fock array
!integer :: nelectrons ! counts the amount of electrons in each state
!integer :: site  !counter to go through all the sites (looking for sites with electrons)

!allocate(fock(max_electrons,6)) 

!nstates = 2**nsites
!position = 1

!do state=0,nstates-1
!  nelectrons = 0
!  do site=1,nsites
!      nelectrons = nelectrons + ibits(state,site-1,1)
!      fock(nelectrons,position(nelectrons)) = state
!      position(nelectrons) = position + 1
!  end do
!end do

!end subroutine fock_state

subroutine site_potentials(nsites,E,delta)              ! randomly creates site potentials and stores them in the array E

real :: random1, random2, temp
real, intent(out) :: E(2,nsites)
real, intent(in) :: delta

do i=1,npairs
   call random_number(random1)
   call random_number(random2)
   random1 = random1 - 0.5
   random2 = random2 - 0.5
   E(i,1) = delta*max(random1,random2)
   E(i,2) = delta*min(random1,random2)
end do
   
end subroutine site_potentials

subroutine make_hamiltonian(E,U,t,iteration,H00,H10,H01,H20,H02,H11,H21,H12,H22)

real, dimension(1,2), intent(in) :: E
real, intent(in) :: t
real, intent(in) :: U
real, intent(out) :: H00, H01(2,2),H10(2,2),H11(4,4),H20,H02,H22,H21(2,2), H12(2,2)
integer :: i,j ! counters
 
H00 = 0
H01 = 0
H10 = 0
H02 = 0
H20 = 0
H11 = 0
H21 = 0
H12 = 0
H22 = 0

! make the H01/H10 matrices

H01(1,1) = E(iteration,1)
H01(2,2) = E(iteration,2)
H01(1,2) = t
H01(2,1) = t

H10 = H01

! make H02/H20 matrices

H02 = E(iteration,1) + E(iteration,2)
H20 = H02

! make H11 matrix
do i=1,4
   do j=1,4
      if(i/=j) then
         H11(i,j) = t
      end if
   end do
end do

H11(1,1) = E(iteration,1) + E(iteration,2)
H11(2,2) = H11(1,1)
H11(3,3) = 2*E(iteration,1) + U
H11(4,4) = 2*E(iteration,2) + U
H11(1,2) = 0  !jump from (up,down) to (down,up)
H11(2,1) = 0  !jump from (down,up) to (up,down)
H11(3,4) = 0
H11(4,3) = 0

!make H21/H12 matrix
H21(1,2) = -t
H21(2,1) = -t
H21(1,1) = 2*E(iteration,1) + E(iteration,2) + U
H21(2,2) = E(iteration,1) + 2*E(iteration,2) + U
H12 = H21
!make H22 

H22 = 2*E(iteration,1) + 2*E(iteration,2) + 2*U
end subroutine make_hamiltonian

end module
