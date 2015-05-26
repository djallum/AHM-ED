module hamiltonian

use math_setup

implicit none 

real, allocatable, dimension(:,:) :: H01,H10,H11,H12,H21
real :: H00,H20,H02,H22
integer,allocatable, dimension(:) :: fockstates

contains

subroutine make_fockstates(nsites,fockstates)

integer, allocatable, dimension(:), intent(out) :: fockstates  !array that stores the fock states. Fock states are in binary and this stores the integers. ex. 3 is (0011)
integer, intent(in) :: nsites
integer :: nstates   !number of states when considering only up electrons (would be the same as if only considering down electrons)
integer :: site, state,i    !counter for going through sites and states
integer :: nelectrons   !counts the number of electrons
integer, allocatable, dimension(:) :: location !location(i) indicates the place in the fockstate_half array when they start having i electrons

nstates = 2**nsites         
allocate(fockstates(nstates))
allocate(location(0:nsites))    !the number of sites is the same as the maximum number of electrons

fockstates(1) = 0  !first state is the vacuum 
location(0) = 1

do nelectrons=1,nsites           !the number of sites is the same as the maximum number of electrons
   location(i) = location(i-1) + choose(nsites,nelectrons-1)
end do

do state=1, nstates-1
   nelectrons = 0
   do site=0,nsites-1                              !loop through all the sites looking for electrons
      nelectrons = nelectrons + ibits(state,site-1,1)  !checks if electron is on that site and adds it to total
   end do
   fockstates(location(nelectrons)) = state
   location(nelectrons) = location(nelectrons) + 1
end do 

deallocate(location)
end subroutine make_fockstates

subroutine make_hamiltonian(nsites,t,U,fockstates)

integer, allocatable, dimension(:) :: location !location(i) indicates the place in the fockstate_half array when they start having i electrons
integer :: nstates,up,down
integer :: nelectrons_up,nelectrons_down
logical :: hoop,place
integer :: neighbour(2), newstate
integer :: hooper(2)
integer :: state,site,i, nelectrons
integer, intent(in) :: nsites
real, intent(in) :: t, U
integer, intent(in) :: fockstates(2**nsites)

allocate(location(0:nsites))    !the number of sites is the same as the maximum number of electrons
nstates = 2**nsites

do nelectrons=1,nsites           !the number of sites is the same as the maximum number of electrons
   location(i) = location(i-1) + choose(nsites,nelectrons-1)
end do

allocate(H01(nsites,nsites)) 
allocate(H10(nsites,nsites))
allocate(H11(nsites*2,nsites*2))
allocate(H21(nsites,nsites))
allocate(H12(nsites,nsites))

H00 = 0

!calculate H01

up = 1
down = 0

do state=1,2
   hooper(fockstates(location(up-1) + state)) = state
end do


neighbour(1) =2
neighbour(2) =1
up = 1
down = 0
do state=1,2
   do site=1,2
      if ( btest(fockstates(location(up-1) + state),site-1)) then
          newstate = ibclr(fockstates(location(up-1) + i),site-1)
          newstate = ibset(fockstates(location(up-1) + i), neighbour(site))
          H01(state,hooper(newstate)) = t
       end if
    end do
end do
 
   
end subroutine make_hamiltonian

end module hamiltonian

