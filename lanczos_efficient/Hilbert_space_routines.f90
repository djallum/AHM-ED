module Hilbert_space_routines
  !
  ! This module contains subroutines that require knowledge of how the Hilbert space is 
  ! structured.  This does not include, for example, routines related to the Lanczos or
  ! recursion algorithms.  
  !
  use math_setup
  use parameters

  implicit none
  
  integer, private :: nsite_, max_electrons_, total_states_  ! variables end with underscore to indicate they are in subroutine.
  integer, allocatable, private :: nstates_(:), states_(:), block_(:)
  integer, allocatable, private :: h_index_(:,:)   ! index of neighbours of lattice sites.
  integer, private :: n_up_,n_dn_  ! number of up/down electrons
  integer, private :: n_           ! total number of states for a given n_up_, n_dn_
  integer, private :: nx_,ny_,lattice_    ! x and y dimension of the lattice (measured in sites) and the integer that specifies lattice shape.
  integer, private, parameter :: nbrs_ = 4  ! short for neighbourhs. Each site can have up to four neighbourhs...
  integer, private, allocatable :: state_lookup_up_(:,:), phase_lookup_up_(:,:)
  integer, private, allocatable :: state_lookup_dn_(:,:), phase_lookup_dn_(:,:)
  real (kind=real_kind), private, allocatable :: n_up_local_(:),n_dn_local_(:),spin_correlations_(:)
  real :: time1_,time2_,time3_

contains

!*******************************
! The subroutines here are meant to replace those in qnummap.f90.
! The main difference is that this is for a single spin only.
! This means that a Fock state is specified by two numbers:
! one for the spin-up part, and one for the spin-down part.
! This saves us a lot of memory.
!*******************************

  subroutine make_hilbert_space(nx,ny,lattice)
    integer, intent (in) :: nx,ny,lattice
    integer, allocatable :: icount(:)
    integer :: ne,istate,isite

    lattice_ = lattice            ! indicates the type of lattice be made (different integers represent different arrangments)
    nx_ = nx                      ! number of sites in x direction for crystal lattice
    ny_ = ny                      ! number of sites in x direction for crystal lattice
    if (lattice < 3) then
       nsite_ = nx * ny           ! number of sites in lattice for 1-band Hubbard model
    else if (lattice==3) then
       nsite_ = nx * ny * 3       ! number of sites in the 3-band Hubbard model
    end if
    max_electrons_ = nsite_       ! maximum number of spin-up electrons
    total_states_ = 2**nsite_     ! number of spin up states
    if (huge(istate) < total_states_) then
       print *, 'Largest integer smaller than total states.'
       stop
    end if
    allocate(nstates_(0:max_electrons_))  ! Number of Fock states for each amount of spin-up electrons.
    allocate(states_(total_states_))      ! Stores all the spin-up Fock states
    allocate(block_(0:max_electrons_))    ! block_(j) contains the position in the array states_ of the 
                                          ! first Fock state with j electrons.  
                                          ! Thus states_(block_(j):block_(j+1)-1) is the set of Fock states
                                          ! with j spin-up electrons.

    block_ = 0  
    nstates_ = 0
    block_(0) = 1
    nstates_(0) = 1
    do ne = 1, max_electrons_
       nstates_(ne) = choose(nsite_,ne)          ! number of spin-up fock states with ne electrons
       block_(ne) = block_(ne-1) + nstates_(ne-1) ! block index updating
    end do

    ! Check that we have constructed block_ correctly
    if (sum(nstates_).ne.total_states_.or.sum(nstates_).eq.0) then
       print *, 'Integer number overflow.'
       stop
    end if

    ! Now make a list of all fock states.  These are ordered such that:
    ! (i) states with more electrons appear later in the list 
    ! (ii) within a block of states with ne electrons, the values of the Fock states (the integer's value thats representing them) are increasing.
    ! example for three cite: 0(000),1(001),2(010),4(100),3(011),5(101),6(110),7(111).

    allocate(icount(0:max_electrons_))
    icount = block_
    do istate = 0, total_states_ - 1
       ne = 0
       do isite = 0, nsite_ - 1           !states are labeled in binary. 0 means empty and 1 means upspin. Ex. for 3-cite 5(101) is up electron 
          ne = ne + ibits(istate,isite,1) !on the first and third site. This loop counts how many 'ones' there are (electrons) for a particular integer
       end do                             !(when it's in binary).
       states_(icount(ne)) = istate       !This now puts that integer (that's representing a state) into the proper location (first place with that
       icount(ne) = icount(ne) + 1        !many electrons). It bumps up the index one so that the next integer with that many electrons goes in the 
    end do                                !next place in the array.
    deallocate(icount)
  end subroutine make_hilbert_space

!*******************************

  subroutine print_hilbert_space()
    integer :: ne
    character(11) :: x
    x = '(50(i6,1x))'
    do ne = 0,max_electrons_-1
       print x, ne,block_(ne),nstates_(ne),states_(block_(ne):block_(ne+1)-1)
    end do
    print x, ne,block_(ne),nstates_(ne),states_(block_(ne))
  end subroutine print_hilbert_space

!*******************************

  subroutine cleanup_hilbert_space()
    if (allocated(nstates_)) deallocate(nstates_)
    if (allocated(states_)) deallocate(states_)
    if (allocated(h_index_)) deallocate(h_index_)
    if (allocated(n_up_local_)) deallocate(n_up_local_,n_dn_local_)
  end subroutine cleanup_hilbert_space

!*******************************

  subroutine num_states(n_up,n_dn,n)
    integer :: n_up,n_dn,n
    n_up_ = n_up
    n_dn_ = n_dn
    n = nstates_(n_up)*nstates_(n_dn)
    n_ = n
  end subroutine num_states

!***********************************

  subroutine make_h_index(lattice)
    integer :: lattice
  
    ! This, and the subroutines which are called from here, are the only places where
    ! the details of the lattice matter (assuming nearest-neighbour hopping only).  
    ! All the details are hidden in the array h_index
    ! in h_index_ there is a row for each cite and four columns for the cite index of each of it's possible four nearest neighbouhrs
  
    if (lattice == 0) then
       call make_h_index_rectangular()
    else if (lattice == -1) then
       call make_h_index_linear()
    else if (lattice == 1) then
       call make_h_index_12h3()       
    else if (lattice == 2) then
       print *,"10 site lattice not yet implemented"
       stop
    else if (lattice == 3) then
       call make_h_index_3band()
    else
       print *,"Error in make_h_index"
       print *,"Lattice type not implemented."
       print *,"You chose: ",lattice
    end if
  end subroutine make_h_index

!***********************************

  subroutine make_h_index_rectangular()
    integer :: ix,iy,i,jx,jy,j

    allocate(h_index_(nsite_,nbrs_)) ! the value for neighbours is defined when initialized
    do ix = 1,nx_
       do iy = 1,ny_
          i = ix + (iy-1)*nx_

          jx = mod(ix,nx_) + 1
          jy = iy
          j = jx + (jy-1)*nx_
          h_index_(i,1) = j
          h_index_(j,2) = i

          jy = mod(iy,ny_) + 1
          jx = ix
          j = jx + (jy-1)*nx_
          h_index_(i,3) = j
          h_index_(j,4) = i
       end do
    end do
  end subroutine make_h_index_rectangular

!***********************************

  subroutine make_h_index_linear()
    integer :: ix,i,jx,j

    allocate(h_index_(nsite_,nbrs_))
    h_index_ = 0
    if (ny_==1) then
       do ix = 1,nx_
          i = ix 
          jx = mod(ix,nx_) + 1
          j = jx 
          h_index_(i,1) = j
          h_index_(j,2) = i
       end do
    else 
       print *,"Error in make_h_index_linear"
       print *,"unsupported lattice shape"
    end if

  end subroutine make_h_index_linear
  

!***********************************

  subroutine make_h_index_12h3()
    integer :: ix,iy,i,jx,jy,j

    allocate(h_index_(nsite_,nbrs_))
    do ix = 1,nx_
       do iy = 1,ny_
          i = ix + (iy-1)*nx_
	
	  jx = mod(ix,nx_) + 1
	  jy = iy
	  j = jx + (jy-1)*nx_
	  h_index_(i,1) = j	! nn right
	  h_index_(j,2) = i	! nn left

	  if (iy.eq.1) then
	     jy = iy + 1
	     if (ix.ge.2) then
	        jx = ix - 1
	     else ! ix.eq.1
	        jx = 4
	     end if
	  else
	     jy = mod(iy,ny_) + 1
	     jx = ix
	  end if
          j = jx + (jy-1)*nx_
	  h_index_(i,3) = j	! nn above
	  h_index_(j,4) = i	! nn below
       end do
    end do
  end subroutine make_h_index_12h3

!***********************************

  subroutine make_h_index_3band()
    integer :: ix,iy,i,jx,jy,j

!!$    The way this is set up, 
!!$    h_index_(X,1) and h_index_(X,3) have +tpd matrix element
!!$    h_index_(X,2) and h_index_(X,4) have -tpd matrix element

    allocate(h_index_(nsite_,nbrs_))
    do ix = 1,nx_
       do iy = 1,ny_
          i = 1 + (ix-1)*3 + (iy-1)*nx_*3  ! index of Cu site in (ix,iy) unit cell

          jx = mod(ix,nx_) + 1
          jy = iy
          j = 1 + (jx-1)*3 + (jy-1)*nx_*3  ! index of Cu site in (jx,jy) unit cell

          h_index_(i,1) = i+1   ! Cu(i) -> px_(i)
          h_index_(i+1,1) = i   ! px(i) -> Cu(i)
          h_index_(i+1,2) = j   ! px(i) -> Cu(j)
          h_index_(j,2) = i+1   ! Cu(j) -> px(i)
          h_index_(i+1,3:4) = 0

          jy = mod(iy,ny_) + 1
          jx = ix
          j = 1 + (jx-1)*3 + (jy-1)*nx_*3  ! index of Cu site in (jx,jy) unit cell

          h_index_(i,3) = i+2   ! Cu(i) -> py(i)
          h_index_(i+2,1) = i   ! py(i) -> Cu(i)
          h_index_(i+2,2) = j   ! py(i) -> Cu(j)
          h_index_(j,4) = i+2   ! Cu(j) -> py(i)
          h_index_(i+2,3:4) = 0

       end do
    end do

  end subroutine make_h_index_3band

!***********************************

  subroutine initialise_lookup_tables()    
    if (.not.allocated(state_lookup_up_)) then
       allocate(state_lookup_up_(0:nbrs_*nsite_,nstates_(n_up_)))      ! nstates_(n_up_) gives the number of fock states with n_up_ electrons
       allocate(phase_lookup_up_(nbrs_*nsite_,nstates_(n_up_)))
    else
       print *,"Error initializing lookup table for up spin"
       stop
    end if

    if (.not.allocated(state_lookup_dn_)) then
       allocate(state_lookup_dn_(0:nbrs_*nsite_,nstates_(n_dn_)))
       allocate(phase_lookup_dn_(nbrs_*nsite_,nstates_(n_dn_)))
    else
       print *,"Error initializing lookup table for down spin"
       stop
    end if
  end subroutine initialise_lookup_tables

!***********************************

  subroutine deallocate_lookup_tables()
    deallocate(state_lookup_up_,state_lookup_dn_)
    deallocate(phase_lookup_up_,phase_lookup_dn_)
  end subroutine deallocate_lookup_tables

!***********************************

  subroutine make_state_lookup_table()
    integer :: istate,state_up,state_dn,i,j,k,new_state,phase
    integer :: icount
    !
    ! Dimensions:  state_lookup_up( nsites_*nbrs_, nstates(n_up_) ). Note that most elements are zero.
    ! state_lookup_XX_(i,j) is a list of all possible states connected to state j by a single electron hop.
    ! If XX = up, then we are looking at the spin-up part
    ! If XX = dn, then we are looking at the spin-down part
    ! state_lookup_XX_(0,j) contains the total number of possible hops for the state j. 
    ! Typically, this is much less than nsites_*nbrs_ (which is what would happen if all electrons could hop
    ! to all neighbours.
    ! Note that state_lookup_XX_ returns the quantum number (ie. the index) for the destination Fock state
    ! and takes the quantum number j as input.
    !
    ! This is not an especially big matrix (~ 7MB for 16 site lattice)
    !
    ! Make lookup table for up spins
    states_up: do istate = 1,nstates_(n_up_)  ! loop through all Fock states in spin-up basis
       state_up = states_(block_(n_up_)+istate-1)  ! take the next state from the table of Fock states
       icount = 0
       sites: do i = 1,nsite_               ! loop through all sites.
          nbrs: do k = 1,nbrs_              ! loop through all hopping directions
             j = h_index_(i,k)              ! j is the new site after hopping in direction k from site i
             hop: if (j/=0) then            ! make sure there is a hop associated with this value of k.
                ! Now check whether a hop is possible.  Requires initial site i be occupied & final site j be empty:
                if ((btest(state_up,i-1)).and.(.not.btest(state_up,j-1))) then
                   icount = icount + 1
                   call new_state_and_phase(state_up,i,j,new_state,phase)     
                   if (lattice_==3) phase = phase*(-1)**(k-1)  ! extra phase for the 3-band Hubbard model
                   ! statesindex uses a bisection search to find the index or quantum number for new_state.
                   ! This is slow, so we'll do this once, and store the results in a lookup table.
                   state_lookup_up_(icount,istate) = statesindex(n_up_,new_state) 
                   phase_lookup_up_(icount,istate) = phase
!!$                   print '(a16," -> ",a16)',convert_to_binary(state_up),convert_to_binary(new_state)
                end if
             end if hop
          end do nbrs
       end do sites
       state_lookup_up_(0,istate) = icount ! store the number of nonzero hops here.
    end do states_up

!!$    print *,1000,state_lookup_up_(1,1,:)
!!$    print *,states_(block_(n_up_)) ,states_(block_(n_up_)+state_lookup_up_(1,1,:)-1)  
!!$    stop

    ! Make lookup table for down spins
    states_dn: do istate = 1,nstates_(n_dn_)  ! loop through all Fock states in spin-down basis
       state_dn = states_(block_(n_dn_)+istate-1)  ! take the next state from the table of Fock states
       icount = 0
       sites2: do i = 1,nsite_              ! loop through all sites.
          nbrs2: do k = 1,nbrs_              ! loop through all hopping directions
             j = h_index_(i,k)              ! j is the new site after hopping in direction k from site i
             hop2: if (j/=0) then                 ! make sure there is a hop associated with this value of k.
                ! Now check whether a hop is possible.  Requires initial site i be occupied & final site j be empty:
                if ((btest(state_dn,i-1)).and.(.not.btest(state_dn,j-1))) then
                   icount = icount + 1
                   call new_state_and_phase(state_dn,i,j,new_state,phase)     
                   if (lattice_==3) phase = phase*(-1)**(k-1)  ! extra phase for the 3-band Hubbard model
                   ! statesindex uses a bisection search to find the index or quantum number for new_state.
                   ! This is slow, so we'll do this once, and store the results in a lookup table.
                   state_lookup_dn_(icount,istate) = statesindex(n_dn_,new_state) 
                   phase_lookup_dn_(icount,istate) = phase
                end if
             end if hop2
          end do nbrs2
       end do sites2
       state_lookup_dn_(0,istate) = icount ! store the number of nonzero hops here
    end do states_dn

  end subroutine make_state_lookup_table


!***********************************

  subroutine new_state_and_phase(state_up,i,j,new_state,phase)
    integer, intent(in) :: i,j,state_up
    integer, intent(out)  :: new_state,phase
    
    integer :: k

    if (.not.btest(state_up,j-1)) then      !this if statement is not neccesary since the very same if statement is used before calling the subroutine.
       new_state = ibset(state_up,j-1)      !add an electron ('one' in binary) to site j 
       new_state = ibclr(new_state,i-1)     !remove an electron from site i
       phase = 1
       if ((j-i==1).or.(j-i==-1)) then      !could be made clearer using the absolute value function.
          phase = 1
       else if (i>j) then
          do k=j,i-2
             phase = phase * (-1)**(ibits(new_state,k,1))  !anticommunitation to give the right sign on the hopping term
          end do
       else if (j>i) then
          do k=i,j-2
             phase = phase * (-1)**(ibits(new_state,k,1))
          end do
       end if
    else
       ! no hopping possible
       new_state = -1
       phase = 0
    end if
  end subroutine new_state_and_phase

!******************************************

  integer function statesindex(ne,selstate)
    integer :: ne
    integer :: selstate
    integer :: i,j,j_lo,j_hi
    integer :: n2site
    logical :: converged
    
!!$    Use a bisection method to determine the quantum number "statesindex"
!!$    associated with a particular state "selstate". 

    j_lo = block_(ne)
    j_hi = block_(ne)+nstates_(ne)-1

    if (selstate==states_(j_lo)) then
       statesindex = j_lo - block_(ne) + 1
       return
    else if (selstate==states_(j_hi)) then
       statesindex = j_hi - block_(ne) + 1
       return
    end if

    converged = .false.
    do while (.not.converged)
       j = (j_lo+j_hi)/2
       if (states_(j)>selstate) then
          j_hi = j
       else if (states_(j)<selstate) then
          j_lo = j
       else if (states_(j)==selstate) then
          statesindex= j - block_(ne) + 1
          converged = .true.
       else
          print *,"ERROR in statesindex"
          stop
       end if
    end do
  end function statesindex



!***********************************

  subroutine multiply_by_hamiltonian(v_in,v_out)
    !
    ! Multiply the vector v_in by the Hamiltonian.
    ! To save space, we don't actually store the Hamiltonian matrix;
    ! we have to work out the hopping on the fly.
    !
    real (real_kind), intent (in) :: v_in(n_)
    real (real_kind), intent (out) :: v_out(n_)
    
    integer :: istate_up,istate_dn,istate,jstate
    integer :: i,j
    integer :: state_up,state_dn
    integer :: nstates_up,nstates_dn
    real (real_kind) :: tmp1

    integer :: phase,new_state
    real :: t1,t2,t3

    ! We assume that the first half of each Fock state is spin-up, and the second half is spin down.

    nstates_up = nstates_(n_up_)
    nstates_dn = nstates_(n_dn_)
    call cpu_time(t1)
    ! Interaction & site energy terms.
    ! Loop over all Fock states making up v_in
    up1: do istate_up=1,nstates_up   ! Spin up Fock states
       state_up = states_(block_(n_up_)+istate_up-1)
       down1: do istate_dn = 1,nstates_dn   ! Spin down Fock states
          istate = istate_dn + (istate_up-1)*nstates_dn   ! This is the index (quantum #) of the state
          state_dn = states_(block_(n_dn_)+istate_dn-1)
          tmp1 = 0d0
          do i = 1,nsite_
             tmp1 = tmp1 + pars%U(i)*ibits(state_up,i-1,1)*ibits(state_dn,i-1,1)  &
                  + pars%Ei(i)*( ibits(state_up,i-1,1) + ibits(state_dn,i-1,1) )
          end do
          v_out(istate) = tmp1*v_in(istate)
       end do down1
    end do up1
    call cpu_time(t2)
    !
    ! SLOWEST PART OF THE CODE. Anything that can be done to speed this up is good.
    !
    ! Hopping term:
    up2: do istate_up=1,nstates_up   ! Spin up Fock states
       state_up = states_(block_(n_up_)+istate_up-1)
       istate = (istate_up-1)*nstates_dn   ! This is the index (quantum #) of the state
       down2: do istate_dn = 1,nstates_dn   ! Spin down Fock states
          istate = istate + 1
          state_dn = states_(block_(n_dn_)+istate_dn-1)
          ! The 0th element of state_lookup_XX_ contains the number of possible hops for this state
          ! Hop the up-spin part:
          hops2: do i = 1,state_lookup_up_(0,istate_up)  
             jstate = istate_dn + ( state_lookup_up_(i,istate_up) - 1 )*nstates_dn
             v_out(jstate) = v_out(jstate) - pars%t * phase_lookup_up_(i,istate_up) * v_in(istate)
          end do hops2
          ! Hop the down-spin part:
          hops3: do i = 1,state_lookup_dn_(0,istate_dn)
             jstate = state_lookup_dn_(i,istate_dn) + (istate_up-1)*nstates_dn
             v_out(jstate) = v_out(jstate) - pars%t * phase_lookup_dn_(i,istate_dn) * v_in(istate)
          end do hops3
       end do down2
    end do up2
    call cpu_time(t3)
    time1_ = time1_ + t2-t1
    time2_ = time2_ + t3-t2
  end subroutine multiply_by_hamiltonian
  
!***********************************

  subroutine get_static_correlations(gvec)
    real (real_kind), intent (in) :: gvec(n_)
    integer :: istate_up,istate_dn,state_up,state_dn,istate,isite
    integer :: nnn

    allocate(n_up_local_(nsite_),n_dn_local_(nsite_))
    allocate(spin_correlations_(nsite_))
    n_up_local_ = 0
    n_dn_local_ = 0
    spin_correlations_ = 0

    ! loop through all states
    do istate_up = 1,nstates_(n_up_)
       state_up = states_(block_(n_up_)+istate_up-1)
       do istate_dn = 1,nstates_(n_dn_)
          state_dn = states_(block_(n_dn_)+istate_dn-1)
          istate = istate_dn + (istate_up-1)*nstates_(n_dn_)

          ! loop through all sites
          do isite = 1,nsite_
             ! if isite has a spin up electron:
             if (btest(state_up,isite-1)) then
                n_up_local_(isite) = n_up_local_(isite) + gvec(istate)**2
             endif
             ! if isite has a spin down electron:
             if (btest(state_dn,isite-1)) then
                n_dn_local_(isite) = n_dn_local_(isite) + gvec(istate)**2
             endif
          end do
          ! Double occupancy of Cu-site at (ix,iy)=(1,1)
          spin_correlations_(1) = spin_correlations_(1) + gvec(istate)**2*ibits(state_up,0,1)*ibits(state_dn,0,1)  
          ! Nearest-neighbour spin correlation <s_z s_z> (Cu-px for 3band; Cu-Cu for 1band)
          spin_correlations_(2) = spin_correlations_(2) + 0.25d0*gvec(istate)**2&
               *(ibits(state_up,0,1)-ibits(state_dn,0,1))*(ibits(state_up,1,1)-ibits(state_dn,1,1))
          ! Next-Nearest-neighbour spin correlation <s_z s_z> (Cu-Cu for 3band and 1band)
          if (lattice_==3) then
             nnn=h_index_(h_index_(1,1),2) ! nearest neighbour Cu to the first Cu, along x
          else
             nnn=h_index_(h_index_(1,1),3) ! nearest neighbour Cu to the first Cu, along diagonal
          end if
          spin_correlations_(3) = spin_correlations_(3) + 0.25d0*gvec(istate)**2&
               *(ibits(state_up,0,1)-ibits(state_dn,0,1))*(ibits(state_up,nnn-1,1)-ibits(state_dn,nnn-1,1))
          
       end do
    end do
    print *,"n_up_local_",n_up_local_
    print *,"n_dn_local_",n_dn_local_
    print *,"Double Occupancy",spin_correlations_(1)
    print *,"nearest neighbour <sz sz>",spin_correlations_(2)
    print *,"next nearest neighbour <sz sz>",spin_correlations_(3)
    print *,sum(gvec**2)

  end subroutine get_static_correlations

!***********************************

  subroutine make_excited_state(isite,icase,n_up_ground,n_dn_ground,GroundStateVec,ExcitedStateVec,Norm)
    integer, intent (in) :: isite,icase,n_up_ground,n_dn_ground
    real (real_kind), intent (in) :: GroundStateVec(:)
    real (real_kind), intent (out) :: ExcitedStateVec(:),Norm

    integer :: istate_up,istate_dn,istate_up_final,istate_dn_final,istate,istate_final
    integer :: state_up,state_dn,state_up_final,state_dn_final
    integer :: phase
    integer :: jsite
    
    print *,"Number of electrons in excited state: ",n_up_,n_dn_

    ExcitedStateVec = 0d0
    SELECT CASE (icase)
       CASE (1) ! Spin-up photoemission: remove spin-up electron          
          do istate_up = 1,nstates_(n_up_ground)
             state_up = states_(block_(n_up_ground)+istate_up-1)
             if (btest(state_up,isite-1)) then  ! remove an electron
                state_up_final = ibclr(state_up,isite-1)
                istate_up_final = statesindex(n_up_,state_up_final)
                phase = 1
                do jsite = isite+1,nsite_
                   if (btest(state_up,jsite-1)) phase = -phase
                end do
                do istate_dn = 1,nstates_(n_dn_ground)
                   istate = istate_dn + (istate_up-1)*nstates_(n_dn_ground)
                   istate_final = istate_dn + (istate_up_final-1)*nstates_(n_dn_)
                   ExcitedStateVec(istate_final) = phase*GroundStateVec(istate)
                   if (istate_final>size(ExcitedStateVec)) then
                      print *,"error",istate_final,size(ExcitedStateVec)
                   end if
                end do
             end if
          end do

       CASE (2) ! Spin-down photoemission: remove a spin-down electron          
          do istate_dn = 1,nstates_(n_dn_ground)
             state_dn = states_(block_(n_dn_ground)+istate_dn-1)
             if (btest(state_dn,isite-1)) then
                state_dn_final = ibclr(state_dn,isite-1)
                istate_dn_final = statesindex(n_dn_,state_dn_final)
                phase = (-1)**n_up_ground
                do jsite = isite+1,nsite_
                   if (btest(state_dn,jsite-1)) phase = -phase
                end do
                do istate_up = 1,nstates_(n_up_ground)
                   istate = istate_dn + (istate_up-1)*nstates_(n_dn_ground)
                   istate_final = istate_dn_final + (istate_up-1)*nstates_(n_dn_)
                   ExcitedStateVec(istate_final) = phase*GroundStateVec(istate)
                   if (istate_final>size(ExcitedStateVec)) then
                      print *,"error",istate_final,size(ExcitedStateVec)
                   end if
                end do
             end if
          end do

       CASE (3) ! Spin-up inverse photoemission: add spin-up electron          
          do istate_up = 1,nstates_(n_up_ground)
             state_up = states_(block_(n_up_ground)+istate_up-1)
             if (.not.btest(state_up,isite-1)) then  ! remove an electron
                state_up_final = ibset(state_up,isite-1)
                istate_up_final = statesindex(n_up_,state_up_final)
                phase = 1
                do jsite = isite+1,nsite_
                   if (btest(state_up,jsite-1)) phase = -phase
                end do
                do istate_dn = 1,nstates_(n_dn_ground)
                   istate = istate_dn + (istate_up-1)*nstates_(n_dn_ground)
                   istate_final = istate_dn + (istate_up_final-1)*nstates_(n_dn_)
                   ExcitedStateVec(istate_final) = phase*GroundStateVec(istate)
                   if (istate_final>size(ExcitedStateVec)) then
                      print *,"error",istate_final,size(ExcitedStateVec)
                   end if
                end do
             end if
          end do

       CASE (4) ! Spin-down inverse photoemission: add a spin-down electron          
          do istate_dn = 1,nstates_(n_dn_ground)
             state_dn = states_(block_(n_dn_ground)+istate_dn-1)
             if (.not.btest(state_dn,isite-1)) then
                state_dn_final = ibset(state_dn,isite-1)
                istate_dn_final = statesindex(n_dn_,state_dn_final)
                phase = (-1)**n_up_ground
                do jsite = isite+1,nsite_
                   if (btest(state_dn,jsite-1)) phase = -phase
                end do
                do istate_up = 1,nstates_(n_up_ground)
                   istate = istate_dn + (istate_up-1)*nstates_(n_dn_ground)
                   istate_final = istate_dn_final + (istate_up-1)*nstates_(n_dn_)
                   ExcitedStateVec(istate_final) = phase*GroundStateVec(istate)
                   if (istate_final>size(ExcitedStateVec)) then
                      print *,"error",istate_final,size(ExcitedStateVec)
                   end if
                end do
             end if
          end do
       END SELECT

       Norm = sqrt(sum(ExcitedStateVec**2))
       ExcitedStateVec = ExcitedStateVec/Norm

  end subroutine make_excited_state

  !**************************

  character(16) function convert_to_binary(a)
    implicit none
    integer :: a
    character (16) :: ac
    integer :: i
    
    ac=''
    do i=0,15
       if (btest(a,i)) then
          ac = '1'//trim(ac)
       else
          ac = '0'//trim(ac)
       end if
    end do
    
    convert_to_binary=ac
  end function convert_to_binary
  

end module Hilbert_space_routines

