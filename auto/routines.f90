module routines

	use lapack

	implicit none

	integer, parameter :: sp = kind(1.0)      !single precison kind
	integer, parameter :: dp = kind(1.0d0)    !double precision kind
	real :: grand_potential_ground=0.0                 ! the lowest grand ensemble energy
	integer, parameter :: nsites = 8
	integer, parameter :: total_states = 65536
	integer, parameter :: int_kind = 4
	integer, parameter :: STDOUT = 10

	integer, dimension(2,total_states) :: fock_states
	integer, dimension(0:2**nsites) :: states_order
	integer :: ne
	integer :: total_states_up
	integer, dimension(nsites,total_states) :: PES_down=0, PES_up=0, IPES_down=0, IPES_up=0  !matrices for PES and IPES 
	integer, dimension(nsites,total_states) :: phase_PES_down=0, phase_PES_up=0, phase_IPES_down=0, phase_IPES_up=0  !to get anticommutation sign right
	
	real :: H00(1,1),H01(8,8),H02(28,28),H03(56,56),H04(70,70),H05(56,56),H06(28,28),H07(8,8),H08(1,1)
	real :: H10(8,8),H11(64,64),H12(224,224),H13(448,448),H14(560,560),H15(448,448),H16(224,224),H17(64,64),H18(8,8)
	real :: H20(28,28),H21(224,224),H22(784,784),H23(1568,1568),H24(1960,1960),H25(1568,1568),H26(784,784),H27(224,224),H28(28,28)
	real :: H30(56,56),H31(448,448),H32(1568,1568),H33(3136,3136),H34(3920,3920),H35(3136,3136),H36(1568,1568),H37(448,448),H38(56,56)
	real :: H40(70,70),H41(560,560),H42(1960,1960),H43(3920,3920),H44(4900,4900),H45(3920,3920),H46(1960,1960),H47(560,560),H48(70,70)
	real :: H50(56,56),H51(448,448),H52(1568,1568),H53(3136,3136),H54(3920,3920),H55(3136,3136),H56(1568,1568),H57(448,448),H58(56,56)
	real :: H60(28,28),H61(224,224),H62(784,784),H63(1568,1568),H64(1960,1960),H65(1568,1568),H66(784,784),H67(224,224),H68(28,28)
	real :: H70(8,8),H71(64,64),H72(224,224),H73(448,448),H74(560,560),H75(448,448),H76(224,224),H77(64,64),H78(8,8)
	real :: H80(1,1),H81(8,8),H82(28,28),H83(56,56),H84(70,70),H85(56,56),H86(28,28),H87(8,8),H88(1,1)

	real, dimension(0:nsites,0:nsites) :: e_ground
	real, dimension(total_states) :: grand_potential      ! grand potentials (eigenenergies - mu*number electrons)
	real, dimension(total_states,4900) :: eigenvectors                  ! the eigenvectors
	integer, dimension(0:nsites) :: block, temp_block
	integer, dimension(0:nsites) :: nstates_up
	integer, allocatable, dimension(:,:) :: neighbours
	integer, dimension(0:nsites,0:nsites) :: msize, mblock

contains

	!------------------------------------------------------

	subroutine make_filename(filename,t,U,mu,delta)

		implicit none

		character(len=50), intent(out) :: filename
		character(len=10) :: mu_str, t_str, U_str, W_str, s_str
		real, intent(in) :: t,U,mu,delta

		write(mu_str,'(F4.1)') mu
		write(t_str,'(I2)') nint(t)
		write(W_str,'(I2)') nint(delta) 
		write(U_str,'(I2)') nint(U)
		write(s_str,'(I1)') nsites

		write(filename,'(A)') trim(adjustl('auto')) // trim(adjustl(s_str)) 
		write(filename,'(A)') trim(adjustl(filename)) // '_t' // trim(adjustl(t_str)) 
		write(filename,'(A)') trim(adjustl(filename)) // 'U' // trim(adjustl(U_str)) 
		write(filename,'(A)') trim(adjustl(filename)) // 'W' // trim(adjustl(W_str)) 
		write(filename,'(A)') trim(adjustl(filename)) // 'mu' // trim(adjustl(mu_str)) // '.dat'

	end subroutine make_filename

	!------------------------------------------------------

	subroutine random_gen_seed()

	  implicit none

	  integer :: seed_size, clock, i
	  integer, allocatable, dimension(:) :: seed

	  call random_seed(size=seed_size)
	  allocate(seed(seed_size))
	  call system_clock(count = clock)
	  seed=clock + 37*(/(i-1,i=1,seed_size)/)
	  seed=3
	  call random_seed(put=seed)
	  deallocate(seed)

	end subroutine random_gen_seed

	!------------------------------------------------------

	subroutine site_potentials(delta,E)

	  implicit none

	  real, intent(in) :: delta
	  real, intent(out) :: E(nsites)
	  real :: random
	  integer :: i ! counter

	  do i=1,nsites
	     call random_number(random)              ! gives a random number between 0 and 1
	     E(i) = delta*(random - 0.5)          ! centers the random numbers about 0
	  end do

	end subroutine site_potentials

	!------------------------------------------------------

	subroutine num_sites()

		implicit none

		integer :: i,j, istate, isite

		integer :: max_electrons

		total_states_up = 2**nsites
		max_electrons = nsites

		block = 0  
	    nstates_up = 0
	    block(0) = 1
	    nstates_up(0) = 1

		do ne = 1, max_electrons
       		nstates_up(ne) = choose(nsites,ne)          ! number of spin-up fock states with ne electrons
       		block(ne) = block(ne-1) + nstates_up(ne-1) ! block index updating
    	end do

    	temp_block = block

		do istate=0,total_states_up-1
			ne = 0
			do isite=1,nsites
				ne = ne + ibits(istate,isite-1,1)     ! count the number of electrons in that state
			end do 
			states_order(temp_block(ne)) = istate
			temp_block(ne) = temp_block(ne) + 1
		end do

		istate = 1
		do ne=0,max_electrons
			do j=1,total_states_up
				do i=block(ne), block(ne) + nstates_up(ne) - 1
					fock_states(1,istate) = states_order(i)
					fock_states(2,istate) = states_order(j)
					istate = istate + 1
				end do
			end do
		end do

	end subroutine num_sites

	!------------------------------------------------------

	subroutine transformations()

		implicit none
		
		integer :: i,j,new_state(2), new_index, position, isite

		PES_up = 0; PES_down = 0
	    IPES_up = 0; IPES_down = 0
	    phase_PES_up = 0; phase_PES_down = 0
	    phase_IPES_up = 0; phase_IPES_down = 0

		! make the PES_up matrices

		do position = 1,nsites
			do i=1,total_states
				ne = 0
				if (ibits(fock_states(1,i),position-1,1) == 1) then
					PES_up(position,i) = 0
					phase_PES_up(position,i) = 0 
				else
					do isite=position,nsites
						ne = ne + ibits(fock_states(1,i),isite-1,1)     ! count the number of up electrons in that state
						end do
					do isite=position,nsites
						ne = ne + ibits(fock_states(2,i),isite-1,1)     ! count the number of down electrons in that state
					end do
					if (MOD(ne,2) == 0) then 
						phase_PES_up(position,i) = 1
					else
						phase_PES_up(position,i) = -1
					end if
					new_state(1) = ibset(fock_states(1,i),position-1)
					new_state(2) = fock_states(2,i)
					do j=1,total_states
						if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
							new_index = j
						end if
					end do
					PES_up(position,i) = new_index
				end if
			end do
		end do

		do position = 1,nsites
			do i=1,total_states
				ne = 0
				if (ibits(fock_states(2,i),position-1,1) == 1) then
					PES_down(position,i) = 0
					phase_PES_down(position,i) = 0 
				else
					do isite=position+1,nsites
						ne = ne + ibits(fock_states(1,i),isite-1,1)     ! count the number of up electrons in that state
					end do
					do isite=position,nsites
						ne = ne + ibits(fock_states(2,i),isite-1,1)     ! count the number of down electrons in that state
					end do
					if (MOD(ne,2) == 0) then 
						phase_PES_down(position,i) = 1
					else
						phase_PES_down(position,i) = -1
					end if
					new_state(2) = ibset(fock_states(2,i),position-1)
					new_state(1) = fock_states(1,i)
					do j=1,total_states
						if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
							new_index = j
						end if
					end do
					PES_down(position,i) = new_index
				end if
			end do
		end do

		sites: do j=1,nsites  ! calculating the IPES matrices by making them the opposite of the PES
         do i=1,total_states
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

	!-------------------------------------------------------

	subroutine make_neighbours()

	! makes matrix that containes information about which sites are nearest neighbours
	! neighbours(i,:) is a list of all the neighbours of site i. Each site has 4 nearest neighbours normally.

		if (nsites == 8) then                  ! betts lattice for 8 sites
			allocate(neighbours(nsites,4))
			neighbours(1,1) = 8; neighbours(1,2) = 7; neighbours(1,3) = 5; neighbours(1,4) = 3
			neighbours(2,1) = 7; neighbours(2,2) = 8; neighbours(2,3) = 3; neighbours(2,4) = 5
			neighbours(3,1) = 1; neighbours(3,2) = 2; neighbours(3,3) = 4; neighbours(3,4) = 6
			neighbours(4,1) = 5; neighbours(4,2) = 3; neighbours(4,3) = 8; neighbours(4,4) = 7
			neighbours(5,1) = 2; neighbours(5,2) = 1; neighbours(5,3) = 6; neighbours(5,4) = 4
			neighbours(6,1) = 3; neighbours(6,2) = 5; neighbours(6,3) = 7; neighbours(6,4) = 8
			neighbours(7,1) = 4; neighbours(7,2) = 6; neighbours(7,3) = 1; neighbours(7,4) = 2
			neighbours(8,1) = 6; neighbours(8,2) = 4; neighbours(8,3) = 2; neighbours(8,4) = 1
		end if

		if (nsites == 4) then                  ! linear 4 site lattice
			allocate(neighbours(nsites,2))
			neighbours(1,1) = 4; neighbours(1,2) = 2
			neighbours(2,1) = 1; neighbours(2,2) = 3
			neighbours(3,1) = 2; neighbours(3,2) = 4
			neighbours(4,1) = 3; neighbours(4,2) = 1
		end if

		if (nsites == 2) then                 ! linear 2 site lattice
			allocate(neighbours(nsites,1))
			neighbours(1,1) = 2; neighbours(2,1) = 1
		end if

	end subroutine make_neighbours

	!-------------------------------------------------------

	subroutine matrix_sizes()

		implicit none

		integer :: n_up,n_dn

		do n_up=0,nsites
			do n_dn=0,nsites
				msize(n_up,n_dn) = choose(nsites,n_up)*choose(nsites,n_dn)
				if (n_dn == 0 .and. n_up == 0) then
					mblock(n_up,n_dn) = 1
				else if (n_dn == 0) then
					mblock(n_up,n_dn) = mblock(n_up-1,nsites) + msize(n_up-1,nsites)
				else 
					mblock(n_up,n_dn) = mblock(n_up,n_dn-1) + msize(n_up,n_dn-1)
				end if
			end do
		end do

	end subroutine matrix_sizes

	!-------------------------------------------------------

	subroutine make_hamiltonian2(t)

		! this program makes the hamiltonians off diagonal terms. On diagonal terms are added during each loop

		implicit none

		real, intent(in) :: t

		integer :: istate,isite, inbr,new_state(2),phase, ne, i, trans_site(2), new_index,j,state_index,y,n_up,n_dn
		
		H00=0.0; H01=0.0; H02=0.0; H03=0.0; H04=0.0; H05=0.0; H06=0.0; H07=0.0; H08=0.0; 
		H10=0.0; H11=0.0; H12=0.0; H13=0.0; H14=0.0; H15=0.0; H16=0.0; H17=0.0; H18=0.0; 
		H20=0.0; H21=0.0; H22=0.0; H23=0.0; H24=0.0; H25=0.0; H26=0.0; H27=0.0; H28=0.0; 
		H30=0.0; H31=0.0; H32=0.0; H33=0.0; H34=0.0; H35=0.0; H36=0.0; H37=0.0; H38=0.0; 
		H40=0.0; H41=0.0; H42=0.0; H43=0.0; H44=0.0; H45=0.0; H46=0.0; H47=0.0; H48=0.0; 
		H50=0.0; H51=0.0; H52=0.0; H53=0.0; H54=0.0; H55=0.0; H56=0.0; H57=0.0; H58=0.0; 
		H60=0.0; H61=0.0; H62=0.0; H63=0.0; H64=0.0; H65=0.0; H66=0.0; H67=0.0; H68=0.0; 
		H70=0.0; H71=0.0; H72=0.0; H73=0.0; H74=0.0; H75=0.0; H76=0.0; H77=0.0; H78=0.0; 
		H80=0.0; H81=0.0; H82=0.0; H83=0.0; H84=0.0; H85=0.0; H86=0.0; H87=0.0; H88=0.0; 

		call matrix_sizes()

		do n_up = 0,nsites
		do n_dn = 0,nsites
		do istate = mblock(n_up,n_dn),mblock(n_up,n_dn) + msize(n_up,n_dn)-1
			do isite = 1,nsites
				if (ibits(fock_states(1,istate),isite-1,1) == 1) then
					do y=1,size(neighbours,2)
						new_state(1) = IBCLR(fock_states(1,istate),isite-1)
						inbr = neighbours(isite,y)
						if (ibits(new_state(1),inbr-1,1) == 0) then
							new_state(2) = fock_states(2,istate)
							ne = 0
							trans_site(1) = inbr; trans_site(2) = isite
							do i=MINVAL(trans_site),MAXVAL(trans_site)-1
								ne = ne + ibits(new_state(1),i-1,1)     ! count the number of up electrons in that state
								ne = ne + ibits(new_state(2),i-1,1)     ! count the number of down electrons in that state
							end do
							if (MOD(ne,2) == 0) then 
								phase = 1
							else
								phase = -1
							end if
							new_state(1) = ibset(new_state(1),inbr-1)
							do j=1,total_states
								if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
								new_index = j
								end if
							end do
							state_index = istate + 1 - mblock(n_up,n_dn)
							new_index = new_index + 1 - mblock(n_up,n_dn)
							select case (n_up)
								case(0)
									select case (n_dn)
										case(0)
											H00(state_index,new_index) = t*phase
										case(1)
											H01(state_index,new_index) = t*phase
										case(2)
											H02(state_index,new_index) = t*phase
										case(3)
											H03(state_index,new_index) = t*phase
										case(4)
											H04(state_index,new_index) = t*phase
										case(5)
											H05(state_index,new_index) = t*phase
										case(6)
											H06(state_index,new_index) = t*phase
										case(7)
											H07(state_index,new_index) = t*phase
										case(8)
											H08(state_index,new_index) = t*phase
									end select
								case(1)
									select case (n_dn)
										case(0)
											H10(state_index,new_index) = t*phase
										case(1)
											H11(state_index,new_index) = t*phase
										case(2)
											H12(state_index,new_index) = t*phase
										case(3)
											H13(state_index,new_index) = t*phase
										case(4)
											H14(state_index,new_index) = t*phase
										case(5)
											H15(state_index,new_index) = t*phase
										case(6)
											H16(state_index,new_index) = t*phase
										case(7)
											H17(state_index,new_index) = t*phase
										case(8)
											H18(state_index,new_index) = t*phase
									end select
								case(2)
									select case (n_dn)
										case(0)
											H20(state_index,new_index) = t*phase
										case(1)
											H21(state_index,new_index) = t*phase
										case(2)
											H22(state_index,new_index) = t*phase
										case(3)
											H23(state_index,new_index) = t*phase
										case(4)
											H24(state_index,new_index) = t*phase
										case(5)
											H25(state_index,new_index) = t*phase
										case(6)
											H26(state_index,new_index) = t*phase
										case(7)
											H27(state_index,new_index) = t*phase
										case(8)
											H28(state_index,new_index) = t*phase
									end select
								case(3)
									select case (n_dn)
										case(0)
											H30(state_index,new_index) = t*phase
										case(1)
											H31(state_index,new_index) = t*phase
										case(2)
											H32(state_index,new_index) = t*phase
										case(3)
											H33(state_index,new_index) = t*phase
										case(4)
											H34(state_index,new_index) = t*phase
										case(5)
											H35(state_index,new_index) = t*phase
										case(6)
											H36(state_index,new_index) = t*phase
										case(7)
											H37(state_index,new_index) = t*phase
										case(8)
											H38(state_index,new_index) = t*phase
									end select
								case(4)
									select case (n_dn)
										case(0)
											H40(state_index,new_index) = t*phase
										case(1)
											H41(state_index,new_index) = t*phase
										case(2)
											H42(state_index,new_index) = t*phase
										case(3)
											H43(state_index,new_index) = t*phase
										case(4)
											H44(state_index,new_index) = t*phase
										case(5)
											H45(state_index,new_index) = t*phase
										case(6)
											H46(state_index,new_index) = t*phase
										case(7)
											H47(state_index,new_index) = t*phase
										case(8)
											H48(state_index,new_index) = t*phase
									end select
								case(5)
									select case (n_dn)
										case(0)
											H50(state_index,new_index) = t*phase
										case(1)
											H51(state_index,new_index) = t*phase
										case(2)
											H52(state_index,new_index) = t*phase
										case(3)
											H53(state_index,new_index) = t*phase
										case(4)
											H54(state_index,new_index) = t*phase
										case(5)
											H55(state_index,new_index) = t*phase
										case(6)
											H56(state_index,new_index) = t*phase
										case(7)
											H57(state_index,new_index) = t*phase
										case(8)
											H58(state_index,new_index) = t*phase
									end select
								case(6)
									select case (n_dn)
										case(0)
											H60(state_index,new_index) = t*phase
										case(1)
											H61(state_index,new_index) = t*phase
										case(2)
											H62(state_index,new_index) = t*phase
										case(3)
											H63(state_index,new_index) = t*phase
										case(4)
											H64(state_index,new_index) = t*phase
										case(5)
											H65(state_index,new_index) = t*phase
										case(6)
											H66(state_index,new_index) = t*phase
										case(7)
											H67(state_index,new_index) = t*phase
										case(8)
											H68(state_index,new_index) = t*phase
									end select
								case(7)
									select case (n_dn)
										case(0)
											H70(state_index,new_index) = t*phase
										case(1)
											H71(state_index,new_index) = t*phase
										case(2)
											H72(state_index,new_index) = t*phase
										case(3)
											H73(state_index,new_index) = t*phase
										case(4)
											H74(state_index,new_index) = t*phase
										case(5)
											H75(state_index,new_index) = t*phase
										case(6)
											H76(state_index,new_index) = t*phase
										case(7)
											H77(state_index,new_index) = t*phase
										case(8)
											H78(state_index,new_index) = t*phase
									end select
								case(8)
									select case (n_dn)
										case(0)
											H80(state_index,new_index) = t*phase
										case(1)
											H81(state_index,new_index) = t*phase
										case(2)
											H82(state_index,new_index) = t*phase
										case(3)
											H83(state_index,new_index) = t*phase
										case(4)
											H84(state_index,new_index) = t*phase
										case(5)
											H85(state_index,new_index) = t*phase
										case(6)
											H86(state_index,new_index) = t*phase
										case(7)
											H87(state_index,new_index) = t*phase
										case(8)
											H88(state_index,new_index) = t*phase
									end select
							end select
						end if
					end do
				end if
				if (ibits(fock_states(2,istate),isite-1,1) == 1) then
					do y=1,size(neighbours,2)
						new_state(2) = IBCLR(fock_states(2,istate),isite-1)
						inbr = neighbours(isite,y)
						if (ibits(new_state(2),inbr-1,1) == 0) then
							new_state(1) = fock_states(1,istate)
							ne = 0
							trans_site(1) = inbr; trans_site(2) = isite
							do i=MINVAL(trans_site)+1,MAXVAL(trans_site)
								ne = ne + ibits(new_state(1),i-1,1)     ! count the number of up electrons in that state
								ne = ne + ibits(new_state(2),i-1,1)     ! count the number of down electrons in that state
							end do
							if (MOD(ne,2) == 0) then 
								phase = 1
							else
								phase = -1
							end if
							new_state(2) = ibset(new_state(2),inbr-1)
							do j=1,total_states
								if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
								new_index = j
								end if
							end do
							state_index = istate + 1 - mblock(n_up,n_dn)
							new_index = new_index + 1 - mblock(n_up,n_dn)
							select case (n_up)
								case(0)
									select case (n_dn)
										case(0)
											H00(state_index,new_index) = t*phase
										case(1)
											H01(state_index,new_index) = t*phase
										case(2)
											H02(state_index,new_index) = t*phase
										case(3)
											H03(state_index,new_index) = t*phase
										case(4)
											H04(state_index,new_index) = t*phase
										case(5)
											H05(state_index,new_index) = t*phase
										case(6)
											H06(state_index,new_index) = t*phase
										case(7)
											H07(state_index,new_index) = t*phase
										case(8)
											H08(state_index,new_index) = t*phase
									end select
								case(1)
									select case (n_dn)
										case(0)
											H10(state_index,new_index) = t*phase
										case(1)
											H11(state_index,new_index) = t*phase
										case(2)
											H12(state_index,new_index) = t*phase
										case(3)
											H13(state_index,new_index) = t*phase
										case(4)
											H14(state_index,new_index) = t*phase
										case(5)
											H15(state_index,new_index) = t*phase
										case(6)
											H16(state_index,new_index) = t*phase
										case(7)
											H17(state_index,new_index) = t*phase
										case(8)
											H18(state_index,new_index) = t*phase
									end select
								case(2)
									select case (n_dn)
										case(0)
											H20(state_index,new_index) = t*phase
										case(1)
											H21(state_index,new_index) = t*phase
										case(2)
											H22(state_index,new_index) = t*phase
										case(3)
											H23(state_index,new_index) = t*phase
										case(4)
											H24(state_index,new_index) = t*phase
										case(5)
											H25(state_index,new_index) = t*phase
										case(6)
											H26(state_index,new_index) = t*phase
										case(7)
											H27(state_index,new_index) = t*phase
										case(8)
											H28(state_index,new_index) = t*phase
									end select
								case(3)
									select case (n_dn)
										case(0)
											H30(state_index,new_index) = t*phase
										case(1)
											H31(state_index,new_index) = t*phase
										case(2)
											H32(state_index,new_index) = t*phase
										case(3)
											H33(state_index,new_index) = t*phase
										case(4)
											H34(state_index,new_index) = t*phase
										case(5)
											H35(state_index,new_index) = t*phase
										case(6)
											H36(state_index,new_index) = t*phase
										case(7)
											H37(state_index,new_index) = t*phase
										case(8)
											H38(state_index,new_index) = t*phase
									end select
								case(4)
									select case (n_dn)
										case(0)
											H40(state_index,new_index) = t*phase
										case(1)
											H41(state_index,new_index) = t*phase
										case(2)
											H42(state_index,new_index) = t*phase
										case(3)
											H43(state_index,new_index) = t*phase
										case(4)
											H44(state_index,new_index) = t*phase
										case(5)
											H45(state_index,new_index) = t*phase
										case(6)
											H46(state_index,new_index) = t*phase
										case(7)
											H47(state_index,new_index) = t*phase
										case(8)
											H48(state_index,new_index) = t*phase
									end select
								case(5)
									select case (n_dn)
										case(0)
											H50(state_index,new_index) = t*phase
										case(1)
											H51(state_index,new_index) = t*phase
										case(2)
											H52(state_index,new_index) = t*phase
										case(3)
											H53(state_index,new_index) = t*phase
										case(4)
											H54(state_index,new_index) = t*phase
										case(5)
											H55(state_index,new_index) = t*phase
										case(6)
											H56(state_index,new_index) = t*phase
										case(7)
											H57(state_index,new_index) = t*phase
										case(8)
											H58(state_index,new_index) = t*phase
									end select
								case(6)
									select case (n_dn)
										case(0)
											H60(state_index,new_index) = t*phase
										case(1)
											H61(state_index,new_index) = t*phase
										case(2)
											H62(state_index,new_index) = t*phase
										case(3)
											H63(state_index,new_index) = t*phase
										case(4)
											H64(state_index,new_index) = t*phase
										case(5)
											H65(state_index,new_index) = t*phase
										case(6)
											H66(state_index,new_index) = t*phase
										case(7)
											H67(state_index,new_index) = t*phase
										case(8)
											H68(state_index,new_index) = t*phase
									end select
								case(7)
									select case (n_dn)
										case(0)
											H70(state_index,new_index) = t*phase
										case(1)
											H71(state_index,new_index) = t*phase
										case(2)
											H72(state_index,new_index) = t*phase
										case(3)
											H73(state_index,new_index) = t*phase
										case(4)
											H74(state_index,new_index) = t*phase
										case(5)
											H75(state_index,new_index) = t*phase
										case(6)
											H76(state_index,new_index) = t*phase
										case(7)
											H77(state_index,new_index) = t*phase
										case(8)
											H78(state_index,new_index) = t*phase
									end select
								case(8)
									select case (n_dn)
										case(0)
											H80(state_index,new_index) = t*phase
										case(1)
											H81(state_index,new_index) = t*phase
										case(2)
											H82(state_index,new_index) = t*phase
										case(3)
											H83(state_index,new_index) = t*phase
										case(4)
											H84(state_index,new_index) = t*phase
										case(5)
											H85(state_index,new_index) = t*phase
										case(6)
											H86(state_index,new_index) = t*phase
										case(7)
											H87(state_index,new_index) = t*phase
										case(8)
											H88(state_index,new_index) = t*phase
									end select
							end select
						end if
					end do
				end if
			end do
		end do
		end do
		end do


	end subroutine make_hamiltonian2

	!-------------------------------------------------------

	subroutine solve_hamiltonian1(E,U,mu)

		! this program makes the on diagonal terms and solves teh eigenenergies

		implicit none

		integer :: n_up, n_dn,i,j,ne,isite,istate

		real, intent(in) :: E(nsites)
  		real, intent(in) :: mu 
  		real, intent(in) :: U
 		real, allocatable, dimension(:) :: wtemp
 		real, allocatable, dimension(:,:) :: htemp, vtemp

		do n_up=0,nsites
			do n_dn=0,nsites
				
				allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))
				allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
				
				select case (n_up)
					case(0)
						select case (n_dn)
							case(0)
								htemp=H00
							case(1)
								htemp=H01
							case(2)
								htemp=H02
							case(3)
								htemp=H03
							case(4)
								htemp=H04
							case(5)
								htemp=H05
							case(6)
								htemp=H06
							case(7)
								htemp=H07
							case(8)
								htemp=H08
						end select
					case(1)
						select case (n_dn)
							case(0)
								htemp=H10
							case(1)
								htemp=H11
							case(2)
								htemp=H12
							case(3)
								htemp=H13
							case(4)
								htemp=H14
							case(5)
								htemp=H15
							case(6)
								htemp=H16
							case(7)
								htemp=H17
							case(8)
								htemp=H18
						end select
					case(2)
						select case (n_dn)
							case(0)
								htemp=H20
							case(1)
								htemp=H21
							case(2)
								htemp=H22
							case(3)
								htemp=H23
							case(4)
								htemp=H24
							case(5)
								htemp=H25
							case(6)
								htemp=H26
							case(7)
								htemp=H27
							case(8)
								htemp=H28
						end select
					case(3)
						select case (n_dn)
							case(0)
								htemp=H30
							case(1)
								htemp=H31
							case(2)
								htemp=H32
							case(3)
								htemp=H33
							case(4)
								htemp=H34
							case(5)
								htemp=H35
							case(6)
								htemp=H36
							case(7)
								htemp=H37
							case(8)
								htemp=H38
						end select
					case(4)
						select case (n_dn)
							case(0)
								htemp=H40
							case(1)
								htemp=H41
							case(2)
								htemp=H42
							case(3)
								htemp=H43
							case(4)
								htemp=H44
							case(5)
								htemp=H45
							case(6)
								htemp=H46
							case(7)
								htemp=H47
							case(8)
								htemp=H48
						end select
					case(5)
						select case (n_dn)
							case(0)
								htemp=H50
							case(1)
								htemp=H51
							case(2)
								htemp=H52
							case(3)
								htemp=H53
							case(4)
								htemp=H54
							case(5)
								htemp=H55
							case(6)
								htemp=H56
							case(7)
								htemp=H57
							case(8)
								htemp=H58
						end select
					case(6)
						select case (n_dn)
							case(0)
								htemp=H60
							case(1)
								htemp=H61
							case(2)
								htemp=H62
							case(3)
								htemp=H63
							case(4)
								htemp=H64
							case(5)
								htemp=H65
							case(6)
								htemp=H66
							case(7)
								htemp=H67
							case(8)
								htemp=H68
						end select
					case(7)
						select case (n_dn)
							case(0)
								htemp=H70
							case(1)
								htemp=H71
							case(2)
								htemp=H72
							case(3)
								htemp=H73
							case(4)
								htemp=H74
							case(5)
								htemp=H75
							case(6)
								htemp=H76
							case(7)
								htemp=H77
							case(8)
								htemp=H78
						end select
					case(8)
						select case (n_dn)
							case(0)
								htemp=H80
							case(1)
								htemp=H81
							case(2)
								htemp=H82
							case(3)
								htemp=H83
							case(4)
								htemp=H84
							case(5)
								htemp=H85
							case(6)
								htemp=H86
							case(7)
								htemp=H87
							case(8)
								htemp=H88
						end select
				end select

				wtemp=0.0

				do istate=1,msize(n_up,n_dn)
					do isite=1,nsites
						ne=0
						ne = ne + IBITS(fock_states(1,istate+mblock(n_up,n_dn)-1),isite-1,1)
						ne = ne + IBITS(fock_states(2,istate+mblock(n_up,n_dn)-1),isite-1,1)
						htemp(istate,istate) = htemp(istate,istate) + ne*E(isite)
						if (ne == 2) then 
							htemp(istate,istate) = htemp(istate,istate) + U
						end if
					end do
				end do

				call ssyevr_lapack1(msize(n_up,n_dn),htemp,wtemp,vtemp)

				e_ground(n_up,n_dn) = wtemp(1) - mu*(n_up+n_dn)  ! grand potentials

				deallocate(htemp,wtemp,vtemp)

			end do
		end do

	end subroutine solve_hamiltonian1

	!-------------------------------------------------------

	subroutine solve_hamiltonian2(E,U,mu,n_up,n_dn)

		! this program makes the on diagonal terms and solves it

		implicit none

		integer :: i,j,ne,isite,istate

		real, intent(in) :: E(nsites)
  		real, intent(in) :: mu 
  		real, intent(in) :: U
  		integer, intent(in) :: n_up,n_dn
 		real, allocatable, dimension(:) :: wtemp
 		real, allocatable, dimension(:,:) :: htemp, vtemp
				
		allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))
		allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
		
		select case (n_up)
			case(0)
				select case (n_dn)
					case(0)
						htemp=H00
					case(1)
						htemp=H01
					case(2)
						htemp=H02
					case(3)
						htemp=H03
					case(4)
						htemp=H04
					case(5)
						htemp=H05
					case(6)
						htemp=H06
					case(7)
						htemp=H07
					case(8)
						htemp=H08
				end select
			case(1)
				select case (n_dn)
					case(0)
						htemp=H10
					case(1)
						htemp=H11
					case(2)
						htemp=H12
					case(3)
						htemp=H13
					case(4)
						htemp=H14
					case(5)
						htemp=H15
					case(6)
						htemp=H16
					case(7)
						htemp=H17
					case(8)
						htemp=H18
				end select
			case(2)
				select case (n_dn)
					case(0)
						htemp=H20
					case(1)
						htemp=H21
					case(2)
						htemp=H22
					case(3)
						htemp=H23
					case(4)
						htemp=H24
					case(5)
						htemp=H25
					case(6)
						htemp=H26
					case(7)
						htemp=H27
					case(8)
						htemp=H28
				end select
			case(3)
				select case (n_dn)
					case(0)
						htemp=H30
					case(1)
						htemp=H31
					case(2)
						htemp=H32
					case(3)
						htemp=H33
					case(4)
						htemp=H34
					case(5)
						htemp=H35
					case(6)
						htemp=H36
					case(7)
						htemp=H37
					case(8)
						htemp=H38
				end select
			case(4)
				select case (n_dn)
					case(0)
						htemp=H40
					case(1)
						htemp=H41
					case(2)
						htemp=H42
					case(3)
						htemp=H43
					case(4)
						htemp=H44
					case(5)
						htemp=H45
					case(6)
						htemp=H46
					case(7)
						htemp=H47
					case(8)
						htemp=H48
				end select
			case(5)
				select case (n_dn)
					case(0)
						htemp=H50
					case(1)
						htemp=H51
					case(2)
						htemp=H52
					case(3)
						htemp=H53
					case(4)
						htemp=H54
					case(5)
						htemp=H55
					case(6)
						htemp=H56
					case(7)
						htemp=H57
					case(8)
						htemp=H58
				end select
			case(6)
				select case (n_dn)
					case(0)
						htemp=H60
					case(1)
						htemp=H61
					case(2)
						htemp=H62
					case(3)
						htemp=H63
					case(4)
						htemp=H64
					case(5)
						htemp=H65
					case(6)
						htemp=H66
					case(7)
						htemp=H67
					case(8)
						htemp=H68
				end select
			case(7)
				select case (n_dn)
					case(0)
						htemp=H70
					case(1)
						htemp=H71
					case(2)
						htemp=H72
					case(3)
						htemp=H73
					case(4)
						htemp=H74
					case(5)
						htemp=H75
					case(6)
						htemp=H76
					case(7)
						htemp=H77
					case(8)
						htemp=H78
				end select
			case(8)
				select case (n_dn)
					case(0)
						htemp=H80
					case(1)
						htemp=H81
					case(2)
						htemp=H82
					case(3)
						htemp=H83
					case(4)
						htemp=H84
					case(5)
						htemp=H85
					case(6)
						htemp=H86
					case(7)
						htemp=H87
					case(8)
						htemp=H88
				end select
		end select

		wtemp=0.0

		do istate=1,msize(n_up,n_dn)
			do isite=1,nsites
				ne=0
				ne = ne + IBITS(fock_states(1,istate+mblock(n_up,n_dn)-1),isite-1,1)
				ne = ne + IBITS(fock_states(2,istate+mblock(n_up,n_dn)-1),isite-1,1)
				htemp(istate,istate) = htemp(istate,istate) + ne*E(isite)
				if (ne == 2) then 
					htemp(istate,istate) = htemp(istate,istate) + U
				end if
			end do
		end do

		call ssyevr_lapack(msize(n_up,n_dn),htemp,wtemp,vtemp)

		do i=1,msize(n_up,n_dn)
	   		grand_potential(i+mblock(n_up,n_dn)-1) = wtemp(i) - mu*(n_up+n_dn)  ! grand potentials
		end do

		do i=1,msize(n_up,n_dn)
	   		eigenvectors(i+mblock(n_up,n_dn)-1,1:msize(n_up,n_dn)) = vtemp(1:msize(n_up,n_dn),i)  ! eigenvectors
		end do

		deallocate(htemp,wtemp,vtemp)

	end subroutine solve_hamiltonian2

	!-------------------------------------------------------

	integer function choose(j,k)

	    implicit none
	    
	    integer :: j,k,i2
	    integer (kind=8) :: tj,tk,tchoose

	    tchoose = 1
	    tj = j
	    tk = k
	    do i2 = tj-tk+1,tj
	       tchoose = tchoose * i2
	    end do
	    do i2 = 2, tk
	       tchoose = tchoose / i2
	    end do
	    choose = tchoose

	end function choose

	!------------------------------------------------------

end module
