program MAIN
  use timemachine
  use Hilbert_space_routines
  use parameters
  use simple_lanczos
  implicit none
  integer :: i,j,iconfig
  integer :: n
  integer :: n_up,n_dn,ne  !  up, down, and total number of electrons
  real (real_kind) :: GroundStateEnergy
  real (real_kind), allocatable :: GroundStateVec(:)
  !---------- Time machine ----------
  integer :: dd,hh,mm,ss,mss


  !---------------------
  ! Initialise everything
  !---------------------
  time1_ = 0.0
  time2_ = 0.0
  time3_ = 0.0
  call random_gen()
  call read_parameters()
  call evaluate_parameters()
  call make_hilbert_space(pars%nx,pars%ny,pars%lattice)
  !call print_hilbert_space()
  call print_parameters()
!!$  n_dn = nint(pars%nsite/3*2.5)
!!$  n_up = nint(pars%nsite/3*2.5)
  n_dn=9; n_up=9
!!$
  call num_states(n_up,n_dn,n)
  allocate(GroundStateVec(n))
  print *,"# number of sites:",pars%nsite
  print *,"# number of up/dn electrons, number of states:",n_up,n_dn,n
  call make_h_index(pars%lattice)  ! Index of the columns of the Hamiltonian, stored in h_index.  
                                   ! This is where the details of the lattice are stored.
                                   ! Everything else in this code is generic.
  call initialise_lookup_tables()
  call make_state_lookup_table()

  !---------------------
  ! Main program
  !---------------------
  call time_starter()
  do iconfig = 1,pars%iconfig
     call filename_parameter(iconfig)
     call pars_Ei()
     print *,"# site energies: ",pars%Ei     

     ! Find the ground state energy:
     call restarted_lanczos(n,GroundStateVec,GroundStateEnergy)  
     call time_elapsed(hh,mm,ss,mss)
     write(*,200) n_up+n_dn,0,n,0,GroundStateEnergy,mm,ss,mss,iteration_lanczos()
     
     open(999,file=pars%file1)
     write(999,'("n_up, n_dn:  ",i3,1x,i3)') n_up,n_dn
     write(999,'("Site_energies: ",20(g16.9,1x))') pars%Ei
     write(999,'("Energy:  ",g16.9)') GroundStateEnergy
     write(999,*) "Eigenvector:"
     write(999,'(g16.9)') GroundStateVec
     close(999)  

     print *,"time doing diagonal multiply",time1_
     print *,"time doing off-diagonal multiply",time2_

     ! Get whatever static correlations we want from the wavefunction
     call get_static_correlations(GroundStateVec)
     call delete_lanczos()
  end do

  !---------------------
  ! Finalize everything
  !---------------------
  call cleanup_parameters()
  call cleanup_hilbert_space()
  call deallocate_lookup_tables()
  deallocate(groundstatevec)
  STOP
200  format(:,"# ne=",i3,1x,"sz=",i4,1x,"n=",i7,1x,"(",i8,")",1x, &
          "E=",f11.6,1x,"time=",i2,":",i2,".",i3,1x,"it="i3)
end program MAIN
