program MAIN
  
  !The following are modules (other files) that are used. The modules themselves contain the subroutines and functions.

  !-------------------------------
  use timemachine
  use Hilbert_space_routines
  use parameters
  use simple_lanczos
  !-------------------------------

  implicit none

  integer :: i,j     ! counters for loops
  integer :: iconfig
  integer :: n
  integer :: n_up,n_dn,ne                  !  number of up, down, and total electrons
  real (real_kind) :: GroundStateEnergy
  real (real_kind), allocatable :: GroundStateVec(:)
  !---------- Time machine ----------
  integer :: dd,hh,mm,ss,mss    !these variables will represent the amount of time elapsed


  !---------------------
  ! Initialise everything
  !---------------------

  time1_ = 0.0
  time2_ = 0.0
  time3_ = 0.0
  call random_gen()                                         ! in module found in random_generator.f90
  call read_parameters()                                    ! in the parameter module found in parameters.f90 will read parameters from file input.dat
  call evaluate_parameters()                                ! in the parameter module found in parameters.f90 will calculate the dependant parameters
  call make_hilbert_space(pars%nx,pars%ny,pars%lattice)     ! makes ordered list of fock states and lookup arrays for when states with n electrons begin. 
  !call print_hilbert_space() 
  call print_parameters()                                   ! in the parameter module found in parameters.f90
  !!$  n_dn = nint(pars%nsite/3*2.5)
  !!$  n_up = nint(pars%nsite/3*2.5)
  n_dn=2; n_up=2                   ! specify that there are 2 up electrons and 2 down electrons
  !!$
  call num_states(n_up,n_dn,n)    ! calculates the total number of fock states possible with n_dn and n_up electrons(in )
  allocate(GroundStateVec(n))      
  print *,"# number of sites:",pars%nsite
  print *,"# number of up/dn electrons, number of states:",n_up,n_dn,n
  call make_h_index(pars%lattice)  ! Index of the columns of the Hamiltonian, stored in h_index.  
                                   ! This is where the details of the lattice are stored.
                                   ! Everything else in this code is generic. h_index(i,1) is site number of first nieghbour site i can hoop to.
  call initialise_lookup_tables()
  call make_state_lookup_table()   ! finds all the states each fock state can jump to (nearest neihgbourgh only)

  !---------------------
  ! Main program
  !---------------------

  call time_starter() ! timer starts here 
  do iconfig = 1,pars%iconfig
     call filename_parameter(iconfig)
     call pars_Ei()                               ! randomly assigns site potentials
     print *,"# site energies: ",pars%Ei     

     ! Find the ground state energy:
     call restarted_lanczos(n,GroundStateVec,GroundStateEnergy)  !n = number of fock states with the asigned amount of electrons
     call time_elapsed(hh,mm,ss,mss) ! timer ends here                    
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
