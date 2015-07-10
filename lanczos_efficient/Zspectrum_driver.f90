module local_routines

  use math_setup
  use parameters
  use Hilbert_space_routines

  implicit none

  real (real_kind) :: GroundStateEnergy
  real (real_kind), allocatable :: GroundStateVec(:),ExcitedStateVec(:)

contains

  !*******************************************************  

  subroutine read_groundstate(input_file,n_up,n_dn,n)
    integer, intent (out) :: n_up,n_dn
    integer, intent (out) :: n
    character (5) :: jnk1,jnk2
    character (32) :: input_file

    integer :: istate

    open(99,file=input_file,status='old')
    read(99,*) jnk1,jnk2,n_up,n_dn
    call num_states(n_up,n_dn,n)
    read(99,*) jnk1,pars%Ei
    print *,"Ei:",pars%Ei
    read(99,*) jnk1,GroundStateEnergy
    read(99,*) jnk1

    allocate(GroundStateVec(n))
    do istate = 1,n
       read(99,*) GroundStateVec(istate)
    end do

    print *, "energy", GroundStateEnergy
    print *, "Vec", GroundStateVec(1:9)
    
    print *,"# number of up/dn electrons in groundstate:",n_up,n_dn

  end subroutine read_groundstate

!*******************************************************  

  subroutine get_final_numbers(spectrum_case,n_up,n_dn,n_up_final,n_dn_final)
    integer :: spectrum_case,n_up,n_dn,n_up_final,n_dn_final    

    SELECT CASE (spectrum_case)
    CASE(1) !PES up
       if (n_up>0) then
          n_up_final = n_up-1
          n_dn_final = n_dn
       else
          print *,"Zero up electrons in ground state: no photoemission spectrum"
          stop
       end if
    CASE(2) !PES down
       if (n_dn>0) then
          n_up_final = n_up
          n_dn_final = n_dn-1
       else
          print *,"Zero down electrons in ground state: no photoemission spectrum"
          stop
       end if
    CASE(3) !IPES up
       if (n_up<pars%nsite) then
          n_up_final = n_up+1
          n_dn_final = n_dn
       else
          print *,"Maximum number of up electrons in ground state: no inverse photoemission spectrum"
          stop
       end if
    CASE(4) !IPES down
       if (n_dn<pars%nsite) then
          n_up_final = n_up
          n_dn_final = n_dn+1
       else
          print *,"Maximum number of down electrons in ground state: no inverse photoemission spectrum"
          stop
       end if
    CASE DEFAULT
       print *,"Error in get_final_numbers: illegal value for spectrum_case.  Must be 1-4."
       stop
    end SELECT
  end subroutine get_final_numbers


end module local_routines

!*******************************************************  

program MAIN

  use local_routines
  use timemachine
  use simple_lanczos
  use simple_recursion_routines

  implicit none

  integer :: i,j,iconfig
  integer :: n_up,n_dn,n  !  up, down, and total number of states (ground state)
  integer :: n_up_final,n_dn_final,n_final  !  up, down, and total number of states (final state)
  real (real_kind) :: norm
  integer :: spectrum_case
  integer, parameter :: maxit=1000
  real (real_kind) :: A(maxit),B(maxit)
  integer :: iterations,isite
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
  call make_h_index(pars%lattice)  ! Index of the columns of the Hamiltonian, stored in h_index.  
  ! This is where the details of the lattice are stored.
  ! Everything else in this code is generic.
  call initialise_spectrum()
  do iconfig = 1,pars%iconfig
     call filename_parameter(iconfig)
     print *,"FILE:",pars%file1
     call read_groundstate(pars%file1,n_up,n_dn,n)


     !---------------------
     ! Main program
     !---------------------  
     do spectrum_case = 1,3,2  ! loop over PES,IPES for up & down 
        do isite = 1,pars%nsite
           print *, "config: ",iconfig,"isite: ",isite,"spectrum_case: ",spectrum_case
           call get_final_numbers(spectrum_case,n_up,n_dn,n_up_final,n_dn_final)
           call num_states(n_up_final,n_dn_final,n_final)   ! in hilbert space routines. Calculates the total amount of electrons from n_up and n_down
           allocate(ExcitedStateVec(n_final))
           call initialise_lookup_tables()
           call make_state_lookup_table()

           call make_excited_state(isite,spectrum_case,n_up,n_dn,GroundStateVec,ExcitedStateVec,Norm)
           call simple_recursion(ExcitedStateVec,maxit,GroundStateEnergy,spectrum_case,A,B,iterations)
           call get_spectrum_alt(iterations,A,B,Norm,spectrum_case,GroundStateEnergy)
           call deallocate_lookup_tables()
           deallocate(ExcitedStateVec)  
        end do
     end do
     deallocate(groundstatevec)
  end do
  !---------------------
  ! Finalize everything
  !---------------------
  call filename_parameter()      ! parameters module in parameters.f90 
  call print_spectrum(pars%file2)  !
  call cleanup_hilbert_space()
  STOP
200 format(:,"# ne=",i3,1x,"sz=",i4,1x,"n=",i7,1x,"(",i8,")",1x, &
       "E=",f11.6,1x,"time=",i2,":",i2,".",i3,1x,"it="i3)
end program MAIN

