program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/06/10            created code 
! P. Daley              15/08/14            change to single precision
! P. Daley              15/08/17            added comments
!
! This code find the density of states (DOS) and generalized inverse participation ratio (GIPR) for an ensemble of 4 site systems. It uses the Anderson-Hubbarb model.
! The hamiltonian matrix is hard coded and the eigenvectors and eigenvalues are found using LAPACK driver 'ssyev'. The lowest grand potential is found and the assosiated
! eigenvector is the lowest many body eigenstate (Psi0). The c^{dagger}_{i,sigma}|Psi0> and c_{i,sigma}|Psi0> vectors are calculated then their inner product is taken with each
! of the other many body eigenstates. This gives the local density of states (LDOS). The LDOS is used the calculate the GIPR and DOS. Final steps normalizes the DOS to 1
! and prints all the data to the output file. Basis and hamiltonians are defined in the excel files.

  use routines

  implicit none

  !-------------------Input Parameters---------------------------------
  integer, parameter :: npairs=100000       ! size of the ensemble
  real :: t = -1.0                          ! nearest neighbour hopping 
  real, parameter :: delta = 12.0           ! width of disorder for the site potentials 
  real, parameter :: U = 4.0                ! on-site interactions
  real, parameter :: mu=U/2                 ! chemical potential (half filling)
  integer, parameter :: nbins = 300         ! number of bins for energy bining to get smooth curves 
  real, parameter :: frequency_max = 12     ! maximum energy considered in energy bining
  real, parameter :: frequency_min = -12    ! lowest energy considered in energy bining

  !-------------------Dependent Variables------------------------------
  integer :: pair, i, j, istate, isite  ! counters for loops
  integer :: error                      ! variable for error message
  integer :: location(1)                ! stores the location in the grand_potential array of the lowest energy 
  integer :: bin                        ! index for the bin number the peak goes in
  real, dimension(nsites) :: E          ! site potentials
  real, dimension(nsites,nstates) :: PES_dn_ground, PES_up_ground   ! ground state vector after an down or up photo emmision respectively c_{i}|Psi_{0}>
  real, dimension(nsites,nstates) :: IPES_dn_ground, IPES_up_ground ! ground state vector after an down or up inverse photo emmision respectively c_{i}^{dagger}|Psi_{0}>
  real, dimension(nsites,2*nstates,2) :: LDOS  ! Local DOS (LDOS(isite,istate,1) is energy of istate's contribution to LDOS of isite and LDOS(i,isite,2) is the weight of the contribution)
  real :: inner_product_up, inner_product_dn ! used in the calculation of the weight of LDOS contributions
  real :: IPR(2*nstates)                       ! contains the Inverse Participation Ratio (IPR) contributions (non weighted)
  real :: frequency_delta                      ! step size between different energy bins
  real, dimension(nbins,2) :: DOS              ! array that stores the DOS peaks and all the frequencies of the energy bins (DOS(i,1) is frequency of bin i and DOS(i,2) is the weight)
  real, dimension(nbins) :: GIPR_num, GIPR_den, GIPR  ! arrays that store the numerator and denominator and the weighted IPR
  real :: dos_sum                              ! sums the entire DOS so that it can be normalized to 1
  real :: half_sum                             ! sums half of the DOS to find the filling 
  character(len=50) :: filename

  !-----------------Preparations for the main loop----------------------

  frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

  call random_gen_seed()
  call transformations()    ! calls subroutines in routines.f90. It makes the lookup tables for PES and IPES.

  call make_filename(filename,t,U,mu,delta)   ! calls subroutine in routines.f90 that assigns a filename based on the input parameters

  open(unit=10,file=filename, status='replace', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
  if (error/=0) then                                                           ! check if there was error when opening the file
     write(*,*) 'error opening output file. Error number:', error              ! if there was an error print the error message
  end if                                                                       ! exit the program 

  ! write informartion about the code to the top of the output file
  write(10,*) "# 4 site DOS and GIPR created by main.f90 with subroutines in routines.f90"
  write(10,*) "#"
  write(10,'(A18,F6.2,2(A5,F6.2))') " # parameters: U =",U,"t =",t,"W =",delta
  write(10,'(A10,I4,A20,I12)') " # nbins =",nbins, "ensemble size =", npairs 
  write(10,*) "#"

  !---------------------------Main loop-----------------------------------

  pairs: do pair=1,npairs
   
    if (npairs < 100) then     ! this section prints the percentage of completion if npairs > 100 and the pair number if npairs < 100
      write(*,*) pair
    else 
      if (MOD(pair,npairs/100) == 0) then
        write(*,'(I3,A1)') nint(real(pair)/npairs*100), "%"
      end if
    end if

    ! zero variables for each loop    
    v_ground = 0.0
    eigenvectors = 0.0
    grand_potential_ground = 0.0
    grand_potential = 0.0
    LDOS = 0.0
    inner_product_up=0.0
    inner_product_dn=0.0
    PES_dn_ground=0.0; PES_up_ground=0.0; IPES_dn_ground=0.0; IPES_up_ground=0.0

    call site_potentials(delta,E)     ! call subroutine to assign the site potentials
    call hamiltonian(E,t,U,mu)        ! make and solve the hamiltonians for their eigenvectors and eigenvalues

    !-----find the ground state grand potential------------------------

    grand_potential_ground = minval(grand_potential)   ! find the lowest grand ensemble energy

    !-----find the corresponding eigenvector----------------

    location = minloc(grand_potential)          ! find the location of the lowest energy  
    v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

    !-----multiply v_ground by the PES and IPES matrices-------------

    do isite=1,nsites
          do istate=1,nstates
              if (PES_up(isite,istate)==0) then
                PES_up_ground(isite,istate) = 0.0
              else 
                PES_up_ground(isite,istate) = v_ground(PES_up(isite,istate))*phase_PES_up(isite,istate)
              end if
              if (PES_dn(isite,istate)==0) then
                PES_dn_ground(isite,istate) = 0.0
              else 
                PES_dn_ground(isite,istate) = v_ground(PES_dn(isite,istate))*phase_PES_dn(isite,istate)
              end if
              if (IPES_up(isite,istate)==0) then
                IPES_up_ground(isite,istate) = 0.0
              else 
                IPES_up_ground(isite,istate) = v_ground(IPES_up(isite,istate))*phase_IPES_up(isite,istate)
              end if
              if (IPES_dn(isite,istate)==0) then
                IPES_dn_ground(isite,istate) = 0.0
              else 
                IPES_dn_ground(isite,istate) = v_ground(IPES_dn(isite,istate))*phase_IPES_dn(isite,istate)
              end if
          end do
      end do 

      ! calculate the LDOS for each site

      do isite=1,nsites
          do istate=1,nstates
              inner_product_up = (dot_product(PES_up_ground(isite,:),eigenvectors(istate,:)))**2   ! dot product the PES up ground vector with each eigenstate
              inner_product_dn =  (dot_product(PES_dn_ground(isite,:),eigenvectors(istate,:)))**2  ! dot product the PES down ground vector with each eigenstate
              LDOS(isite,istate,1) = grand_potential_ground - grand_potential(istate)              ! location of the contribution
              LDOS(isite,istate,2) = (inner_product_up + inner_product_dn)*0.5                     ! weight of the contribution (average up and down spin components)
          end do
      end do

      do isite=1,nsites
          do istate=1,nstates
              inner_product_up = (dot_product(IPES_up_ground(isite,:),eigenvectors(istate,:)))**2
              inner_product_dn =  (dot_product(IPES_dn_ground(isite,:),eigenvectors(istate,:)))**2
              LDOS(isite,istate+nstates,1) = grand_potential(istate) - grand_potential_ground       ! location of the contribution
              LDOS(isite,istate+nstates,2) = (inner_product_up + inner_product_dn)*0.5              ! weight of the contribution (average up and down spin components)
          end do
      end do

      do istate=1,2*nstates
          bin = floor(LDOS(1,istate,1)/frequency_delta) + nbins/2  + 1              ! find the bin number for energy bining
          DOS(bin,2) = DOS(bin,2) + (SUM(LDOS(:,istate,2)))/real(nsites)            ! make the contribution to that bin
          if (SUM(LDOS(:,istate,2)) /= 0) then
              IPR(istate) = SUM(LDOS(:,istate,2)**2)/(SUM(LDOS(:,istate,2))**2)
              GIPR_num(bin) = GIPR_num(bin) + IPR(istate)*(SUM(LDOS(:,istate,2)))/real(nsites)  ! numerator of the weighted GIPR
              GIPR_den(bin) = GIPR_den(bin) + (SUM(LDOS(:,istate,2)))/real(nsites)              ! denominator of the weighted GIPR
          end if
      end do

  end do pairs

  !-----------------Normalize DOS and print to output file--------------------------

  dos_sum = DOS(1,2)                  ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  DOS(1,1) = frequency_min            ! set the first DOS energy to the mininum frequency

  do i=2,nbins                                  ! calculate sum to normalize the area under DOS to 1 and add frequencies to bins
    DOS(i,1) = DOS(i-1,1) + frequency_delta     ! step the energy by frequency_delta
    dos_sum = dos_sum + DOS(i,2)                ! add the bin to the total
  end do

  half_sum = DOS(1,2)                  ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  do i=2,nbins/2                       ! calculate half sum
    half_sum = half_sum + DOS(i,2)
  end do

  ! print additional information to the file
  write(10,*) "# Filling:", half_sum/dos_sum
  write(10,*) "# Filling Error:", DOS(nbins/2+1,2)/dos_sum
  write(10,*) "#"
  write(10,*) "#  Frequency           DOS            GIPR"

  do i=1,nbins
    GIPR(i) = GIPR_num(i)/GIPR_den(i)
    write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)  ! print final data to the output file
  end do

  close(10)

end program main
