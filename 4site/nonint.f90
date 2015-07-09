program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/07/08            created code 
! P. Daley              15/07/09            added comments to code
!
! This program calculates the density of states (DOS) for a 4site system with no on-site interactions (U=0). It finds the single particle transitions by solving the single particle hamiltonian. The contributions to the DOS occur at the eigenvalues of this matrix. The program then uses energy binning to produce a smooth DOS. The DOS is normalized when printed to the output file. The site potentials are chosen from a bounded distribution. The program needs the subroutines found in routines.f90 and the LAPACK library. 

  use routines

  implicit none

  integer, parameter :: npairs = 10000000   ! the number of pairs that will be in the ensemble
  real(dp) :: t = -1.0_dp                   ! hopping term (dp means double precision the declaration is in routines.f90)
  real(dp), parameter :: delta = 6.0_dp     ! width of disorder for the site potentials (W)
  real(dp), parameter :: mu = U/2           ! chemical potential (half filling) 
  integer, parameter :: nbins = 400         ! number of bins for energy bining to get smooth curves
  real, parameter :: frequency_max = 10     ! highest energy considered in energy bining
  real, parameter :: frequency_min = -10    ! lowest energy considered in energy bining

  real(dp) :: frequency_delta = 0.0_dp             ! step size between different energy bins (calculated from the max and min)
  real(dp), dimension(4) :: E = 0.0_dp             ! site potentials  
  real(dp), dimension(nbins,2) :: DOS = 0.0_dp     ! array that stores the DOS peaks and all the energy bin frequencies 
  real(dp), dimension(4,4) :: hamiltonian(4,4)     ! matrix for the single particle hamiltonian
  real(dp), dimension(4) :: eigenvalues = 0.0_dp   ! array outputed from LAPACK containing eigenvalues
  real(dp) :: sum=0.0_dp                           ! the area under the DOS (used to normalize graph to 1)
  integer :: bin=0                                 ! index for the bin number the contribution goes in
  integer :: error=0                               ! variable for error message in file opening
  integer :: pair=0,i=0,j=0, k=0                   ! counters
!--- ----For LAPACK---------------  
  real(dp) :: WORK(20) = 0.0_dp         ! array for LAPACK must be > (size of matrix) * 3
  integer :: LWORK = 20                 ! tells LAPACK the size of the matrix work
  integer :: INFO = 0                   ! holds the error message for the LAPACK program


  frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

  call random_gen_seed()

  open(unit=10,file='nonint_data.dat', status='replace', action='write',IOSTAT = error)  ! open the file that DOS will be printed to
  if (error/=0) then
     write(*,*) 'error opening output file. Error number:', error
  end if

  ! write informartion about the code to the output file
  write(10,*) "#created by nonint.f90 with subroutines in routines.f90"
  write(10,*) "#"
  write(10,*) "#parameters: t=",t,"W=",delta,"nbins=",nbins
  write(10,*) "#"
  write(10,*) "   #Frequency                    DOS"

  pairs: do pair=1,npairs    ! loop over each pair in the ensemble

     call site_potentials(delta,E)  ! assigns the site potentials

     ! contruct the hamiltonian matrix
     hamiltonian(1,1) = E(1);   hamiltonian(2,1) = t 
     hamiltonian(1,2) = t;      hamiltonian(2,2) = E(2)
     hamiltonian(1,3) = 0.0_dp; hamiltonian(2,3) = t
     hamiltonian(1,4) = t;      hamiltonian(2,4) = 0.0_dp
     hamiltonian(3,1) = 0.0_dp; hamiltonian(4,1) = t
     hamiltonian(3,2) = t;      hamiltonian(4,2) = 0.0_dp
     hamiltonian(3,3) = E(3);   hamiltonian(4,3) = t
     hamiltonian(3,4) = t;      hamiltonian(4,4) = E(4)

     call dsyev('v','u',4,hamiltonian,4,eigenvalues,WORK,LWORK,INFO)  ! LAPACK program to find the eigenvalues and eigenvector

     do i=1,4
        bin = floor(eigenvalues(i)/frequency_delta) + nbins/2  + 1  ! finding the bin each contribution should be made to
        DOS(bin,2) = DOS(bin,2) + 1                         ! making contribution to that bin
     end do

  end do pairs

  sum = DOS(1,2)                                  ! set sum equal to the amount of contributions in first bin (loop sums rest in next step)
  DOS(1,1) = frequency_min                        ! set the frequency of the first bin to the minimum

  do i=2,nbins
     DOS(i,1) = DOS(i-1,1) + frequency_delta      ! step each bin frequency value by the step size (frequency_delta)
     sum = sum + DOS(i,2)                         ! calculate sum to normalize the area under DOS to 1
  end do

  do i=1,nbins
     write(10,*), DOS(i,1), DOS(i,2)/sum/frequency_delta  ! write the information to the file (normalize the DOS)
  end do

  close(10)  ! close the output file


end program main
