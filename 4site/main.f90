program main

  use routines

  implicit none

  integer, parameter :: npairs=10000             ! size of the ensemble
  real(dp) :: t = -1.0_dp                      ! hopping term
  real(dp), parameter :: delta = 6.0_dp        ! width of disorder for the site potentials 
  real(dp), parameter :: U = 0.0_dp            ! on-site interactions
  real(dp), parameter :: mu=U/2                ! chemical potential (half filling) 
  real(dp) :: E(4)=0.0_dp                ! site potentials
  integer :: pair=0,i=0,j=0, k=0         ! counters
  integer :: error=0                     ! variable for error message
  integer :: location(1)=0               ! stores the location in the grand_potential array of the lowest energy 
  real(dp), dimension(4,256) :: PES_down_ground=0.0_dp, PES_up_ground=0.0_dp, IPES_down_ground=0.0_dp, IPES_up_ground=0.0_dp
  real(dp), dimension(4,512,2) :: LDOS=0.0_dp
  real(dp) :: inner_product_up=0.0_dp, inner_product_down=0.0_dp
  real(dp) :: IPR(512)=0.0_dp
  integer, parameter :: nbins = 300                  ! number of bins for energy bining to get smooth curves
  real, parameter :: frequency_max = 10              ! maximum energy considered in energy bining
  real, parameter :: frequency_min = -10             ! lowest energy considered in energy bining
  real(dp) :: frequency_delta=0.0_dp                 ! step size between different energy bins
  integer :: bin=0                                   ! index for the bin number the peak goes in
  real(dp), dimension(nbins,2) :: DOS=0.0_dp                                      ! array that stores the DOS peaks and all the frequencies of the energy bins
  real(dp), dimension(nbins) :: GIPR_num=0.0_dp, GIPR_den=0.0_dp, GIPR=0.0_dp     ! arrays that store the numerator and denominator and the GIPR
  real(dp) :: dos_sum=0.0_dp
  real(dp) :: eps(4)=0.0_dp

  frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

  call random_gen_seed()
  call transformations()

  open(unit=10,file='4site_data.dat', status='replace', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
  if (error/=0) then
     write(*,*) 'error opening output file. Error number:', error
  end if


  ! write informartion about the code to the top of the output file
  write(10,*) "# created by main.f90 with subroutines in routines.f90"
  write(10,*) "#"
  write(10,'(A17,F6.2,2(A5,F6.2))') "# parameters: U=",U,"t=",t,"W=",delta
  write(10,'(A9,I4,A20,I12)') "# nbins=",nbins, "ensemble size=", npairs 
  write(10,*) "#"
  write(10,*) "#     Frequency                      DOS                       GIPR" 

  pairs: do pair=1,npairs
   
    v_ground(256)=0.0_dp; eigenvectors(256,256)=0.0_dp; grand_potential_ground=0.0_dp; grand_potential(256)=0.0_dp
    eigenvectors = 0.0_dp
    grand_potential_ground = 0.0_dp
    grand_potential = 0.0_dp
    LDOS = 0.0_dp
    inner_product_up=0.0_dp
    inner_product_down=0.0_dp
    PES_down_ground=0.0_dp; PES_up_ground=0.0_dp; IPES_down_ground=0.0_dp; IPES_up_ground=0.0_dp

    call site_potentials(delta,E)
    call hamiltonian(E,t,U,mu,eps)

    !-----find ground state energy------------------------

    grand_potential_ground = minval(grand_potential)   ! find the lowest grand ensemble energy

    !-----find the corresponding eigenvector----------------

    location = minloc(grand_potential)          ! find the location of the lowest energy  
    v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

    !multiply ground state vector by the matrices

    do j=1,4
       do i=1,256
          if (PES_up(j,i)==0) then
             PES_up_ground(j,i) = 0.0_dp
          else 
             PES_up_ground(j,i) = v_ground(PES_up(j,i))*phase_PES_up(j,i)
          end if
          if (PES_down(j,i)==0) then
             PES_down_ground(j,i) = 0.0_dp
          else 
             PES_down_ground(j,i) = v_ground(PES_down(j,i))*phase_PES_down(j,i)
          end if
          if (IPES_up(j,i)==0) then
             IPES_up_ground(j,i) = 0.0_dp
          else 
             IPES_up_ground(j,i) = v_ground(IPES_up(j,i))*phase_IPES_up(j,i)
          end if
          if (IPES_down(j,i)==0) then
             IPES_down_ground(j,i) = 0.0_dp
          else 
             IPES_down_ground(j,i) = v_ground(IPES_down(j,i))*phase_IPES_down(j,i)
          end if
       end do
    end do 

    ! calculate the LDOS for all the cites

    do j=1,4
       do i=1,256
          inner_product_up = (dot_product(PES_up_ground(j,:),eigenvectors(i,:)))**2
          inner_product_down =  (dot_product(PES_down_ground(j,:),eigenvectors(i,:)))**2
          LDOS(j,i,1) = grand_potential_ground - grand_potential(i)           ! location of the peak
          LDOS(j,i,2) = (inner_product_up + inner_product_down)*0.5_dp           ! weight of the peak (average up and down spin components)
       end do
    end do

    do j=1,4
       do i=1,256
          inner_product_up = (dot_product(IPES_up_ground(j,:),eigenvectors(i,:)))**2
          inner_product_down =  (dot_product(IPES_down_ground(j,:),eigenvectors(i,:)))**2
          LDOS(j,i+256,1) = grand_potential(i) - grand_potential_ground           ! location of the peak
          LDOS(j,i+256,2) = (inner_product_up + inner_product_down)*0.5_dp        ! weight of the peak
       end do
    end do

    do i=1,512
       bin = floor(LDOS(2,i,1)/frequency_delta) + nbins/2  +1              !find the bin number for energy bining
       DOS(bin,2) = DOS(bin,2) + (SUM(LDOS(:,i,2)))/4.0_dp
        if (SUM(LDOS(:,i,2)) /= 0) then
          IPR(i) = SUM(LDOS(:,i,2)**2)/(SUM(LDOS(:,i,2))**2)
          GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(SUM(LDOS(:,i,2)))/4.0_dp  ! numerator of the weighted GIPR
          GIPR_den(bin) = GIPR_den(bin) + (SUM(LDOS(:,i,2)))/4.0_dp         ! denominator of the weighted GIPR
       end if
    end do

  end do pairs

  dos_sum = DOS(1,2)
  DOS(1,1) = frequency_min

  do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
     DOS(i,1) = DOS(i-1,1) + frequency_delta
     dos_sum = dos_sum + DOS(i,2)
  end do

  do i=1,nbins
    GIPR(i) = GIPR_num(i)/GIPR_den(i)
    if(DOS(i,2)/dos_sum/frequency_delta < 0.000000001) then
      GIPR(i) = 0.25
    end if
    write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)
  end do

  close(10)

end program main
