module parameters
  use math_setup
  use random_generator
  implicit none

  type :: list_parameters
     !########## INPUT ##########
     ! Lattice Structure
     integer :: nx, ny, lattice   !nx is number of horizontal sites and ny is number of veriticle sites. Lattice specifies the type of structure. See subroutine make_hindex in hilbert space routines for the different types.
     ! Hamiltonian
     real :: t                   ! nearest neighbor hopping term
     real :: t1                  ! 2nd nearest neighbor hopping term
     real, dimension(:), pointer :: U  ! on-site repulsion
     real :: V                   ! nearest neighbor interaction
     real :: mu                  ! chemical potential
     ! disorder
     real :: delta               ! Box distribution
     integer :: Nimp             ! Number of impurity sites for binary alloy disorder
     real :: Vimp                ! Impurity potential for alloy sites
     ! configuration
     integer :: iconfig          ! Amount of iterations to do (the ensemble)?
     ! Output filenames
     character (len=32) :: file1,file2

     !########## DEPENDENT ##########
     integer :: nsite
     integer :: ensemble
     real, dimension(:), pointer :: Ei
  end type list_parameters

  type :: dependent_parameters
     integer :: n        ! total number of electrons
     integer :: ne
     integer :: sz
     integer :: nonzero
     real :: E
     real, dimension(:), pointer :: vec
     integer :: it
  end type dependent_parameters

  !---------- Variables ----------
  type (list_parameters) :: pars
  type (dependent_parameters) :: ground

contains

  subroutine read_parameters()
    allocate(pars%U(1))                          ! make the variable U be an array of length 1.
    open(100,file="input.dat",status="old")      ! open the input file called input.dat
    read(100,*) pars%nx,pars%ny,pars%lattice     ! read the x and y dimension of the lattice and the shape
    read(100,*) pars%t,pars%t1,pars%U(1),pars%V,pars%mu,pars%delta  ! read all the parameters needed to make the lattice
    read(100,*) pars%iconfig                     ! the amount of configurations you want to do. Used in calculation of the size of the ensemble.
    read(100,*) pars%Nimp,pars%Vimp              ! 
    close(100)
    if ((pars%lattice==1).and.(pars%nx/=4)) then             ! check if it was a valid lattice type
       print *,"Warning: MUST use a 4x3 lattice for 12 sites (cannot use 3x4)"
       STOP
    end if
    ! Lanczos
    ! Filename
    pars%file1 = ""      ! make the file names contain and empty string. There name will be assignend in evaluate_parameters subroutine
    pars%file2 = ""
    ! Lattice type
    !call choose_lattice(pars%lattice)
    !call choose_lattice(-1) ! linear lattice
    !call choose_lattice(0) ! square/rectangular
    !call choose_lattice(2) ! 10-site
    !call choose_lattice(1) ! 12-site
    !call choose_lattice(3) ! 3-band Hubbard
  end subroutine read_parameters

  !---------- Evaluate dependent parameters ----------

  subroutine evaluate_parameters()
    integer :: ict
    real(kind=8) :: U_tmp

    U_tmp=pars%U(1)                          ! % sign means the element called U from the variable pars. (Pars is of type list parameter which is a user defined data type). The type is defined at the top of this file.
    if (pars%lattice<3) then                 ! if lattice < 3 then it is a single band hubbard model
       pars%nsite = pars%nx * pars%ny        ! the number of sites is equal to the x dimension times the y dimension (it's a square)
       deallocate(pars%U)                    ! deallocate the array U so that the size can be respecified. Previous value is stored as U_temp
       allocate(pars%U(pars%nsite),pars%Ei(pars%nsite)) ! allocate the array U (on-site interactions) and Ei(site potentials) to be of size nsite.
       pars%U(:) = U_tmp                     ! set all the values in the array U equal to the value specified in the input file
    else if (pars%lattice==3) then           ! if this is true then it is a three band hubbard model so it will need 3 times as many sites
       pars%nsite = pars%nx * pars%ny * 3    
       deallocate(pars%U)
       allocate(pars%U(pars%nsite),pars%Ei(pars%nsite))  ! reallocate everything to the proper size of the 3 band hubbard model
       do ict=0,pars%nx*pars%ny-1
          pars%U(ict*3+1) = U_tmp      ! every third onsite interaction energy is equal to that specified in the input file
          pars%U(ict*3+2) = 0d0        ! the rest are all equal to zero
          pars%U(ict*3+3) = 0d0
       end do
    end if
    print *,"U: ",U_tmp                ! print the values that you calculated to the screen to make sure the new value is same as the inputed one.
    print *,"U: ",pars%U

    if (pars%delta.eq.0) pars%iconfig = 1  ! if there is no disorder (delta=0) then all configurations will be the same so only need to do 1.
    pars%ensemble = pars%nsite * pars%iconfig  ! the size of the ensemble is the number of sites times the amount of configurations you want to do.
    call filename_parameter()      ! assign file names based on what the input values into the program are.

  end subroutine evaluate_parameters

  subroutine pars_Ei()
    integer :: kk
    integer, allocatable :: isite(:)
    real :: r
  
    if (pars%nsite.eq.2) then
       pars%Ei(1) =  0.5 * pars%delta
       pars%Ei(2) = -0.5 * pars%delta
    else if (pars%lattice < 3) then
       call random_number(pars%Ei)
       pars%Ei = pars%delta * (pars%Ei - 0.5)                 ! randomly assigns site potentials to all sites
    else if (pars%lattice == 3) then
       ! In this case, delta is the charge-transfer gap
       call random_number(pars%Ei)
       do kk=0,pars%nx*pars%ny-1
          pars%Ei(kk*3+1) = pars%delta
          pars%Ei(kk*3+2) = 0d0
          pars%Ei(kk*3+3) = 0d0
       end do
    end if

     if (pars%Nimp/=0) then                                 ! finding site to put the inpurities
        allocate(isite(pars%Nimp))
        do kk=1,pars%Nimp
           call random_number(r)
           isite(kk) = floor(r*pars%nsite) + 1                     ! assign impurity to randoms site between 1 and isite
           do while ( minval(abs(isite(kk)-isite(1:kk-1))) < 0.5 )   !make sure that site hasn't already been assigned an impurity
              call random_number(r)           
              isite(kk) = floor(r*pars%nsite) + 1                   ! try another site if site already has impurity
           end do
           pars%Ei(isite(kk)) = pars%Ei(isite(kk)) + pars%Vimp       ! add the inpurity potential to the impurity sites current potential
        end do
        print *,"impurity locations:  ",isite                   ! print the locations of the impurities
        deallocate(isite)
     end if
  end subroutine pars_Ei

  subroutine print_parameters()
    write(*,100)
    write(*,210) pars%nx,pars%ny,pars%lattice
    write(*,220) pars%t,pars%t1
    write(*,230) pars%U(1),pars%V,pars%delta,pars%mu
    write(*,270) pars%iconfig
    write(*,280) pars%file1
    write(*,280) pars%file2
    write(*,100)

100 format("#")
210 format("# nx=",i2,1x,"ny=",i2,1x,"lattice=",i2)
220 format("# t=",f5.1,1x,"t1=",f5.1)
230 format("# U=",f5.1,1x,"V=",f5.1,1x,"Delta=",f5.1,1x,"mu=",f6.2)
240 format("# wi=",f6.1,1x,"wf=",f5.1,1x,"step=",f6.3,1x,"number of step=",i5)
250 format("# gamma=",f5.2)
260 format("# block_size=",i3)
270 format("# disorder configurations=",i5)
280 format("# output file=",1x,a)
  end subroutine print_parameters

  subroutine cleanup_parameters()
    deallocate(pars%U)
    deallocate(pars%Ei)
  end subroutine cleanup_parameters

  subroutine filename_parameter(iconfig)
    character(len=4) :: config,ui,uf,vi,vf,di,df,mi,mf,nsite,ic
    integer, optional :: iconfig
    integer :: u1,u2,v1,v2,d1,d2,m1,m2
    call r2i(pars%U(1),u1,u2)
    call r2i(pars%V,v1,v2)
    call r2i(pars%Delta,d1,d2)
    call r2i(pars%mu,m1,m2)
    write(config,'(i4)') pars%iconfig
    write(ui,'(i4)') u1
    write(uf,'(i4)') u2
    write(vi,'(i4)') v1
    write(vf,'(i4)') v2
    write(di,'(i4)') d1
    write(df,'(i4)') d2
    write(mi,'(i4)') m1
    write(mf,'(i4)') m2
    write(nsite,'(i4)') pars%nsite
    pars%file1 = "GS"//trim(adjustl(nsite))//"u"//trim(adjustl(ui))
    if (u2.ne.0) pars%file1 = trim(pars%file1)//"."//trim(adjustl(uf))
    pars%file1 = trim(pars%file1)//"v"//trim(adjustl(vi))
    if (v2.ne.0) pars%file1 = trim(pars%file1)//"."//trim(adjustl(vf))
    pars%file1 = trim(pars%file1)//"d"//trim(adjustl(di))
    if (d2.ne.0) pars%file1 = trim(pars%file1)//"."//trim(adjustl(df))
    pars%file1 = trim(pars%file1)//"m"//trim(adjustl(mi))
    if (m2.ne.0) pars%file1 = trim(pars%file1)//"."//trim(adjustl(mf))
    if (present(iconfig)) then
       write(ic,'(i4)') iconfig
       if (iconfig.ne.0) pars%file1 = trim(pars%file1)//"c"//trim(adjustl(ic))
    end if
    pars%file1 = trim(pars%file1)//".dat"

    pars%file2 = "spectrum"//trim(adjustl(nsite))//"u"//trim(adjustl(ui))
    if (u2.ne.0) pars%file2 = trim(pars%file2)//"."//trim(adjustl(uf))
    pars%file2 = trim(pars%file2)//"v"//trim(adjustl(vi))
    if (v2.ne.0) pars%file2 = trim(pars%file2)//"."//trim(adjustl(vf))
    pars%file2 = trim(pars%file2)//"d"//trim(adjustl(di))
    if (d2.ne.0) pars%file2 = trim(pars%file2)//"."//trim(adjustl(df))
    pars%file2 = trim(pars%file2)//"m"//trim(adjustl(mi))
    if (m2.ne.0) pars%file2 = trim(pars%file2)//"."//trim(adjustl(mf))
    if (present(iconfig)) then
       if (iconfig.ne.0) pars%file2 = trim(pars%file1)//"c"//trim(adjustl(ic))
    end if
    pars%file2 = trim(pars%file2)//".dat"
    contains
      subroutine r2i(a,i1,i2)
        real :: a
        integer :: i1,i2
        real :: b,c
        i1 = int(a)
        b = a - i1 + 0.0000005
        c = b
        do while (c.gt.0.001)
           b = 10 * b
           i2 = int(b)
           c = b - i2
        end do
        i2 = int(b)
      end subroutine r2i
  end subroutine filename_parameter

end module parameters
