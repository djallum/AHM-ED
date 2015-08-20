module simple_recursion_routines

  use math_setup
  use Hilbert_space_routines
  use parameters

  implicit none

  real (real_kind), parameter, private :: abstol_ = 1d-4
  complex (real_kind), parameter, private :: igam_ = (0d0,1d-1)
  integer, parameter :: nw_=401
  real (real_kind), parameter :: wmin_=-20.0, wmax_=20.0, dw_=(wmax_-wmin_)/(nw_-1)
  real (real_kind), save :: spectrum_(nw_,4)

contains
  
  !**************************************

  subroutine simple_recursion(q0,maxit,E0,icase,A,B,k)
    integer, intent (in) :: maxit,icase
    integer, intent (out) :: k
    real (real_kind), intent (in) :: E0
    real (real_kind), intent (inout), target :: q0(:)
    real (real_kind), intent (out) :: A(:),B(:)

    integer :: n
    real (real_kind), allocatable :: r(:)
    real (real_kind), allocatable, target :: q1(:)
    real (real_kind), pointer :: q_new(:),q_old(:)
    complex (real_kind) :: g0_old,g0_new,z

    integer :: ik
    logical :: converged

    n = size(q0)
    allocate(r(n),q1(n))
    
    converged = .false.
    q_new => q0
    q_old => q1
    g0_new = 0d0

    k=1
    call multiply_by_hamiltonian(q_new,r)    
    A(k) = dot_product(q_new,r)
    r = r - A(k)*q_new 
    B(k) = sqrt(dot_product(r,r))
    !print *,k,A(k),B(k)
    do while (.not.converged)
       k = k+1
       if (mod(k,2)==0) then
          q_new => q1; q_old => q0
       else
          q_new => q0; q_old => q1
       end if
       q_new = r/B(k-1)
       call multiply_by_hamiltonian(q_new,r)    
       A(k) = dot_product(q_new,r)
       r = r - A(k)*q_new - B(k-1)*q_old
       B(k) = sqrt(dot_product(r,r))
       if (mod(k,5)==0) then
          g0_old = g0_new
          if ((icase==1).or.(icase==2)) then
             z = igam_ - E0 + pars%mu
             g0_new = pars%mu + igam_ + a(k) 
             do ik = k-1,1,-1
                g0_new = pars%mu + igam_ + a(ik) - b(ik)**2/(g0_new)
             end do
          else if ((icase==3).or.(icase==4)) then
             z = igam_ + E0 + pars%mu
             g0_new = pars%mu + igam_ - a(k) 
             do ik = k-1,1,-1
                g0_new = pars%mu + igam_ - a(ik) - b(ik)**2/(g0_new)
             end do
          end if
          g0_new = 1d0/g0_new
          converged = ((abs(g0_new-g0_old)<abstol_) .and. (k>95))
       end if

       if (abs(B(k)) < abstol_) then
          converged = .true.
       else if (k==maxit) then
          converged = .true.
       end if          
       !print *,k,A(k),B(k),abs(g0_new-g0_old)
    end do
    deallocate(r,q1)
    nullify(q_new,q_old)
  end subroutine simple_recursion


  !**************************************
  subroutine initialise_spectrum()
    spectrum_ = 0d0
  end subroutine initialise_spectrum


  !**************************************

  subroutine get_spectrum(n,A,B,norm,icase,E0)
    integer, intent (in) :: n,icase
    real (real_kind), intent (in) :: A(:),B(:),norm,E0

    real (real_kind) :: Z(n,n)
    real (real_kind) :: work(2*n-2)
    integer :: info

    integer :: iw,in
    real (real_kind) :: w

    call DSTEV('v',n,A,B,Z,n, WORK, INFO )
    if (info/=0) then
       print *, "error in get spectrum",info
    end if
    print *,A(1:n)
    do in = 1,n
       if ((icase==1).or.(icase==2)) then
          iw = floor((E0-A(in)-pars%mu-wmin_)/dw_)+1
       else if ((icase==3).or.(icase==4)) then
          iw = floor((A(in)-E0-pars%mu-wmin_)/dw_)+1
       end if
       if ((iw>0).and.(iw<=nw_)) then
          spectrum_(iw,icase) = spectrum_(iw,icase) + Z(1,in)**2*norm**2/dw_
       end if
    end do
  end subroutine get_spectrum

  !**************************************

  subroutine get_spectrum_alt(n,A,B,norm,icase,E0)
    integer, intent (in) :: n,icase
    real (real_kind), intent (in) :: A(:),B(:),norm,E0

    integer :: iw,in
    complex :: z,g

    do iw = 1,nw_
       if ((icase==1).or.(icase==2)) then
          z = wmin_ + (iw-1)*dw_ + igam_ - E0 + pars%mu
          g = z + a(n)
          do in = n-1,1,-1
             g = z + a(in) - b(in)**2/g
          end do
       else if ((icase==3).or.(icase==4)) then
          z = wmin_ + (iw-1)*dw_ + igam_ + E0 + pars%mu
          g = z - a(n)
          do in = n-1,1,-1
             g = z - a(in) - b(in)**2/g
          end do
       end if
       g = 1d0/g
       spectrum_(iw,icase) = spectrum_(iw,icase)-aimag(g)*norm**2/pi_
    end do
  end subroutine get_spectrum_alt

  !************************************

  subroutine print_spectrum(output_file)
    integer :: iw,in         ! counter used for making the different frequency bins
    real (real_kind) :: w    ! the frequency value of each bin will be assigned to this during loop
    character(32) :: output_file  ! the name of the file being outputed. Specified in the parameters module and printed to terminal
    
    open(10,file=output_file)
    do iw = 1,nw_
       w = wmin_ + (iw-1)*dw_
       write(10,100) w,spectrum_(iw,:)/pars%nsite/pars%iconfig    ! this normalizes the values
    end do
    close(10)
100 format(6(1x,g13.6))
  end subroutine print_spectrum

end module simple_recursion_routines
