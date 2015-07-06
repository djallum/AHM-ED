module simple_lanczos
!  step 1: arbitrary choose |1^>
!  step 2: normalized |1^>
!          (1) b1 = <1^|1^>
!          (2) |1> = |1^> / b1
!  step 3:
!          (1) r = H|1>
!          (2) a1 = <1|H|1>
!          (3) r = H|1> - a1*|1>
!          (4) b1 = sqrt(r*r)
!          (5) |2> = r / b1
!  note: As b(j) = 0, it means that the subspace is completed.
  use math_setup
  use Hilbert_space_routines
  use lapack
  !use sparse_matrix
  implicit none
  private

  public :: parameters_lanczos
  public :: iteration_lanczos
  public :: groundVec_lanczos, eigenstates_lanczos
  public :: restarted_lanczos
  public :: delete_lanczos

  !---------- Static variables ----------
  integer :: max_iteration_ = 1001
  integer :: iteration_ = 5
  integer :: inc_step_ = 5
  real (real_kind) :: accuracy_ = 1.0d-6
  real (real_kind) :: cutoff_ = 1.0d-4
  !---------- Run-Time variables ----------
  integer :: n_
  integer :: it_
  logical :: converged_ = .false.
  real (real_kind), dimension(:,:), allocatable :: Q_
  real (real_kind), dimension(:), allocatable :: a_
  real (real_kind), dimension(:), allocatable :: b_
  real (real_kind), dimension(:), allocatable, target :: val_
  real (real_kind), dimension(:,:), allocatable, target :: vec_
  real (real_kind) :: groundE_ = 0
contains

  !########## PUBLIC ##########

!***********************************
  subroutine parameters_lanczos(max,acc,cut)  ! use the default value when zero
    integer :: max
    real (real_kind) :: acc,cut
    if (max.ne.0) max_iteration_ = max
    if (acc.ne.0) accuracy_ = acc
    if (cut.ne.0) cutoff_ = cut
  end subroutine parameters_lanczos


!***********************************
  integer function iteration_lanczos()
    iteration_lanczos = it_
  end function iteration_lanczos

!***********************************
  subroutine groundVec_lanczos(j,vec)
    real (real_kind), dimension(:) :: vec
    integer :: j
    vec = matmul(Q_(:,1:j),vec_(1:j,1))
    vec = vec / sqrt(dot_product(vec,vec))
  end subroutine groundVec_lanczos

!***********************************
  subroutine eigenstates_lanczos(val,vec)
    real (real_kind), dimension(:), pointer :: val
    real (real_kind), dimension(:,:), pointer :: vec
    val => val_
    vec => vec_
  end subroutine eigenstates_lanczos

!***********************************
  subroutine restarted_lanczos(n,r,groundenergy_lanczos)
    implicit none
    integer, intent (in) :: n
    real (real_kind), dimension(n), intent (out) :: r        !ground state vector
    real (real_kind), intent (out) :: groundenergy_lanczos
    integer :: j,k
    real (real_kind) :: b

    if (n.eq.1) return
    call new_lanczos(n)
    allocate(val_(iteration_),vec_(iteration_,iteration_))
    allocate(a_(iteration_),b_(iteration_))
    j = 1
    converged: do while (.not.converged_)
       !print *,it_,j,groundE_
       call single_recursion_step(j,r,Q_,a_,b_)
       if (j.eq.iteration_) then
          ! Restart the lanczos algorithm
          call diagonalize_restarted(j)
          r = matmul(Q_,vec_(:,1))
          b = sqrt(dot_product(r,r))
          Q_(:,1) = r / b
          j = 1
       else
          Q_(:,j+1) = r / b_(j)
          j = j + 1       
       end if
       if (it_ > max_iteration_) then
          write(*,*) '# (restarted) The maximum iteration numbers have been reached.'
          !deallocate(a_,b_)
          !deallocate(val_,vec_)
          !return
          converged_=.true.
       end if
    end do converged
    call groundVec_lanczos(j,r) ! store the ground state vector in r
    groundenergy_lanczos = groundE_
    deallocate(a_,b_)    
    deallocate(val_,vec_)
  end subroutine restarted_lanczos

!***********************************
  subroutine delete_lanczos()
    if (allocated(Q_)) deallocate(Q_)
    if (allocated(val_)) deallocate(val_)
    if (allocated(vec_)) deallocate(vec_)
  end subroutine delete_lanczos

  !########## PRIVATE ##########

!***********************************
  subroutine new_lanczos(n)
    integer :: n
    integer :: i
    real  (real_kind) :: b

    groundE_ = 1000d0
    it_ = 0
    n_ = n
    converged_ = .false.
    if (n.eq.1) return

    allocate(Q_(n_,iteration_))
    !call random_number(Q_(:,1))
    !Q_(:,1) = Q_(:,1)/sqrt(sum(Q_(:,1)*Q_(:,1)))
    Q_(:,1) = 0d0
    Q_(1,1) = 1d0

  end subroutine new_lanczos

!***********************************
  subroutine diagonalize_restarted(step)
    real (real_kind) :: WORK(2*iteration_),tmp(iteration_)
    integer :: info
    integer :: step

    val_ = a_
    tmp = b_
    call dstev('V',step,val_,tmp,vec_,iteration_,WORK,INFO)
    if (INFO/=0) then
       print '("# (dstev) INFO=",i6)', INFO
       stop
    end if

    if (abs(groundE_-a_(1)) < accuracy_) then
       converged_ = .true.
    end if
    groundE_ = a_(1)

  end subroutine diagonalize_restarted


!***********************************
  subroutine single_recursion_step(j,r,Q,a,b)
    integer :: j
    real (real_kind), dimension(:) :: r,a,b
    real (real_kind), dimension(:,:) :: Q

    call multiply_by_hamiltonian(Q(:,j),r)
    a(j) = dot_product(Q(:,j),r)
    if (j.eq.1) then
       r = r - a(1) * Q(:,1)
    else
       r = r - Q(:,j)*a(j) - Q(:,j-1)*b(j-1) 
    end if
    b(j) = sqrt(dot_product(r,r))
    converged_ = (b(j).lt.cutoff_)
    it_ = it_ + 1

  end subroutine single_recursion_step


end module simple_lanczos
