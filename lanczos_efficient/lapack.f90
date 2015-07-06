module lapack
  use math_setup

contains

  !########## REAL SINGLE PRECISION ##########

  !---------- Symmetric ----------

  subroutine ssyevx_lapack(n,H,val,vec)
    integer :: n
    real :: H(n,n),val,vec(n)
    ! for LAPACK
    real :: ABSTOL,W(n),Z(n,n),WORK(8*n)
    integer :: M,LWORK,IWORK(5*n),IFAIL(n),INFO
    LWORK = 8*n
    call ssyevx("V","I","L",n,H,n,0,0,1,1,ABSTOL,M,W,Z,n,WORK,LWORK,IWORK,IFAIL,INFO)
    if (INFO/=0) then
       print '("# (ssyevx) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine ssyevx_lapack

  !########## REAL DOUBLE PRECISION ##########

  !---------- Symmetric Tridiagonal ----------

  subroutine dstev_lapack(n,a,b,val,vec)
    integer :: n
    real (real_kind) :: a(n),b(n-1),val(n),vec(n,n)
    ! for LAPACK
    real (real_kind) :: E(n-1),WORK(2*n-2)
    integer :: INFO
    val = a
    E = b
    call dstev('V',n,val,E,vec,n,WORK,INFO)
    if (INFO/=0) then
       print '("# (dstev) INFO=",i6)', INFO
       stop
    end if
  end subroutine dstev_lapack

  subroutine dstevd_lapack(n,a,b,val,vec)
    integer :: n
    real (real_kind) :: a(n),b(n-1),val(n),vec(n,n)
    ! for LAPACK
    real (real_kind) :: E(n-1)
    real (real_kind) :: WORK(n**2+4*n+1)
    integer :: LWORK,IWORK(5*n+3),LIWORK,INFO
    val = a
    E = b
    LWORK = n**2+4*n+1
    LIWORK = 5*n+3
    call dstevd('V',n,val,E,vec,n,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("(dstevd) INFO=",i6,"use DSTEV instead")', INFO
       call dstev_lapack(n,a,b,val,vec)
       stop
    end if
  end subroutine dstevd_lapack

  subroutine dstevr_ground_lapack(n,a,b,val,vec)
    integer :: n
    real (real_kind) :: a(n),b(n-1),val,vec(n)
    ! for LAPACK
    real (real_kind) :: D(n),E(n-1),ABSTOL,W(n),Z(n,n)
    real (real_kind) :: WORK(20*n)
    integer :: M,ISUPPZ(n),LWORK,IWORK(10*n),LIWORK,INFO
    D = a
    E = b
    LWORK = 20*n
    LIWORK = 10*n
    call dstevr('V','I',n,D,E,0,0,1,1,ABSTOL,M,W,Z,n,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("# (dstevr_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dstevr_ground_lapack

  subroutine dstevx_ground_lapack(n,a,b,val,vec)
    integer :: n
    real (real_kind) :: a(n),b(n-1),val,vec(n)
    ! for LAPACK
    real (real_kind) :: ABSTOL,W(n),Z(n,n),WORK(5*n)
    integer :: M,IWORK(5*n),IFAIL(n),INFO
    call dstevx('V','I',n,A,B,0,0,1,1,ABSTOL,M,W,Z,n,WORK,IWORK,IFAIL,INFO)
    if (INFO/=0) then
       print '("# (dstevx_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dstevx_ground_lapack

  !---------- Symmetric Band ----------

  subroutine dsbev_lapack(n,kd,H,val,vec)
    integer :: n,kd
    real (real_kind) :: H(kd+1,n),val(n),vec(n,n)
    ! for LAPACK
    real (real_kind) :: AB(kd+1,n),WORK(3*n-2)
    integer :: INFO
    AB = H
    call dsbev('V','L',n,kd,AB,kd+1,val,vec,n,WORK,INFO)
    if (INFO/=0) then
       print '("# (dsbev) INFO=",i6)', INFO
       stop
    end if
  end subroutine dsbev_lapack

  subroutine dsbevd_lapack(n,kd,H,val,vec)
    integer :: n,kd
    real (real_kind) :: H(kd+1,n),val(n),vec(n,n)
    ! for LAPACK
    real (real_kind) :: AB(kd+1,n)
    real (real_kind) :: WORK(2*n**2+5*n+1)
    integer :: LWORK,IWORK(5*n+3),LIWORK,INFO
    AB = H
    LWORK = 2*n**2+5*n+1
    LIWORK = 5*n+3
    call dsbevd('V','L',n,kd,AB,kd+1,val,vec,n,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("(dsbevd) INFO=",i6,"use DSBEV instead")', INFO
       call dsbev_lapack(n,kd,H,val,vec)
       stop
    end if
  end subroutine dsbevd_lapack

  subroutine dsbev_ground_lapack(n,kd,H,val,vec)
    integer :: n,kd
    real (real_kind) :: H(kd+1,n),val,vec(n)
    ! for LAPACK
    real (real_kind) :: AB(kd+1,n),W(n),Z(n,n),WORK(3*n-2)
    integer :: INFO
    AB = H
    call dsbev('V','L',n,kd,AB,kd+1,W,Z,n,WORK,INFO)
    if (INFO/=0) then
       print '("# (dsbev_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dsbev_ground_lapack

  subroutine dsbevd_ground_lapack(n,kd,AB,val,vec)
    integer :: n,kd
    real (real_kind) :: AB(kd+1,n),val,vec(n)
    ! for LAPACK
    real (real_kind) :: W(n),Z(n,n),WORK(2*n**2+5*n+1)
    integer :: LWORK,IWORK(5*n+3),LIWORK,INFO
    LWORK = 2*n**2+5*n+1
    LIWORK = 5*n+3
    call dsbevd('V','L',n,kd,ab,kd+1,W,Z,n,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("# (dsbevd_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dsbevd_ground_lapack

  subroutine dsbevx_ground_lapack(n,kd,H,val,vec)
    integer :: n,kd
    real (real_kind) :: H(kd+1,n),val,vec(n)
    ! for LAPACK
    real (real_kind) :: AB(kd+1,n)
    real (real_kind) :: Q(n,n),ABSTOL,W(n),Z(n,n),WORK(7*n)
    integer :: M,IWORK(5*n),IFAIL(n),INFO
    AB = H
    call dsbevx('V','I','L',n,kd,AB,kd+1,Q,n,0,0,1,1,ABSTOL,M,W,Z,n,WORK,IWORK,IFAIL,INFO)
    if (INFO/=0) then
       print '("# (dsbevx_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dsbevx_ground_lapack

  !---------- Symmetric ----------

  subroutine dsyevd_lapack(n,h,val,vec)
    integer :: n
    real (real_kind) :: h(n,n),val(n),vec(n,n)
    ! for LAPACK
    real (real_kind) :: WORK(2*n**2+6*n+1)
    integer :: LWORK,IWORK(5*n+3),LIWORK,INFO
    vec = h
    LWORK = 2*n**2+6*n+1
    LIWORK = 5*n+3
    call dsyevd("V","L",n,vec,n,val,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("# (dsyevd) INFO=",i6)', INFO
       stop
    end if
  end subroutine dsyevd_lapack

  subroutine dsyev_ground_lapack(n,H,val,vec)
    integer :: n
    real (real_kind) :: H(n,n),val,vec(n)
    ! for LAPACK
    real (real_kind) :: T(n,n),W(n),WORK(3*n)
    integer :: LWORK,INFO
    T = H
    LWORK = 3*n-1
    call dsyev('V','U',n,T,n,W,WORK,LWORK,INFO)
    if (INFO/=0) then
       print '("# (dsyev_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = T(:,1)
  end subroutine dsyev_ground_lapack

  subroutine dsyevx_ground_lapack(n,H,val,vec)
    integer :: n
    real (real_kind) :: H(n,n),val,vec(n)
    ! for LAPACK
    real (real_kind) :: T(n,n),ABSTOL,W(n),Z(n,n),WORK(8*n)
    integer :: M,LWORK,IWORK(5*n),IFAIL(n),INFO
    T = H
    LWORK = 8*n
    call dsyevx("V","I","L",n,T,n,0,0,1,1,ABSTOL,M,W,Z,n,WORK,LWORK,IWORK,IFAIL,INFO)
    if (INFO/=0) then
       print '("# (dsyevx_ground) INFO=",i6)', INFO
       stop
    end if
    val = W(1)
    vec = Z(:,1)
  end subroutine dsyevx_ground_lapack

  !########## COMPLEX DOUBLE PRECISION ##########

  !---------- Symmetric Band ----------

  subroutine zhbevd_lapack(n,kd,ab,val,vec)
    integer :: n,kd
    complex (real_kind) :: ab(kd+1,n),vec(n,n)
    real (real_kind) :: val(n)
    ! for LAPACK
    complex (real_kind) :: WORK(2*n**2)
    real (real_kind) :: RWORK(2*n**2+5*n+1)
    integer :: LWORK,LRWORK,IWORK(5*n+3),LIWORK,INFO
    LWORK = 2*n**2
    LRWORK = 2*n**2+5*n+1
    LIWORK = 5*n+3
    call zhbevd('V','L',n,kd,ab,kd+1,val,vec,n,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("# (zhbevd) INFO=",i6)', INFO
       stop
    end if
  end subroutine zhbevd_lapack

  !---------- Symmetric ----------

  subroutine zheevd_lapack(n,A,val)
    integer :: n
    complex (real_kind) :: A(n,n)
    real (real_kind) :: val(n)
    ! for LAPACK
    complex (real_kind) :: WORK(n**2+2*n)
    real (real_kind) :: RWORK(2*n**2+5*n+1)
    integer :: LWORK,LRWORK,IWORK(5*n+3),LIWORK,INFO
    LWORK = n**2+2*n
    LRWORK = 2*n**2+5*n+1
    LIWORK = 5*n+3
    call zheevd("V","L",n,A,n,val,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
    if (INFO/=0) then
       print '("# (zheevd) INFO=",i6)', INFO
       stop
    end if
  end subroutine zheevd_lapack

  subroutine zhesv_lapack(n,A,B)
    integer :: n
    complex (real_kind) :: A(n,n),B(n,n)
    ! for LAPACK
    integer :: IPIV(n),LWORK,INFO
    complex (real_kind) :: WORK(16*n)
    B = (0.0,0.0)
    do INFO = 1, n
       B(INFO,INFO) = 1.d0
    end do
    LWORK = 16*n
    call zhesv("L",n,n,A,n,IPIV,B,n,WORK,LWORK,INFO)
    if (INFO/=0) then
       print '("# (zhesv) INFO=",i6)', INFO
       stop
    end if
  end subroutine zhesv_lapack

end module lapack
