module lapack

contains

subroutine ssyev_lapack(dim,matrix,eigvalues)

	implicit none

	integer, intent(in) :: dim
	real, intent(inout) :: matrix(dim,dim)
	real, intent(out) :: eigvalues(dim)

	integer :: INFO
	integer :: LWORK
	real, allocatable,dimension(:) :: WORK

	LWORK = 3*dim

	allocate(WORK(LWORK))

	call ssyev('v','u',dim,matrix,dim,eigvalues,WORK,LWORK,INFO)
	if (INFO /= 0) then
		write(*,*) "Error with LAPACK:", INFO
	end if

end subroutine ssyev_lapack

subroutine ssyev_lapack1(dim,matrix,eigvalues)

	implicit none

	integer, intent(in) :: dim
	real, intent(inout) :: matrix(dim,dim)
	real, intent(out) :: eigvalues(dim)

	integer :: INFO
	integer :: LWORK
	real, allocatable,dimension(:) :: WORK

	LWORK = 3*dim

	allocate(WORK(LWORK))

	call ssyev('n','u',dim,matrix,dim,eigvalues,WORK,LWORK,INFO)
	if (INFO /= 0) then
		write(*,*) "Error with LAPACK:", INFO
	end if

end subroutine ssyev_lapack1

subroutine ssyevr_lapack(dim,matrix,eigvalues,eigvectors)

	implicit none

	integer, intent(in) :: dim
	real, intent(in) :: matrix(dim,dim)
	real, intent(out) :: eigvalues(dim)
	real, intent(out) :: eigvectors(dim,dim)

	real :: VL=0,VU=0
	real :: ABSTOL = -1
	integer :: IL=0,IU=0

	integer :: M
	integer :: INFO
	integer :: LWORK, LIWORK
	integer, allocatable, dimension(:) :: ISUPPZ, IWORK
	real, allocatable,dimension(:) :: WORK

	if (dim == 1) then
		eigvectors = 1
		eigvalues = matrix(1,1)
		return
	end if

	LWORK = -1
	LIWORK = -1
	
	allocate(ISUPPZ(2*dim))
	allocate(WORK(1),IWORK(1))

	call ssyevr('V','A','U',dim,matrix,dim,VL,VU,IL,IU,ABSTOL,dim,eigvalues,eigvectors,dim,WORK,LWORK,IWORK,LIWORK,INFO)

	LWORK= int(WORK(1))
	LIWORK = IWORK(1)

	deallocate(WORK,IWORK)
	allocate(WORK(LWORK),IWORK(LIWORK))

	call ssyevr('V','A','U',dim,matrix,dim,VL,VU,IL,IU,ABSTOL,dim,eigvalues,eigvectors,dim,WORK,LWORK,IWORK,LIWORK,INFO)

end subroutine ssyevr_lapack

end module lapack