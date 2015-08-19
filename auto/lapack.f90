module lapack

! This module contains subroutines to simplify the calling of LAPACK drivers for solving the eigenvectors and eigenvalues of real
! symmetric matrices. To find out more about each individual subroutines read the comments at the top of the subroutines.

contains

!************************************************************************************

subroutine ssyev_lapack(dim,matrix,eigvalues)

	!  %-----------------------------------------------------------------%
	!  |  This subroutine calls the driver SSYEV.                        |
	!  |                                                                 |
	!  |  Solves the eigenvalues and eigenvectors of the matrix.         |
	!  |      - Eigenvalues returned in eigvalues                        |
	!  |      - Eigenvectors returned in matrix (the input matrix)       |
	!  %-----------------------------------------------------------------%

	implicit none

	integer, intent(in) :: dim               ! dimension of the matrix
	real, intent(inout) :: matrix(dim,dim)   ! the matrix to be solved
	real, intent(out) :: eigvalues(dim)      ! the outputed eigenvalues

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

!************************************************************************************

subroutine ssyev_lapack1(dim,matrix,eigvalues)

	!  %-----------------------------------------------------------------%
	!  |  This subroutine calls the driver SSYEV.                        |
	!  |                                                                 |
	!  |  Solves the only the eigenvalues of the matrix.                 |
	!  |      - Eigenvalues returned in eigvalues                        |
	!  %-----------------------------------------------------------------%

	implicit none

	integer, intent(in) :: dim                 ! dimension of the matrix
	real, intent(inout) :: matrix(dim,dim)     ! the matrix to be solved
	real, intent(out) :: eigvalues(dim)        ! the outputed eigenvalues

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

!************************************************************************************

subroutine ssyevr_lapack1(dim,matrix,eigvalues,eigvectors)

	!  %-----------------------------------------------------------------%
	!  |  This subroutine calls the driver SSYEVR.                       |
	!  |                                                                 |
	!  |  Solves the only the lowest eigenvalues of the matrix.          |
	!  |      - Eigenvalues returned in eigvalues                        |
	!  %-----------------------------------------------------------------%


	implicit none

	integer, intent(in) :: dim                  ! dimension of the matrix
	real, intent(in) :: matrix(dim,dim)         ! the matrix to be solved
	real, intent(out) :: eigvalues(dim)         ! the outputed eigenvalues
	real, intent(out) :: eigvectors(dim,dim)    ! array contains no useful information this case

	integer, parameter :: NSELECT=1

	integer :: LDA, LDZ
	real :: VL,VU
	real :: ABSTOL = -1
	integer :: IL=1,IU=NSELECT

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


	LDA = dim; LDZ = dim
	LWORK = -1
	LIWORK = -1
	
	allocate(ISUPPZ(2*dim))
	allocate(WORK(1),IWORK(1))

	call ssyevr('N','I','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)

	LWORK= int(WORK(1))
	LIWORK = IWORK(1)
	
	deallocate(WORK,IWORK)
	allocate(WORK(LWORK),IWORK(LIWORK))

	call ssyevr('N','I','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)

end subroutine ssyevr_lapack1

!************************************************************************************

subroutine ssyevr_lapack(dim,matrix,eigvalues,eigvectors)

	!  %-----------------------------------------------------------------%
	!  |  This subroutine calls the driver SSYEVR.                       |
	!  |                                                                 |
	!  |  Solves the eigenvalues and eigenvectors of the matrix.         |
	!  |      - Eigenvalues returned in eigvalues                        |
	!  |      - Eigenvectors returned in eigvectors                      |
	!  %-----------------------------------------------------------------%

	implicit none

	integer, intent(in) :: dim
	real, intent(in) :: matrix(dim,dim)
	real, intent(out) :: eigvalues(dim)
	real, intent(out) :: eigvectors(dim,dim)

	integer, parameter :: NSELECT=1

	integer :: LDA, LDZ
	real :: VL,VU
	real :: ABSTOL = -1
	integer :: IL=1,IU=1

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


	LDA = dim; LDZ = dim
	LWORK = -1
	LIWORK = -1
	
	allocate(ISUPPZ(2*dim))
	allocate(WORK(1),IWORK(1))

	call ssyevr('V','A','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)

	LWORK= int(WORK(1)) * 2
	LIWORK = IWORK(1) * 2
	
	deallocate(WORK,IWORK)
	allocate(WORK(LWORK),IWORK(LIWORK))

	call ssyevr('V','A','U',dim,matrix,LDA,VL,VU,IL,IU,ABSTOL,M,eigvalues,eigvectors,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)

end subroutine ssyevr_lapack

!************************************************************************************

subroutine ssyevd_lapack(dim,matrix,eigvalues)

	!  %-----------------------------------------------------------------%
	!  |  This subroutine calls the driver SSYEVD.                       |
	!  |                                                                 |
	!  |  Solves the eigenvalues and eigenvectors of the matrix.         |
	!  |      - Eigenvalues returned in eigvalues                        |
	!  |      - Eigenvectors returned in matrix (the input matrix)       |
	!  |                                                                 |
	!  |  * not extensively tested but most likely works                 |
	!  %-----------------------------------------------------------------%

	implicit none

	integer, intent(in) :: dim
	real, intent(inout) :: matrix(dim,dim)
	real, intent(out) :: eigvalues(dim)

	integer :: INFO
	integer :: LWORK, LIWORK
	real, allocatable,dimension(:) :: WORK
	integer, allocatable,dimension(:) :: IWORK

	LWORK = -1
	LIWORK = -1

	allocate(WORK(1))
	allocate(IWORK(1))

	call ssyevd('v','u',dim,matrix,dim,eigvalues,WORK,LWORK,IWORK,LIWORK,INFO)
	if (INFO /= 0) then
		write(*,*) "Error with LAPACK:", INFO
	end if

	LWORK = WORK(1)
	LIWORK = IWORK(1)

	deallocate(WORK,IWORK)
	allocate(WORK(LWORK))
	allocate(IWORK(LIWORK))

	call ssyevd('v','u',dim,matrix,dim,eigvalues,WORK,LWORK,IWORK,LIWORK,INFO)
	if (INFO /= 0) then
		write(*,*) "Error with LAPACK:", INFO
	end if

end subroutine ssyevd_lapack

!************************************************************************************

end module lapack