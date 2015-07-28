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

end module lapack