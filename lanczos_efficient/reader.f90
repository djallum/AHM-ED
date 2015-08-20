program main

	implicit none

	integer :: i,j
	integer, parameter :: sizes = 3
	character(len=4) :: ichar,jchar
	character(len=10) :: junk1, junk2, junk3, junk4
	character(len=20) :: junk5, junk6, junk7, junk8, junk9, junk10
	character(len=100) :: filename
	integer, parameter :: real_kind = kind(1.d0)
	real (real_kind), allocatable, dimension(:,:) :: GroundStateEnergy

	allocate(GroundStateEnergy(sizes,sizes))
	GroundStateEnergy = 0

	do i=1,sizes
		do j=1,sizes
			write(ichar,'(I4)') i
			write(jchar,'(I4)') j
			filename = "GS4u0v0d20m0n_up"//trim(adjustl(ichar))//"n_dn"//trim(adjustl(jchar))//".dat"
			!write(*,*), filename
			open(unit=1,file=filename,status="old")
			read(1,*) junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8, junk9, junk10, GroundStateEnergy(i,j)
		end do
	end do
	write(*,*) minloc(GroundStateEnergy)

end program main