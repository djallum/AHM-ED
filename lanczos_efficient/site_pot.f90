program main

	implicit none
	
	integer, parameter :: dp = kind(1.0d0)  !double precision kind
	integer :: seed_size, clock, i
	integer, allocatable, dimension(:) :: seed
	real :: delta
  	real :: E(8)
  	real :: random

  	E = 0
  	delta = 12

	call random_seed(size=seed_size)
	allocate(seed(seed_size))
	call system_clock(count = clock)
	seed=clock + 37*(/(i-1,i=1,seed_size)/)
	call random_seed(put=seed)
	deallocate(seed)

	do i=1,8
     	call random_number(random)              ! gives a random number between 0 and 1
     	E(i) = delta*(random - 0.5_dp)          ! centers the random numbers about 0
  	end do

  write(*,*) E

end program main