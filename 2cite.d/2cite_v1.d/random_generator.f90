module random_gen

contains

subroutine gen_seed()

integer :: seed_size, clock, i
integer, allocatable, dimension(:) :: seed


call random_seed(size=seed_size)
allocate(seed(seed_size))
call system_clock(count = clock)
seed=clock + 37*(/(i-1,i=1,seed_size)/)

call random_seed(put=seed)

deallocate(seed)

end subroutine gen_seed

end module random_gen
