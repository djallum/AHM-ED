module random_generator
contains

  subroutine random_gen()
    integer :: nsd,clock
    integer, allocatable :: sd(:)
    integer :: i

    call random_seed(size  = nsd)
    allocate(sd(nsd))
    call system_clock(count=clock)
    sd = clock * 37 * (/ (i, i = 1, nsd) /)
!    sd = (/ 0,-1535653932,1223659432,-311994500,-1847648432,911664932,-623989000,2135324364 /)

    call random_seed(put=sd)
    !print 100,sd
    deallocate(sd)
100 format("# random seed ",40(i12,1x))
  end subroutine random_gen

  subroutine random_gen2(my_id)
    integer :: my_id
    integer :: nsd
    integer :: values(8)
    integer :: clock
    integer, allocatable :: sd(:)
    integer :: i

    call random_seed(size  = nsd)
    allocate(sd(nsd))
    call date_and_time(values=values)
!    call system_clock(count=clock)
    clock = values(7)*1000 + values(8)
    sd = (my_id+1) * clock * 37 * (/ (i - 1, i = 1, nsd) /)
    call random_seed(put=sd)
    deallocate(sd)
  end subroutine random_gen2

end module random_generator
