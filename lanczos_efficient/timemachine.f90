module timemachine
  implicit none
  !---------- Run-Time variables ----------
  integer, private :: values(8)
  integer, private :: days=0,delta_days=0,total_days=0
  integer, private :: hours=0,delta_hours=0,total_hours=0
  integer, private :: minutes=0,delta_mins=0,total_mins=0
  integer, private :: seconds=0,delta_secs=0,total_secs=0
  integer, private :: milliseconds=0,delta_msecs=0,total_msecs=0

contains

  subroutine time_starter()
    call date_and_time(values=values)
    hours = values(5)
    minutes = values(6)
    seconds = values(7)
    milliseconds = values(8)
  end subroutine time_starter

  subroutine time_elapsed(dh,dm,ds,dms)
    integer :: dh,dm,ds,dms
    call date_and_time(values=values)
    delta_hours = values(5) - hours
    delta_mins = values(6) - minutes
    delta_secs = values(7) - seconds
    delta_msecs = values(8) - milliseconds
    if (delta_msecs.lt.0) then          ! everything from here down is just tidying up the time outputs
       delta_msecs = delta_msecs + 1000
       delta_secs = delta_secs - 1
    end if
    if (delta_secs.lt.0) then          ! same as "if (delta_secs < 0) then"
       delta_secs = delta_secs + 60
       delta_mins = delta_mins - 1
    end if
    if (delta_mins.lt.0) then
       delta_mins = delta_mins + 60
       delta_hours = delta_hours - 1
    end if
    if (delta_hours.lt.0) then
       delta_hours = delta_hours + 24
    end if
    hours = values(5)
    minutes = values(6)
    seconds = values(7)
    milliseconds = values(8)

    dh = delta_hours
    dm = delta_mins
    ds = delta_secs
    dms = delta_msecs

    total_msecs = total_msecs + delta_msecs
    total_secs = total_secs + delta_secs
    total_mins = total_mins + delta_mins
    total_hours = total_hours + delta_hours
    if (total_msecs.ge.1000) then          ! same as "if (delta_secs > 0) then"
       total_msecs = total_msecs - 1000
       total_secs = total_secs + 1
    end if
    if (total_secs.ge.60) then
       total_secs = total_secs - 60
       total_mins = total_mins + 1
    end if
    if (total_mins.ge.60) then
       total_mins = total_mins - 60
       total_hours = total_hours + 1
    end if
    if (total_hours.ge.24) then
       total_hours = total_hours - 24
       total_days = total_days + 1
    end if
  end subroutine time_elapsed

  subroutine time_total(td,th,tm,ts,tms)
    integer :: td,th,tm,ts,tms
    td = total_days
    th = total_hours
    tm = total_mins
    ts = total_secs
    tms = total_msecs
  end subroutine time_total

end module timemachine
