program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/07/07            created code 
! P. Daley              15/07/09            added comments to code
!
! This program tests if the constructed photo emmision spectrum (PES) and inverse photo emmision spectrum (IPES) lookup tables for the 4 site code are working properly. You can assign the ground vector to be at any state any perform PES and IPES transitions from every site and see what states it is transitioned to. Every states should transition to 8 other states. To check all the states to see if there is 8 transtions uncomment the do loop called states.Uses subroutines in routines.f90.

  use routines

  implicit none
  
  integer :: state             ! the state that the that the ground vector started in 
  integer :: i,j               ! counters
  integer :: count=0           ! sum of how many transtions have occured (should be 8)
  real, dimension(4,256) :: PES_down_ground=0, PES_up_ground=0, IPES_down_ground=0, IPES_up_ground=0 ! resulting vectors from each transition

  75 continue
  write(*,*) "what state is the ground vector?"
  read(*,'(I4)') state                                     ! take the state number as imput
  if (state < 0 .or. state > 256) then                     ! check if it's a valid state number
     write(*,*) "Invalid state number please try again"    ! print error message if not valid state
     write(*,*) ""                                         ! add blank line
     goto 75                                               ! make them try again
  end if
  

  call transformations() ! subroutines in routines.f90 that makes the IPES and PES lookup tables

  !states: do state=1,256

     write(*,*) "--------------------------------------------"

     v_ground = 0
     v_ground(state) = 1 

     ! do each of the transformations on the ground vector
     do j=1,4
        do i=1,256
           if (PES_up(j,i)==0) then
              PES_up_ground(j,i) = 0   ! zero if there is no transition
           else 
              PES_up_ground(j,i) = v_ground(PES_up(j,i))*phase_PES_up(j,i)
           end if
           if (PES_down(j,i)==0) then
              PES_down_ground(j,i) = 0
           else 
              PES_down_ground(j,i) = v_ground(PES_down(j,i))*phase_PES_down(j,i)
           end if
           if (IPES_up(j,i)==0) then
              IPES_up_ground(j,i) = 0
           else 
              IPES_up_ground(j,i) = v_ground(IPES_up(j,i))*phase_IPES_up(j,i)
           end if
           if (IPES_down(j,i)==0) then
              IPES_down_ground(j,i) = 0
           else 
              IPES_down_ground(j,i) = v_ground(IPES_down(j,i))*phase_IPES_down(j,i)
           end if
        end do
     end do

     do j=1,4
        do i=1,256
           if (PES_up_ground(j,i) /= 0) then
              write(*,100)"PES_up site",j, "location", i ! prints the site of the transformation and the new state
              count = count + 1                        ! add one to the counter for how many transformation there are
           end if
           if (PES_down_ground(j,i) /= 0) then
              write(*,100)"PES_down site",j, "location", i
              count = count + 1
           end if
           if (IPES_up_ground(j,i) /= 0) then
              write(*,100)"IPES_up site",j, "location", i
              count = count + 1
           end if
           if (IPES_down_ground(j,i) /= 0) then
              write(*,100)"IPES_down site",j, "location", i
              count = count + 1
           end if
        end do
     end do

     100 format (A15,I4,A10,I4)
     write(*,*) "--------------------------------------------"

     if (count /= 8) then
        write(*,*) "error there was ot 8 transitions for state", state  ! check if right amount of transformations occured
     end if

     write(*,*) "count:", count  ! prints how many transformations there was

     !count=0           !reset counter for the next loop
  !end do states

end program main
