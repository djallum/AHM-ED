module routines

implicit none

contains

subroutine site_potentials(delta,E)

real, intent(in) :: delta
real, intent(out) :: E(2)
real :: random1, random2

call random_number(random1)
call random_number(random2)
random1 = random1 - 0.5          !centering the random numbers about 0
random2 = random2 - 0.5
E(1) = delta*max(random1,random2) !multiplying the random numbers by the width of disorder
E(2) = delta*min(random1,random2)

end subroutine site_potentials

subroutine transformations(PES_up,PES_down,IPES_up,IPES_down,phase_PES_up,phase_PES_down,phase_IPES_up,phase_IPES_down)

integer, dimension(2,16), intent(out) :: PES_up,PES_down,IPES_up,IPES_down
integer, dimension(2,16), intent(out) :: phase_PES_up,phase_PES_down,phase_IPES_up,phase_IPES_down

PES_up = 0; PES_down = 0
IPES_up = 0; IPES_down = 0
phase_PES_up = 0; phase_PES_down = 0
phase_IPES_up = 0; phase_IPES_down = 0

PES_up(1,1) = 2; phase_PES_up(1,1) = 1
PES_up(1,3) = 6; phase_PES_up(1,3) = -1 
PES_up(1,4) = 10; phase_PES_up(1,4) = -1
PES_up(1,5) = 8;  phase_PES_up(1,5) = -1 
PES_up(1,7) = 14;  phase_PES_up(1,7) = 1
PES_up(1,9) = 12; phase_PES_up(1,9) = 1
PES_up(1,11) = 13; phase_PES_up(1,11) = 1
PES_up(1,15) = 16; phase_PES_up(1,15) = -1

PES_down(1,1) = 4; phase_PES_down(1,1) = 1
PES_down(1,2) = 10; phase_PES_down(1,2) = 1
PES_down(1,3) = 9; phase_PES_down(1,3) = -1 
PES_down(1,5) = 7; phase_PES_down(1,4) = -1
PES_down(1,6) = 12;  phase_PES_down(1,5) = -1
PES_down(1,8) = 14; phase_PES_down(1,8) = -1
PES_down(1,11) = 15; phase_PES_down(1,11) = 1
PES_down(1,13) = 16; phase_PES_down(1,13) = 1

IPES_up(1,2) = 1; phase_IPES_up(1,2) = 1
IPES_up(1,6) = 3; phase_IPES_up(1,6) = -1 
IPES_up(1,8) = 5; phase_IPES_up(1,8) = -1
IPES_up(1,10) = 4; phase_IPES_up(1,10) = -1 
IPES_up(1,12) = 9;  phase_IPES_up(1,12) = 1
IPES_up(1,13) = 11; phase_IPES_up(1,13) = 1
IPES_up(1,14) = 7; phase_IPES_up(1,14) = 1
IPES_up(1,16) = 15; phase_IPES_up(1,16) = -1           !found by hand

IPES_down(1,4) = 1; phase_IPES_down(1,4) = 1
IPES_down(1,7) = 5; phase_IPES_down(1,7) = -1
IPES_down(1,9) = 3; phase_IPES_down(1,9) = -1 
IPES_down(1,10) = 2; phase_IPES_down(1,10) = 1
IPES_down(1,12) = 6; phase_IPES_down(1,12) = -1
IPES_down(1,14) = 8; phase_IPES_down(1,14) = -1
IPES_down(1,15) = 11; phase_IPES_down(1,15) = 1
IPES_down(1,16) = 13; phase_IPES_down(1,16) = 1

PES_up(2,1) = 3; phase_PES_up(2,1) = 1
PES_up(2,2) = 6; phase_PES_up(2,2) = 1 
PES_up(2,4) = 9; phase_PES_up(2,4) = 1
PES_up(2,5) = 11; phase_PES_up(2,5) = -1 
PES_up(2,7) = 15; phase_PES_up(2,7) = -1
PES_up(2,8) = 13; phase_PES_up(2,8) = -1
PES_up(2,10) = 12; phase_PES_up(2,10) = 1
PES_up(2,14) = 16; phase_PES_up(2,14) = -1

PES_down(2,1) = 5; phase_PES_down(2,1) = 1
PES_down(2,2) = 8; phase_PES_down(2,2) = 1
PES_down(2,3) = 11; phase_PES_down(2,3) = 1 
PES_down(2,4) = 7; phase_PES_down(2,4) = 1
PES_down(2,6) = 13; phase_PES_down(2,6) = 1
PES_down(2,9) = 15; phase_PES_down(2,9) = 1
PES_down(2,10) = 14; phase_PES_down(2,10) = 1
PES_down(2,12) = 16; phase_PES_down(2,12) = 1

IPES_up(2,3) = 1; phase_IPES_up(2,3) = 1
IPES_up(2,6) = 2; phase_IPES_up(2,6) = 1 
IPES_up(2,9) = 4; phase_IPES_up(2,9) = 1
IPES_up(2,11) = 5; phase_IPES_up(2,11) = -1 
IPES_up(2,12) = 10; phase_IPES_up(2,12) = 1
IPES_up(2,13) = 8; phase_IPES_up(2,13) = -1
IPES_up(2,15) = 7; phase_IPES_up(2,15) = -1
IPES_up(2,16) = 14; phase_IPES_up(2,16) = -1           !found by hand

IPES_down(2,5) = 1; phase_IPES_down(2,5) = 1
IPES_down(2,7) = 4; phase_IPES_down(2,7) = 1
IPES_down(2,8) = 2; phase_IPES_down(2,8) = 1 
IPES_down(2,11) = 3; phase_IPES_down(2,11) = 1
IPES_down(2,13) = 6; phase_IPES_down(2,13) = 1
IPES_down(2,14) = 10; phase_IPES_down(2,14) = 1
IPES_down(2,15) = 9; phase_IPES_down(2,15) = 1
IPES_down(2,16) = 12; phase_IPES_down(2,16) = 1
end subroutine transformations

end module routines
