program main

  ! Programmer         Date(yy/mm/dd)         Modification
  ! -----------        --------------         ------------
  ! P. Daley              15/08/19            created code
  ! P. Daley              15/10/02            added filling recalculation 
  !
  ! %------------------------------------------------------------%
  ! |  This code takes two data files with DOS and GIPR with     |
  ! |  the same paramter and number of bins and superimposes     |
  ! |  them to create a larger ensemble. The DOS and GIPR        |
  ! |  are weighted according the ensemble size of the data      |
  ! |  file they came from.                                      |
  ! |                                                            |
  ! |  Output file name is the first input filename with a 'c'   |
  ! |  appended to the end.                                      |
  ! |                                                            |
  ! |  To combine more then two simply run the program with      |
  ! |  the already joined data as one input and third data file  | 
  ! |  to be joined as the other input. Repeat as many times as  |
  ! |  wanted.                                                   |
  ! %------------------------------------------------------------%

  implicit none

  character(len=50) :: filename1, filename2      ! the two input file names
  character(len=50) :: outfilename               ! the outputed data file name 
  character (len=100) :: preamble1(9)            ! array to contain the preamble of the input data file 1
  character (len=100) :: preamble2(9)            ! array to contain the preamble of the input data file 2
  real :: sum, half_sum                          ! for calculating density of states filling
  integer :: npairs1, npairs2                    ! ensemble sizes of the input files
  integer :: npairs_out                          ! ensemble size of the output file
  character(len=15) :: str_npairs_out            ! ensemble size of the output file as string
  character(len=15) :: str_nbins                 ! number of bins as a string
  integer :: nbins1,nbins2,nbins                 ! number of energy bins for each data file and the output (must be the same)
  integer :: i                                   ! counter for loops
  integer :: istart 
  real, allocatable, dimension(:,:) :: indata1, indata2     ! array to hold the frequency of bins, DOS and GIPR of input data
  real, allocatable, dimension(:,:) :: outdata              ! array to hold the frequency of bins, DOS and GIPR of output data


  !---------------------Input from user----------------------------------
  write(*,*) "What is the first data file name?"
  read(*,*) filename1

  write(*,*) "What is the second data file name?"
  read(*,*) filename2

  !------------------Make the new filename-------------------------------
  outfilename = filename1(1:LEN_TRIM(filename1)-4)                  ! remove the ".dat"
  write(outfilename,'(A)') trim(adjustl(outfilename)) // 'c.dat'    ! add the "c.dat"

  open(25,file=filename1,status='old',action='read')
  open(26,file=filename2,status='old',action='read')

  open(35,file=outfilename,status='replace',action='write')

  !-----------------Read preamble of data files--------------------------
  do i=1,9
     read(25,'(A)') preamble1(i)
     read(26,'(A)') preamble2(i)
  end do

  read(preamble1(4)(10:13),'(I4)') nbins1
  read(preamble2(4)(10:13),'(I4)') nbins2
  if(nbins1 == nbins2) then
     nbins=nbins1                 ! data files can be combined
  else
     write(*,*) "input files has different amount of bins" ! data files cannot be combined
     stop
  end if
 

  !--------------------Find ensemble sizes-------------------------------
  do i=16,50
     if (preamble1(4)(i:i) == '=')  then              ! find the column that the ensemble size integer starts at
        read(preamble1(4)(i+1:),'(I15)') npairs1     ! print the ensemble size to npairs1
        istart = i                                   ! remember column it start
     end if
     if (preamble2(4)(i:i) == '=')  then              ! find the column that the ensemble size integer starts at
        read(preamble2(4)(i+1:),'(I15)') npairs2     ! print the ensemble size to npairs2
     end if
  end do

  npairs_out = npairs1 + npairs2                       ! calculate new ensemble size
  write(str_npairs_out,'(I15)') npairs_out             ! convert it to a string
  
  !---------------------Read Old DOS and GIPR----------------------------
  allocate(indata1(nbins,3))
  allocate(indata2(nbins,3))
  allocate(outdata(nbins,3))
  
  do i=1,nbins
     read(25,*) indata1(i,1), indata1(i,2), indata1(i,3)
     read(26,*) indata2(i,1), indata2(i,2), indata2(i,3)
  end do
  
  !--------------------Calculate output and filling------------------------
  do i=1,nbins 
     outdata(i,1) = (indata1(i,1)*npairs1 + indata2(i,1)*npairs2)/npairs_out    ! weighting each data file
     outdata(i,2) = (indata1(i,2)*npairs1 + indata2(i,2)*npairs2)/npairs_out    ! weighting each data file
     outdata(i,3) = (indata1(i,3)*npairs1 + indata2(i,3)*npairs2)/npairs_out    ! weighting each data file
  end do
  
  do i=1,nbins
     sum = sum + outdata(i,2)
  end do

  do i=1,nbins/2
     half_sum = half_sum + outdata(i,2)
  end do
  
  !-----------------Write preamble to output file------------------------
  write(str_nbins,'(I15)') nbins
  
  write(35,*) "# created by main.f90 with subroutines in routines.f90"
  write(35,*) "#"
  write(35,*) trim(adjustl(preamble1(3)))
  write(35,'(A)',advance='no') " # nbins= " // trim(adjustl(str_nbins)) 
  write(35,*) " ensemble size=" // adjustl(trim(str_npairs_out))  
  write(35,*) "#"
  write(35,*) "# Filling:", half_sum/sum
  write(35,*) "# Filling Error:", outdata(nbins/2+1,2)/sum
  write(35,*) "#"
  write(35,*) "#  Frequency           DOS             GIPR"
  
  !---------------------Write new DOS and GIPR----------------------------
  
  do i=1,nbins 
     write(35,*) outdata(i,1), outdata(i,2), outdata(i,3)                       ! printing the data to output file
  end do

  close(25)
  close(26)
  close(35)

end program main
