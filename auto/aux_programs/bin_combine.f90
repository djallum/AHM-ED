program main

	! This program takes a dos+ipr input file and will make the energy binning width larger by combining 2 bins
	! into one larger bin. Will out put it into new file will same name but with a 'b' appended to it. Can apply
	! it multiple times to get 4 bins into one or even 8.

	implicit none

	integer :: i                                      ! counters for loops
	integer :: nbins                                  ! number of bins in the input file 
	integer :: nbinsc                                 ! number of bins in the output file (nbins/2)
	character(len=10) :: str_nbins                    ! string that contains the number of bins (ex. "240")
	character (len=100) :: preamble(9)                ! array to contain the preamble of the input data file
	character (len=100) :: infilename, outfilename    ! input and output data files
	real, allocatable, dimension(:,:) :: old_data     ! array to hold the frequency of bins, DOS and GIPR of input data

	!---------------------Input from user----------------------------------
	write(*,*) "What is the data file name?"
	read(*,*) infilename

	!------------------Make the new filename-------------------------------
	outfilename = infilename(1:LEN_TRIM(infilename)-4)
	write(outfilename,'(A)') trim(adjustl(outfilename)) // 'b.dat'

	open(25,file=infilename,status='old',action='read')

	open(35,file=outfilename,status='replace',action='write')

	!------------Read and Modify preamble of data file--------------------
	do i=1,9
		read(25,'(A)') preamble(i)
		if (i == 4) then
			read(preamble(i)(10:13),'(I4)') nbins
			write(str_nbins,'(I4)') nbins/2
			write(35,*) trim(adjustl(preamble(i)(1:9))) // ' ' // trim(adjustl(str_nbins)) // ' ' // trim(adjustl(preamble(i)(14:)))
		else
			write(35,*) trim(adjustl(preamble(i)))
		end if
	end do

	allocate(old_data(nbins,3))

	!---------------------Read Old DOS and GIPR----------------------------
	do i=1,nbins
		read(25,*) old_data(i,1), old_data(i,2), old_data(i,3)
	end do

	!---------------------Write new DOS and GIPR---------------------------
	do i=1,nbins,2
		write(35,*) old_data(i,1), (old_data(i,2) + old_data(i+1,2))/2, (old_data(i,3) + old_data(i+1,3))/2
	end do

	close(25)
	close(35)

end program main