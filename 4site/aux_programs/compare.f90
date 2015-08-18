program main
  
  ! Programmer         Date(yy/mm/dd)         Modification
  ! -----------        --------------         ------------
  ! P. Daley              15/07/08            created code 
  ! P. Daley              15/08/17            added comments to code
  !
  ! This program takes too files that contain a matrix and will compare the two. It will print out the indices for which the two matrices are not equal
  ! and optionally the matrices themselves.

  implicit none
  
  real, allocatable, dimension(:,:) :: matrix1, matrix2   ! arrays that store the matrices
  character(len=20) :: file_name1, file_name2             ! names of the files the matrices are found in
  character(len=20) :: FMstring                           ! string that contains the information for the formatted output of the matrices
  character(len=1)  :: mprint              ! contains either 'y' or 'n' to specify if matrices should be printed to terminal
  integer :: matrix_size                   ! size of the matrices   
  integer :: error                         ! variable to store error messages when opening files
  integer :: i, j                          ! counters for loops

  !------------Get info from user---------------------
  write(*,*) "What is the size of the matrices?"
  read(*,*) matrix_size
  write(*,*) "Enter the filename of the matrix"
  read(*,*) file_name1
  write(*,*) "Enter the filename of the matrix"
  read(*,*) file_name2
  write(*,*) "Print matrices (y or n)"
  read(*,*) mprint

  !------------Read in the matrices-------------------
  allocate(matrix1(matrix_size,matrix_size))
  allocate(matrix2(matrix_size,matrix_size))

  open(unit=111,file=file_name1,status='old',action='read',IOSTAT=error)
  if (error/=0) then
     write(*,*) 'error opening matrix 1 file. Error number:', error
  end if
  open(unit=222,file=file_name2,status='old',action='read',IOSTAT=error)
  if (error/=0) then
     write(*,*) 'error opening matrix 2 file. Error number:', error
  end if

  read(111,*) matrix1
  read(222,*) matrix2

  !------------Compare matrices---------------------
  do i=1,matrix_size
     do j=1,matrix_size
        if(matrix1(i,j) /= matrix2(i,j) .and. i /= j) then
           write(*,*) "i:", i, "j:", j
        end if
     end do
  end do

  !------------Print matrices---------------------

  if (mprint == 'n') go to 15        ! if matrices aren't to be printed skip to the end

  write(FMstring,'("(I2,2X,", I0,"F8.4)")') matrix_size

  write(*,*) "------------------------------------------------------------"
  do i=1,matrix_size
     write(*,FMstring) i, matrix1(i,(/(j,j=1,matrix_size)/))
  end do
  write(*,*) "-----------------------------------------------------------"
   do i=1,matrix_size
     write(*,FMstring) i, matrix2(i,(/(j,j=1,matrix_size)/))
  end do
  
  15 continue
  
  close(111)
  close(222)

end program main
