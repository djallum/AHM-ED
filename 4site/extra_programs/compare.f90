program main

  implicit none
  
  real, allocatable, dimension(:,:) :: matrix1, matrix2
  character(len=20) :: file_name1, file_name2, FMstring
  integer :: matrix_size, error, i, j

  write(*,*) "What is the size of the matrices?"
  read(*,*) matrix_size
  write(*,*) "Enter the filename of the matrix"
  read(*,*) file_name1
  write(*,*) "Enter the filename of the matrix"
  read(*,*) file_name2

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


  do i=1,matrix_size
     do j=1,matrix_size
        if(matrix1(i,j) /= matrix2(i,j) .and. i /= j) then
           write(*,*) "i:", i, "j:", j
        end if
     end do
  end do

  write(FMstring,'("(I2,2X,", I0,"F8.4)")') matrix_size

  write(*,*) "------------------------------------------------------------"
  do i=1,matrix_size
     write(*,FMstring) i, matrix1(i,(/(j,j=1,matrix_size)/))
  end do
  write(*,*) "-----------------------------------------------------------"
   do i=1,matrix_size
     write(*,FMstring) i, matrix2(i,(/(j,j=1,matrix_size)/))
  end do
  
  close(111)
  close(222)

end program main
