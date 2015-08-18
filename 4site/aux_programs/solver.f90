program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/07/06            created code 
! P. Daley              15/07/09            added comments to code
!
! This program will solve for the eigenvalues and eigenvectors of any matrix. The program will ask for the input file name and the size of the matrix. Use the LAPACK library to do this. It is useful for verifying if code is functioning properly by making the program output it's matrices to a file then solving them with this program and then comparing.

  implicit none
  
  double precision, allocatable,  dimension(:,:) :: x
  integer :: i,j, error
  double precision, allocatable, dimension(:) :: y, work
  integer :: lwork, matrix_size
  integer :: info
  character(len=20) :: FM1, FM2, FM3, file_name
  character(len=1) :: check1, check2, check3


  write(*,*) "What is the size of the matrix?"
  read(*,*) matrix_size
  write(*,*) "What is the file name of the matrix?"
  read(*,*) file_name
  write(*,*) "Would you like to print the matrix? (y/n)"
  read(*,*) check3
  write(*,*) "Do you want the eigenvalues outputed? (y/n)"
  read(*,*) check1
  write(*,*) "Do you want the eigenvectors outputed? (y/n)"
  read(*,*) check2

  allocate(x(matrix_size,matrix_size))
  allocate(y(matrix_size))
  lwork = matrix_size*5
  allocate(work(lwork))
  
  open(unit=99,file=file_name,status='old', action='read', IOSTAT  = error)

  read(99,*) x

  close(99)
  
  write(FM1, '("(I2,2X", I0,"F7.4)")') matrix_size
  write(FM2, '("(A12,2X", I0,"F10.4)")') matrix_size
  write(FM3,'("(A14,", I0,"F8.4)")') matrix_size

  if (check3 == "y") then
     do i=1,matrix_size
        write(*,FM1) i, x(i,(/(j,j=1,matrix_size)/))
     end do
  end if

  call dsyev('v','u',matrix_size,x,matrix_size,y,work,lwork,info)
  
  if (check1 == "y") then
     write(*,FM2) "eigenvalues:", y((/(i,i=1,matrix_size)/))
  end if

  if (check2 == "y") then 
     do i=1,matrix_size
        write(*,FM3) "eigenvectors:", x((/(j,j=1,matrix_size)/),i)
     end do
  end if

end program main
