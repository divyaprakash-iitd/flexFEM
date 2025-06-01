module matrix_writer
  use, intrinsic :: iso_c_binding
  implicit none
contains

  subroutine write_vector(vec, prefix, iter)
    ! Writes a 1D vector to a file with optional iterator
    real(8), intent(in) :: vec(:)
    character(len=*), intent(in) :: prefix
    integer, intent(in), optional :: iter
    character(len=100) :: filename
    integer :: i, ios

    if (present(iter)) then
      write(filename, '(A, "_", I2.2, ".txt")') trim(prefix), iter
    else
      write(filename, '(A, ".txt")') trim(prefix)
    end if

    open(unit=10, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file: ", filename
      return
    end if

    do i = 1, size(vec)
      write(10, *) vec(i)
    end do

    close(10)
  end subroutine write_vector

  subroutine write_matrix(mat, prefix, iter)
    ! Writes a 2D matrix to a file with optional iterator
    real(8), intent(in) :: mat(:,:)
    character(len=*), intent(in) :: prefix
    integer, intent(in), optional :: iter
    character(len=100) :: filename
    integer :: i, j, ios
    integer :: rows, cols

    rows = size(mat,1)
    cols = size(mat,2)

    if (present(iter)) then
      write(filename, '(A, "_", I2.2, ".txt")') trim(prefix), iter
    else
      write(filename, '(A, ".txt")') trim(prefix)
    end if

    open(unit=11, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file: ", filename
      return
    end if

    do i = 1, rows
      do j = 1, cols
        write(11, "(F10.5)", advance='no') mat(i, j)
        if (j < cols) write(11, "(A)", advance='no') " "
      end do
      write(11, *)
    end do

    close(11)
  end subroutine write_matrix

end module matrix_writer

