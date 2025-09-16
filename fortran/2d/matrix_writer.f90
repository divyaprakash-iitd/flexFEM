module matrix_writer
  implicit none
  private
  public :: write_to_file

  interface write_to_file
    module procedure write_real_matrix
    module procedure write_int_matrix
    module procedure write_real_vector
    module procedure write_int_vector
  end interface write_to_file

contains

  subroutine write_real_matrix(filename, A)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: A(:,:)
    integer :: i, j, nrows, ncols
    open(unit=10, file=filename, status='replace')
    nrows = size(A, 1)
    ncols = size(A, 2)
    do i = 1, nrows
      do j = 1, ncols
        write(10, '(F20.12)', advance='no') A(i,j)
        if (j < ncols) write(10, '(A)', advance='no') ' '
      end do
      write(10, *)
    end do
    close(10)
  end subroutine write_real_matrix

  subroutine write_int_matrix(filename, A)
    character(len=*), intent(in) :: filename
    integer(kind=8), intent(in) :: A(:,:)
    integer :: i, j, nrows, ncols
    open(unit=10, file=filename, status='replace')
    nrows = size(A, 1)
    ncols = size(A, 2)
    do i = 1, nrows
      do j = 1, ncols
        write(10, '(I12)', advance='no') A(i,j)
        if (j < ncols) write(10, '(A)', advance='no') ' '
      end do
      write(10, *)
    end do
    close(10)
  end subroutine write_int_matrix

  subroutine write_real_vector(filename, A)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: A(:)
    integer :: i, n
    open(unit=10, file=filename, status='replace')
    n = size(A)
    do i = 1, n
      write(10, '(F20.12)') A(i)
    end do
    close(10)
  end subroutine write_real_vector

  subroutine write_int_vector(filename, A)
    character(len=*), intent(in) :: filename
    integer(kind=8), intent(in) :: A(:)
    integer :: i, n
    open(unit=10, file=filename, status='replace')
    n = size(A)
    do i = 1, n
      write(10, '(I12)') A(i)
    end do
    close(10)
  end subroutine write_int_vector

end module matrix_writer
