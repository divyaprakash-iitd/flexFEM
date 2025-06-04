module mod_io
    use iso_fortran_env, only: int32, real64
    implicit none
    
    private
    public :: write_field, write_field_lim, read_column, read_fem_data
contains

     subroutine write_field_lim(F,fieldname,timestep,xlb, xub, ylb, yub)
        real(real64), intent(in) :: F(xlb:xub,ylb:yub)
        character(len=1), intent(in) :: fieldname
        integer(int32), intent(in) :: timestep, xlb, xub, ylb, yub

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        character(len=8) :: itnumber
        integer(int32) :: i,j

        write(itnumber,"(I8.8)") timestep
        filename = fieldname // '_' // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", STATUS="unknown", position="append")
        do j = 1,size(F,2)
            write(fileunit, '(*(1p1e20.11))')( F(i,j) , i = 1,size(F,1))
        end do
        close(fileunit)
    end subroutine write_field_lim


    subroutine write_field(F,fieldname,timestep)
        real(real64), intent(in) :: F(:,:)
        character(len=1), intent(in) :: fieldname
        integer(int32), intent(in) :: timestep

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        character(len=8) :: itnumber
        integer(int32) :: i,j

        write(itnumber,"(I8.8)") timestep
        filename = fieldname // '_' // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", STATUS="unknown", position="append")
        do j = 1,size(F,2)
            write(fileunit, '(*(1p1e20.11))')( F(i,j) , i = 1,size(F,1))
        end do
        close(fileunit)
    end subroutine write_field

    subroutine read_column(file_name, column_array, num_rows)
        implicit none
        character(len=*), intent(in) :: file_name
        real(real64), dimension(:), allocatable, intent(out) :: column_array
        integer, intent(out) :: num_rows
        integer :: unit_number, i

        ! Open the file
        open(newunit=unit_number, file=file_name, status='old', action='read')

        ! Count the number of rows in the file
        num_rows = 0
        do while (.true.)
        read(unit_number, *, iostat=i)
        if (i /= 0) exit
        num_rows = num_rows + 1
        end do

        ! Allocate memory for the array
        allocate(column_array(num_rows))

        ! Rewind the file to the beginning
        rewind(unit_number)

        ! Read the column of numbers into the array
        do i = 1, num_rows
        read(unit_number, *) column_array(i)
        end do

        ! Close the file
        close(unit_number)
  end subroutine read_column

  subroutine read_fem_data(mp,paelem,pb,pp)
    implicit none

    real(real64), intent(inout) :: pb(:,:,:)
    integer, intent(inout) :: mp(:,:)
    real(real64), intent(inout) :: paelem(:)
    real(real64), intent(inout) :: pp(:,:)

    !------- pb ----------!
    OPEN(33, FILE="pb.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) pb
    close(33)
    !---------------------!
    !------- paelem----------!
    OPEN(33, FILE="paelem.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) paelem
    close(33)
    !---------------------!

    !------- mp ----------!
    OPEN(33, FILE="mp.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) mp
    close(33)
    !---------------------!

    !------- pp ----------!
    OPEN(33, FILE="pp.bin",&
     FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
    READ(33) pp
    close(33)
    !---------------------!
    
    ! print *, shape(a)
    ! print *, shape(b)
    ! print *, b

  end subroutine read_fem_data

end module mod_io
