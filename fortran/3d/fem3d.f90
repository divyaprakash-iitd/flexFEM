module fem3d
    use iso_fortran_env, only: int32, real64
    implicit none

    ! integer(int32), parameter :: ncor = 3
    ! integer(int32), parameter :: nnodes = 4

    type :: festruct
        integer(int32), allocatable :: M(:,:) ! Indices providing connectivity information
        real(real64), allocatable :: XE(:,:) ! Coordinates of points
        real(real64), allocatable :: fden(:,:)
        real(real64), allocatable :: U(:,:)
        real(real64), allocatable :: bc(:,:,:) ! Shape function coefficients
        real(real64), allocatable :: velem(:)
        real(real64) :: co
        real(real64) :: kval
        real(real64) :: dl

        contains
            procedure :: calculate_forces
            procedure :: update_position
            procedure :: set_force
            procedure :: set_velocity
            procedure :: shapecoefficients

    end type festruct

    interface festruct
        module procedure :: festruct_constructor
    end interface festruct

    contains

    ! type(festruct) function festruct_constructor(M, XE, bc, boundary, aelem, fden, U, co, kval, dl) result(self)
    type(festruct) function festruct_constructor(M, XE, fden, U, co, kval, dl) result(self)
        implicit none
        integer(int32), intent(in):: M(:,:) ! Indices providing connectivity information
        real(real64), intent(in) :: XE(:,:) ! Coordinates of points
        real(real64), intent(in) :: fden(:,:)
        real(real64), intent(in) :: U(:,:)
        real(real64), intent(in) :: co
        real(real64), intent(in) :: kval
        real(real64), intent(in) :: dl

        self%M = M
        self%XE = XE
        self%fden = fden
        self%U = U
        self%co = co
        self%kval = kval
        self%dl = dl

        call self%shapecoefficients()
        
    end function festruct_constructor
    
    subroutine shapecoefficients(self)
        implicit none
        class(festruct), intent(inout) :: self

        integer(int32) :: ielem, nelems, npoints, inode
        real(real64), allocatable :: bcoff(:,:,:), acoff(:,:), velem(:)
        real(real64) :: A(4,4), C(4,4) ! Co-factor matrix
        
        npoints = size(self%XE,1)
        nelems = size(self%M,1)

        allocate(bcoff(nelems,4,3), acoff(nelems,4), velem(nelems))
        bcoff   = 0.0d0
        acoff   = 0.0d0
        A       = 0.0d0
        C       = 0.0d0
        velem   = 0.0d0

        do ielem = 1,nelems
            A(:,1) = 1.0d0
            do inode = 1,4
                A(inode,2) = self%XE(self%M(ielem,inode),1)
                A(inode,3) = self%XE(self%M(ielem,inode),2)
                A(inode,4) = self%XE(self%M(ielem,inode),3)
            end do 

        
            ! Find the volume of the tetrahedrals
            velem(ielem) = (1.0d0/6.0d0)*abs(determinant(A))

            ! Go over each node of the element and
            ! calculate the coefficients for each node (a and b)
            C = cofactor(A)
            bcoff(ielem,:,:) = (1.0d0/6.0d0)*C(:,2:4)
            acoff(ielem,:) = (1.0d0/6.0d0)*C(:,1)

            ! Divide by the volume of the element
            bcoff(ielem,:,:) = bcoff(ielem,:,:)/velem(ielem)
            acoff(ielem,:) = acoff(ielem,:)/velem(ielem)

        end do
        self%bc = bcoff
        self%velem = velem

    end subroutine shapecoefficients

    elemental subroutine set_force(self,fden)
        implicit none
        class(festruct), intent(inout) :: self
        real(real64), intent(in) :: fden
        
        integer(int32) :: ipoints
        ! ! Fix the bottom nodes (Make force zero)
        ! forall (ipoints = 1:size(self%XE,1), self%boundary(ipoints,2).eqv..True.) &
        !     self%fden(ipoints,:) = fden
    end subroutine set_force
    
    elemental subroutine set_velocity(self,UN)
        implicit none
        class(festruct), intent(inout) :: self
        real(real64), intent(in) :: UN
        self%U = UN 
    end subroutine set_velocity

    elemental impure subroutine calculate_forces(self)
        implicit none
        class(festruct), intent(inout) :: self
        real(real64) :: D(size(self%M,1),3,3)
        real(real64) :: didx(size(self%M,1),3,4)
        real(real64) :: djdx(size(self%M,1),3,4)
        real(real64) :: jac(size(self%M,1))
        real(real64) :: F(size(self%M,1),3,4)
        integer :: j, k, nelem, nidx, ipoints

        nelem = size(self%M,1)

        ! Calculate deformation gradient
        D = deformationgradient(self%M,self%XE,self%bc)
        ! Calculate the jacobian
        jac = jacobian(D)
        ! Calculate the invariant derivative term
        didx = invariantderivative(D,self%bc)
        ! Calculate the jacobian derivative term
        djdx = jacobianderivative(D,self%bc)

        ! Calculate forces at each node per element 
        ! (4 nodes, 3 components of force per element)
        F = 0.0d0
        do k = 1, nelem
            F(k,:,:) = -0.5d0 * self%co * self%velem(k) * &
                        ( didx(k,:,:) - 2.0d0/jac(k) * djdx(k,:,:) * &
                        (1.0d0 - self%kval/self%co * log(jac(k)))) 
        end do

        ! Adding forces for all nodes
        self%fden = 0.0d0
        do k = 1, nelem
            do j = 1, 4
               nidx = self%M(k,j)
               self%fden(nidx,:) = self%fden(nidx,:) + F(k,:,j)            
            end do
        end do 
        
        ! ! ! Fix the bottom nodes (Make force zero)
        ! ! forall (ipoints = 1:size(self%XE,1), self%boundary(ipoints,2).eqv..True.) &
        ! !     self%fden(ipoints,:) = 0.0d0

    end subroutine calculate_forces
    
    elemental subroutine update_position(self, dt)
        implicit none
        class(festruct), intent(inout) :: self
        real(real64), intent(in) :: dt

        integer(int32) :: ipoint, npoints

        npoints = size(self%XE,1)

        do ipoint = 1,npoints
            self%XE(ipoint,:) = self%XE(ipoint,:) + dt * self%U(ipoint,:)
        end do
            

    end subroutine update_position

    subroutine read_data(file_path,M)
          implicit none
          character(len=*), intent(in) :: file_path
          real*8, intent(inout), allocatable :: M(:,:)
          integer :: i, j, dims(2)

          ! Open the binary file for reading
          open(unit=10, file=file_path, access='stream', form='unformatted', status='old')

          ! Read the size of the matrix
          read(10) dims

          allocate(M(dims(1),dims(2)))
          ! Read the data
            do j = 1, dims(2)
                do i = 1, dims(1)
                    read(10) M(i, j)
                end do
            end do

          ! Close the file
          close(10)

          ! Now you have the data in the 'data' array
          ! Do whatever processing you need here
          print*, dims(1), dims(2)
    end subroutine read_data

    ! Helper functions
    ! Jacobian derivative term
    function jacobianderivative(D,b) result(djdx)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: djdx(size(b,1),3,4)
        real(dp) :: CM(3,3) ! Co-factor matrix

        integer :: k, nelem

        nelem = size(b,1)

        do concurrent (k = 1:nelem)
            ! Calculate the co-factor matrix of the deformation tensor of the current element
            CM = cofactor(D(k,:,:))
            
            djdx(k,:,:) = matmul(CM, transpose(b(k,:,:)))
        end do

    end function jacobianderivative 

    ! Invariant derivative term
    function invariantderivative(D,b) result(didx)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: didx(size(b,1),3,4)

        integer :: k, nelem

        nelem = size(b,1)
        
        do concurrent (k = 1:nelem)
            didx(k,:,:) = 2 * matmul(D(k,:,:), transpose(b(k,:,:))) 
        end do

    end function invariantderivative

    function deformationgradient(M,P,b) result(D)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        integer, intent(in) :: M(:,:)  ! Triangulation
        real(dp), intent(in) :: P(:,:) ! Points
        real(dp), intent(in) :: b(:,:,:) ! Shape function coefficient
        real(dp) :: D(size(M,1),3,3) ! Deformation gradient
        real(dp), allocatable :: X(:,:) ! Coordinates matrix of nodes of a single tetrahedron

        integer :: i, j, nelem, idxelem, nvertex, ncor

        ncor = 3 ! In the case of 3D
        nelem = size(M,1) ! No of elements 
        nvertex = 4 ! No of nodes in an element\triangle 
        allocate(X(ncor,nvertex))

        do i = 1, nelem
            do concurrent (j = 1:nvertex)
                idxelem = M(i,j)
                X(:,j) = P(idxelem,:)
            end do
            D(i,:,:) = matmul(X, b(i,:,:))
        end do
    end function deformationgradient

    function jacobian(F) result(J)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: F(:,:,:)
        real(dp)             :: J(size(F,1))
        ! real(dp), intent(out)             :: J(size(F,1))

        integer             :: k, nelem

        nelem = size(F,1)
        do k = 1, nelem
            J(k) = determinant(F(k,:,:))
        end do

    end function jacobian
    ! end subroutine jacobian

    pure function determinant(M) result(D)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: M(:,:)
        real(dp) :: D

        integer :: n
        n = size(M,1)

        if (n.eq.2) then
            D = M(1,1) * M(2,2) - M(1,2) * M(2,1)
        else if (n.eq.3) then
            D = M(1,1) * (M(2,2) * M(3,3) - M(3,2) * M(2,3)) - &
                M(1,2) * (M(2,1) * M(3,3) - M(3,1) * M(2,3)) + &
                M(1,3) * (M(2,1) * M(3,2) - M(2,2) * M(3,1))
        else if (n.eq.4) then
            D = M(1,1) * (M(2,2) * (M(3,3) * M(4,4) - M(4,3) * M(3,4)) - &
                        M(2,3) * (M(3,2) * M(4,4) - M(4,2) * M(3,4)) + &
                        M(2,4) * (M(3,2) * M(4,3) - M(3,3) * M(4,2))) - &
                M(1,2) * (M(2,1) * (M(3,3) * M(4,4) - M(4,3) * M(3,4)) - &
                        M(2,3) * (M(3,1) * M(4,4) - M(4,1) * M(3,4)) + &
                        M(2,4) * (M(3,1) * M(4,3) - M(3,3) * M(4,1))) + &
                M(1,3) * (M(2,1) * (M(3,2) * M(4,4) - M(4,2) * M(3,4)) - &
                        M(2,2) * (M(3,1) * M(4,4) - M(4,1) * M(3,4)) + &
                        M(2,4) * (M(3,1) * M(4,2) - M(3,2) * M(4,1))) - &
                M(1,4) * (M(2,1) * (M(3,2) * M(4,3) - M(4,2) * M(3,3)) - &
                        M(2,2) * (M(3,1) * M(4,3) - M(4,1) * M(3,3)) + &
                        M(2,3) * (M(3,1) * M(4,2) - M(3,2) * M(4,1)))
        end if


    end function determinant

    pure function cofactor(A) result(CM)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: A(:,:)
        real(dp), allocatable :: CM(:,:) ! result
        real(dp), allocatable :: C(:,:) ! result

        integer :: n,i,j,l,m,ielem
        real(dp), allocatable :: Minor(:)

        n = size(A,1)

        ! Allocate and initialize
        allocate(Minor((n-1)*(n-1)), CM(n,n), C((n-1),(n-1)))
        Minor   = 0.0d0
        CM  = 0.0d0

        do i = 1,n
            do j = 1,n
                ielem = 1
                do l = 1,n
                    do m = 1,n
                        if ((l.NE.i).AND.(m.NE.j)) then
                            Minor(ielem) = A(l,m)
                            ielem = ielem + 1
                        end if
                    end do
                end do
                C = transpose(reshape(Minor,[n-1,n-1]))
                ! print *, C
                CM(i,j) = (-1)**(i+j)*determinant(C)
            end do
        end do

    end function cofactor

    subroutine print_matrix(matrix)
        implicit none
        real(8), dimension(:,:) :: matrix
        integer :: nrows, ncols
        integer :: i, j
    
        ! Determine the number of rows and columns
        nrows = size(matrix, 1)
        ncols = size(matrix, 2)
    
        ! Print the matrix
        do i = 1, nrows
            print *, (matrix(i, j), j = 1, ncols)
        end do
    
    end subroutine print_matrix


    ! function determinant(A) result(D)
    !     implicit none
    !     real(real64), intent(in) :: A(:,:)
    !     real(real64) :: D

    !     integer(int32) :: n   !Size of the matrix
    !     integer, allocatable :: ipiv(:)
    !     integer::  info
    !     real(real64), allocatable :: work(:)
    !     integer :: i, j

    !     ! External LAPACK routine for LU factorization
    !     external dgetrf

    !     ! Allocate
    !     n = size(A,1)
    !     allocate(ipiv(n),work(n))

    !     ! Compute the LU factorization of the matrix A
    !     call dgetrf(n, n, A, n, ipiv, info)

    !     ! Check for successful factorization
    !     if (info /= 0) then
    !         print *, "Error: Matrix is singular or other error occurred"
    !         stop
    !     endif

    !     ! Compute the determinant using the LU factors
    !     D = 1.0
    !     do i = 1, n
    !         D = D * A(i, i)
    !     end do

    !     ! Apply the permutation from the factorization
    !     do i = 1, n
    !         if (ipiv(i) /= i) then
    !             D = -D
    !         endif
    !     end do

    !     ! Output the determinant
    !     ! print *, "Determinant:", determinant
    ! end function determinant

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
 
end module fem3d
