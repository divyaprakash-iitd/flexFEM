module fem2d
    use iso_fortran_env, only: int32, real64
    implicit none

    type :: festruct
        integer(int32), allocatable :: M(:,:) ! Indices providing connectivity information
        real(real64), allocatable :: XE(:,:) ! Coordinates of points
        real(real64), allocatable :: bc(:,:,:) ! Shape function coefficients
        logical, allocatable :: boundary(:,:) ! Top, bottom, left and right boundaries
        real(real64), allocatable :: aelem(:)
        real(real64), allocatable :: fden(:,:)
        real(real64), allocatable :: U(:,:)
        real(real64) :: co
        real(real64) :: kval = 10.0d0
        real(real64) :: dl

        contains
            procedure :: calculate_forces
            procedure :: update_position
            procedure :: shapecoefficients

    end type festruct

    interface festruct
        module procedure :: festruct_constructor
    end interface festruct

    contains

    type(festruct) function festruct_constructor(M, XE, fden, U, co, kval, dl) result(self)
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

        call shapecoefficients(self)
    end function festruct_constructor
    
    elemental impure subroutine calculate_forces(self)
        class(festruct), intent(inout) :: self
        real(real64) :: D(size(self%M,1),2,2)
        real(real64) :: didx(size(self%M,1),2,3)
        real(real64) :: djdx(size(self%M,1),2,3)
        real(real64) :: jac(size(self%M,1))
        real(real64) :: F(size(self%M,1),2,3)
        integer :: j, k, nelem, nidx, ipoints

        nelem = size(self%M,1)

        ! Calculate deformation gradient
        D = deformationgradient(self%M,self%XE,self%bc)

        ! Calculate the invariant derivative term
        didx = invariantderivative(D,self%bc)

        ! Calculate the jacobian derivative term
        djdx = jacobianderivative(D,self%bc)

        ! Calculate the jacobian
        jac = jacobian(D)

        ! Calculate forces at each node per element 
        ! (3 nodes, 2 components of force per element)
        F = 0.0d0
        do k = 1, nelem
            F(k,:,:) = -0.5d0 * self%co * self%aelem(k) * &
                        ( didx(k,:,:) - 2.0d0/jac(k) * djdx(k,:,:) * &
                        (1.0d0 - self%kval/self%co * log(jac(k)))) 
        end do

        ! Adding forces for all nodes
        ! self%fden = 0.0d0
        do k = 1, nelem
            do j = 1, 3
               nidx = self%M(k,j)
               self%fden(nidx,1) = self%fden(nidx,1) + F(k,1,j)            
               self%fden(nidx,2) = self%fden(nidx,2) + F(k,2,j)            
            end do
        end do 
        
    end subroutine calculate_forces
    
    elemental subroutine update_position(self, dt)
        class(festruct), intent(inout) :: self
        real(real64), intent(in) :: dt

        integer(int32) :: ipoint, npoints

        npoints = size(self%XE,1)

        do ipoint = 1,npoints
            self%XE(ipoint,1) = self%XE(ipoint,1) + dt * self%U(ipoint,1)
            self%XE(ipoint,2) = self%XE(ipoint,2) + dt * self%U(ipoint,2)
        end do
            
    end subroutine update_position

    ! Helper functions
    ! Jacobian derivative term
    pure function jacobianderivative(D,b) result(djdx)
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: djdx(size(b,1),2,3)

        integer :: k, nelem
        real(dp) :: MD(2,2)

        nelem = size(b,1)

        do concurrent (k = 1:nelem)
            ! Rearrange the elements of the deformation matrix
            MD(1,1) = D(k,2,2)
            MD(1,2) = -D(k,2,1)
            MD(2,1) = -D(k,1,2)
            MD(2,2) = D(k,1,1)

            djdx(k,:,:) = matmul(MD, b(k,:,:))
        end do

    end function jacobianderivative 

    ! Invariant derivative term
    pure function invariantderivative(D,b) result(didx)
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: didx(size(b,1),2,3)

        integer :: k, nelem

        nelem = size(b,1)
        
        do concurrent (k = 1:nelem)
            didx(k,:,:) = 2 * matmul(D(k,:,:), b(k,:,:)) 
        end do

    end function invariantderivative

    pure function deformationgradient(M,P,b) result(D)
        integer, parameter :: dp = kind(0.d0)
        integer, intent(in) :: M(:,:)  ! Triangulation
        real(dp), intent(in) :: P(:,:) ! Points
        real(dp), intent(in) :: b(:,:,:) ! Shape function coefficient
        real(dp) :: D(size(M,1),2,2) ! Deformation gradient

        integer :: i, j, nelem, idxelem, nvertex

        real(dp) :: X(2,3) ! Coordinates matrix of nodes of a triangle

        nelem = size(M,1) ! No of elements 
        nvertex = 3 ! No of nodes in an element\triangle 

        do i = 1, nelem
            do concurrent (j = 1:nvertex)
                idxelem = M(i,j)
                X(:,j) = P(idxelem,:)
            end do
            D(i,:,:) = matmul(X, transpose(b(i,:,:)))
        end do
    end function deformationgradient

    pure function jacobian(F) result(J)
    ! subroutine jacobian(F,J)
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

    subroutine shapecoefficients(self)
        implicit none
        class(festruct), intent(inout) :: self

        integer(int32) :: ielem, nelems, npoints, inode
        real(real64), allocatable :: bcoff(:,:,:), acoff(:,:), velem(:)
        real(real64) :: A(3,3), C(3,3) ! Co-factor matrix
        
        npoints = size(self%XE,1)
        nelems = size(self%M,1)

        allocate(bcoff(nelems,2,3), acoff(nelems,3), velem(nelems))
        bcoff   = 0.0d0
        acoff   = 0.0d0
        A       = 0.0d0
        C       = 0.0d0
        velem   = 0.0d0

        do ielem = 1,nelems
            A(3,:) = 1.0d0
            do inode = 1,3
                A(1,inode) = self%XE(self%M(ielem,inode),1)
                A(2,inode) = self%XE(self%M(ielem,inode),2)
            end do 
        
            ! Find the volume of the tetrahedrals
            velem(ielem) = abs(determinant(A))

            ! Go over each node of the element and
            ! calculate the coefficients for each node (a and b)
            C = cofactor(A)
            bcoff(ielem,:,:) = C(1:2,:)
            acoff(ielem,:) = C(3,:)

            ! Divide by the volume of the element
            bcoff(ielem,:,:) = bcoff(ielem,:,:)/velem(ielem)
            acoff(ielem,:) = acoff(ielem,:)/velem(ielem)

        end do
        self%bc = bcoff
        self%aelem = velem

    end subroutine shapecoefficients

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

end module fem2d
