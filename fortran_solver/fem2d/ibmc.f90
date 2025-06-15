program ibmc
    use soft_particles
    use iso_c_binding, only: c_int, c_double, c_loc
    implicit none

    integer(C_INT) :: n, niter, iter
    logical, allocatable :: isOnPerimeter(:), isOnLeftSector(:), isOnRightSector(:), & 
                                isOnLeft(:), isOnRight(:), isLeftPatch(:), isRightPatch(:)
    real(c_double) :: a, b, eps, dt, radius
    real(c_double), allocatable :: XP(:), YP(:), X1(:), Y1(:)
    real(c_double), allocatable :: FXC(:),FYC(:), Fxleft(:), Fxright(:), fxboundary(:), fyboundary(:)
    real(c_double), allocatable :: F1XC(:),F1YC(:), F1ZC(:), U(:), V(:), W(:)
    real(c_double), allocatable :: X(:),Y(:),Z(:), nangle(:), xb(:), yb(:)
    ! real(c_double), parameter :: pi = acos(-1.0)
    real(c_double) :: angle45
    logical, allocatable :: rmask(:), lmask(:)
    
    call sayhello()
   
    ! Particle (Ellipse) parameters
    ! To-DO: Add functionality to read them using namelist
    a     = 2.5e-4
    b     = 1.25e-4
    eps    = 1.0e-4
    radius = 1.0d0


    ! Generate ellipse
    ! print *, "n = ", n
    call generateellipse(n)
    ! print *, "n = ", n
    allocate(FXC(n), FYC(n),X1(n),Y1(n))
    allocate(F1XC(n), F1YC(n), F1ZC(n))    
    allocate(X(n), Y(n), Z(n),nangle(n))

    ! Use get positions to get the cooridinates of the particles
    ! Then create masks based on a criteria that gets you the points 
    ! on a sector spanning -30 to 30 degrees and -150 to 150 degrees
    call getpositions(X,Y,Z,n)

    ! Allocate mask
    allocate(isOnPerimeter(n))
    
    ! Mask based on ellipse condition
    ! isOnPerimeter = abs((X/a)**2 + (Y/b)**2 - 1.0d0) < eps
    isOnPerimeter = abs(X**2 + Y**2 - radius**2) < eps
    
    ! Mask based on x-location
    isOnLeft = X < 0.0d0
    isOnRight = X > 0.0d0

    ! Mask based on angle
    nangle = atan2(Y,X)
    angle45 = 45.0 * pi / 180.0
    isOnLeftSector = abs(nangle) > angle45
    isOnRightSector = abs(nangle) < angle45

    ! Final masks
    isLeftPatch = isOnPerimeter .and. isOnLeft .and. isOnLeftSector
    isRightPatch = isOnPerimeter .and. isOnRight .and. isOnRightSector

    ! XP = merge(X, 0.0d0, isOnPerimeter)
    ! YP = merge(Y, 0.0d0, isOnPerimeter)

    XP = pack(X, isOnPerimeter)
    YP = pack(Y, isOnPerimeter)
    xb = pack(X, isRightPatch)
    yb = pack(Y, isRightPatch)

    ! print *, xp
    ! print *, yp
    ! print *, xb
    ! print *, yb

    !Apply boundary forces
    FXC = 500.0d0 !0.250d0
    FYC = 0.0d0
    fxleft = merge(-1.0D0*fxc,0.0d0,isLeftPatch)
    fxright = merge(1.0D0*fxc,0.0d0,isRightPatch)
    fxboundary = fxleft + fxright
    fyboundary = 0.0d0*fxboundary

    print *, fxboundary
    niter = 10000
    dt = 0.001
    do iter = 1, niter
        ! print *, fxboundary
        ! Call applyboundaryforces inside of calculateforces and make the fx/fyboundary optional
        ! Make sure that fden is initialised as zero at the start of every iteration
        call applyboundaryforces(fxboundary,fyboundary,n)
        ! call applyboundaryforces(Fright,FYC,n)
    
        call calculateforces()
        call updatepositions(dt)
    end do

end program ibmc
