program main
    use soft_particles, only: generatefestructures, getpositions, &
        applyboundaryforces, calculateforces, getforces, updatepositions 
    use iso_c_binding, only: c_int, c_double, c_loc
    use matrix_writer, only: write_to_file
    implicit none

    real(8), parameter :: PI = 3.141592653589793
    integer(C_INT) :: n, niter, iter
    logical, allocatable :: isOnPerimeter(:), isBetweenAngleRange(:), & 
                                isOnLeft(:), isOnRight(:), isLeftPatch(:), isRightPatch(:)
    real(c_double) :: eps, dt, radius
    real(c_double), allocatable :: XE(:,:), FE(:,:), X1(:), Y1(:), XP(:), YP(:)
    real(c_double), allocatable :: FXC(:),FYC(:), Fxleft(:), Fxright(:), fxboundary(:), fyboundary(:)
    real(c_double), allocatable :: F1XC(:),F1YC(:), F1ZC(:)
    real(c_double), allocatable :: X(:),Y(:),Z(:), nangle(:), xb(:), yb(:)
    ! real(c_double), parameter :: pi = acos(-1.0)
    real(c_double) :: angle45
 
    ! Measure time
    real :: t_start, t_end, t_elapsed
    character(len=20)           :: filename

    ! Particle geometrical parameters
    ! To-DO: Add functionality to read them using namelist
    eps    = 1.0e-4
    radius = 1.0d0


    ! Generate ellipse
    call generatefestructures(n)
    allocate(FXC(n), FYC(n),X1(n),Y1(n))
    allocate(F1XC(n), F1YC(n), F1ZC(n))    
    allocate(X(n), Y(n), Z(n),nangle(n))
    allocate(XE(n,3), FE(n,3))

    ! Use get positions to get the cooridinates of the particles
    ! Then create masks based on a criteria that gets you the points 
    ! on a sector spanning -30 to 30 degrees and -150 to 150 degrees
    call getpositions(XE,n)
    X = XE(:,1)
    Y = XE(:,2)


    ! Allocate mask
    allocate(isOnPerimeter(n))
    
    ! Mask based on circle condition
    isOnPerimeter = abs(X**2 + Y**2 - radius**2) < eps
    
    ! Mask based on x-location
    isOnLeft = X < 0.0d0
    isOnRight = X > 0.0d0

    ! Mask based on angle
    nangle = atan(Y/X)
    angle45 = 45.0 * pi / 180.0
    isBetweenAngleRange = abs(nangle) <= angle45

    ! Final masks
    isLeftPatch = isOnPerimeter .and. isOnLeft .and. isBetweenAngleRange
    isRightPatch = isOnPerimeter .and. isOnRight .and. isBetweenAngleRange
    
    ! call write_to_file('isLeftPatch.txt', int(merge(1, 0, isLeftPatch),8))
    ! call write_to_file('isRightPatch.txt', int(merge(1, 0, isRightPatch),8))

    XP = pack(X, isOnPerimeter)
    YP = pack(Y, isOnPerimeter)
    xb = pack(X, isRightPatch)
    yb = pack(Y, isRightPatch)

    !Apply boundary forces
    FXC = 500.0d0 !0.250d0
    FYC = 0.0d0
    fxleft = merge(-1.0D0*fxc,0.0d0,isLeftPatch)
    fxright = merge(1.0D0*fxc,0.0d0,isRightPatch)
    fxboundary = fxleft + fxright
    fyboundary = 0.0d0*fxboundary

    niter = 50000
    dt = 0.001
    call cpu_time(t_start)
    do iter = 1, niter
        if (iter .gt. niter/2) then
            fxboundary = 0.0d0
        end if
        ! Call applyboundaryforces inside of calculateforces and make the fx/fyboundary optional
        ! Make sure that fden is initialised as zero at the start of every iteration
        
        call applyboundaryforces(fxboundary,fyboundary,n)
        call calculateforces()
        call getforces(FE,n)
        call updatepositions(dt)
        call getpositions(XE,n)
        
        if (mod(iter,200).eq.0) then
            write(filename, '(A,I8.8,A)') 'F_', iter, '.txt'
            call write_to_file(filename, FE)
            write(filename, '(A,I8.8,A)') 'P_', iter, '.txt'
            call write_to_file(filename, XE)

        end if

        ! Print progress bar every 500 iterations
        if (mod(iter, 1000) == 0 .or. iter == niter) then
            write(*,'(A,I6,A,I6,A,F6.2,A)', advance='no') achar(13)//"Simulation Progress: Iter " &
                , iter, " /", niter, " (", 100.0*iter/niter, "%)"
            if (iter == niter) write(*,*) ! Move to next line at the end
        end if
    end do
    call cpu_time(t_end)
    t_elapsed = t_end - t_start
    print *, "Simulation loop time (seconds): ", t_elapsed

end program main
