program main
    use fem_interface, only: generatefestructures, getpositions, &
        calculateforces, getforces, updatepositions 
    use iso_c_binding, only: c_int, c_double, c_loc
    use matrix_writer, only: write_to_file
    implicit none

    real(8), parameter :: PI = 3.141592653589793
    integer(C_INT) :: n, niter, iter
    logical, allocatable :: isOnPerimeter(:), isBetweenAngleRange(:), & 
                                isOnLeft(:), isOnRight(:), isLeftPatch(:), isRightPatch(:)
    real(c_double) :: eps, dt, radius
    real(c_double), allocatable :: XN(:,:), FN(:,:)
    real(c_double), allocatable :: fxmag(:),fymag(:), fxleft(:), fxright(:), fxboundary(:), fyboundary(:)
    real(c_double), allocatable :: X(:),Y(:),Z(:), nangle(:)
    real(c_double) :: angle45
 
    ! Measure time
    real :: t_start, t_end, t_elapsed
    character(len=20) :: filename

    ! Particle geometrical parameters
    eps    = 1.0e-4
    radius = 1.0d0

    ! Generate fe structures
    call generatefestructures(n)
    allocate(fxmag(n), fymag(n), X(n), Y(n), Z(n), nangle(n), XN(n,3), FN(n,3))

    !------------------------- Boundary conditions -------------------------!
    ! Use get positions to get the cooridinates of the particles
    ! Then create masks based on a criteria that gets you the points 
    ! on a sector spanning -30 to 30 degrees and -150 to 150 degrees
    call getpositions(XN,n)
    X = XN(:,1)
    Y = XN(:,2)
    
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
    
    !Apply boundary forces
    fxmag = 500.0d0
    fymag = 0.0d0
    fxleft = merge(-1.0D0*fxmag,0.0d0,isLeftPatch)
    fxright = merge(1.0D0*fxmag,0.0d0,isRightPatch)
    fxboundary = fxleft + fxright
    fyboundary = 0.0d0*fxboundary
    !-----------------------------------------------------------------------------!

    niter = 50000
    dt = 0.001
    call cpu_time(t_start)
    do iter = 1, niter
        if (iter .gt. niter/2) then
            fxboundary = 0.0d0
        end if
        
        call calculateforces(fxboundary,fyboundary,n)
        call getforces(FN,n)
        call updatepositions(dt)
        call getpositions(XN,n)
        
        if (mod(iter,200).eq.0) then
            write(filename, '(A,I8.8,A)') 'F_', iter, '.txt'
            call write_to_file(filename, FN)
            write(filename, '(A,I8.8,A)') 'P_', iter, '.txt'
            call write_to_file(filename, XN)

        end if

        ! Print progress bar every nth iterations
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
