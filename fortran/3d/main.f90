program main
    use fem3d_interface, only: generatefestructures, getpositions, &
        calculateforces, getforces, updatepositions 
    use iso_c_binding, only: c_int, c_double, c_loc
    use matrix_writer, only: write_to_file
    implicit none

    real(8), parameter :: PI = 3.141592653589793
    integer(C_INT) :: n, niter, iter
    logical, allocatable :: isOnSurface(:), isBetweenAngleRange(:), & 
                                isOnTop(:), isOnBottom(:), isTopPatch(:), isBottomPatch(:)
    real(c_double) :: eps, dt, a, b, c
    real(c_double), allocatable :: XN(:,:), FN(:,:)
    real(c_double), allocatable :: fzmag(:), fztop(:), fzbottom(:), fxboundary(:), fyboundary(:), fzboundary(:), &
                                    fboundary(:,:)
    real(c_double), allocatable :: X(:),Y(:),Z(:), nangle(:)
    real(c_double) :: angle45
 
    ! Measure time
    real :: t_start, t_end, t_elapsed
    character(len=20) :: filename

    ! Particle geometrical parameters
    eps = 1.0e-4
    a = 1.0d0
    b = 1.0d0
    c = 1.0d0

    ! Generate fe structures
    call generatefestructures(n)
    allocate(fzmag(n), X(n), Y(n), Z(n), nangle(n), XN(n,3), FN(n,3), fboundary(n,3))

    !------------------------- Boundary conditions -------------------------!
    ! Use get positions to get the cooridinates of the particles
    ! Then create masks based on a criteria that gets you the points 
    ! on a sector spanning -30 to 30 degrees and -150 to 150 degrees
    call getpositions(XN,n)
    X = XN(:,1)
    Y = XN(:,2)
    Z = XN(:,3)
    iter = 0
    write(filename, '(A,I8.8,A)') 'F_', iter, '.txt'
    call write_to_file(filename, FN)
    write(filename, '(A,I8.8,A)') 'P_', iter, '.txt'
    call write_to_file(filename, XN)
    
    ! Mask based on ellipsoid condition
    isOnSurface = abs((X**2)/(a**2) + (Y**2)/(b**2) + (Z**2)/(c**2) - 1.0) < eps
    ! Mask based on z-location
    isOnTop = Z > 0.40d0
    isOnBottom = Z < -0.40d0

    ! Final masks
    isTopPatch = isOnSurface .and. isOnTop
    isBottomPatch = isOnSurface .and. isOnBottom
    !Apply boundary forces
    fzmag = 200.0d0
    fztop = merge(1.0D0*fzmag,0.0d0,isTopPatch)
    fzbottom = merge(-1.0D0*fzmag,0.0d0,isBottomPatch)
    fzboundary = fztop + fzbottom
    fxboundary = 0.0d0*fzboundary
    fyboundary = 0.0d0*fzboundary

    fboundary(:,1) = fxboundary
    fboundary(:,2) = fyboundary
    fboundary(:,3) = fzboundary

    !-----------------------------------------------------------------------------!

    niter = 25000
    dt = 0.001
    call cpu_time(t_start)
    do iter = 1, niter
        if (iter .gt. niter/3) then
            fzboundary = 0.0d0
        end if

        call calculateforces(fboundary,n)
        call getforces(FN,n)
        call updatepositions(dt)
        call getpositions(XN,n)
        ! print *, maxval(XN(:,1)), maxval(XN(:,2)), maxval(XN(:,3))
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
