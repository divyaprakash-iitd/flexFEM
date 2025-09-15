program main
    use fem3d_interface, only: generatefestructures, getpositions, &
        calculateforces, getforces, updatepositions 
    use iso_c_binding, only: c_int, c_double, c_loc
    use matrix_writer, only: write_to_file
    implicit none

    real(8), parameter :: PI = 3.141592653589793
    integer(C_INT) :: n, niter, iter, i
    logical, allocatable :: isOnSurface(:), isBetweenAngleRange(:), & 
                                isOnTop(:), isOnBottom(:), isTopPatch(:), isBottomPatch(:)
    real(c_double) :: eps, dt, a, b, c, radius, height, dvec(3)
    real(c_double), allocatable :: XN(:,:), FN(:,:), fboundary(:,:), XORG(:,:), & 
                                    topcircle(:,:)
    real(c_double), allocatable :: fzmag(:), fztop(:), fzbottom(:), fxboundary(:), fyboundary(:), & 
                                    fzboundary(:), fpenalty(:,:), fbending(:,:)
    real(c_double), allocatable :: X(:),Y(:),Z(:), nangle(:)
    real(c_double) :: angle45, kpenalty
 
    ! Measure time
    real :: t_start, t_end, t_elapsed
    character(len=20) :: filename

    ! Particle geometrical parameters
    eps = 1.0e-4
    a = 1.0d0
    b = 1.0d0
    c = 1.0d0
    radius = 0.1d0
    height = 2.0d0

    ! Anchor properties
    kpenalty = 10.0d0 ! Penalty force constant

    ! Generate fe structures
    call generatefestructures(n)
    allocate(fzmag(n), X(n), Y(n), Z(n), nangle(n), XN(n,3), & 
                FN(n,3), fboundary(n,3), XORG(n,3), fpenalty(n,3), fbending(n,3))

    !------------------------- Boundary conditions -------------------------!
    ! Use get positions to get the cooridinates of the particles
    ! Then create masks based on a criteria that gets you the points 
    ! on a sector spanning -30 to 30 degrees and -150 to 150 degrees
    call getpositions(XN,n)
    XORG = XN ! Store original positions for later use
    X = XN(:,1)
    Y = XN(:,2)
    Z = XN(:,3)
    iter = 0
    write(filename, '(A,I8.8,A)') 'F_', iter, '.txt'
    call write_to_file(filename, FN)
    write(filename, '(A,I8.8,A)') 'P_', iter, '.txt'
    call write_to_file(filename, XN)
    
    ! Mask for top and bottom patches
    isOnTop = abs(Z-height) < eps
    isOnBottom = abs(Z) < eps

    topcircle(:,1) = pack(X,isOnTop)
    topcircle(:,2) = pack(Y,isOnTop)
    topcircle(:,3) = pack(X,isOnTop)

    !Apply boundary forces
    fzmag = 100.0d0
    fztop = merge(1.0D0*fzmag,0.0d0,isOnTop)
    fzbottom = merge(-1.0D0*fzmag,0.0d0,isOnBottom)
    ! To-Do: Apply anchoring force on the bottom
    fboundary = 0.0d0
    fboundary(:,3) = fztop + fzbottom*0.0d0
    !-----------------------------------------------------------------------------!

    niter = 25000
    dt = 0.001
    call cpu_time(t_start)
    do iter = 1, niter
        if (iter .gt. niter/2) then
            fboundary = 0.0d0
        end if
        ! Calculate forces here which uses the difference between the original and the displaced
        ! position of the nodes along the x, y and z direction.
        ! PRIORITY: Refactor the code to use force vector instead of individual comments
        fpenalty = 0.0d0
        where(isOnBottom) 
          fpenalty(:,1) = - kpenalty*abs((xorg(:,1) - xn(:,1)))
          fpenalty(:,2) = - kpenalty*abs((xorg(:,2) - xn(:,2))) 
          fpenalty(:,3) = - kpenalty*abs((xorg(:,3) - xn(:,3)))
        end where

        ! Apply the bending forces
        fbending = 0.0d0
        dvec = calculate_direction_vector(topcircle)
        ! where(isOnTop)
        ! end where

        ! fbend = getparallelforces(x,y,z)
        ! end where
        ! Apply the twisting forces
        ! where(isOnTop)
        ! fbend = getorthogonalforces(x,y,z)
        ! end where
        ! To-Do: Implement a C++ API
        call calculateforces(fboundary+fpenalty,n)
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

contains
    function calculate_direction_vector(patch) result(direction_vector)
    implicit none
    real(8), intent(in) :: patch(:,:)
    real(8) :: point1(3), point2(3), direction_vector(3)
    real(8) :: norm_factor
   
   
    point1 = [patch(1,1), patch(1,2), patch(1,3)]
    point2 = [patch(5,1), patch(5,2), patch(5,3)]

    ! Calculate direction vector (point2 - point1)
    direction_vector(1) = point2(1) - point1(1)  ! dx
    direction_vector(2) = point2(2) - point1(2)  ! dy  
    direction_vector(3) = point2(3) - point1(3)  ! dz
    
    ! Calculate magnitude for normalization
    norm_factor = sqrt(direction_vector(1)**2 + direction_vector(2)**2 + direction_vector(3)**2)
    
    direction_vector = direction_vector / norm_factor
    
end function calculate_direction_vector
end program main
