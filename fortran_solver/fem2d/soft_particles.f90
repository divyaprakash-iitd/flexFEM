module soft_particles
    use fem2d
    use mesh_module, only: get_nodes_connectivity, get_physical_group_nodes
    use matrix_writer, only: write_to_file
    use, intrinsic :: iso_c_binding

    implicit none

    ! Parameters
    real(8), parameter :: PI = 3.141592653589793
    
    ! Particle information
    real(8)         :: Kp, Bp

    ! FEM data
    type(festruct), allocatable :: structures(:)
    integer(int32)              :: ntri, npp, nparticle
    real(real64), allocatable   :: pb(:,:,:)
    integer,      allocatable   :: mp(:,:)
    real(real64), allocatable   :: paelem(:)
    real(real64), allocatable   :: pp(:,:)
    real(real64), allocatable   :: FN(:,:)
    real(real64), allocatable   :: UN(:,:)
    logical, allocatable        :: pboundary(:,:)
    integer(int32)              :: femdata(2)

    ! Namelists for input
    namelist /particleprops/ Kp, Bp

    !---------------------- Begin Calculations ------------------------------------!

    contains 

    subroutine generatefestructures(nnodes) bind(C)
        use iso_c_binding, only: C_INT, C_CHAR
        implicit none

        integer(C_INT), intent(inout) :: nnodes
        integer(c_size_t), allocatable :: nodeTagsAll(:)
        integer(c_size_t), allocatable :: connectivity(:)
        real(c_double), allocatable :: coordAll(:)
        integer :: ierr


        call get_nodes_connectivity("donut2d_mesh.msh", connectivity, nodeTagsAll, coordAll, ierr)
        pp = transpose(reshape(coordAll,[3,size(coordAll)/3]))
        mp = transpose(reshape(connectivity,[3,size(connectivity)/3]))
        
        FN = pp*0.0d0
        UN = pp*0.0d0

        ! Read input data from file
        open(1004,file="input_params.dat",form='formatted')
        READ(unit=1004,nml=particleprops,iostat=ierr)
        close(1004)
    
        ! Construct the structure
        nparticle   = 1
        allocate(structures(nparticle))
        structures(1) = festruct(MP,PP,FN,UN,bp,kp,1.0d0) ! kp = kval, co = bp , dl = 1.0d0

        nnodes = size(structures(1)%XE,1)

        call write_to_file('connectivity.txt', reshape(connectivity,[3,size(connectivity)/3]))
    
    end subroutine generatefestructures

    subroutine getforces(FE,nn) bind(C)
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)      :: nn
        real(c_double), intent(inout)   :: FE(nn,3)

        integer(int32) :: nparticles, npoints, i

        nparticles = 1
        npoints = size(structures(1)%XE,1)

        do i = 1,npoints
            FE(i,1)   = structures(nparticles)%fden(i,1)
            FE(i,2)   = structures(nparticles)%fden(i,2)
            FE(i,3)   = 0.0d0 ! 2D problem, so Z is always 0
        end do
    end subroutine getforces

    subroutine getpositions(XE, nn) bind(C)
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)      :: nn
        real(c_double), intent(inout)   :: XE(nn,3)

        integer(int32) :: nparticles, npoints, i

        nparticles = 1

        npoints = size(structures(1)%XE,1)

        do i = 1,npoints
            XE(i,1)   = structures(nparticles)%XE(i,1)
            XE(i,2)   = structures(nparticles)%XE(i,2)
            XE(i,3)   = 0.0d0 ! 2D problem, so Z is always 0
        end do
    end subroutine getpositions
    
    subroutine applyboundaryforces(FXC,FYC,nn) bind(C)
        ! Applies the boundary forces to the particle
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        integer(c_int), intent(in)   :: nn
        real(c_double), intent(in)   :: FXC(nn),FYC(nn)

        integer(int32) :: i, npoints

        npoints = size(structures(1)%XE,1)

        ! Initialize the force arrays
        structures(1)%fden = 0.0d0

        ! Apply boundary forces
        do i = 1, nn
            structures(1)%fden(i,1) = FXC(i)
            structures(1)%fden(i,2) = FYC(i)
        end do

    end subroutine applyboundaryforces

   
    subroutine calculateforces() bind(C)
        ! Calculates the forces in the particle
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        call structures(1)%calculate_forces()

    end subroutine calculateforces


    subroutine updatepositions(dt) bind(C)
        use iso_c_binding, only: c_int, c_double, c_loc
        implicit none

        real(c_double), intent(in)      :: dt

        integer(int32) :: i, nparticles, npoints
        real(c_double) :: mu

        mu = 1000

        nparticles = 1

        npoints = size(structures(1)%XE,1)

        do i = 1,npoints
            structures(1)%U(i,1) = structures(1)%fden(i,1) / mu
            structures(1)%U(i,2) = structures(1)%fden(i,2) / mu
        end do


        do i = 1,nparticles
            call structures(i)%update_position(dt)
        end do

    end subroutine updatepositions

end module soft_particles
