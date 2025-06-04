module mod_solid
    use iso_fortran_env, only: int32, real64
    implicit none

    ! real(real64), parameter :: PI = 3.141592653589793
    real(real64), parameter :: PI=4.D0*DATAN(1.0D0)

    private :: PI

    type, abstract :: solid
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
    end type solid 

end module mod_solid
