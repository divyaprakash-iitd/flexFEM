program ibmc
    use soft_particles
    use iso_c_binding, only: C_INT
    implicit none

    integer(C_INT)  :: n
    call sayhello()
    call generateellipse(n)
    print *, "n = ", n

end program ibmc
