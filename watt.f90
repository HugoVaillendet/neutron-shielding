program watt_spectrum
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64         ! Allow float or double precision real
    implicit none

    integer :: n, i
    real(dp) :: a, b, K, L, M, x, y, xi1, xi2
    logical :: accept

    real(dp), allocatable :: E(:)

    n = 10000000

    a = 0.988
    b = 2.249

    K = 1.0_dp + (a * b) / 8.0_dp
    L = a * (K + sqrt(K * K - 1.0_dp))
    M = (L / a) - 1.0_dp

    print *, "K : ", K
    print *, "L : ", L
    print *, "M : ", M

    allocate(E(n))

    E = 0.0_dp

    do i = 1, n
        accept = .false.
        do while (accept .eqv. .false.)

            call random_number(xi1); call random_number(xi2)

            x = -log(xi1)
            y = -log(xi2)

            if ((y - M * (x + 1.0_dp))**2 <= b * L * x) then
                E(i) = L * x
                accept = .true.
            else
                accept = .false.
            end if
        end do
    end do

    open(1, file = 'watt.dat', status = 'replace')
    do i = 1, n
        write(1, *) E(i)
    end do

end program watt_spectrum