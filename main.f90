program monte_carlo
    use neutron_module
    use omp_lib
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64         ! Allow float or double precision real
    implicit none

    !--------------------- Declarations ---------------------
    type(neutron) :: ntr
    integer :: nseed, n, i, transmitted, scattered, absorbed, eval
    integer, allocatable :: seed(:)
    real(dp) :: R, mass, thickness, E_0, A_, a, b, K, L, M

    real(dp) :: sigma_a(3), sigma_t(3)

    !--------------------- Execution ---------------------
    n = 1000000
    transmitted = 0
    absorbed = 0
    scattered = 0

    R = 10.0_dp
    mass = 1.0_dp
    thickness = 10.0_dp
    E_0 = 2.3e6_dp
    A_ = 1.0_dp

    a = 0.988
    b = 2.249

    K = 1.0_dp + (a * b) / 8.0_dp
    L = a * (K + sqrt(K * K - 1.0_dp))
    M = (L / a) - 1.0_dp

    ! - Light Water - 
    sigma_a = [3.0e-4_dp, 3.0e-3_dp, 2.2e-2_dp]
    sigma_t = [6.6e-1_dp, 1.2_dp,1.5_dp]

    ! - Concrete -
    !sigma_a = 0.003_dp
    !sigma_t = 0.10_dp

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 12345

    call ntr%init(n, b, L, M, seed)

    do i = 1, n
        eval = ntr%evaluate_step(thickness, sigma_a, sigma_t, A_, i)
        select case (eval)
            case (0)
                transmitted = transmitted + 1
            case (1)
                absorbed = absorbed + 1
        end select
    end do

    print *, "Transmitted : ", transmitted, "hits", real(transmitted, dp) / real(n, dp) * 100.0_dp, "%"
    print *, "Absorbed : ", absorbed, "hits", real(absorbed, dp) / real(n, dp) * 100.0_dp, "%"
    print *, "Total scattering events : ", sum(ntr%scatter_tally), "hits"

    open(1, file = 'scatter_histogram.dat', status = 'replace')
    open(2, file = 'energy_histogram.dat', status = 'replace')    
    do i = 1, n  
        write(1,*) ntr%scatter_tally(i)  
        write(2,*) ntr%E(i)  
    end do
    close(1)
    close(2)  

end program monte_carlo