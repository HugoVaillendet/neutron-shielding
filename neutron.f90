module neutron_module
    use vector_module
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64         ! Allow float or double precision real
    implicit none
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)

    type :: neutron
        !SoA structure of position, velocity and energy
        real(dp), allocatable :: x(:), y(:), z(:), vx(:), vy(:), vz(:), E(:)
        !Tally of scatters for each neutron
        integer, allocatable :: scatter_tally(:)
        type(vector_dp) :: fx, fy, fz
        !real(dp), allocatable :: fx(:), fy(:), fz(:)

        !Unit vector of speed and mass
        real(dp) :: m
    contains
        !Initializes the arrays random seed.
        procedure :: init
        ! Samples energy.
        procedure :: sample_energy
        ! Checks if neutron has escaped.
        procedure :: boundary_check
        ! Updates energy after scatter.
        procedure :: update_energy
        ! Samples a new direction.
        procedure :: sample_direction
        ! Samples free path.
        procedure :: sample_free_path
        ! Ealuates one random step.
        procedure :: evaluate_step
    end type neutron

    contains

        subroutine init(this, n, b, L, M, seed)
            class(neutron), intent(inout) :: this

            integer, intent(in) :: n, seed(:)
            real(dp), intent(in) :: b, L, M

            integer :: i

            ! Deallocates arrays if they are aleady allocated.
            if (allocated(this%x)) deallocate(this%x)
            if (allocated(this%y)) deallocate(this%y)
            if (allocated(this%z)) deallocate(this%z)

            if (allocated(this%vx)) deallocate(this%vx)
            if (allocated(this%vy)) deallocate(this%vy)
            if (allocated(this%vz)) deallocate(this%vz)

            if (allocated(this%E)) deallocate(this%E)
            if (allocated(this%scatter_tally)) deallocate(this%scatter_tally)
            
            call this%fx%init_vector()
            call this%fy%init_vector()
            call this%fz%init_vector()

            ! Allocates size.
            allocate(this%x(n), this%y(n), this%z(n), this%vx(n), this%vy(n), this%vz(n), this%E(n))
            allocate(this%scatter_tally(n))
            
            ! Initialises a random seed.
            call random_seed(put=seed)

            ! Assignes all neutrons to (0,0).
            this%x = 0.0_dp
            this%y = 0.0_dp
            this%z = 0.0_dp

            ! Fills energy and velocity arrays.
            do i = 1, n
                call this%sample_direction(i)
                this%E(i) = this%sample_energy(b, L, M)
            end do
            
            ! Set scatter tally to 0.
            this%scatter_tally = 0

        end subroutine init

        real(dp) function sample_energy(this, b, L, M)

            ! We apply here the R12 algorithm for Watt spectrum energy smapling
            ! proposed by C. J. EVERETT and E. CASHWELL.
            ! This is an implementation of their algorithm published in "3rd Monte Carlo Sampler".
            ! b, L and M parameters are all set for U-235 emission spectrum.

            class(neutron), intent(inout) :: this
            real(dp), intent(in) :: b, L, M

            real(dp) :: x, y, xi1, xi2
            logical :: accept
            
            accept = .false.
            ! While accept is false we repeat the random samplign attempt.
            do while (accept .eqv. .false.)

                call random_number(xi1); call random_number(xi2)

                x = -log(xi1)
                y = -log(xi2)

                if ((y - M * (x + 1.0_dp))**2 <= b * L * x) then
                    sample_energy = L * x
                    accept = .true.
                else
                    accept = .false.
                end if
            end do
        end function sample_energy

        logical function boundary_check(this, R, h, i)
            class(neutron), intent(inout) :: this
            integer, intent(in) :: i
            real(dp), intent(in) :: R, h          ! Note h is the half height here.
            real(dp) :: xy2, z2, R2, H2

            ! We test the boudnary against the square of the thickness to avoid sqrt()
            ! and reduce calculation costs.

            R2 = R * R
            H2 = h * h

            xy2 = this%x(i) * this%x(i) + this%y(i) * this%y(i)
            z2 = this%z(i) * this%z(i)

            if (xy2 > R2 .or. (z2) > H2) then
                boundary_check = .true.
                return
            end if
            
            boundary_check = .false.
            return
        end function boundary_check

        subroutine sample_direction(this, i)
            class(neutron), intent(inout) :: this
            integer, intent(in) :: i

            real(dp) :: theta, phi
           
            ! Sample random speed unit vector.
            call random_number(theta); theta = theta * 2 * PI
            call random_number(phi); phi = acos(1.0_dp - 2.0_dp * phi)
            this%vx(i) = cos(theta) * sin(phi)
            this%vy(i) = sin(theta) * sin(phi)
            this%vz(i) = cos(phi)
            
        end subroutine sample_direction

        real(dp) function sample_free_path(this, sigma_t)
            class(neutron), intent(inout) :: this
            real(dp), intent(in) :: sigma_t 
            
            real(dp) :: xi
            
            call random_number(xi)
            xi = max(xi, tiny(1.0_dp))

            ! Sample the free path following the CDF : -log(xi) / sigma_t formula.
            sample_free_path = -log(xi) / sigma_t
            
        end function sample_free_path

        real(dp) function update_energy(this, A, i)
            class(neutron), intent(inout) :: this
            real(dp), intent(in) :: A 
            integer, intent(in) :: i
            
            real(dp) :: u, v, w
            real(dp) :: mu_cm, mu_lab, phi, alpha, E_scatter, C_1, C_2
            real(dp) :: xi1, xi2

            u = this%vx(i)
            v = this%vy(i)
            w = this%vz(i)
            
            alpha = ((A - 1) / (A + 1)) * ((A - 1) / (A + 1))

            call random_number(xi1)
            call random_number(xi2)
            xi1 = max(xi1, tiny(1.0_dp))
            xi2 = max(xi2, tiny(1.0_dp))

            mu_cm = 2.0_dp * xi1 - 1.0_dp
            phi = 2.0_dp * PI * xi2

            ! We update the neutron energy upon scatter.
            E_scatter = 0.5_dp * this%E(i) * ((1.0_dp - alpha) * mu_cm + alpha + 1.0_dp)

            mu_lab = (1.0_dp + A * mu_cm) / sqrt(1.0_dp + A * A + 2.0_dp * A * mu_cm)

            C_1 = sqrt(1 - mu_lab * mu_lab)
            C_2 = sqrt(1 - w * w)
            
            if ((0.9999_dp - abs(w)) > 0 ) then
                this%vx(i) = u * mu_lab + (C_1/C_2) * (u * w * cos(phi) - v * sin(phi))
                this%vy(i) = v * mu_lab + (C_1/C_2) * (v * w * cos(phi) + u * sin(phi))
                this%vz(i) = w * mu_lab - (C_1 * C_2) * cos(phi)
            else
                this%vx(i) = C_1 * cos(phi)
                this%vy(i) = C_1 * sin(phi)
                this%vz(i) = mu_lab * (w / abs(w))
            end if
        
            update_energy = E_scatter
            
        end function update_energy

        integer function evaluate_step(this, R, h, sigma_a, sigma_t, A, i)
            class(neutron), intent(inout) :: this
            integer, intent(in) :: i
            real(dp), intent(in) :: R, h, A, sigma_a(3), sigma_t(3)

            logical :: active
            real(dp) :: p, path, s_a, s_t
            integer :: regime

            active = .true.

            do while (active)

                if (this%E(i) > 100.0e3_dp) then
                    regime = 1
                else if (this%E(i) <= 100.0e3_dp .and. this%E(i) > 0.625_dp) then
                    regime = 2
                else
                    regime = 3
                end if

                select case (regime)
                    case (1)
                        s_a = sigma_a(1)
                        s_t = sigma_t(1)
                    case (2)
                        s_a = sigma_a(2)
                        s_t = sigma_t(2)
                    case (3)
                        s_a = sigma_a(3)
                        s_t = sigma_t(3)
                end select

                path = this%sample_free_path(s_t)
                this%x(i) = this%x(i) + path * this%vx(i)
                this%y(i) = this%y(i) + path * this%vy(i)
                this%z(i) = this%z(i) + path * this%vz(i)
                
                if (this%boundary_check(R, h, i)) then
                    evaluate_step = 0
                    !$omp critical
                    call this%fx%push(this%x(i))
                    call this%fy%push(this%y(i))
                    call this%fz%push(this%z(i))
                    !$omp end critical
                    return
                end if

                call random_number(p)
                if (p < s_a / s_t) then
                    evaluate_step = 1
                    return
                else
                    this%scatter_tally(i) = this%scatter_tally(i) + 1
                    this%E(i) = this%update_energy(A, i)
                end if
            end do

        end function evaluate_step

end module neutron_module