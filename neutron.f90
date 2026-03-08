module neutron_module
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64         ! Allow float or double precision real
    implicit none
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)

    type :: neutron
        !SoA structure of position and energy
        real(dp), allocatable :: x(:), y(:), E(:)
        !Tally of scatters for each neutron
        integer, allocatable :: scatter_tally(:)

        !Unit vector of speed and mass
        real(dp) :: vx, vy, m
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
            if (allocated(this%E)) deallocate(this%E)
            if (allocated(this%scatter_tally)) deallocate(this%scatter_tally)

            ! Allocates size.
            allocate(this%x(n), this%y(n), this%E(n))
            allocate(this%scatter_tally(n))

            ! Assignes all neutrons to (0, 0).
            this%x = 0.0_dp
            this%y = 0.0_dp

            ! Fills energy array.
            do i = 1, n
                this%E(i) = this%sample_energy(b, L, M)
            end do
            
            ! Set scatter tally to 0.
            this%scatter_tally = 0

            ! Initialises a random seed.
            call random_seed(put=seed)

            !Samples a random direction.
            call this%sample_direction()

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

        logical function boundary_check(this, thickness, i)
            class(neutron), intent(inout) :: this
            integer, intent(in) :: i
            real(dp), intent(in) :: thickness
            real(dp) :: R2, thickness2

            ! We test the boudnary against the sqaure of the thickness to avoid sqrt()
            ! and reduce calculation costs.

            thickness2 = thickness * thickness
            R2 = this%x(i) * this%x(i) + this%y(i) * this%y(i)

            if (R2 > thickness2) then
                boundary_check = .true.
                return
            end if
            
            boundary_check = .false.
            return
        end function boundary_check

        subroutine sample_direction(this)
            class(neutron), intent(inout) :: this

            real(dp) :: theta
           
            ! Sample random speed unit vector.
            call random_number(theta); theta = theta * 2 * PI
            this%vx = cos(theta)
            this%vy = sin(theta)
            
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
            
            real(dp) :: xi, alpha, r
            
            alpha = ((A - 1) / (A + 1)) * ((A - 1) / (A + 1))

            call random_number(xi)
            xi = max(xi, tiny(1.0_dp))

            ! We update the neutron energy upon scatter.
            r = alpha + (1.0_dp - alpha) * xi
            
            update_energy = r * this%E(i)
            
        end function update_energy

        integer function evaluate_step(this, thickness, sigma_a, sigma_t, A, i)
            class(neutron), intent(inout) :: this
            integer, intent(in) :: i
            real(dp), intent(in) :: thickness, A, sigma_a(3), sigma_t(3)

            logical :: active
            real(dp) :: p, path, s_a, s_t
            integer :: regime

            active = .true.

            call this%sample_direction()

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
                this%x(i) = this%x(i) + path * this%vx
                this%y(i) = this%y(i) + path * this%vy
                
                if (this%boundary_check(thickness, i)) then
                    evaluate_step = 0
                    return
                end if

                call random_number(p)
                if (p < s_a / s_t) then
                    evaluate_step = 1
                    return
                else
                    call this%sample_direction()
                    this%scatter_tally(i) = this%scatter_tally(i) + 1
                    this%E(i) = this%update_energy(A, i)
                end if
            end do

        end function evaluate_step

end module neutron_module