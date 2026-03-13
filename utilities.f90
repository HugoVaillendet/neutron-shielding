module vector_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    type :: vector_dp
        real(dp), allocatable :: x(:)
        integer :: size
        integer :: capacity
    contains
        procedure :: init_vector
        procedure :: push
        procedure :: grow
        procedure :: at        
    end type vector_dp

    contains

        subroutine init_vector(this)
            class(vector_dp), intent(inout) :: this
            if (allocated(this%x)) deallocate(this%x)
            this%capacity = 16
            allocate(this%x(this%capacity))
            this%x = 0.0_dp
            this%size = 0
        end subroutine init_vector

        subroutine grow(this)
            class(vector_dp), intent(inout) :: this
            real(dp), allocatable :: nx(:)
            allocate(nx(2 * this%capacity))
            nx = 0.0_dp
            nx(1:this%capacity) = this%x
            this%capacity = 2 * this%capacity

            call move_alloc(nx, this%x)
        end subroutine grow

        subroutine push(this, x)
            class(vector_dp), intent(inout) :: this
            real(dp), intent(in) :: x
            if (this%size == this%capacity) call grow(this)
            this%size = this%size + 1
            this%x(this%size) = x
        end subroutine push

        real(dp) function at(this, i) result(val)
            class(vector_dp), intent(in) :: this
            integer, intent(in) :: i

            if(i <= this%size .and. i >= 1) then
                val = this%x(i)
            else
                val = 0.0_dp
            end if
        end function at

end module vector_module