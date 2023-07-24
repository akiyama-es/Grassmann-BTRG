module module_phase

    contains

    subroutine set_phase(size_total,size_even,phase)

        implicit none

        integer, intent(in) :: size_total, size_even
        integer, intent(out) :: phase(size_total)
        integer i

        phase = 1

        do i = 1, size_even

            phase(i) = 0

        enddo

    end subroutine

end module module_phase