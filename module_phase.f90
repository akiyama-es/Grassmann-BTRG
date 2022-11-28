module module_phase

    contains

    subroutine set_basic_phase(phase_1,phase_2)

        use module_declaration

        implicit none

        integer i
        integer, intent(inout) :: phase_1(bond(1)), phase_2(bond(2))

        phase_1 = 1
        phase_2 = 1

        do i = 1, bond_even(1)

            phase_1(i) = 0

        enddo

        do i = 1, bond_even(2)

            phase_2(i) = 0

        enddo

    end subroutine

    subroutine set_basic_phase_new(phase_1,phase_2)

        use module_declaration

        implicit none

        integer i
        integer, intent(inout) :: phase_1(bond_new(1)), phase_2(bond_new(2))

        phase_1 = 1
        phase_2 = 1

        do i = 1, bond_even_new(1)

            phase_1(i) = 0

        enddo

        do i = 1, bond_even_new(2)

            phase_2(i) = 0

        enddo

    end subroutine

end module module_phase