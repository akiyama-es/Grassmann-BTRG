module module_time

    contains

    subroutine start_clock()

        use module_declaration

        implicit none

        call system_clock(time_1)

    end subroutine

    subroutine end_clock()

        use module_declaration

        implicit none

        call system_clock(time_2,t_rate,t_max)
        
        if ( time_2 < time_1 ) then

            diff = (t_max - time_1)+time_1+1

        else

            diff = time_2-time_1

        endif

        print*, ' COMPUTATIONAL TIME = ', diff/dble(t_rate)

    end subroutine

end module