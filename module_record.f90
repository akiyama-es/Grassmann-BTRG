module module_record

    contains

    subroutine print_result()

        use module_declaration

        implicit none

        if ( iter == 0 ) then

            if ( coupling == 0d0 ) then

                if ( flag_exact == 'y' ) then

                    write(*,"( a5, a1, a16, a1, a16 )") 'ITER', ' ', 'lnZ/V', ' ', 'Relative Error'

                elseif ( flag_exact == 'n' ) then

                    write(*,"( a5, a1, a16 )") 'ITER', ' ', 'lnZ/V'

                endif

            elseif ( coupling /= 0d0 ) then

                write(*,"( a5, a1, a16 )") 'ITER', ' ', 'lnZ/V'

            endif

        endif

        if ( coupling /= 0d0 ) then

            write(*,"( i5, a1, f16.13 )") iter, ' ', real(lnz(iter))

        elseif ( coupling == 0d0 ) then

            if ( mod(iter,2) == 0 ) then

                if ( flag_exact == 'y' ) then

                    write(*,"( i5, a1, f16.13, a1, f16.13 )") iter, ' ', real(lnz(iter)), ' ', relative_error(iter)

                elseif ( flag_exact == 'n' ) then

                    write(*,"( i5, a1, f16.13 )") iter, ' ', real(lnz(iter))
    
                endif

            elseif ( mod(iter,2) /= 0 ) then

                write(*,"( i5, a1, f16.13 )") iter, ' ', real(lnz(iter))

            endif

        endif

    end subroutine 

    subroutine record_svd_result(s,size_s,s_trunc,size_s_trunc,size_even,index,iteration)

        !use module_declaration

        implicit none

        character(len=256) file
        character(len=5) index
        integer, intent(in) :: size_s, size_s_trunc, size_even, iteration
        integer i
        double precision tsum, psum, max_singular
        double precision, intent(in) :: s(size_s), s_trunc(size_s_trunc)

        tsum = 0d0

        do i = 1, size_s
            
            tsum = tsum + s(i)*s(i)
            
        enddo

        psum = 0d0

        do i = 1, size_s_trunc
            
            psum = psum + s_trunc(i)*s_trunc(i)
            
        enddo

        write(file,"(a,a,a,i2.2)")  "svd_", index, "_", iteration

        if ( size(s) /= size_even ) then

            max_singular = max(s(1),s(size_even+1))

        else 

            max_singular = s(1)

        endif

        open(90,file=file)

        if ( max_singular /= 0d0 ) then

            do i = 1, size_s
            
                write(90,*) i, s(i), s(i)/max_singular
            
            enddo

            write(90,*) ' '

            !write(20,*) ' ORIGINAL NORM = ', tsum
            write(90,*) ' % ', sqrt(psum/tsum)*100d0

            write(90,*) ' '

        else

            write(90,*) ' zero singular value '

        endif

        close(90)

    end subroutine

end module module_record
