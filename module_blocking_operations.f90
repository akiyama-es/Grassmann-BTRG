module module_blocking_operations

    contains

    subroutine blocked_svd(size_matrix,size_even_matrix,size_even_result,matrix,s,u,v)

        use module_lapack

        implicit none

        integer, intent(in)  :: size_matrix(4)
        integer, intent(in)  :: size_even_matrix(4)
        complex(kind(0d0)), intent(in)  :: matrix(size_matrix(1)*size_matrix(2)*size_matrix(3)*size_matrix(4))

        integer, intent(out) :: size_even_result
        double precision, intent(out)   :: s(min(size_matrix(1)*size_matrix(2),size_matrix(3)*size_matrix(4)))
        complex(kind(0d0)), intent(out) :: u(size_matrix(1)*size_matrix(2)*size(s))
        complex(kind(0d0)), intent(out) :: v(size_matrix(3)*size_matrix(4)*size(s)) !!! NOT V_TRANSPOSED

        integer size_odd_matrix(4)
        integer size_even_row, size_odd_row
        integer size_even_column, size_odd_column
        integer size_1(2), size_2(2)
        integer size_even_1(2), size_even_2(2)
        integer size_odd_result

        complex(kind(0d0)), allocatable :: matrix_even(:), matrix_odd(:)
        double precision, allocatable   :: s_even(:), s_odd(:)
        complex(kind(0d0)), allocatable :: u_even(:), u_odd(:)
        complex(kind(0d0)), allocatable :: v_even(:), v_odd(:)

        integer i

        do i = 1, 4

            size_odd_matrix(i) = size_matrix(i) - size_even_matrix(i)

        enddo

        do i = 1, 2

            size_1(i) = size_matrix(i)

            size_2(i) = size_matrix(i+2)

            size_even_1(i) = size_even_matrix(i)

            size_even_2(i) = size_even_matrix(i+2)

        enddo

        size_even_row = size_even_matrix(1)*size_even_matrix(2)+size_odd_matrix(1)*size_odd_matrix(2)

        size_odd_row  = size_even_matrix(1)*size_odd_matrix(2)+size_odd_matrix(1)*size_even_matrix(2)

        size_even_column = size_even_matrix(3)*size_even_matrix(4)+size_odd_matrix(3)*size_odd_matrix(4)

        size_odd_column  = size_even_matrix(3)*size_odd_matrix(4)+size_odd_matrix(3)*size_even_matrix(4)

        size_even_result = min(size_even_row,size_even_column)

        size_odd_result  = min(size_odd_row,size_odd_column)

        !!!!!!!!!!

        if ( size(matrix) /= (size_even_row+size_odd_row)*(size_even_column+size_odd_column) ) then

            print *, ' ERROR IN SIZE OF THE MATRIX IN BLOCKED SVD '
            stop

        endif

        allocate(matrix_even(size_even_row*size_even_column))
        allocate(s_even(min(size_even_row,size_even_column)))
        allocate(u_even(size_even_row*min(size_even_row,size_even_column)))
        allocate(v_even(size_even_column*min(size_even_row,size_even_column)))

        matrix_even = (0d0,0d0)
        s_even = 0d0
        u_even = (0d0,0d0)
        v_even = (0d0,0d0)

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            allocate(matrix_odd(size_odd_row*size_odd_column))
            allocate(s_odd(min(size_odd_row,size_odd_column)))
            allocate(u_odd(size_odd_row*min(size_odd_row,size_odd_column)))
            allocate(v_odd(size_odd_column*min(size_odd_row,size_odd_column)))

        else

            allocate(matrix_odd(0:size_odd_row*size_odd_column))
            allocate(s_odd(0:min(size_odd_row,size_odd_column)))
            allocate(u_odd(0:size_odd_row*min(size_odd_row,size_odd_column)))
            allocate(v_odd(0:size_odd_column*min(size_odd_row,size_odd_column)))
            
        endif

        matrix_odd = (0d0,0d0)
        s_odd = 0d0
        u_odd = (0d0,0d0)
        v_odd = (0d0,0d0)

        call blocking_matrix(size_matrix,size_even_matrix,size_even_row,size_odd_row,size_even_column,size_odd_column,matrix,matrix_even,matrix_odd)

        call zsvd(matrix_even,size_even_row,size_even_column,s_even,u_even,v_even)

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            call zsvd(matrix_odd,size_odd_row,size_odd_column,s_odd,u_odd,v_odd)

        endif

        call unblocking_singular(size_even_row,size_odd_row,size_even_column,size_odd_column,s_even,s_odd,s)

        call unblocking_unitary(size_1,size_even_1,size_even_row,size_odd_row,size_even_column,size_odd_column,u_even,u_odd,u)

        call unblocking_unitary(size_2,size_even_2,size_even_row,size_odd_row,size_even_column,size_odd_column,v_even,v_odd,v)

        deallocate(matrix_even,s_even,u_even,v_even)
        deallocate(matrix_odd,s_odd,u_odd,v_odd)

    end subroutine 

    subroutine blocking_matrix(size_matrix,size_even_matrix,size_even_row,size_odd_row,size_even_column,size_odd_column,matrix,matrix_even,matrix_odd)

        implicit none

        integer, intent(in) :: size_matrix(4)
        integer, intent(in) :: size_even_matrix(4)
        integer, intent(in) :: size_even_row, size_odd_row, size_even_column, size_odd_column
        complex(kind(0d0)), intent(in) :: matrix(size_matrix(1)*size_matrix(2)*size_matrix(3)*size_matrix(4))

        complex(kind(0d0)), intent(out) :: matrix_even(size_even_row*size_even_column)
        complex(kind(0d0)), intent(out) :: matrix_odd(size_odd_row*size_odd_column)

        integer i, j, k
        integer size_odd_matrix(4)
        integer size_double(2), size_even_double(2), size_odd_double(2)
        
        complex(kind(0d0)), allocatable :: column_vect(:), even_vect(:), odd_vect(:)

        do i = 1, 4

            size_odd_matrix(i) = size_matrix(i) - size_even_matrix(i)

        enddo

        do i = 1, 2

            size_double(i) = size_matrix(i)

            size_even_double(i) = size_even_matrix(i)

            size_odd_double(i) = size_double(i) - size_even_double(i)

        enddo


        allocate(column_vect(size_even_row+size_odd_row))
        allocate(even_vect(size_even_row))

        column_vect = (0d0,0d0)
        even_vect = (0d0,0d0)

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            allocate(odd_vect(size_odd_row))

            odd_vect = (0d0,0d0)

        endif

        
        k = 1

        do while ( k <= size_even_matrix(4) )

            do j = 1, size_even_matrix(3)

                do i = 1, size_matrix(1)*size_matrix(2)

                    column_vect(i) = matrix(i+(j-1)*size_matrix(1)*size_matrix(2)+(k-1)*size_matrix(1)*size_matrix(2)*size_matrix(3))

                enddo

                call gather_even_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,size_odd_column,column_vect,even_vect)

                do i = 1, size_even_row

                    matrix_even(i+(j-1)*size_even_row+(k-1)*size_even_row*size_even_matrix(3)) = even_vect(i)

                enddo

                column_vect = (0d0,0d0)
                even_vect = (0d0,0d0)

            enddo

            if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

                do j = 1, size_odd_matrix(3)

                    do i = 1, size_matrix(1)*size_matrix(2)

                        column_vect(i) = matrix(i+((j+size_even_matrix(3))-1)*size_matrix(1)*size_matrix(2)+(k-1)*size_matrix(1)*size_matrix(2)*size_matrix(3))

                    enddo

                    call gather_odd_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,column_vect,odd_vect)

                    do i = 1, size_odd_row

                        matrix_odd(i+(j-1)*size_odd_row+(k-1)*size_odd_row*size_odd_matrix(3)) = odd_vect(i)

                    enddo

                    column_vect = (0d0,0d0)
                    odd_vect = (0d0,0d0)

                enddo

            endif

            k = k + 1

        enddo

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            k = 1

            do while ( k <= size_odd_matrix(4) )

                do j = 1, size_even_matrix(3)

                    do i = 1, size_matrix(1)*size_matrix(2)

                        column_vect(i) = matrix(i+(j-1)*size_matrix(1)*size_matrix(2)+((k+size_even_matrix(4))-1)*size_matrix(1)*size_matrix(2)*size_matrix(3))

                    enddo

                    call gather_odd_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,column_vect,odd_vect)

                    do i = 1, size_odd_row

                        matrix_odd(i+(j-1)*size_odd_row+(k-1)*size_odd_row*size_even_matrix(3)+size_odd_row*size_odd_matrix(3)*size_even_matrix(4)) = odd_vect(i)

                    enddo

                    column_vect = (0d0,0d0)
                    odd_vect = (0d0,0d0)

                enddo

                do j = 1, size_odd_matrix(3)

                    do i = 1, size_matrix(1)*size_matrix(2)

                        column_vect(i) = matrix(i+((j+size_even_matrix(3))-1)*size_matrix(1)*size_matrix(2)+((k+size_even_matrix(4))-1)*size_matrix(1)*size_matrix(2)*size_matrix(3))

                    enddo

                    call gather_even_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,size_odd_column,column_vect,even_vect)

                    do i = 1, size_even_row

                        matrix_even(i+(j-1)*size_even_row+(k-1)*size_even_row*size_odd_matrix(3)+size_even_row*size_even_matrix(3)*size_even_matrix(4)) = even_vect(i)

                    enddo

                    column_vect = (0d0,0d0)
                    even_vect = (0d0,0d0)

                enddo

                k = k + 1

            enddo

        endif



        deallocate(column_vect,even_vect)

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            deallocate(odd_vect)

        endif

    end subroutine 

    subroutine gather_even_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,size_odd_column,input_column,output_column)

        implicit none

        integer, intent(in) :: size_double(2), size_even_double(2), size_odd_double(2)
        integer, intent(in) :: size_even_row, size_odd_row, size_odd_column
        complex(kind(0d0)), intent(in)  :: input_column(size_even_row+size_odd_row)
        complex(kind(0d0)), intent(out) :: output_column(size_even_row)

        integer i, j
        
        do j = 1, size_even_double(2)

            do i = 1, size_even_double(1)

                output_column(i+(j-1)*size_even_double(1)) = input_column(i+(j-1)*size_double(1))

            enddo

        enddo

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            do j = 1, size_odd_double(2)

                do i = 1, size_odd_double(1)

                    output_column(size_even_double(1)*size_even_double(2)+i+(j-1)*size_odd_double(1)) = input_column((i+size_even_double(1))+((j+size_even_double(2))-1)*size_double(1))

                enddo

            enddo

        endif

    end subroutine 

    subroutine gather_odd_block(size_double,size_even_double,size_odd_double,size_even_row,size_odd_row,input_column,output_column)

        implicit none

        integer, intent(in) :: size_double(2), size_even_double(2), size_odd_double(2)
        integer, intent(in) :: size_even_row, size_odd_row
        complex(kind(0d0)), intent(in)  :: input_column(size_even_row+size_odd_row)
        complex(kind(0d0)), intent(out) :: output_column(size_odd_row)

        integer i, j

        do j = 1, size_even_double(2)

            do i = 1, size_odd_double(1)

                output_column(i+(j-1)*size_odd_double(1)) = input_column((i+size_even_double(1))+(j-1)*size_double(1))

            enddo

        enddo

        do j = 1, size_odd_double(2)

            do i = 1, size_even_double(1)

                output_column(size_odd_double(1)*size_even_double(2)+i+(j-1)*size_even_double(1)) = input_column(i+((j+size_even_double(2))-1)*size_double(1))

            enddo

        enddo

    end subroutine 

    subroutine unblocking_singular(size_even_row,size_odd_row,size_even_column,size_odd_column,s_even,s_odd,s)

        implicit none

        integer, intent(in) :: size_even_row, size_odd_row, size_even_column, size_odd_column
        double precision, intent(in)  :: s_even(min(size_even_row,size_even_column))
        double precision, intent(in)  :: s_odd(min(size_odd_row,size_odd_column))
        double precision, intent(out) :: s(min(size_even_row+size_odd_row,size_even_column+size_odd_column))

        integer i

        do i = 1, min(size_even_row,size_even_column)

            s(i) = s_even(i)

        enddo

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            do i = 1, min(size_odd_row,size_odd_column)

                s(i+min(size_even_row,size_even_column)) = s_odd(i)

            enddo

        endif

    end subroutine 

    subroutine unblocking_unitary(size_double,size_even_double,size_even_row,size_odd_row,size_even_column,size_odd_column,u_even,u_odd,u)

        implicit none

        integer, intent(in) :: size_double(2)
        integer, intent(in) :: size_even_double(2)
        integer, intent(in) :: size_even_row, size_odd_row, size_even_column, size_odd_column
        complex(kind(0d0)), intent(in)  :: u_even(size_even_row*min(size_even_row,size_even_column))
        complex(kind(0d0)), intent(in)  :: u_odd(size_odd_row*min(size_odd_row,size_odd_column))
        complex(kind(0d0)), intent(out) :: u((size_even_row+size_odd_row)*min(size_even_row+size_odd_row,size_even_column+size_odd_column))

        integer i, j, k
        integer size_odd_double(2)
        
        complex(kind(0d0)), allocatable :: column_vect(:), even_vect(:), odd_vect(:)


        do i = 1, 2

            size_odd_double(i) = size_double(i) - size_even_double(i)

        enddo

        allocate(column_vect(size_even_row+size_odd_row))
        allocate(even_vect(size_even_row))

        column_vect = (0d0,0d0)
        even_vect = (0d0,0d0)

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            allocate(odd_vect(size_odd_row))

            odd_vect = (0d0,0d0)

        endif


        do k = 1, min(size_even_row,size_even_column)

            do i = 1, size_even_row

                even_vect(i) = u_even(i+(k-1)*size_even_row)

            enddo

            do j = 1, size_even_double(2)

                do i = 1, size_even_double(1)
    
                    column_vect(i+(j-1)*size_double(1)) = even_vect(i+(j-1)*size_even_double(1)) 
    
                enddo
    
            enddo

            if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

                do j = 1, size_odd_double(2)

                    do i = 1, size_odd_double(1)

                        column_vect((i+size_even_double(1))+((j+size_even_double(2))-1)*size_double(1))= even_vect(size_even_double(1)*size_even_double(2)+i+(j-1)*size_odd_double(1)) 

                    enddo

                enddo

            endif

            do i = 1, size_double(1)*size_double(2)

                u(i+(k-1)*size_double(1)*size_double(2)) = column_vect(i)

            enddo

            even_vect = (0d0,0d0)
            column_vect = (0d0,0d0)

        enddo

        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            do k = 1, min(size_odd_row,size_odd_column)

                do i = 1, size_odd_row
    
                    odd_vect(i) = u_odd(i+(k-1)*size_odd_row)
    
                enddo

                do j = 1, size_even_double(2)

                    do i = 1, size_odd_double(1)
        
                        column_vect((i+size_even_double(1))+(j-1)*size_double(1)) = odd_vect(i+(j-1)*size_odd_double(1)) 
        
                    enddo
        
                enddo
        
                do j = 1, size_odd_double(2)
        
                    do i = 1, size_even_double(1)
        
                        column_vect(i+((j+size_even_double(2))-1)*size_double(1)) = odd_vect(size_odd_double(1)*size_even_double(2)+i+(j-1)*size_even_double(1))
        
                    enddo
        
                enddo
    
                do i = 1, size_double(1)*size_double(2)
    
                    u(i+((k+min(size_even_row,size_even_column))-1)*size_double(1)*size_double(2)) = column_vect(i)
    
                enddo

                odd_vect = (0d0,0d0)
                column_vect = (0d0,0d0)
    
            enddo

        endif



        deallocate(column_vect,even_vect)


        if ( size_odd_row /= 0 .and. size_odd_column /= 0 ) then

            deallocate(odd_vect)

        endif

    end subroutine 

end module module_blocking_operations