module module_trg

    contains

    subroutine svd_even_site(size_input,size_input_even,input,n_iter,bond_dim,s_trunc,u_trunc,v_trunc,size_even_trunc)

        use module_blocking_operations
        use module_record

        implicit none
        
        integer, intent(in) :: size_input(2), size_input_even(2), n_iter, bond_dim
        integer, intent(out) :: size_even_trunc
        integer i
        integer size_tensor(4), size_tensor_even(4)
        integer row, col, size_even_result

        double precision, intent(out) :: s_trunc(min(size_input(1)*size_input(2),bond_dim))
        double precision, allocatable :: s(:)

        complex(kind(0d0)), intent(in) :: input(bond_dim**4)
        complex(kind(0d0)), intent(out) :: u_trunc(size_input(1)*size_input(2)*min(size_input(1)*size_input(2),bond_dim))
        complex(kind(0d0)), intent(out) :: v_trunc(size_input(1)*size_input(2)*min(size_input(1)*size_input(2),bond_dim))
        complex(kind(0d0)), allocatable :: temporary(:)
        complex(kind(0d0)), allocatable :: u(:), v(:)

        row = size_input(1)*size_input(2)
        col = row

        do i = 1, 2

            size_tensor(i) = size_input(i)
            size_tensor(i+2) = size_tensor(i)

            size_tensor_even(i) = size_input_even(i)
            size_tensor_even(i+2) = size_tensor_even(i)

        enddo

        allocate(temporary(row*col))

        allocate(s(min(row,col)))
        allocate(u(size(temporary)))
        allocate(v(size(temporary)))
        
        temporary = (0d0,0d0)

        s = 0d0
        u = (0d0,0d0)
        v = (0d0,0d0)

        do i = 1, size(temporary)

            temporary(i) = input(i)

        enddo

        call blocked_svd(size_tensor,size_tensor_even,size_even_result,temporary,s,u,v)

        call truncation(row,col,min(row,col),size_even_result,s,u,v,bond_dim,size_even_trunc,s_trunc,u_trunc,v_trunc)

        call record_svd_result(s,size(s),s_trunc,size(s_trunc),size_even_result,'esite',n_iter)

        !!! Set new even index

        !bond_even_new(2) = size_even_trunc

        

        deallocate(temporary,s,u,v)

    end subroutine

    subroutine svd_odd_site(size_input,size_input_even,input,n_iter,bond_dim,s_trunc,u_trunc,v_trunc,size_even_trunc)

        use module_phase
        use module_blocking_operations
        use module_record

        implicit none
        
        integer, intent(in) :: size_input(2), size_input_even(2), n_iter, bond_dim
        integer, intent(out) :: size_even_trunc
        integer, allocatable :: phase_1(:), phase_2(:)
        integer i, j, k, l
        integer size_tensor(4), size_tensor_even(4)
        integer row, col, size_even_result

        double precision, intent(out) :: s_trunc(min(size_input(1)*size_input(2),bond_dim))
        double precision, allocatable :: s(:)

        complex(kind(0d0)), intent(in) :: input(bond_dim**4)
        complex(kind(0d0)), intent(out) :: u_trunc(size_input(1)*size_input(2)*min(size_input(1)*size_input(2),bond_dim))
        complex(kind(0d0)), intent(out) :: v_trunc(size_input(1)*size_input(2)*min(size_input(1)*size_input(2),bond_dim))
        complex(kind(0d0)), allocatable :: temporary(:)
        complex(kind(0d0)), allocatable :: u(:), v(:)

        row = size_input(1)*size_input(2)
        col = row

        do i = 1, 2

            size_tensor(i) = size_input(i)
            size_tensor(i+2) = size_tensor(i)

            size_tensor_even(i) = size_input_even(i)
            size_tensor_even(i+2) = size_tensor_even(i)

        enddo

        allocate(phase_1(size_input(1)))
        allocate(phase_2(size_input(2)))

        allocate(temporary(row*col))

        allocate(s(min(row,col)))
        allocate(u(size(temporary)))
        allocate(v(size(temporary)))
        
        temporary = (0d0,0d0)

        s = 0d0
        u = (0d0,0d0)
        v = (0d0,0d0)

        !call set_basic_phase(phase_1,phase_2)
        call set_phase(size_input(1),size_input_even(1),phase_1)
        call set_phase(size_input(2),size_input_even(2),phase_2)

        do l = 1, size_tensor(2)

            do k = 1, size_tensor(1)

                do j = 1, size_tensor(2)

                    do i = 1, size_tensor(1)

                        temporary(i+(l-1)*size_tensor(1)+(k-1)*row+(j-1)*row*size_tensor(1)) &
                        = input(i+(j-1)*size_tensor(1)+(k-1)*row+(l-1)*row*size_tensor(1)) &
                        * (-1d0)**(phase_2(l)*(phase_2(j)+phase_1(k))) &
                        * (-1d0)**(phase_1(k)*phase_2(j))

                    enddo

                enddo

            enddo

        enddo

        call blocked_svd(size_tensor,size_tensor_even,size_even_result,temporary,s,u,v)

        call truncation(row,col,min(row,col),size_even_result,s,u,v,bond_dim,size_even_trunc,s_trunc,u_trunc,v_trunc)

        call record_svd_result(s,size(s),s_trunc,size(s_trunc),size_even_result,'osite',n_iter)

        !!! Set new even index

        !bond_even_new(1) = size_even_trunc

        deallocate(phase_1,phase_2)
        deallocate(temporary,s,u,v)

    end subroutine

    subroutine truncation(m,n,k,size_even_k,s,u,v,trunc,size_even_trunc,s_trunc,u_trunc,v_trunc)

        implicit none

        integer, intent(in) :: m, n, k
        integer, intent(in) :: size_even_k
        integer, intent(in) :: trunc
        double precision, intent(in) :: s(k)
        complex(kind(0d0)), intent(in) :: u(m*k), v(n*k)

        integer, intent(out) :: size_even_trunc
        double precision, intent(out) :: s_trunc(min(k,trunc))
        complex(kind(0d0)), intent(out) :: u_trunc(m*min(k,trunc)), v_trunc(n*min(k,trunc))

        integer i, j
        integer flag, size_odd_k, size_odd_trunc
        double precision s_even_cand, s_odd_cand
        double precision new_s_even_cand, new_s_odd_cand

        size_odd_k = k - size_even_k

        size_even_trunc = 0
        size_odd_trunc = 0

        if ( k >= trunc ) then

            if ( size_odd_k /= 0 ) then

                flag = 1
    
                s_even_cand = s(1)
                s_odd_cand = s(size_even_k+1)
    
                do while ( flag <= trunc )
    
                    if ( size_even_trunc < size_even_k .and. size_odd_trunc < size_odd_k ) then
    
                        if ( s_even_cand >= s_odd_cand ) then
    
                            size_even_trunc = size_even_trunc + 1

                            if ( size_even_trunc < size_even_k ) then
    
                                new_s_even_cand = s(size_even_trunc+1)

                            else

                                new_s_even_cand = 0d0

                            endif

                            new_s_odd_cand = s_odd_cand
    
                        elseif ( s_even_cand < s_odd_cand ) then
    
                            size_odd_trunc = size_odd_trunc + 1

                            if ( size_odd_trunc < size_odd_k ) then
    
                                new_s_odd_cand = s(size_even_k+size_odd_trunc+1)

                            else

                                new_s_odd_cand = 0d0

                            endif

                            new_s_even_cand = s_even_cand

                        endif
    
                    elseif ( size_even_trunc == size_even_k ) then
    
                        size_odd_trunc = size_odd_trunc + 1
    
                    elseif ( size_odd_trunc == size_odd_k ) then
    
                        size_even_trunc = size_even_trunc + 1
    
                    endif
    
                    flag = flag + 1

                    s_even_cand = new_s_even_cand
                    s_odd_cand = new_s_odd_cand
    
                enddo
    
                do i = 1, size_even_trunc
    
                    s_trunc(i) = s(i)
    
                enddo
    
                do i = 1, size_odd_trunc
    
                    s_trunc(size_even_trunc+i) = s(size_even_k+i)
    
                enddo
    
                do j = 1, size_even_trunc
    
                    do i = 1, m
    
                        u_trunc(i+(j-1)*m) = u(i+(j-1)*m)
    
                    enddo
    
                enddo
    
                do j = 1, size_odd_trunc
    
                    do i = 1, m
    
                        u_trunc(i+(size_even_trunc+j-1)*m) = u(i+(size_even_k+j-1)*m)
    
                    enddo
    
                enddo
    
                do j = 1, size_even_trunc
    
                    do i = 1, n
    
                        v_trunc(i+(j-1)*n) = v(i+(j-1)*n)
    
                    enddo
    
                enddo
    
                do j = 1, size_odd_trunc
    
                    do i = 1, n
    
                        v_trunc(i+(size_even_trunc+j-1)*n) = v(i+(size_even_k+j-1)*n)
    
                    enddo
    
                enddo
    
            elseif ( size_odd_k == 0 ) then
    
                size_even_trunc = trunc
    
                size_odd_trunc = 0
    
                do i = 1, size_even_trunc
    
                    s_trunc(i) = s(i)
    
                enddo
    
                do j = 1, size_even_trunc
    
                    do i = 1, m
    
                        u_trunc(i+(j-1)*m) = u(i+(j-1)*m)
    
                    enddo
    
                enddo
    
                do j = 1, size_even_trunc
    
                    do i = 1, n
    
                        v_trunc(i+(j-1)*n) = v(i+(j-1)*n)
    
                    enddo
    
                enddo
    
            endif

        elseif ( k < trunc ) then

            size_even_trunc = size_even_k !+ int((trunc-k)*0.5d0+0.5d0)

            size_odd_trunc = size_odd_k !+ int((trunc-k)*0.5d0)

            do i = 1, k

                s_trunc(i) = s(i)

            enddo

            do j = 1, k

                do i = 1, m

                    u_trunc(i+(j-1)*m) = u(i+(j-1)*m)

                enddo

            enddo

            do j = 1, k

                do i = 1, n

                    v_trunc(i+(j-1)*n) = v(i+(j-1)*n)

                enddo

            enddo

        endif

    end subroutine 

    subroutine contraction(size_old,size_old_even,size_new,size_new_even,s_0,u_0,v_0,s_1,u_1,v_1,input_bw_1,input_bw_2,bond_dim,h_param,output)

        use module_phase
        use module_lapack

        implicit none

        integer, intent(in) :: size_old(2), size_old_even(2)
        integer, intent(in) :: size_new(2), size_new_even(2)
        integer, intent(in) :: bond_dim

        double precision, intent(in) :: h_param, input_bw_1(bond_dim**2), input_bw_2(bond_dim**2)
        double precision, intent(in) :: s_0(size_new(2))
        complex(kind(0d0)), intent(in) :: u_0(size_old(1)*size_old(2)*size_new(2))
        complex(kind(0d0)), intent(in) :: v_0(size_old(1)*size_old(2)*size_new(2))
        double precision, intent(in) :: s_1(size_new(1))
        complex(kind(0d0)), intent(in) :: u_1(size_old(1)*size_old(2)*size_new(1))
        complex(kind(0d0)), intent(in) :: v_1(size_old(1)*size_old(2)*size_new(1))

        complex(kind(0d0)), intent(out) :: output(bond_dim**4)

        integer i, j, k, l
        integer, allocatable :: phase_new_1(:), phase_new_2(:)
        integer, allocatable :: phase_old_1(:), phase_old_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:), temporary_3(:), temporary_4(:)
        complex(kind(0d0)), allocatable :: temporary_5(:), temporary_6(:), temporary_7(:), temporary_8(:)
        complex(kind(0d0)), allocatable :: temporary_9(:)

        allocate(phase_new_1(size_new(1)))
        allocate(phase_new_2(size_new(2)))
        allocate(phase_old_1(size_old(1)))
        allocate(phase_old_2(size_old(2)))

        allocate(temporary_1(size(u_0)))
        allocate(temporary_2(size(u_0)))
        allocate(temporary_3(size(u_1)))
        allocate(temporary_4(size(u_1)))
        allocate(temporary_5(size_old(2)*size_new(1)*size_old(2)*size_new(2)))
        allocate(temporary_6(size_old(2)*size_new(2)*size_old(2)*size_new(1)))
        allocate(temporary_7(size_old(2)*size_new(1)*size_old(2)*size_new(2)))
        allocate(temporary_8(size_old(2)*size_new(2)*size_old(2)*size_new(1)))
        allocate(temporary_9(size_new(2)*size_new(1)*size_new(1)*size_new(2)))

        temporary_1 = (0d0,0d0)
        temporary_2 = (0d0,0d0)
        temporary_3 = (0d0,0d0)
        temporary_4 = (0d0,0d0)
        temporary_5 = (0d0,0d0)
        temporary_6 = (0d0,0d0)
        temporary_7 = (0d0,0d0)
        temporary_8 = (0d0,0d0)
        temporary_9 = (0d0,0d0)



        call set_phase(size_new(1),size_new_even(1),phase_new_1)
        call set_phase(size_new(2),size_new_even(2),phase_new_2)
        call set_phase(size_old(1),size_old_even(1),phase_old_1)
        call set_phase(size_old(2),size_old_even(2),phase_old_2)

        

        do k = 1, size_new(2)

            do j = 1, size_old(2)

                do i = 1, size_old(1)

                    temporary_1(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = u_0(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) * sqrt(s_0(k))**(1d0-h_param) 

                    temporary_2(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = v_0(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) * sqrt(s_0(k))**(1d0-h_param) 

                enddo

            enddo

        enddo

        do k = 1, size_new(1)

            do j = 1, size_old(2)

                do i = 1, size_old(1)

                    temporary_3(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = u_1(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) * sqrt(s_1(k))**(1d0-h_param) 

                    temporary_4(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = v_1(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) * sqrt(s_1(k))**(1d0-h_param) 

                enddo

            enddo

        enddo

        call contract_bond_weights(size_old,size_new(2),input_bw_1,input_bw_2,bond_dim,temporary_1)
        call contract_bond_weights(size_old,size_new(2),input_bw_1,input_bw_2,bond_dim,temporary_2)

        !!! Pretreatment to restore the Grassmann algebra

        do k = 1, size_new(2)

            do j = 1, size_old(2)

                do i = 1, size_old(1)

                    temporary_1(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = temporary_1(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    * (-1d0)**(phase_old_1(i)*(phase_old_2(j)+phase_new_2(k)))

                    temporary_2(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = temporary_2(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    * (-1d0)**(phase_old_1(i)*phase_new_2(k))

                enddo

            enddo

        enddo

        do k = 1, size_new(1)

            do j = 1, size_old(2)

                do i = 1, size_old(1)

                    temporary_3(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = temporary_3(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    * (-1d0)**(phase_old_1(i)*(phase_old_2(j)+phase_new_1(k)))

                    temporary_4(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    = temporary_4(i+(j-1)*size_old(1)+(k-1)*size_old(1)*size_old(2)) &
                    * (-1d0)**(phase_old_1(i)*phase_new_1(k))

                enddo

            enddo

        enddo

        call general_zgemm('T','N',size_old(2)*size_new(1),size_old(2)*size_new(2),size_old(1),temporary_3,temporary_2,temporary_5)
        call general_zgemm('T','N',size_old(2)*size_new(2),size_old(2)*size_new(1),size_old(1),temporary_1,temporary_4,temporary_6)

        do l = 1, size_new(2)

            do k = 1, size_old(2)

                do j = 1, size_new(1)

                    do i = 1, size_old(2)

                        temporary_7(i+(k-1)*size_old(2)+(j-1)*size_old(2)*size_old(2)+(l-1)*size_old(2)*size_old(2)*size_new(1)) &
                        = temporary_5(i+(j-1)*size_old(2)+(k-1)*size_old(2)*size_new(1)+(l-1)*size_old(2)*size_new(1)*size_old(2)) &
                        * (-1d0)**(phase_old_2(k)*(phase_new_1(j)+phase_new_2(l)))

                    enddo

                enddo

            enddo

        enddo

        do l = 1, size_new(1)

            do k = 1, size_old(2)

                do j = 1, size_new(2)

                    do i = 1, size_old(2)

                        temporary_8(i+(k-1)*size_old(2)+(j-1)*size_old(2)*size_old(2)+(l-1)*size_old(2)*size_old(2)*size_new(2)) &
                        = temporary_6(i+(j-1)*size_old(2)+(k-1)*size_old(2)*size_new(2)+(l-1)*size_old(2)*size_new(2)*size_old(2)) &
                        * (-1d0)**(phase_old_2(i)*(phase_new_2(j)+phase_new_1(l)+phase_old_2(k)))

                    enddo

                enddo

            enddo

        enddo

        call general_zgemm('T','N',size_new(2)*size_new(1),size_new(1)*size_new(2),size_old(2)*size_old(2),temporary_8,temporary_7,temporary_9)

        output = (0d0,0d0)

        do l = 1, size_new(2)

            do k = 1, size_new(1)

                do j = 1, size_new(1)

                    do i = 1, size_new(2)

                        output(k+(i-1)*size_new(1)+(j-1)*size_new(1)*size_new(2)+(l-1)*size_new(1)*size_new(2)*size_new(1)) &
                        = temporary_9(i+(j-1)*size_new(2)+(k-1)*size_new(2)*size_new(1)+(l-1)*size_new(2)*size_new(1)*size_new(1)) &
                        * (-1d0)**(phase_new_1(k)*(phase_new_2(i)+phase_new_1(j)))

                    enddo

                enddo

            enddo

        enddo

        deallocate(phase_new_1,phase_new_2)
        deallocate(phase_old_1,phase_old_2)

        deallocate(temporary_1)
        deallocate(temporary_2)
        deallocate(temporary_3)
        deallocate(temporary_4)
        deallocate(temporary_5)
        deallocate(temporary_6)
        deallocate(temporary_7)
        deallocate(temporary_8)
        deallocate(temporary_9)


    end subroutine

    subroutine contract_bond_weights(size_old,size_new,input_bw_1,input_bw_2,bond_dim,array)

        use module_lapack

        implicit none

        integer i, j
        integer, intent(in) :: size_old(2), size_new, bond_dim
        double precision, intent(in) :: input_bw_1(bond_dim**2), input_bw_2(bond_dim**2)
        complex(kind(0d0)), intent(inout) :: array(size_old(1)*size_old(2)*size_new)


        complex(kind(0d0)), allocatable :: bw_1(:), bw_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:)

        allocate(bw_1(size_old(1)*size_old(1)))
        allocate(bw_2(size_old(2)*size_old(2)))

        allocate(temporary_1(size(array)))
        allocate(temporary_2(size(array)))

        bw_1 = (0d0,0d0)
        bw_2 = (0d0,0d0)

        temporary_1 = (0d0,0d0)
        temporary_2 = (0d0,0d0)

        do i = 1, size(bw_1)

            bw_1(i) = dcmplx(input_bw_1(i))

        enddo

        do i = 1, size(bw_2)

            bw_2(i) = dcmplx(input_bw_2(i))

        enddo

        call general_zgemm('T','N',size_old(2)*size_new,size_old(1),size_old(1),array,bw_1,temporary_1)
        call general_zgemm('T','N',size_new*size_old(1),size_old(2),size_old(2),temporary_1,bw_2,temporary_2)

        array = (0d0,0d0)

        do j = 1, size_old(1)*size_old(2)

            do i = 1, size_new

                array(j+(i-1)*size_old(1)*size_old(2)) &
                = temporary_2(i+(j-1)*size_new)

            enddo

        enddo


        deallocate(bw_1,bw_2)
        deallocate(temporary_1)
        deallocate(temporary_2)


    end subroutine

    subroutine update_bond_weights(size_s,bond_dim,h_param,s,output)

        use module_record

        implicit none

        integer i, j
        integer, intent(in) :: size_s, bond_dim
        double precision, intent(in) :: h_param, s(size_s)
        double precision, intent(out) :: output(bond_dim**2)

        output = 0d0

        do j = 1, size_s

            do i = 1, size_s

                if ( i == j ) then

                    if ( s(i) /= 0d0 ) then

                        output(i+(j-1)*size_s) = s(i)**h_param

                    endif

                endif

            enddo

        enddo

    end subroutine

end module module_trg