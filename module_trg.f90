module module_trg

    contains

    subroutine btrg_loop()

        use module_declaration
        use module_setup
        use module_record
        
        implicit none

        do while ( iter < maxiter )

            iter = iter + 1

            call set_new_index()

            call get_new_tensor()

            call measurement()

            call update_index()

            call print_result()

        enddo

    end subroutine

    subroutine get_new_tensor()

        use module_declaration

        implicit none
        
        double precision, allocatable :: s_0(:), s_1(:)
        complex(kind(0d0)), allocatable :: u_0(:), v_0(:), u_1(:), v_1(:)

        allocate(s_0(min(bond(1)*bond(2),D1)))
        allocate(u_0(bond(1)*bond(2)*min(bond(1)*bond(2),D1)))
        allocate(v_0(bond(1)*bond(2)*min(bond(1)*bond(2),D1)))
        allocate(s_1(size(s_0)))
        allocate(u_1(size(u_0)))
        allocate(v_1(size(v_0)))

        s_0 = 0d0
        u_0 = (0d0,0d0)
        v_0 = (0d0,0d0)
        s_1 = 0d0
        u_1 = (0d0,0d0)
        v_1 = (0d0,0d0)


        call svd_even_site(s_0,u_0,v_0)
        call svd_odd_site(s_1,u_1,v_1)

        call contraction(s_0,u_0,v_0,s_1,u_1,v_1)


        deallocate(s_0,u_0,v_0)
        deallocate(s_1,u_1,v_1)

    end subroutine

    subroutine svd_even_site(s_trunc,u_trunc,v_trunc)

        use module_declaration
        use module_blocking_operations
        use module_record

        implicit none

        integer i
        integer size_tensor(4), size_tensor_even(4), size_even_result, size_even_trunc
        double precision, intent(out) :: s_trunc(min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(out) :: u_trunc(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(out) :: v_trunc(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), allocatable :: temporary(:)
        double precision, allocatable :: s(:)
        complex(kind(0d0)), allocatable :: u(:), v(:)

        do i = 1, 2

            size_tensor(i) = bond(i)
            size_tensor(i+2) = size_tensor(i)

            size_tensor_even(i) = bond_even(i)
            size_tensor_even(i+2) = size_tensor_even(i)

        enddo

        allocate(temporary(bond(1)*bond(2)*bond(1)*bond(2)))

        allocate(s(bond(1)*bond(2)))
        allocate(u(size(temporary)))
        allocate(v(size(temporary)))
        
        temporary = (0d0,0d0)

        s = 0d0
        u = (0d0,0d0)
        v = (0d0,0d0)

        do i = 1, size(temporary)

            temporary(i) = tensor(i)

        enddo

        call blocked_svd(size_tensor,size_tensor_even,size_even_result,temporary,s,u,v)

        call truncation(bond(1)*bond(2),bond(1)*bond(2),bond(1)*bond(2),size_even_result,s,u,v,D1,size_even_trunc,s_trunc,u_trunc,v_trunc)

        call record_svd_result(s,size(s),s_trunc,size(s_trunc),size_even_result,'esite',iter)

        !!! Set new even index

        bond_even_new(2) = size_even_trunc

        

        deallocate(temporary,s,u,v)

    end subroutine

    subroutine svd_odd_site(s_trunc,u_trunc,v_trunc)

        use module_declaration
        use module_phase
        use module_blocking_operations
        use module_record

        implicit none

        integer i, j, k, l
        integer size_tensor(4), size_tensor_even(4), size_even_result, size_even_trunc
        integer, allocatable :: phase_1(:), phase_2(:)
        double precision, intent(out) :: s_trunc(min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(out) :: u_trunc(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(out) :: v_trunc(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), allocatable :: temporary(:)
        double precision, allocatable :: s(:)
        complex(kind(0d0)), allocatable :: u(:), v(:)

        do i = 1, 2

            size_tensor(i) = bond(i)
            size_tensor(i+2) = size_tensor(i)

            size_tensor_even(i) = bond_even(i)
            size_tensor_even(i+2) = size_tensor_even(i)

        enddo

        allocate(phase_1(bond(1)))
        allocate(phase_2(bond(2)))

        allocate(temporary(bond(1)*bond(2)*bond(1)*bond(2)))

        allocate(s(bond(1)*bond(2)))
        allocate(u(size(temporary)))
        allocate(v(size(temporary)))
        
        temporary = (0d0,0d0)

        s = 0d0
        u = (0d0,0d0)
        v = (0d0,0d0)

        call set_basic_phase(phase_1,phase_2)

        do l = 1, bond(2)

            do k = 1, bond(1)

                do j = 1, bond(2)

                    do i = 1, bond(1)

                        temporary(i+(l-1)*bond(1)+(k-1)*bond(1)*bond(2)+(j-1)*bond(1)*bond(2)*bond(1)) &
                        = tensor(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)+(l-1)*bond(1)*bond(2)*bond(1)) &
                        * (-1d0)**(phase_2(l)*(phase_2(j)+phase_1(k))) &
                        * (-1d0)**(phase_1(k)*phase_2(j))

                    enddo

                enddo

            enddo

        enddo

        call blocked_svd(size_tensor,size_tensor_even,size_even_result,temporary,s,u,v)

        call truncation(bond(1)*bond(2),bond(1)*bond(2),bond(1)*bond(2),size_even_result,s,u,v,D1,size_even_trunc,s_trunc,u_trunc,v_trunc)

        call record_svd_result(s,size(s),s_trunc,size(s_trunc),size_even_result,'osite',iter)

        !!! Set new even index

        bond_even_new(1) = size_even_trunc

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

    subroutine contraction(s_0,u_0,v_0,s_1,u_1,v_1)

        use module_declaration
        use module_phase
        use module_lapack

        implicit none

        double precision, intent(in) :: s_0(min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(in) :: u_0(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(in) :: v_0(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        double precision, intent(in) :: s_1(min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(in) :: u_1(bond(1)*bond(2)*min(bond(1)*bond(2),D1))
        complex(kind(0d0)), intent(in) :: v_1(bond(1)*bond(2)*min(bond(1)*bond(2),D1))


        integer i, j, k, l
        integer, allocatable :: phase_new_1(:), phase_new_2(:)
        integer, allocatable :: phase_1(:), phase_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:), temporary_3(:), temporary_4(:)
        complex(kind(0d0)), allocatable :: temporary_5(:), temporary_6(:), temporary_7(:), temporary_8(:)
        complex(kind(0d0)), allocatable :: temporary_9(:)

        allocate(phase_new_1(bond_new(1)))
        allocate(phase_new_2(bond_new(2)))
        allocate(phase_1(bond(1)))
        allocate(phase_2(bond(2)))

        allocate(temporary_1(bond(1)*bond(2)*bond_new(2)))
        allocate(temporary_2(bond(1)*bond(2)*bond_new(2)))
        allocate(temporary_3(bond(1)*bond(2)*bond_new(1)))
        allocate(temporary_4(bond(1)*bond(2)*bond_new(1)))
        allocate(temporary_5(bond(2)*bond_new(1)*bond(2)*bond_new(2)))
        allocate(temporary_6(bond(2)*bond_new(2)*bond(2)*bond_new(1)))
        allocate(temporary_7(bond(2)*bond_new(1)*bond(2)*bond_new(2)))
        allocate(temporary_8(bond(2)*bond_new(2)*bond(2)*bond_new(1)))
        allocate(temporary_9(bond_new(2)*bond_new(1)*bond_new(1)*bond_new(2)))

        temporary_1 = (0d0,0d0)
        temporary_2 = (0d0,0d0)
        temporary_3 = (0d0,0d0)
        temporary_4 = (0d0,0d0)
        temporary_5 = (0d0,0d0)
        temporary_6 = (0d0,0d0)
        temporary_7 = (0d0,0d0)
        temporary_8 = (0d0,0d0)
        temporary_9 = (0d0,0d0)



        call set_basic_phase_new(phase_new_1,phase_new_2)
        call set_basic_phase(phase_1,phase_2)

        

        do k = 1, bond_new(2)

            do j = 1, bond(2)

                do i = 1, bond(1)

                    temporary_1(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = u_0(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) * sqrt(s_0(k))**(1d0-hyper) 

                    temporary_2(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = v_0(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) * sqrt(s_0(k))**(1d0-hyper) 

                enddo

            enddo

        enddo

        do k = 1, bond_new(1)

            do j = 1, bond(2)

                do i = 1, bond(1)

                    temporary_3(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = u_1(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) * sqrt(s_1(k))**(1d0-hyper) 

                    temporary_4(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = v_1(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) * sqrt(s_1(k))**(1d0-hyper) 

                enddo

            enddo

        enddo

        call contract_bond_weights(temporary_1)
        call contract_bond_weights(temporary_2)

        !!! Pretreatment to restore the Grassmann algebra

        do k = 1, bond_new(2)

            do j = 1, bond(2)

                do i = 1, bond(1)

                    temporary_1(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = temporary_1(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    * (-1d0)**(phase_1(i)*(phase_2(j)+phase_new_2(k)))

                    temporary_2(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = temporary_2(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    * (-1d0)**(phase_1(i)*phase_new_2(k))

                enddo

            enddo

        enddo

        do k = 1, bond_new(1)

            do j = 1, bond(2)

                do i = 1, bond(1)

                    temporary_3(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = temporary_3(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    * (-1d0)**(phase_1(i)*(phase_2(j)+phase_new_1(k)))

                    temporary_4(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    = temporary_4(i+(j-1)*bond(1)+(k-1)*bond(1)*bond(2)) &
                    * (-1d0)**(phase_1(i)*phase_new_1(k))

                enddo

            enddo

        enddo

        call general_zgemm('T','N',bond(2)*bond_new(1),bond(2)*bond_new(2),bond(1),temporary_3,temporary_2,temporary_5)
        call general_zgemm('T','N',bond(2)*bond_new(2),bond(2)*bond_new(1),bond(1),temporary_1,temporary_4,temporary_6)

        do l = 1, bond_new(2)

            do k = 1, bond(2)

                do j = 1, bond_new(1)

                    do i = 1, bond(2)

                        temporary_7(i+(k-1)*bond(2)+(j-1)*bond(2)*bond(2)+(l-1)*bond(2)*bond(2)*bond_new(1)) &
                        = temporary_5(i+(j-1)*bond(2)+(k-1)*bond(2)*bond_new(1)+(l-1)*bond(2)*bond_new(1)*bond(2)) &
                        * (-1d0)**(phase_2(k)*(phase_new_1(j)+phase_new_2(l)))

                    enddo

                enddo

            enddo

        enddo

        do l = 1, bond_new(1)

            do k = 1, bond(2)

                do j = 1, bond_new(2)

                    do i = 1, bond(2)

                        temporary_8(i+(k-1)*bond(2)+(j-1)*bond(2)*bond(2)+(l-1)*bond(2)*bond(2)*bond_new(2)) &
                        = temporary_6(i+(j-1)*bond(2)+(k-1)*bond(2)*bond_new(2)+(l-1)*bond(2)*bond_new(2)*bond(2)) &
                        * (-1d0)**(phase_2(i)*(phase_new_2(j)+phase_new_1(l)+phase_2(k)))

                    enddo

                enddo

            enddo

        enddo

        call general_zgemm('T','N',bond_new(2)*bond_new(1),bond_new(1)*bond_new(2),bond(2)*bond(2),temporary_8,temporary_7,temporary_9)

        tensor = 0d0

        do l = 1, bond_new(2)

            do k = 1, bond_new(1)

                do j = 1, bond_new(1)

                    do i = 1, bond_new(2)

                        tensor(k+(i-1)*bond_new(1)+(j-1)*bond_new(1)*bond_new(2)+(l-1)*bond_new(1)*bond_new(2)*bond_new(1)) &
                        = temporary_9(i+(j-1)*bond_new(2)+(k-1)*bond_new(2)*bond_new(1)+(l-1)*bond_new(2)*bond_new(1)*bond_new(1)) &
                        * (-1d0)**(phase_new_1(k)*(phase_new_2(i)+phase_new_1(j)))

                    enddo

                enddo

            enddo

        enddo

        call update_bond_weights(s_0,s_1)

        deallocate(phase_new_1,phase_new_2)
        deallocate(phase_1,phase_2)

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

    subroutine contract_bond_weights(array)

        use module_declaration
        use module_lapack

        implicit none

        integer i, j
        complex(kind(0d0)), intent(inout) :: array(bond(1)*bond(2)*min(bond(1)*bond(2),D1))


        complex(kind(0d0)), allocatable :: bw_1(:), bw_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:)

        allocate(bw_1(bond(1)*bond(1)))
        allocate(bw_2(bond(2)*bond(2)))

        allocate(temporary_1(size(array)))
        allocate(temporary_2(size(array)))

        bw_1 = (0d0,0d0)
        bw_2 = (0d0,0d0)

        temporary_1 = (0d0,0d0)
        temporary_2 = (0d0,0d0)

        do i = 1, size(bw_1)

            bw_1(i) = dcmplx(bond_weight_1(i))

        enddo

        do i = 1, size(bw_2)

            bw_2(i) = dcmplx(bond_weight_2(i))

        enddo

        call general_zgemm('T','N',bond(2)*min(bond(1)*bond(2),D1),bond(1),bond(1),array,bw_1,temporary_1)
        call general_zgemm('T','N',min(bond(1)*bond(2),D1)*bond(1),bond(2),bond(2),temporary_1,bw_2,temporary_2)

        array = (0d0,0d0)

        do j = 1, bond(1)*bond(2)

            do i = 1, min(bond(1)*bond(2),D1)

                array(j+(i-1)*bond(1)*bond(2)) &
                = temporary_2(i+(j-1)*min(bond(1)*bond(2),D1))

            enddo

        enddo


        deallocate(bw_1,bw_2)
        deallocate(temporary_1)
        deallocate(temporary_2)


    end subroutine

    subroutine update_bond_weights(s_0,s_1)

        use module_declaration
        use module_record

        implicit none

        integer i, j
        double precision, intent(in) :: s_0(min(bond(1)*bond(2),D1))
        double precision, intent(in) :: s_1(min(bond(1)*bond(2),D1))

        !!! Set new bond_weight_2

        bond_weight_2 = 0d0

        do j = 1, min(bond(1)*bond(2),D1)

            do i = 1, min(bond(1)*bond(2),D1)

                if ( i == j ) then

                    if ( s_0(i) /= 0d0 ) then

                        bond_weight_2(i+(j-1)*min(bond(1)*bond(2),D1)) = s_0(i)**hyper

                    endif

                endif

            enddo

        enddo


        !!! Set new bond_weight_1

        bond_weight_1 = 0d0

        do j = 1, min(bond(1)*bond(2),D1)

            do i = 1, min(bond(1)*bond(2),D1)

                if ( i == j ) then

                    if ( s_1(i) /= 0d0 ) then

                        bond_weight_1(i+(j-1)*min(bond(1)*bond(2),D1)) = s_1(i)**hyper

                    endif

                endif

            enddo

        enddo

    end subroutine

end module module_trg