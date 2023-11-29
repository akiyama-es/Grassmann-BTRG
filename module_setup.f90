module module_setup

    contains

    subroutine load_input_trg()

        use module_declaration
        use module_initial_tensor

        implicit none

        write(*,*) ' Set the parameters in Gross-Neveu-Wilson model : '
        write(*,*) ' Mass ? '
        read(*,*) mass
        write(*,*) ' Chemical potential ? '
        read(*,*) mu
        write(*,*) ' Four-point coupling ? '
        read(*,*) coupling

        if ( coupling == 0d0 ) then

            write(*,*) ' The free Wilson fermion theory will be computed. '
            write(*,*) ' Do you want to calculate the relative error ? [y/n] '
            write(*,*) ' REMARK: '
            write(*,*) ' The computational time can be much longer to obtain the exact solution. '
            read(*,*) flag_exact

        endif

        write(*,*) ' Set the parameters in Grassmann BTRG : '
        write(*,*) ' Bond dimension D ? '
        write(*,*) ' REMARK: '
        write(*,*) ' The initial tensor network is written by D=4.  '
        write(*,*) ' D must be equal to or greater than 4. '
        read(*,*) D1

        if ( D1 < 4 ) then

            write(*,*) ' ERROR: '
            write(*,*) ' D must be equal to or greater than 4. '
            stop

        endif

        write(*,*) ' Iteration number n ? ( n iterations give us the lattice volume of 2^n ) '
        read(*,*) maxiter
        write(*,*) ' hyperparameter k in BTRG ? '
        write(*,*) ' Optimal BTRG >>> k = -0.5 ( type -0.5 ) '
        write(*,*) ' Levin-Nave TRG >>> k = 0.0 ( type 0.0 ) '
        read(*,*) hyper

        iter = 0

        call load_input_init_size(bond,bond_even,bond_new,bond_even_new)

    end subroutine

    subroutine load_input_init_tensor()

        use module_declaration
        use module_initial_tensor

        implicit none

        integer i
        complex(kind(0d0)), allocatable :: temporary(:)

        allocate(temporary(bond(1)*bond(2)*bond(1)*bond(2)))

        temporary = (0d0,0d0)

        call get_init_tensor(mu,mass,coupling,temporary)

        do i = 1, size(temporary)

            tensor(i) = temporary(i)

        enddo

        deallocate(temporary)
        
    end subroutine

    subroutine btrg_loop()

        use module_declaration
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
        use module_trg

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

        call svd_even_site(bond,bond_even,tensor,iter,D1,s_0,u_0,v_0,bond_even_new(2))
        call svd_odd_site(bond,bond_even,tensor,iter,D1,s_1,u_1,v_1,bond_even_new(1))

        call contraction(bond,bond_even,bond_new,bond_even_new,s_0,u_0,v_0,s_1,u_1,v_1,bond_weight_1,bond_weight_2,D1,hyper,tensor)

        call update_bond_weights(size(s_0),D1,hyper,s_0,bond_weight_2)
        call update_bond_weights(size(s_1),D1,hyper,s_1,bond_weight_1)

        deallocate(s_0,u_0,v_0)
        deallocate(s_1,u_1,v_1)

    end subroutine

    subroutine measurement()

        use module_declaration
        use module_measurement

        implicit none

        call get_trace(bond_new,bond_even_new,tensor,bond_weight_1,bond_weight_2,D1,boundary_condition,part(iter))

        call get_lnz(iter,maxiter,part,lnz(iter))

        if ( flag_exact == 'y' .and. mod(iter,2) == 0 ) then

            call get_relative_error(iter,mass,mu,lnz(iter),relative_error(iter))

        endif 

        call normalize_tensor()

    end subroutine

    subroutine normalize_tensor()

        use module_declaration

        implicit none

        f_norm_1(iter) = 0.5d0*(maxval(bond_weight_1)+maxval(bond_weight_2))
        f_norm_2(iter) = part(iter)/(f_norm_1(iter)*f_norm_1(iter))

        tensor = tensor/f_norm_2(iter)

        bond_weight_1 = bond_weight_1/f_norm_1(iter)
        bond_weight_2 = bond_weight_2/f_norm_1(iter)

    end subroutine

    subroutine set_new_index()

        use module_declaration

        implicit none

        integer i

        do i = 1, 2

            bond_new(i) = min(bond(1)*bond(2),D1)

        enddo

    end subroutine

    subroutine update_index()

        use module_declaration

        implicit none

        integer i

        do i = 1, 2

            bond(i) = bond_new(i)

            bond_even(i) = bond_even_new(i)

        enddo

    end subroutine

    subroutine allocate_memory()

        use module_declaration

        implicit none

        allocate(tensor(D1*D1*D1*D1))
        allocate(bond_weight_1(D1*D1))
        allocate(bond_weight_2(D1*D1))
        allocate(part(0:maxiter))
        allocate(lnz(0:maxiter))
        allocate(relative_error(0:maxiter))
        allocate(f_norm_1(0:maxiter))
        allocate(f_norm_2(0:maxiter))

        tensor = (0d0,0d0)
        bond_weight_1 = 0d0
        bond_weight_2 = 0d0
        part = (0d0,0d0)
        lnz = (0d0,0d0)
        relative_error = 0d0
        f_norm_1 = 1d0
        f_norm_2 = 1d0

    end subroutine

    subroutine deallocate_memory()

        use module_declaration

        implicit none

        deallocate(tensor,part,lnz)
        deallocate(bond_weight_1,bond_weight_2)
        deallocate(relative_error)
        deallocate(f_norm_1,f_norm_2)

    end subroutine

    subroutine initialize_bond_weights()

        use module_declaration

        implicit none

        integer i, j
        double precision, allocatable :: temporary_1(:), temporary_2(:)

        allocate(temporary_1(bond(1)*bond(1)))
        allocate(temporary_2(bond(2)*bond(2)))

        temporary_1 = 0d0
        temporary_2 = 0d0

        do j = 1, bond(1)

            do i = 1, bond(1)

                if ( i == j ) then

                    temporary_1(i+(j-1)*bond(1)) = 1d0

                endif

            enddo

        enddo

        do j = 1, bond(2)

            do i = 1, bond(2)

                if ( i == j ) then

                    temporary_2(i+(j-1)*bond(2)) = 1d0

                endif

            enddo

        enddo

        do i = 1, size(temporary_1)

            bond_weight_1(i) = temporary_1(i)

        enddo

        do i = 1, size(temporary_2)

            bond_weight_2(i) = temporary_2(i)

        enddo

        deallocate(temporary_1)
        deallocate(temporary_2)

    end subroutine

end module module_setup