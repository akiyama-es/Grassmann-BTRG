module module_measurement

    contains

    subroutine get_trace(size_input,size_input_even,input_tensor,input_bw_1,input_bw_2,bond_dim,b_cond,output)

        use module_phase
    
        implicit none

        integer, intent(in) :: size_input(2), size_input_even(2), bond_dim
        double precision, intent(in) :: b_cond
        double precision, intent(in) :: input_bw_1(bond_dim**2), input_bw_2(bond_dim**2)
        complex(kind(0d0)), intent(in) :: input_tensor(bond_dim**4)
        complex(kind(0d0)), intent(out) :: output

        integer i, j
        integer, allocatable :: phase_1(:), phase_2(:)
        complex(kind(0d0)), allocatable :: temporary(:)
        complex(kind(0d0)) summ

        allocate(phase_1(size_input(1)))
        allocate(phase_2(size_input(2)))

        allocate(temporary(size_input(1)*size_input(2)*size_input(1)*size_input(2)))

        temporary = (0d0,0d0)

        call set_phase(size_input(1),size_input_even(1),phase_1)
        call set_phase(size_input(2),size_input_even(2),phase_2)

        do i = 1, size(temporary)

            temporary(i) = input_tensor(i)

        enddo

        call multiply_bond_weights(size_input,input_bw_1,input_bw_2,bond_dim,temporary)  

        !!! Phase from reordering

        do j = 1, size_input(2)

            do i = 1, size_input(1)

                temporary(i+(j-1)*size_input(1)+(i-1)*size_input(1)*size_input(2)+(j-1)*size_input(1)*size_input(2)*size_input(1)) &
                = temporary(i+(j-1)*size_input(1)+(i-1)*size_input(1)*size_input(2)+(j-1)*size_input(1)*size_input(2)*size_input(1)) &
                * (-1d0)**(phase_2(j)*phase_1(i)) &
                * (-1d0)**(phase_2(j)*b_cond)

            enddo

        enddo

        summ = (0d0,0d0)

        do i = 1, size_input(1)*size_input(2)

            summ = summ + temporary(i+(i-1)*size_input(1)*size_input(2))

        enddo

        !!! Recording

        output = summ

        deallocate(phase_1,phase_2)
        deallocate(temporary)

    end subroutine

    subroutine get_lnz(n_iter,n_iter_max,normalization_factor,output)
    
        implicit none
    
        integer, intent(in) :: n_iter, n_iter_max
        complex(kind(0d0)), intent(in) :: normalization_factor(0:n_iter_max)
        complex(kind(0d0)), intent(out) :: output

        integer i
        complex(kind(0d0)) summ
        complex(kind(0d0)), allocatable :: temporary(:)
        

        allocate(temporary(0:n_iter))

        temporary = (0d0,0d0)

    
        do i = 0, n_iter
    
            temporary(i) = log(normalization_factor(i))
    
        enddo
    
        summ = (0d0,0d0)
    
        do i = 0, n_iter
    
            summ = summ + dble(temporary(i)/(2d0**i))
    
        enddo
    
        !!! Recording

        output = summ
    
    end subroutine 

    subroutine get_relative_error(n_iter,input_mass,input_mu,result_trg,rel_error)

        implicit none

        integer, intent(in) :: n_iter
        double precision, intent(in) :: input_mass, input_mu
        complex(kind(0d0)), intent(in) :: result_trg
        double precision, intent(out) :: rel_error
        
        integer i, j, system_size
        double precision pi
        double precision, allocatable :: series_1(:), series_2(:), p_1(:), p_2(:)
        complex(kind(0d0)) im, exact
        complex(kind(0d0)), allocatable :: a_11(:,:), a_12(:,:), a_21(:,:), a_22(:,:)
        complex(kind(0d0)), allocatable :: det_D(:,:)

        pi = 4d0*atan(1d0)
        im = (0d0,1d0)

        if ( n_iter == 0 ) then
            
            system_size = 1

        elseif ( n_iter /= 0 ) then

            system_size = 2**(n_iter-1)

        endif

        allocate(series_1(system_size))
        allocate(series_2(system_size))
        allocate(p_1(system_size))
        allocate(p_2(system_size))
        allocate(a_11(system_size,system_size))
        allocate(a_12(system_size,system_size))
        allocate(a_21(system_size,system_size))
        allocate(a_22(system_size,system_size))
        allocate(det_D(system_size,system_size))
    
        series_1 = 1d0
        series_2 = 1d0
        p_1 = 0d0
        p_2 = 0d0

        a_11 = (0d0,0d0)
        a_12 = (0d0,0d0)
        a_21 = (0d0,0d0)
        a_22 = (0d0,0d0)
        det_D = (0d0,0d0)

        do i = 1, system_size-1

            series_1(i+1) = series_1(i) + 1d0
            series_2(i+1) = series_2(i) + 1d0

        enddo

        do i = 1, system_size

            p_1(i) = 2d0*pi*series_1(i)/dble(system_size)

            p_2(i) = 2d0*pi*series_2(i)/dble(system_size)
            p_2(i) = p_2(i)+pi/dble(system_size) !!! ANTI-PERIODIC
    
        enddo

        do j = 1, system_size

            do i = 1, system_size
    
                a_11(i,j) = -cos(p_1(i))-cos(p_2(j))*cosh(input_mu)-im*sin(p_2(j))*sinh(input_mu)+input_mass+2d0
                a_22(i,j) = a_11(i,j)
    
                a_12(i,j) = im*sin(p_1(i)) + sin(p_2(j))*cosh(input_mu)-im*cos(p_2(j))*sinh(input_mu)
                a_21(i,j) = im*sin(p_1(i)) - sin(p_2(j))*cosh(input_mu)+im*cos(p_2(j))*sinh(input_mu)
    
                det_D(i,j) = a_11(i,j)*a_22(i,j)-a_12(i,j)*a_21(i,j)
    
            enddo
    
        enddo

        do j = 1, system_size

            do i = 1, system_size
    
                exact = exact + log(det_D(i,j))/dble(system_size*system_size)
    
            enddo
    
        enddo

        rel_error = abs(exact-result_trg)/abs(exact)

        deallocate(series_1,series_2)
        deallocate(p_1,p_2)
        deallocate(a_11,a_12,a_21,a_22)
        deallocate(det_D)

    end subroutine

    subroutine multiply_bond_weights(size_input,input_bw_1,input_bw_2,bond_dim,array)

        use module_lapack

        implicit none

        integer, intent(in) :: size_input(2), bond_dim
        double precision, intent(in) :: input_bw_1(bond_dim**2), input_bw_2(bond_dim**2)
        complex(kind(0d0)), intent(inout) :: array(size_input(1)*size_input(2)*size_input(1)*size_input(2))

        integer i, j
        complex(kind(0d0)), allocatable :: bw_1(:), bw_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:), temporary_3(:)

        allocate(bw_1(size_input(1)*size_input(1)))
        allocate(bw_2(size_input(2)*size_input(2)))
        allocate(temporary_1(size(array)))
        allocate(temporary_2(size(array)))
        allocate(temporary_3(size(array)))

        do i = 1, size(bw_1)

            bw_1(i) = dcmplx(input_bw_1(i))

        enddo

        do i = 1, size(bw_2)

            bw_2(i) = dcmplx(input_bw_2(i))

        enddo

        call general_zgemm('T','N',size_input(2)*size_input(1)*size_input(2),size_input(1),size_input(1),array,bw_1,temporary_1)
        call general_zgemm('T','N',size_input(1)*size_input(2)*size_input(1),size_input(2),size_input(2),temporary_1,bw_2,temporary_2)

        do j = 1, size_input(1)*size_input(2)

            do i = 1, size_input(1)*size_input(2)

                array(j+(i-1)*size_input(1)*size_input(2)) &
                = temporary_2(i+(j-1)*size_input(1)*size_input(2))

            enddo

        enddo

        deallocate(bw_1,bw_2)
        deallocate(temporary_1)
        deallocate(temporary_2)
        deallocate(temporary_3)

    end subroutine

end module module_measurement