module module_measurement

    contains

    subroutine get_trace()

        use module_declaration
        use module_phase
    
        implicit none

        integer i, j
        integer, allocatable :: phase_1(:), phase_2(:)
        complex(kind(0d0)), allocatable :: temporary(:)
        complex(kind(0d0)) summ

        allocate(phase_1(bond_new(1)))
        allocate(phase_2(bond_new(2)))

        allocate(temporary(size(tensor)))

        temporary = (0d0,0d0)

        call set_basic_phase_new(phase_1,phase_2)

        do i = 1, size(temporary)

            temporary(i) = tensor(i)

        enddo

        call multiply_bond_weights(temporary)  

        !!! Phase from reordering

        do j = 1, bond_new(2)

            do i = 1, bond_new(1)

                temporary(i+(j-1)*bond_new(1)+(i-1)*bond_new(1)*bond_new(2)+(j-1)*bond_new(1)*bond_new(2)*bond_new(1)) &
                = temporary(i+(j-1)*bond_new(1)+(i-1)*bond_new(1)*bond_new(2)+(j-1)*bond_new(1)*bond_new(2)*bond_new(1)) &
                * (-1d0)**(phase_2(j)*phase_1(i))

            enddo

        enddo

        !!! Phase from boundary condition

        do j = 1, bond_new(2)

            do i = 1, bond_new(1)

                temporary(i+(j-1)*bond_new(1)+(i-1)*bond_new(1)*bond_new(2)+(j-1)*bond_new(1)*bond_new(2)*bond_new(1)) &
                = temporary(i+(j-1)*bond_new(1)+(i-1)*bond_new(1)*bond_new(2)+(j-1)*bond_new(1)*bond_new(2)*bond_new(1)) &
                * (-1d0)**(phase_2(j)*boundary_condition)

            enddo

        enddo
        
        summ = (0d0,0d0)

        do i = 1, bond_new(1)*bond_new(2)

            summ = summ + temporary(i+(i-1)*bond_new(1)*bond_new(2))

        enddo

        !!! Recording

        part(iter) = summ

        deallocate(phase_1,phase_2)
        deallocate(temporary)

    end subroutine

    subroutine get_lnz()

        use module_declaration
    
        implicit none
    
        integer i
        complex(kind(0d0)) summ
        complex(kind(0d0)), allocatable :: temporary(:)

        allocate(temporary(0:iter))

        temporary = (0d0,0d0)

    
        do i = 0, iter
    
            temporary(i) = log(part(i))
    
        enddo
    
        summ = (0d0,0d0)
    
        do i = 0, iter
    
            summ = summ + dble(temporary(i)/(2d0**i))
    
        enddo
    
        !!! Recording

        lnz(iter) = summ
    
    end subroutine 

    subroutine get_relative_error()

        use module_declaration

        implicit none

        integer i, j, system_size
        double precision pi
        double precision, allocatable :: series_1(:), series_2(:), p_1(:), p_2(:)
        complex(kind(0d0)) im
        complex(kind(0d0)), allocatable :: a_11(:,:), a_12(:,:), a_21(:,:), a_22(:,:)
        complex(kind(0d0)), allocatable :: det_D(:,:)

        pi = 4d0*atan(1d0)
        im = (0d0,1d0)

        if ( iter == 0 ) then
            
            system_size = 1

        elseif ( iter /= 0 ) then

            system_size = 2**(iter-1)

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
    
                a_11(i,j) = -cos(p_1(i))-cos(p_2(j))*cosh(mu)-im*sin(p_2(j))*sinh(mu)+mass+2d0
                a_22(i,j) = a_11(i,j)
    
                a_12(i,j) = im*sin(p_1(i)) + sin(p_2(j))*cosh(mu)-im*cos(p_2(j))*sinh(mu)
                a_21(i,j) = im*sin(p_1(i)) - sin(p_2(j))*cosh(mu)+im*cos(p_2(j))*sinh(mu)
    
                det_D(i,j) = a_11(i,j)*a_22(i,j)-a_12(i,j)*a_21(i,j)
    
            enddo
    
        enddo

        do j = 1, system_size

            do i = 1, system_size
    
                exact(iter) = exact(iter) + log(det_D(i,j))/dble(system_size*system_size)
    
            enddo
    
        enddo

        relative_error(iter) = abs(exact(iter)-lnz(iter))/abs(exact(iter))

        deallocate(series_1,series_2)
        deallocate(p_1,p_2)
        deallocate(a_11,a_12,a_21,a_22)
        deallocate(det_D)

    end subroutine

    subroutine multiply_bond_weights(array)

        use module_declaration
        use module_lapack

        implicit none

    
        complex(kind(0d0)), intent(inout) :: array(bond_new(1)*bond_new(2)*bond_new(1)*bond_new(2))

        integer i, j
        complex(kind(0d0)), allocatable :: bw_1(:), bw_2(:)
        complex(kind(0d0)), allocatable :: temporary_1(:), temporary_2(:), temporary_3(:)

        allocate(bw_1(bond_new(1)*bond_new(1)))
        allocate(bw_2(bond_new(2)*bond_new(2)))
        allocate(temporary_1(size(array)))
        allocate(temporary_2(size(array)))
        allocate(temporary_3(size(array)))

        do i = 1, size(bw_1)

            bw_1(i) = dcmplx(bond_weight_1(i))

        enddo

        do i = 1, size(bw_2)

            bw_2(i) = dcmplx(bond_weight_2(i))

        enddo

        call general_zgemm('T','N',bond_new(2)*bond_new(1)*bond_new(2),bond_new(1),bond_new(1),array,bw_1,temporary_1)
        call general_zgemm('T','N',bond_new(1)*bond_new(2)*bond_new(1),bond_new(2),bond_new(2),temporary_1,bw_2,temporary_2)
        !call general_dgemm('T','N',bond_new(1)*bond_new(2)*bond_new(1),bond_new(1),bond_new(1),temporary_2,bw_1,temporary_3)
        !call general_dgemm('T','N',bond_new(1)*bond_new(2)*bond_new(1),bond_new(2),bond_new(2),temporary_3,bw_2,array)

        do j = 1, bond_new(1)*bond_new(2)

            do i = 1, bond_new(1)*bond_new(2)

                array(j+(i-1)*bond_new(1)*bond_new(2)) &
                = temporary_2(i+(j-1)*bond_new(1)*bond_new(2))

            enddo

        enddo

        deallocate(bw_1,bw_2)
        deallocate(temporary_1)
        deallocate(temporary_2)
        deallocate(temporary_3)

    end subroutine

end module module_measurement