module module_initial_tensor

    contains

    subroutine get_init_tensor(input_mu,input_mass,input_coupling,output)

        implicit none

        integer i_1, j_1, i_2, j_2
        integer pi_1, pj_1, pi_2, pj_2
        integer nx, nt, px, pt
        integer i, pow(6)
        integer convert(0:1,0:1)
        double precision, intent(in) :: input_mass, input_coupling, input_mu
        complex(kind(0d0)), allocatable :: temporary_1(:,:,:,:,:,:,:,:), temporary_2(:,:,:,:,:,:,:,:)
        complex(kind(0d0)), intent(out) :: output(4*4*4*4)

        pow(1) = 2

        do i = 1, 5

            pow(i+1) = pow(i) * pow(1)

        enddo

        convert(0,0) = 1
        convert(1,1) = 2
        convert(1,0) = 3
        convert(0,1) = 4

        allocate(temporary_1(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1))
        allocate(temporary_2(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1))

        call get_phase_tensor(input_mu,temporary_1)

        call get_coeff_tensor(input_mass,input_coupling,temporary_2)

        do pj_2 = 0, 1

            do pi_2 = 0, 1

                pt = convert(pi_2,pj_2)

                do pj_1 = 0, 1

                    do pi_1 = 0, 1

                        px = convert(pi_1,pj_1)

                        do j_2 = 0, 1
                            
                            do i_2 = 0, 1

                                nt = convert(i_2,j_2)

                                do j_1 = 0, 1

                                    do i_1 = 0, 1

                                        nx = convert(i_1,j_1)

                                        output(nx+(nt-1)*pow(2)+(px-1)*pow(4)+(pt-1)*pow(6)) &
                                        = temporary_2(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2) &
                                        * temporary_1(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2)
                                        
                                    enddo

                                enddo

                            enddo

                        enddo

                    enddo

                enddo

            enddo

        enddo

        deallocate(temporary_1,temporary_2)

    end subroutine

    subroutine get_phase_tensor(input_mu,output)

        implicit none

        integer i_1, j_1, i_2, j_2
        integer pi_1, pj_1, pi_2, pj_2
        integer, allocatable :: phase(:,:,:,:,:,:,:,:)
        double precision, intent(in) :: input_mu
        complex(kind(0d0)), intent(out) :: output(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1)

        allocate(phase(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1))

        phase = 0

        output = (0d0,0d0)

        do pj_2 = 0, 1

            do pi_2 = 0, 1

                do pj_1 = 0, 1

                    do pi_1 = 0, 1

                        do j_2 = 0, 1
                            
                            do i_2 = 0, 1

                                do j_1 = 0, 1

                                    do i_1 = 0, 1

                                        phase(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2) &
                                        = (-1d0)**(i_1  * (j_1+j_2+pi_1+pi_2)) &
                                        * (-1d0)**(i_2  * (j_2+pi_1+pi_2)) &
                                        * (-1d0)**(pj_1 * (pi_1+pi_2)) &
                                        * (-1d0)**(pj_2 * pi_2)

                                        output(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2) &
                                        = phase(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2) &
                                        * (-1d0)**(pi_1+pi_2) &
                                        * exp(0.5d0*input_mu*dble(i_2-j_2+pi_2-pj_2)) &
                                        * (1d0/sqrt(2d0))**(i_1+j_1+i_2+j_2+pi_1+pj_1+pi_2+pj_2)

                                    enddo

                                enddo

                            enddo

                        enddo

                    enddo

                enddo

            enddo

        enddo

        deallocate(phase)

    end subroutine

    subroutine get_coeff_tensor(input_mass,input_coupling,output)

        implicit none

        integer i_1, j_1, i_2, j_2
        integer pi_1, pj_1, pi_2, pj_2

        double precision, intent(in) :: input_mass, input_coupling
        double precision delta(0:4,0:1)
        complex(kind(0d0)) temporary_1(0:1,0:1,0:1,0:1), temporary_2(0:1,0:1,0:1,0:1)
        complex(kind(0d0)), intent(out) :: output(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1)

        delta = 0d0
        delta(0,0) = 1d0
        delta(1,1) = 1d0

        call get_fundamental_tensor_1(temporary_1)
        call get_fundamental_tensor_2(temporary_2)

        output = (0d0,0d0)

        do pj_2 = 0, 1

            do pi_2 = 0, 1

                do pj_1 = 0, 1

                    do pi_1 = 0, 1

                        do j_2 = 0, 1
                            
                            do i_2 = 0, 1

                                do j_1 = 0, 1

                                    do i_1 = 0, 1

                                        output(i_1,j_1,i_2,j_2,pi_1,pj_1,pi_2,pj_2) &
                                        = ((input_mass+2d0)*(input_mass+2d0)+input_coupling) &
                                        *delta(pj_2+pj_1+i_2+i_1,0)*delta(pi_2+pi_1+j_2+j_1,0) &
                                        -(input_mass+2d0)*delta(pj_2+pj_1+i_2+i_1,1)*delta(pi_2+pi_1+j_2+j_1,1) &
                                        -(input_mass+2d0)*delta(pj_2+pj_1+i_2+i_1,1)*delta(pi_2+pi_1+j_2+j_1,1) &
                                        *(-1d0)**(i_2+i_1+pi_1+j_2)*(0d0,1d0)**(pj_2+i_2+pi_2+j_2) &
                                        -temporary_2(pj_2,pj_1,i_2,i_1)*temporary_1(pi_2,pi_1,j_2,j_1)

                                    enddo

                                enddo

                            enddo

                        enddo

                    enddo

                enddo

            enddo

        enddo

    end subroutine

    subroutine get_fundamental_tensor_1(output)

        implicit none

        integer j_2, j_1
        integer pi_2, pi_1
        integer i
        complex(kind(0d0)) array_1(2), array_2(2), array_3(2), array_4(2)
        complex(kind(0d0)) matrix(4)
        complex(kind(0d0)), intent(out) :: output(0:1,0:1,0:1,0:1)

        array_1(1) = (1d0,0d0)
        array_1(2) = (0d0,1d0)

        array_2(1) = (1d0,0d0)
        array_2(2) = (-1d0,0d0)

        array_3(1) = (1d0,0d0)
        array_3(2) = (0d0,-1d0)

        array_4(1) = (1d0,0d0)
        array_4(2) = (1d0,0d0)


        do j_1 = 0, 1

            do j_2 = 0, 1

                do pi_1 = 0, 1

                    do pi_2 = 0, 1

                        if ( pi_2+pi_1+j_2+j_1 == 2 ) then

                            if ( pi_2 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_1(i)

                                enddo

                                if ( pi_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_2(i)
    
                                    enddo

                                endif

                                if ( j_2 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_3(i)
    
                                    enddo

                                endif

                                if ( j_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            elseif ( pi_1 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_2(i)

                                enddo

                                if ( j_2 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_3(i)
    
                                    enddo

                                endif

                                if ( j_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            elseif ( j_2 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_3(i)

                                enddo

                                if ( j_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            endif

                            output(pi_2,pi_1,j_2,j_1) = matrix(1)*matrix(4)-matrix(2)*matrix(3)

                        endif

                    enddo

                enddo

            enddo

        enddo

    end subroutine

    subroutine get_fundamental_tensor_2(output)

        implicit none

        integer i_2, i_1
        integer pj_2, pj_1
        integer i
        complex(kind(0d0)) array_1(2), array_2(2), array_3(2), array_4(2)
        complex(kind(0d0)) matrix(4)
        complex(kind(0d0)), intent(out) :: output(0:1,0:1,0:1,0:1)

        array_1(1) = (1d0,0d0)
        array_1(2) = (0d0,1d0)

        array_2(1) = (1d0,0d0)
        array_2(2) = (1d0,0d0)

        array_3(1) = (1d0,0d0)
        array_3(2) = (0d0,-1d0)

        array_4(1) = (1d0,0d0)
        array_4(2) = (-1d0,0d0)


        do i_1 = 0, 1

            do i_2 = 0, 1

                do pj_1 = 0, 1

                    do pj_2 = 0, 1

                        if ( pj_2+pj_1+i_2+i_1 == 2 ) then

                            if ( pj_2 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_1(i)

                                enddo

                                if ( pj_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_2(i)
    
                                    enddo

                                endif

                                if ( i_2 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_3(i)
    
                                    enddo

                                endif

                                if ( i_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            elseif ( pj_1 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_2(i)

                                enddo

                                if ( i_2 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_3(i)
    
                                    enddo

                                endif

                                if ( i_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            elseif ( i_2 == 1 ) then

                                do i = 1, 2

                                    matrix(i) = array_3(i)

                                enddo

                                if ( i_1 == 1 ) then

                                    do i = 1, 2

                                        matrix(i+2) = array_4(i)
    
                                    enddo

                                endif

                            endif

                            output(pj_2,pj_1,i_2,i_1) = matrix(1)*matrix(4)-matrix(2)*matrix(3)

                        endif

                    enddo

                enddo

            enddo

        enddo

    end subroutine
    
end module module_initial_tensor