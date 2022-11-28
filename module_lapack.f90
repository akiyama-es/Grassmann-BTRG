module module_lapack

    contains

    subroutine general_zgemm(transa,transb,size_result_row,size_result_column,size_sum,a,b,c)

        implicit none

        character(1), intent(in) :: transa, transb
        integer, intent(in) :: size_result_row, size_result_column, size_sum
        double complex c1, c2
        double complex, dimension(:), intent(in) :: a, b
        double complex, dimension(:), intent(out) :: c

        c1 = dcmplx(1,0)
        c2 = dcmplx(0,0)

        c = (0d0,0d0)

        if ( transa == 'N' .and. transb == 'N' ) then

            call zgemm('N','N',size_result_row,size_result_column,size_sum,c1,a,size_result_row,b,size_sum,c2,c,size_result_row)

        elseif ( transa == 'N' .and. transb == 'T' ) then

            call zgemm('N','T',size_result_row,size_result_column,size_sum,c1,a,size_result_row,b,size_result_column,c2,c,size_result_row)

        elseif ( transa == 'T' .and. transb == 'N' ) then

            call zgemm('T','N',size_result_row,size_result_column,size_sum,c1,a,size_sum,b,size_sum,c2,c,size_result_row)

        elseif ( transa == 'T' .and. transb == 'T' ) then

            call zgemm('T','T',size_result_row,size_result_column,size_sum,c1,a,size_sum,b,size_result_column,c2,c,size_result_row)

        else 

            print*, " ERROR IN DGEMM "
            stop

        endif

    end subroutine 

    subroutine zsvd(a,m,n,s,u,v)
    
        implicit none
    
        character jobu, jobvt
        integer, intent(in) :: m, n
        integer lda, ldu, ldvt, lwork, info
        integer i, j
        complex(kind(0d0)), intent(in) :: a(m*n)
        double precision, allocatable :: rwork(:)
        double precision, intent(out) :: s(min(m,n))
        complex(kind(0d0)), intent(out) :: u(m*min(m,n)), v(n*min(m,n))
        complex(kind(0d0)), allocatable :: vt(:), work(:)

        allocate(vt(min(m,n)*n))
        allocate(rwork(4*max(1,5*min(m,n))))

        vt = (0d0,0d0)
        rwork = 0d0
    
        s = 0d0
        u = (0d0,0d0)
        v = (0d0,0d0)

        jobu = 'S' 
        jobvt = 'S'
        lda = max(1,m)
        ldu = m 
        ldvt = min(m,n)
    
        lwork = -1
    
        allocate(work(max(1,lwork)))
        work = (0d0,0d0)

        call zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)

        if ( int(work(1)) >= max(1,2*min(m,n)+max(m,n)) ) then
    
           lwork = int(work(1))

        else

           !print*, ' QUERY DOES NOT WORK '
           lwork = 2*max(1,2*min(m,n)+max(m,n))

        endif
    
        deallocate(work)
    
        allocate(work(max(1,lwork)))
        work = (0d0,0d0)
    
        call zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
    
        if ( info /= 0 ) then
            
            print *, ' ERRORS IN SVD PROCESS : INFO = ', info
            print *, ' THE PROGRAM IS STOPPED ' 
    
            stop
    
        endif

        do i = 1, min(m,n)

            if ( s(1) /= 0d0 ) then

                if ( s(i)/s(1) <= 1e-16 ) then
 
                    s(i) = 0d0
 
                    do j = 1, m
 
                        u(j+(i-1)*m) = (0d0,0d0)
 
                    enddo
 
                    do j = 1, n
 
                        vt(i+(j-1)*min(m,n)) = (0d0,0d0)
 
                    enddo

                endif
 
            endif
 
        enddo

        do j = 1, n

            do i = 1, min(m,n)

                v(j+(i-1)*n) = vt(i+(j-1)*min(m,n))

            enddo

        enddo

        deallocate(vt,work)
    
    end subroutine 

end module module_lapack
