module module_declaration

    !!! Bond dimension

    integer D1

    !!! Variables to control the TRG iteration

    integer maxiter, iter

    !!! Variables to control the size of local tensor

    integer bond(2), bond_even(2)
    integer bond_new(2), bond_even_new(2)
    
    !!! Variables to compute the execution time

    integer time_1, time_2, t_rate, t_max, diff

    !!! Parameters in Gross-Neveu-Wilson model

    double precision mass, mu, ext, coupling

    !!! Hyperparameter in the BTRG
    
    double precision hyper

    !!! Variable to control the boundary condition (0 -> PBC, 1 -> APBC)

    double precision, parameter :: boundary_condition = 1d0

    !!! Bond weights

    double precision, allocatable :: bond_weight_1(:), bond_weight_2(:)

    !!! Rormalization factors

    double precision, allocatable :: f_norm_1(:), f_norm_2(:)

    !!! Local tensor

    complex(kind(0d0)), allocatable :: tensor(:)

    !!! Partition function and free energy

    complex(kind(0d0)), allocatable :: part(:), lnz(:)

    !!! Exact solution for free theory and the relative error

    character(len=1) flag_exact

    complex(kind(0d0)), allocatable :: exact(:)
    double precision, allocatable :: relative_error(:)

end module module_declaration