program main

    use module_setup
    use module_record
    use module_time

    implicit none

    call start_clock()

    call load_input_trg()

    call allocate_memory()

    call load_input_init_tensor()

    call initialize_bond_weights()

    call measurement()

    call print_result()

    call btrg_loop()

    call deallocate_memory()

    call end_clock()

    print*, ' NORMAL END '

end program main