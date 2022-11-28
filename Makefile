  f90       = ifort -mkl
  MYFLAGS   = 

  FFLAGS      = $(MYFLAGS)
  LDFLAGS=


.SUFFIXES:
.SUFFIXES: .a .o .f90 .mod

.f90.o:
	$(f90) $(FFLAGS) -c $*.f90

OBJS = \
	module_lapack.o\
	module_blocking_operations.o\
	module_declaration.o\
	module_time.o\
	module_record.o\
	module_phase.o\
	module_measurement.o\
	module_initial_tensor.o\
	module_setup.o\
	module_trg.o\
	main.o \

GBTRG: $(OBJS)
	$(f90) -o $@ $(OBJS)

clean:
	rm *.o *.mod *.lst GBTRG
