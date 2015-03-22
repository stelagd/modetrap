FC = gfortran

%.o :		%.f90
		$(FC) -c $(FFLAGS) -o $@ $<

% :		%.o
		$(FC) $(FFLAGS) -o $@ $<

clean:
	rm *.o

modetrap_sub:	modetrap_sub.f
		f2py --opt=-O2 -c -m modetrap_sub modetrap_sub.f

