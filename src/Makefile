.PHONY: default library executable all clean
FC = gfortran
F2PY = f2py

default: all

all: library executable

library: radial.so

executable: hydrogenic.e

radial.so: frontend.for radialf.for
	f2py -c $^ -m radial
	cp radial.so ..

hydrogenic.e: hydrogenic.for radialf.for
	gfortran $^ -o $@

clean:
	rm -rfv *.so *.e ../radial.so
