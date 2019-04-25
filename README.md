# Simulation of hydrogenic atoms

This software is used in Quantum Chemistry course in University of
Buenos Aires. For more reference, go to the
[webpage](http://materias.df.uba.ar/e3a2016c1/) [in Spanish!]. It
ports a numerical solution for central potentials as described in
[Salvat and Mayol, Comp. Phys. Comm. 62 (65) 1991][1]

[1]: http://www.sciencedirect.com/science/article/pii/0010465591901222

DISCLAIMER: The files here do not follow any guideline respecting to
the usual python module distribution. It is done so in purpose, so
students *know* what the sources are.


## Installation
After downloading the package, go to `src/` folder and run `make
library`

```
cd src/
make library
```

## Running
In the root folder of the project, there is an example called
`hydrogen.py` that calculates some bound states for the hydrogen atom
in the un-screened case. You can (*and should*) take a look at its
source code and modify it properly.

## Help
The file `atomic.py` has the help embedded in the docstrings

## Another way
In order to keep the previous format working, you can also build the old
(and less user friendly) executable file with:

```
cd src/
make executable
```

After that, you can run the file `./hydrogenic.e` and work with the
executable file directly. Keep in mind that if you do that, you'll
have to manually extract all info like the wavefunction and the
energy. We insist on *not using* this file, except it is strictly
necessary.
