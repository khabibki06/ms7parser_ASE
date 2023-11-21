# ms7parser_ASE
It is a parser for ASE (Atomic Simulation Environment). It can read the atomic coordinate, lattice parameter (if any), total energy (if any), and forces (if any) from the material studio 7.0 outmol file. 

How to install:
1. Copy the ms7out.py to $HOME/.local/lib/python3.X/site-packages/ase/io
2. Edit the $HOME/.local/lib/python3.X/site-packages/ase/io/format.py 
   a. Find "F = define_io_format"
   b. Add the ms7-out module, by adding this line after F('xyz', 'XYZ-file', '+F') line
   F('ms7-out', 'Material Studio Outmol', '1S', module='ms7out', ext='outmol')

Note: 
1. replace python3.X with the python version; It may be different depending on the OS system and ASE installation 
