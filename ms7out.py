from ase import Atom, Atoms
from ase.units import Bohr, Ha
from ase.cell import Cell
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import reader

@reader
def read_ms7_out(fd):
    lines = fd.readlines()
    
    atoms = Atoms()
    searchmode = False
    mode = ""
    readcoor = False
    readforce = False
    readlistatom = True
    readpbc = False
    listatom = []
    listcoor = []
    listforce = []
    Epot = 0.0
    index = 0
    pbc = False
    lattice = []
    
    ### Read the information
    for line in lines:
        if "$cell vectors" in line:
            readpbc = True
            pbc = True
        if readpbc:
            data = list(filter(None, line.split()))
            if len(data) == 3:
                lattice.append([float(data[0])*Bohr, float(data[1])*Bohr, float(data[2])*Bohr])
        if "$coordinates" in line:
            readpbc = False
            if not pbc:
                lattice.append([0, 0, 0], [0, 0, 0], [0, 0, 0])
        if "# Task parameters" in line:
            searchmode = True
        if "Calculate" in line and searchmode:
            mode = list(filter(None, line.split()))[1]
            searchmode = False
        if "ATOM               X                   Y                   Z" in line and not readcoor:
            readcoor = True
        if readcoor:
            data = list(filter(None, line.split()))
            if len(data) == 5:
                if readlistatom:
                    listatom.append(data[1])
                    listcoor.append([float(data[2]), float(data[3]), float(data[4])])
                else:
                    listcoor[index] = [float(data[2]), float(data[3]), float(data[4])]
                    index += 1
        if "-------------------------------------------" in line and readcoor:
            readcoor = False
            readlistatom = False
            index = 0
        if line.startswith("df "):
            if readforce:
                if "binding energy" in line:
                    readforce = False
                    index = 0
                data = list(filter(None, line.split()))
                if len(data) == 8:
                    if len(listforce) == len(listatom):
                        listforce[index] = [float(data[5])*Ha/Bohr, float(data[6])*Ha/Bohr, float(data[7])*Ha/Bohr]
                        index += 1
                    else:
                        listforce.append([float(data[5])*Ha/Bohr, float(data[6])*Ha/Bohr, float(data[7])*Ha/Bohr])              
            else:
                if "ATOMIC  COORDINATES" in line:
                    readforce = True
        if "Total DFT-D energy" in line:
            Epot = float(list(filter(None, line.split()))[4].split("Ha")[0]) * Ha
    
    ##construct the Atoms
    for i in range(0,len(listatom), 1):
        atom = Atom(listatom[i], listcoor[i])
        atoms.append(atom)
        if len(listforce) < len(listatom):
            listforce.append([0, 0, 0])
    cell = Cell.new(lattice)
    atoms.set_cell(cell)
    atoms.set_pbc(pbc)
    atoms.calc = SinglePointDFTCalculator(atoms, energy = Epot, forces=listforce)
    atoms.calc.name = 'dmol'
    return atoms
