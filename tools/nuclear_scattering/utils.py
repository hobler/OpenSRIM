import os


atom = [ None,
          'H', 'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca',
         'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
         'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr',
         'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
         'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
         'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
         'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
         'Pa',  'U' ]



def ask_if_save(fname=None):
    """Ask whether to save the plot. Return filename if yes, None if no."""
    if fname is None:
        answer = input(f"Save plot ? (y/n) ")
        answer = 'o' if answer.lower() == 'y' else 'n'
    else:
        answer = input(f"Save plot as '{fname}'? (y/n/o) ")
    
    if answer.lower() == 'y':
        return fname
    elif answer.lower() == 'o':
        return input("Enter filename to save plot as: ")
    else:
        return None


def get_mass(Z):
    """Get the mass of an atom given its atomic number.

    Parameters:
        Z (int): atomic number of the atom

    Returns:
        (float): mass of the atom in amu
    """
    fname = os.path.join(os.path.dirname(__file__), 
                         '../../data/atom_data/ATOMDATA')
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.split()
            Z_val = int(items[0])
            M_val = float(items[4])
            if Z_val == Z:
                return M_val
            
    raise ValueError(f"Atomic number {Z} not found in atomic masses data.")


def get_density(Z):
    """Get the density of an atom given its atomic number.

    Parameters:
        Z (int): atomic number of the atom

    Returns:
        (float): density of the atom in 1/Å^3
    """
    fname = os.path.join(os.path.dirname(__file__), 
                         '../../data/atom_data/ATOMDATA')
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.split()
            Z_val = int(items[0])
            dens_val = float(items[7]) * 1e-24
            if Z_val == Z:
                return dens_val
            
    raise ValueError(f"Atomic number {Z} not found in atomic densities data.")


