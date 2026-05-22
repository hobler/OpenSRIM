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