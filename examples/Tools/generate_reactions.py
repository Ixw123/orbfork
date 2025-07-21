from utils import generate_configuration
from mace.calculators import mace_mp

# Also add the number of atoms in the system or a full breackdown of the types of atoms and numbers of each

mp = mace_mp(model='large', default_dtype='float64', device='cuda')

smiles_strings = ["[CH4]", "O=O"] # ["[HH]", "O=O"]
N = 50
counts = [N, 2*N] # [2*N, N]
# use string to get density of mmoleucles in the unit cell moles/vol = density (pV=n*r*t)

atomCountsStr = ""

for i, smiles in enumerate(smiles_strings):
    atomCountsStr += f"{counts[i]}_{smiles}_"
    
atomCountsStr = atomCountsStr[:-1]  # Remove the trailing underscore
print("Atom counts string:", atomCountsStr)

dimensions = [30, 30, 30]
conformers = 25
temperature = 800
reaction_name = "CH4-O2"
# Assumes the dimensions are hemogenous 
generate_configuration(smiles_strings, counts, mp, dimensions, conformers, temperature, outfile=f"./examples/{reaction_name}_{dimensions[0]}_{temperature}_{atomCountsStr}.xyz")
