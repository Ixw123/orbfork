from utils import generate_configuration
from mace.calculators import mace_mp

mp = mace_mp(model='large', default_dtype='float64', device='cuda')

smiles_strings = ["[CH4]", "O=O"] # ["[HH]", "O=O"]
N = 100
counts = [N, 2*N] # [2*N, N]
# use string to get density of mmoleucles in the unit cell moles/vol = density (pV=n*r*t)
dimensions = [25, 25, 25]
conformers = 25
temperature = 873.15 # 800
reaction_name = "CH4-O2"
generate_configuration(smiles_strings, counts, mp, dimensions, conformers, temperature, outfile=f"{reaction_name}_{temperature}.xyz")
