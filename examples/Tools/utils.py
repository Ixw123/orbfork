from dataclasses import dataclass
import numpy as np
import logging
from tqdm import tqdm
import numpy.typing as npt
import urllib.request
import os
from typing import Optional, Union
from os.path import join
from rdkit import Chem
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem.rdDetermineBonds import DetermineConnectivity
from ase import Atoms, optimize
from ase.io import write
from ase.io.extxyz import write_extxyz
from sklearn.linear_model import LinearRegression
from scipy.spatial.transform import Rotation
from scipy.stats import maxwell


@dataclass
class Molecule3D:
    smiles_or_inchi: str
    atoms: npt.NDArray[np.int32]
    coordinates: npt.NDArray[np.float64]


def XYZBlockFromMolecule3D(mol: Molecule3D) -> str:
    xyz_block = ""
    for atom, coord in zip(mol.atoms, mol.coordinates):
        xyz_block += f"{atom} {coord[0]} {coord[1]} {coord[2]}\n"
    return xyz_block


def calculate_conformers(
        mol: Chem.Mol,
        conformers: int = 1
    ) -> Optional[Chem.Mol]:
    """Calculate a conformer for an RDKit molecule."""
    mol.UpdatePropertyCache()
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)

    try:
        AllChem.EmbedMultipleConfs(mol, conformers)
        # AllChem.MMFFOptimizeMolecule(mol)
        return mol
    except Exception as e:
        print(e)
        return None


def optimize_conformer(
        configuration: Atoms,
        calculator = None,
        min_fmax: float = 0.01,
        max_steps: int = 1000,
        time_step: float = 0.1,
        ase_optimize: bool = True,
        gd_optimize: bool = False,
        outfile: Optional[str] = None,
        logfile: str = '/dev/null',
    ) -> Atoms:
    """
    Optimize an ASE configuration with a given ASE MMFF calculator.
    """
    if calculator is not None:
        configuration.calc = calculator

    if ase_optimize:
        try:
            dyn = optimize.sciopt.SciPyFminCG(configuration, logfile=logfile)
            dyn.run(fmax=min_fmax, steps=max_steps)
        except Exception as e:
            pass
    if gd_optimize:
        step = 0
        fmax_prev = np.Inf

        prev_config = Atoms()
        if calculator is not None:
            prev_config.calc = calculator
        while step < max_steps:
            # F = dp/dt => p_f = p_i + F*dt
            forces = configuration.get_forces()
            momenta = configuration.get_momenta()
            momenta += forces * time_step
            configuration.set_momenta(momenta)

            # v = dx/dt => x_f = x_i + v*dt
            velocities = configuration.get_velocities()
            positions = configuration.get_positions()
            positions += velocities * time_step
            configuration.set_positions(positions)
            positions = configuration.get_positions()

            fmax = np.max(np.sqrt(np.sum(forces ** 2, axis=1)))
            if abs(fmax - min_fmax) < 1e-5:
                break
            if fmax - fmax_prev > 1e-5:
                time_step /= 2
                configuration = prev_config
            fmax_prev = fmax
            min_fmax = min(fmax, min_fmax)
            step += 1
            prev_config = Atoms(configuration)
            if calculator is not None:
                prev_config.calc = calculator
    if outfile is not None:
        write(outfile, configuration, append=True)
    return configuration


def optimize_smiles(
        smiles: str,
        calculator,
        min_fmax: float = 0.01,
        max_steps: int = 1000,
        time_step: float = 0.1,
        ase_optimize: bool = True,
        gd_optimize: bool = False,
        conformers: int = 1,
        outfile: Optional[str] = None,
        logfile: str = '/dev/null'
    ) -> list[Atoms]:
    """
    Optimize a SMILES string with a given ASE MMFF calculator, possibly generating multiple conformers. Returns an ase.Atoms.
    """
    mol: Optional[Chem.Mol] = None
    try:
        mol = Chem.MolFromSmiles(smiles)         # RDKit.Mol
        mol = calculate_conformers(mol, conformers)
    except ValueError:
        print(f'Can\'t generate conformer for {smiles}')
        return None

    configurations = []
    for i in range(mol.GetNumConformers()):
        m = Chem.Mol(mol, confId=i)
        configuration = optimize_conformer(
            Atoms_from_Mol(Chem.Mol(mol, confId=i)),
            calculator,
            min_fmax,
            max_steps,
            time_step,
            ase_optimize,
            gd_optimize,
            None,
            logfile,
        )
        configurations.append(configuration)
        if outfile is not None:
            configuration.info = {
                'smiles': smiles
            }
            write(outfile, configuration, append=True)
    return configurations


def get_Molecule3D(mol):
    """Get a Molecule3D object from an RDKit molecule.
    If conformer fails random coordinates are generated."""
    atom_positions = []
    atoms = []
    if mol.GetNumConformers() != 0:
        for i, atom in enumerate(mol.GetAtoms()):
            positions = mol.GetConformer().GetAtomPosition(i)
            atom_positions.append(positions)
            atoms.append(atom.GetAtomicNum())

        atom_positions = np.array(atom_positions)
        atoms = np.array(atoms)
        return Molecule3D(
            smiles_or_inchi=Chem.MolToSmiles(mol),
            atoms=atoms,
            coordinates=atom_positions,
        )


    mol_with_conformer = calculate_conformers(mol)

    if mol_with_conformer is not None:
        for i, atom in enumerate(mol.GetAtoms()):
            positions = mol.GetConformer().GetAtomPosition(i)
            atom_positions.append(positions)
            atoms.append(atom.GetAtomicNum())

        atom_positions = np.array(atom_positions)
        atoms = np.array(atoms)
        return Molecule3D(
            smiles_or_inchi=Chem.MolToSmiles(mol),
            atoms=atoms,
            coordinates=atom_positions,
        )
    else:
        atoms = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
        num_atoms = len(atoms)
        random_coords = np.random.rand(num_atoms, 3)
        return Molecule3D(
            smiles_or_inchi=Chem.MolToSmiles(mol),
            atoms=atoms,
            coordinates=random_coords,
        )


def XYZToMolecule3D(xyz_file, file_parser):
    """
    Convert an XYZ file to a Molecule3D object.
    """
    lines = [line.decode("UTF-8") for line in xyz_file.readlines()]
    return file_parser(lines)


def Atoms_from_Molecule3D(mol: Molecule3D) -> Atoms:
    return Atoms(numbers=mol.atoms, positions=mol.coordinates)

def Atoms_from_Mol(mol):
    return Atoms_from_Molecule3D(get_Molecule3D(mol))

def Molecule3D_from_Atoms(mol: Atoms) -> Molecule3D:
    return Molecule3D(atoms=mol.numbers, coordinates=mol.positions)


def Mol_from_Atoms(atoms: Atoms):
    columns = ['symbols', 'positions']
    write_extxyz('/tmp/tmpmolfromatoms.xyz', atoms, columns=columns)
    mol = Chem.rdmolfiles.MolFromXYZFile('/tmp/tmpmolfromatoms.xyz')
    mol.UpdatePropertyCache()
    DetermineConnectivity(mol)
    return mol


def XYZ_Writer(num_atoms, element_list, coordinates, inchi_or_smiles):
    """
    Write an XYZ String
    """
    xyz_string = f"{num_atoms}\n"
    xyz_string += f"{inchi_or_smiles}\n"
    for element, coord in zip(element_list, coordinates):
        xyz_string += f"{element}\t{coord[0]:.8f}\t{coord[1]:.8f}\t{coord[2]:.8f}\n"
    return xyz_string


def download_dataset(
    datadir,
    dataname,
    fname="dsgdb9nsd.xyz.tar.bz2",
):
    """
    Download and prepare the QM9 (GDB9) dataset.
    """
    # Define directory for which data will be output.
    gdb9dir = join(*[datadir, dataname])

    gdb9_url_data = "https://springernature.figshare.com/ndownloader/files/3195389"
    gdb9_tar_data = join(gdb9dir, fname)

    if os.path.exists(gdb9_tar_data):
        logging.info("GDB9 dataset already downloaded!")
        return
    # Important to avoid a race condition
    os.makedirs(gdb9dir, exist_ok=True)
    logging.info(
        "Downloading and processing GDB9 dataset. Output will be in directory: {}.".format(
            gdb9dir
        )
    )

    logging.info("Beginning download of GDB9 dataset!")

    urllib.request.urlretrieve(gdb9_url_data, filename=gdb9_tar_data)
    logging.info("GDB9 dataset downloaded successfully!")


def extract_tarfile(fname, outputdir, file_parser):
    """
    Extract a tarfile to a specified output directory.
    """
    import tarfile

    logging.info(
        "Extracting tarfile: {} to output directory: {}".format(tarfile, outputdir)
    )

    mols = []
    if tarfile.is_tarfile(fname):
        logging.info("File is a valid tarfile.")
        tardata = tarfile.open(fname, "r")
        file = tardata.getmembers()
        count = 0
        for f in tqdm(file):
            mol_data = XYZToMolecule3D(tardata.extractfile(f), file_parser)
            mols.append(mol_data)
            count += 1
        tardata.close()
    else:
        logging.error("File is not a valid tarfile. Exiting extraction.")

    logging.info("Extraction complete!")
    return mols


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def ensure_xyz_format(lines):
    """
    Ensure that the XYZ format is correct.
    """
    n_atoms = int(lines[0])


def get_bde(
        mol: Chem.Mol,
        z1: int,
        z2: int,
        energy_calc,
        enthalpy_calc,
        outfile: Optional[str] = None,
        logfile: str = '/dev/null',
        conversion_factor: float = 1,
        enthalpy_shift: float = 0
    ) -> dict[str, Union[float, list[float]]]:
    """
    Gets the bond dissociation energies when removing z1 (bonded to z2)
    """
    atoms = Atoms_from_Mol(mol)
    atoms.calc = enthalpy_calc
    enthalpy_initial = atoms.get_potential_energy() * conversion_factor

    test_confs = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == z1 and bond.GetEndAtom().GetAtomicNum() == z2:
            test_mol = Chem.RWMol(mol)
            test_mol.RemoveAtom(bond.GetBeginAtomIdx())
            atoms = Atoms_from_Mol(test_mol)
            test_confs.append(
                optimize_conformer(
                    atoms,
                    energy_calc,
                    outfile=None,
                    logfile=logfile
                )
            )
        elif bond.GetBeginAtom().GetAtomicNum() == z2 and bond.GetEndAtom().GetAtomicNum() == z1:
            test_mol = Chem.RWMol(mol)
            test_mol.RemoveAtom(bond.GetEndAtomIdx())
            test_confs.append(
                optimize_conformer(
                    Atoms_from_Mol(test_mol),
                    energy_calc,
                    outfile=None,
                    logfile=logfile
                )
            )

    enthalpies_final = []
    for test_conf in test_confs:
        test_conf.calc = enthalpy_calc
        enthalpy_final = test_conf.get_potential_energy() * conversion_factor
        enthalpies_final.append(enthalpy_final)
        #print(f'h_f = {energies_final[-1]}')
        if outfile is not None:
            test_conf.info['h_i'] = enthalpy_initial
            test_conf.info['h_f'] = enthalpy_final
            test_conf.info['smiles'] = atoms.info['smiles']
            test_conf.calc = energy_calc
            test_conf.write(outfile, append=True)
    return {
            'bdes': [enthalpy_final - enthalpy_initial + enthalpy_shift for enthalpy_final in enthalpies_final],
            'h_i': enthalpy_initial,
            'h_f': enthalpies_final
            }


def max_magnitude(l: list[list[float]]) -> float:
    return np.max(np.sqrt(np.sum(l**2, axis=-1)))


def linear_fit(x: list[float], y: list[float]) -> dict[str, float]:
    x = np.array(x).reshape(-1, 1)
    reg = LinearRegression().fit(x.reshape(-1, 1), y.reshape(-1, 1))
    return {
            'm': reg.coef_.item(),
            'b': reg.intercept_.item(),
            'R^2': reg.score(x.reshape(-1, 1), y.reshape(-1, 1))
            }

def generate_configuration(
        smiles_strings: list[str],
        counts: list[int],
        calc,
        dimensions: list[float],
        conformers: int = 1,
        temperature: float = 293.15,
        outfile: Optional[str] = None,
        seed: Optional[int] = None
        ) -> Atoms:
    """
    smiles_strings: string[MOLECULE_TYPES], containing the smiles representations of the molecules
    counts: int[MOLECULE_TYPES], where counts[i] = the number of smiles_strings[i]s to generate
    calc: ase.Calculator, for generating the conformations and checking the energy after generation
    dimensions: float[3], the dimensions of the bounding box. eg: (0, 0, 0) to (a, b, c)
    conformers: int, number of conformers to generate for each molecule
    temperature: float, temperature in Kelvin
    outfile: str, path to write the generated configuration to, if any
    returns the generated configuration
    """
    configuration = Atoms(cell=dimensions, pbc=True)
    configuration.info["T_K"] = temperature
    mols = {}
    
    # generate conformations
    for smiles_string in smiles_strings:
        if smiles_string in mols:
            mols[smiles_string].append(optimize_smiles(smiles_string, calc, conformers=conformers))
        else:
            mols[smiles_string] = optimize_smiles(smiles_string, calc, conformers=conformers)

    rng = np.random.default_rng(seed=seed)

    for (smiles_string, count) in zip(smiles_strings, counts):
        for _ in range(count):
            atoms = Atoms()
            # select conformation
            for atom in rng.choice(mols[smiles_string]):
                atoms.extend(atom)

            # rotate atoms around Center Of Positions
            (phi, theta, psi) = Rotation.random(rng=rng).as_euler('zxy', degrees=True)
            atoms.euler_rotate(phi, theta, psi, center='COP')

            # translate atoms
            translation = [0, 0, 0]
            radius = max_magnitude(atoms.get_positions())
            for i in range(3):
                distance = dimensions[i] - radius
                translation[i] = rng.uniform(low=0, high=distance)
            atoms.translate(translation)

            # set velocities according to temperature
            # Use center of mass to have this and calc vels from center of mass translation of vectors
            velocity_direction = Rotation.random(rng=rng)
            speed = maxwell.rvs()
            velocity = velocity_direction.apply([speed, 0, 0])
            atoms.set_velocities(velocity)

            # add atoms to simulation
            configuration.extend(atoms)

    if outfile is not None:
        configuration.write(outfile)
    return configuration
