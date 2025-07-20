import sys

def getRatio(xyz_file):
    print(f"Processing file: {xyz_file}")
    """
    Reads an XYZ file and calculates the ratio of the number of atoms in the inner box and the number of atoms in the outer box.
    """
    with open(xyz_file, 'r') as file:
        lines = file.readlines()
    
    numInnerAtoms = getNumAtoms(lines[2:], 10.0)  # Assuming a distance of 10.0 for inner box
    print(f"Number of inner atoms: {numInnerAtoms}")
    numOutterAtoms = getNumAtoms(lines[2:], 20.0)  # Assuming a distance of 20.0 for outer box
    
    return numOutterAtoms / numInnerAtoms

def getDistance(atom):
    """
    Extracts the distance of an atom from the center.
    """
    # print(f"Calculating distance for atom: {atom}")
    # Assuming the atom is in the format: 'atom_type x y z'
    parts = atom.split()
    if len(parts) != 4:
        print(f"Invalid atom format: {atom}")
        return float('inf')  # Return infinity for invalid atoms
    # print(f"Parts: {parts}")
    x, y, z = map(float, parts[1:])
    
    return (x**2 + y**2 + z**2)**0.5  # Euclidean distance from origin

def getNumAtoms(atoms, distance):
    """
    Counts the number of atoms within a certain distance from the center.
    """
    count = 0
    for atom in atoms:
        if getDistance(atom) <= distance:
            count += 1
            
    return count

def main():
    xyz1 = sys.argv[1]
    xyz2 = sys.argv[2]

    r1 = getRatio(xyz1)
    r2 = getRatio(xyz2)
    
    print(f"Ratio for {xyz1}: {r1}")
    print(f"Ratio for {xyz2}: {r2}")

if __name__ == "__main__":
    # Entry point for the script
    print("This script compares XYZ ratios.")
    # Add your main logic here
    main()