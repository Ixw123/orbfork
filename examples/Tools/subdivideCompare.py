import sys

def correctCoordinate(a, unitCell):
    a = unitCell if float(a%unitCell) == 0.0 else float(a % unitCell)
    
    return a
    
def correctXyz(xyz_data, unitCell=25.25):
    
    correctedXyz = []
    
    for atom in xyz_data:
        x, y, z = atom
        # Correct the coordinates to fit within the unit cell
        x = correctCoordinate(x, unitCell)
        y = correctCoordinate(y, unitCell)
        z = correctCoordinate(z, unitCell)
        yield (x, y, z)
        correctedXyz.append((x, y, z))
        
    return correctedXyz

def inSubCell(a, subCell):
    """ 
    Check if a coordinate is within the unit cell.
    Args:        
        a (float): Coordinate to check.
        unitCell (float): Size of the unit cell
    """
    
    if a >= 0 and a <= subCell:
        return 1
    elif a > subCell and a <= 2 * subCell:
        return 2
    elif a > 2 * subCell and a <= 3 * subCell:
        return 3
    elif a > 3 * subCell and a <= 4 * subCell:
        return 4
    else:
        return 0

def determineSubCell(atom, unitCell, subCellSize=.25):
    '''
    Determine the subcell that the atom belongs in and return it as a string for a key
    
    possible subcells are 1, 2, 3, 4
    
    111, 112, 113, 114
    121, 122, 123, 124
    131, 132, 133, 134
    141, 142, 143, 144
    
    211, 212, 213, 214
    221, 222, 223, 224
    231, 232, 233, 234
    241, 242, 243, 244
    
    311. 312, 313, 314
    321, 322, 323, 324
    331, 332, 333, 334
    341, 342, 343, 344
    
    411, 412, 413, 414
    421, 422, 423, 424
    431, 432, 433, 434
    441, 442, 443, 444
    
    '''
    x, y, z = atom
    
    subCell = (unitCell * subCellSize)
    
    xSubCell = inSubCell(x, subCell)
    ySubCell = inSubCell(y, subCell)
    zSubCell = inSubCell(z, subCell)
    if xSubCell == 0 or ySubCell == 0 or zSubCell == 0:
        return None  # Atom is outside the sub-cell range
    
    return str(xSubCell) + str(ySubCell) + str(zSubCell)

def getSubCellCounts(atoms, unitCell, subCellSize=0.25):
    """
    Count the number of atoms in each sub-cell.
    
    Args:
        atoms (list): List of tuples containing XYZ coordinates of atoms.
        unitCell (float): Size of the unit cell.
        subCellSize (float): Size of each sub-cell.
    
    Returns:
        dict: Dictionary with sub-cell coordinates as keys and atom counts as values.
    """
    
    subCellCounts = {}
    
    for atom in atoms:
        
        key = determineSubCell(atom, unitCell, subCellSize)
        
        if key is None:
            print(f"Atom {atom} is outside the sub-cell range.")
            exit(1)  # Exit if an atom is outside the sub-cell range
        
        # Increment the count for this sub-cell
        if key not in subCellCounts:
            subCellCounts[key] = 0
        subCellCounts[key] += 1
        
    subCellCounts = sorted(subCellCounts.items())
        
    return subCellCounts    

def verifyCounts(counts):
    """
    Verify that the counts of atoms in each sub-cell are correct.
    
    Args:
        counts (list): List of tuples containing sub-cell keys and their respective counts.
    
    Returns:
        bool: True if all counts are correct, False otherwise.
    """
    total_atoms = sum(count for _, count in counts)
    print(f"Total atoms counted: {total_atoms}")
    
# TODO make this save to the correct subdir
def writeOutput(outStr, outFile="subCountsCheck.txt"):
    """
    Write the output string to a file.
    
    Args:
        outStr (str): The output string to write.
    """
    with open(outFile, "w") as f:
        f.write(outStr)
    print("Output written to ", outFile)
    
def processFile(file_path):
    """
    Load XYZ data from a file.
    
    Args:
        file_path (str): Path to the XYZ file.
    
    Returns:
        list: List of tuples containing XYZ coordinates.
    """
    xyzData = []
    counts = {}
    fsStr = ""
    outStr = ""
    
    with open(file_path, 'r') as file:
        for line in file:
            # print(line)
            parts = line.strip().split()
            if len(parts) == 1:
                if xyzData:
                    # Process the previous fsStr and reset the dictionary
                    counts = getSubCellCounts(xyzData, unitCell)
                    # print(next(file))
                    fsStr = next(file).strip().split()[1]
                    outStr += "Time step: " + fsStr + " fs\n"
                    # verifyCounts(counts)
                    outStr = printCounts(counts, outStr)
                    xyzData = []
                    counts = {}
                    fsStr = ""
            if len(parts) == 4:
                x, y, z = map(float, parts[1:])
                x = correctCoordinate(x, unitCell)
                y = correctCoordinate(y, unitCell)
                z = correctCoordinate(z, unitCell)
                xyzData.append((x, y, z))
            # line that has some meta data but skips the number of atoms
            elif len(parts) >= 5:
                unitCell = float(line.strip().split()[-1][:-2])
                
            else:
                continue
                # print(f"Skipping invalid line: {line.strip()}")
                
    writeOutput(outStr)

def printCounts(counts, outStr):
    """
    Print the counts of atoms in each sub-cell.
    
    Args:
        counts (list): List of tuples containing sub-cell keys and their respective counts.
    """
    outStr += "Sub-cell counts:\n"
    for key, count in counts:
        outStr += "Sub-cell " + key + ": " + str(count) + " atoms\n"
    
    return outStr

def main():
    xyz = sys.argv[1]
    processFile(xyz)
    # atoms = correctXyz(atoms, unitCell)
    # counts = getSubCellCounts(atoms, unitCell)
    
    # printCounts(counts)

if __name__ == "__main__":
    # Entry point for the script
    print("This script compares XYZ ratios.")
    # Add your main logic here
    main()