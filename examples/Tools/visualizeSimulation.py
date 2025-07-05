import os
import sys
import cv2
import numpy as np
from PIL import Image

# Dubug why there seems to be no change between 90 fs and 200 fs on dir output_2025_6_28_23_0001-01-01 00:00:00-1000.0-873.15-25

# Extract the number after 'reaction_traj_' and sort the files
def extract_number(file_path):
    base_name = os.path.basename(file_path)
    try:
        # Extract the number after 'reaction_traj_' and before the next '_'
        start_idx = base_name.index('reaction_traj_') + len('reaction_traj_')
        end_idx = base_name.index('.', start_idx)
        number = int(base_name[start_idx:end_idx])
        # print(f"Extracted number: {number} from file: {base_name}")
        return number
    except (ValueError, IndexError):
        return float(0)  # Place files without a valid number at the beginning
    
def find_xyz_files(directory):
    """
    Walk through the given directory and its subdirectories to find .xyz files
    that contain the substring 'reaction_traj' and order them in ascending order
    based on the number after 'reaction_traj_'.

    Additionally, extract all unique atom types from the .xyz files. The atom type
    is the first character (or string) in each line of atom data and assign a unique
    color to each atom type.

    Args:
        directory (str): The root directory to start the search.

    Returns:
        tuple: A sorted list of paths to .xyz files, a set of unique atom types,
               and a dictionary mapping atom types to unique colors.
    """
    xyz_files = []
    atom_types = set()

    for root, _, files in os.walk(directory):
        for file in files:
            print(f"Checking file: {file} in {root}")
            if file.endswith('.xyz') and 'reaction_traj' in file:
                print(f"Found .xyz file: {file} in {root}")
                file_path = os.path.join(root, file)
                xyz_files.append(file_path)

                # Extract atom types from the file
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    num_atoms = int(lines[0].strip())
                    atom_data = lines[2:2 + num_atoms]
                    for atom in atom_data:
                        parts = atom.split()
                        if parts:
                            atom_types.add(parts[0])  # First part is the atom type

    # Assign unique colors to each atom type
    atom_types = sorted(atom_types)  # Sort for consistent color assignment
    color_map = {}
    np.random.seed(42)  # Seed for reproducibility
    for atom_type in atom_types:
        color_map[atom_type] = tuple(np.random.randint(0, 256, size=3).tolist())

    '''
    def extract_number(file_path):
        base_name = os.path.basename(file_path)
        try:
            # Extract the number after 'reaction_traj_'
            start_idx = base_name.index('reaction_traj_') + len('reaction_traj_')
            number = int(''.join(filter(str.isdigit, base_name[start_idx:])))
            return number
        except (ValueError, IndexError):
            return float('inf')  # Place files without a valid number at the end
    '''
    
    
    # Include the .xyz file with 'md.xyz' in its name at the beginning of the list
    for root, _, files in os.walk(directory):
        for file in files:
            if 'md.xyz' in file:
                file_path = os.path.join(root, file)
                xyz_files.insert(0, file_path)
                break
    
    xyz_files.sort(key=extract_number)
    
    return xyz_files, atom_types, color_map

if __name__ == "__main__":
    directory_to_search = sys.argv[1]  # Replace with your directory path
    # print(f"Searching for .xyz files in directory: {directory_to_search}")
    xyz_files, atom_types, color_map = find_xyz_files(directory_to_search)
    # print("Found .xyz files:")
    
    # for file in xyz_files:
    #     print(file)

    def visualize_xyz(file_path, color_map):
        """
        Visualize the atoms in the .xyz file using OpenCV, centering around the centroid
        and scaling the image to fit all atoms based on the cell size.

        Args:
            file_path (str): Path to the .xyz file.
            color_map (dict): Dictionary mapping atom types to unique colors.
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()

            # Parse the .xyz file
            num_atoms = int(lines[0].strip())
            atom_data = lines[2:2 + num_atoms]

            # Extract atom coordinates and types
            coordinates = []
            atom_types_in_file = []
            for atom in atom_data:
                parts = atom.split()
                if len(parts) < 4:
                    continue
                atom_type = parts[0]
                x, y, z = map(float, parts[1:4])
                coordinates.append((x, y, z))
                atom_types_in_file.append(atom_type)

        # Calculate centroid
        coordinates = np.array(coordinates)
        centroid = np.mean(coordinates, axis=0)

        # Determine the bounding box
        min_coords = np.min(coordinates, axis=0)
        max_coords = np.max(coordinates, axis=0)
        cell_size = max(max_coords - min_coords)

        # Create a blank image with padding
        img_size = 500
        circle_radius = max(1, int(10 / num_atoms))  # Ensure a minimum radius of 1
        scale = img_size / cell_size
        padding = 50
        img = np.zeros((img_size + 2 * padding, img_size + 2 * padding, 3), dtype=np.uint8)

        # Draw atoms centered around the centroid
        for (x, y, z), atom_type in zip(coordinates, atom_types_in_file):
            # Normalize coordinates to fit in the image
            x = int((x - centroid[0]) * scale + img_size / 2 + padding)
            y = int((y - centroid[1]) * scale + img_size / 2 + padding)

            # Get the color for the atom type
            color = color_map.get(atom_type, (255, 255, 255))  # Default to white if not found

            # Draw the atom as a circle
            cv2.circle(img, (x, y), circle_radius, color, -1)

        # Add a key to the image
        key_x, key_y = 10, 10
        for atom_type, color in color_map.items():
            cv2.circle(img, (key_x, key_y), circle_radius, color, -1)
            cv2.putText(img, atom_type, (key_x + 20, key_y + 5), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1)
            key_y += 30

        return img

    # Create a list to store frames for the GIF
    frames = []

    # Generate frames for each .xyz file
    # Change this so the putText includes the bonding info to some degree
    for file in xyz_files:
        # print(f"Visualizing {file}")
        img = visualize_xyz(file, color_map)
        
        # Extract the timestamp from the file name
        base_name = os.path.basename(file)
        number = extract_number(file)

        # Add the timestamp as text to the image
        cv2.putText(img, f"Time: {number} fs", (10, img.shape[0] - 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1)

        # Convert BGR to RGB for GIF and append to frames
        frames.append(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))

    # Save the frames as a GIF
    gif_path = os.path.join(directory_to_search, "visualization.gif")
    Image.fromarray(frames[0]).save(
    gif_path, save_all=True, append_images=[Image.fromarray(frame) for frame in frames[1:]], duration=500, loop=0
    )
    print(f"GIF saved at {gif_path}")