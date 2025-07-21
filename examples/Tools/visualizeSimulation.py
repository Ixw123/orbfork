import os
import ffmpeg
import sys
import cv2
import numpy as np
from PIL import Image
import colorsys
import subprocess
import json


# Dubug why there seems to be no change between 90 fs and 200 fs on dir output_2025_6_28_23_0001-01-01 00:00:00-1000.0-873.15-25

# TODO make different colored circles around various bonds based on a threshold distance if provided one
# TODO make sure that this properly parses the xyz files with reaction_traj_ in the name since right now it just looks at the first entry and is missing things that the reaction_traj_ files have picks up outside the steps 
# TODO make color map better contrast
# TODO make the gif determine the frames based on fs steps 
# TODO remove the reaction_traj_ files from the gif and only keep the every 10 fs xyz file 
# TODO make gif mpeg

def correctCoordinate(x, y, z, cellSize):
    """
    Correct the coordinates to fit within the unit cell.
    
    Args:
        x, y, z (float): Coordinates of the atom.
        cellSize (float): Size of the unit cell.
    
    Returns:
        tuple: Corrected coordinates (x, y, z).
    """
    x = cellSize if float(x % cellSize) == 0.0 else float(x % cellSize)
    y = cellSize if float(y % cellSize) == 0.0 else float(y % cellSize)
    z = cellSize if float(z % cellSize) == 0.0 else float(z % cellSize)
    
    return x, y, z



# Extract the number after 'frame_' and sort the files
def extract_number(file_path):
    base_name = os.path.basename(file_path)
    # print(base_name)
    try:
        # Extract the number after 'reaction_traj_' and before the next '_'
        start_idx = base_name.index('frame_') + len('frame_')
        # print(f"Start index: {start_idx} in file: {base_name}")
        numberStr = ''
        for c in base_name[start_idx:]:
            if c.isdigit():
                numberStr += c
            else:
                break
        # end_idx = base_name.index('.', start_idx) - 2
        # print(f"End index: {end_idx} in file: {base_name}")
        number = int(numberStr)
        # print(f"Extracted number: {number} from file: {base_name}")
        # print(f"Extracted number: {number} from file: {base_name}")
        return number
    except (ValueError, IndexError):
        print(f"Failed to extract number from file: {base_name}")
        return float(0)  # Place files without a valid number at the beginning
    
"""
def findBonds(atoms, threshold):
    
    # Find bonds between atoms based on a distance threshold.
    
    # Args:
    #     atoms (ase.Atoms): Atoms object containing atomic positions.
    #     threshold (float): Distance threshold for bond detection.
    
    # Returns:
    #     list: List of tuples representing bonded atom pairs.
    
    from ase.neighborlist import NeighborList
    nl = NeighborList([threshold] * len(atoms), self_interaction=False)
    nl.update(atoms)
    bonds = {}
    
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            if i < j:  # Avoid double counting
                key = str(i) + ',' + str(j)
                bonds[key] = 
    
    return bonds

def drawBonds(atoms, bonds, color_map, img_size=500, padding=50):
    
    # Draw bonds between atoms on a blank image.
    
    # Args:
    #     atoms (ase.Atoms): Atoms object containing atomic positions.
    #     bonds (list): List of tuples representing bonded atom pairs.
    #     color_map (dict): Dictionary mapping atom types to unique colors.
    #     img_size (int): Size of the output image.
    #     padding (int): Padding around the image.
    
    # Returns:
    #     np.ndarray: Image with drawn bonds.
    
    coordinates = atoms.get_positions()
    atom_types = atoms.get_chemical_symbols()
    
    # Create a blank image
    img = np.zeros((img_size + 2 * padding, img_size + 2 * padding, 3), dtype=np.uint8)
    
    # Draw atoms as circles
    for i, (x, y, z) in enumerate(coordinates):
        x = int(x + img_size / 2 + padding)
        y = int(y + img_size / 2 + padding)
        color = color_map.get(atom_types[i], (255, 255, 255))  # Default to white if not found
        cv2.circle(img, (x, y), 5, color, -1)  # Draw atom as a circle
    
    # Draw bonds as lines
    for i, j in bonds:
        x1, y1 = coordinates[i][:2]
        x2, y2 = coordinates[j][:2]
        x1 = int(x1 + img_size / 2 + padding)
        y1 = int(y1 + img_size / 2 + padding)
        x2 = int(x2 + img_size / 2 + padding)
        y2 = int(y2 + img_size / 2 + padding)
        
        cv2.line(img, (x1, y1), (x2, y2), (200, 200, 200), 1)  # Draw bond as a line
    
    return img
"""
    
def distinct_colors(n):
    hsv_colors = [(i / n, 0.8, 0.9) for i in range(n)]  # hue-spaced
    rgb_colors = [colorsys.hsv_to_rgb(*hsv) for hsv in hsv_colors]
    return [tuple(int(255 * c) for c in rgb) for rgb in rgb_colors]
    
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
            # print(f"Checking file: {file} in {root}")
            if file.endswith('.xyz') and 'frame' in file:
                # print(f"Found .xyz file: {file} in {root}")
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
    colors = distinct_colors(len(atom_types))
    for atom_type, color in zip(atom_types, colors):
        color_map[atom_type] = color
    
    # color_map = {}
    # np.random.seed(42)  # Seed for reproducibility
    # for atom_type in atom_types:
    #     color_map[atom_type] = tuple(np.random.randint(0, 256, size=3).tolist())

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
    # for root, _, files in os.walk(directory):
    #     for file in files:
    #         if 'md.xyz' in file:
    #             file_path = os.path.join(root, file)
    #             xyz_files.insert(0, file_path)
    #             break
    
    xyz_files.sort(key=extract_number)
    
    return xyz_files, atom_types, color_map

def makeVideo(xyz_files, color_map, video_path):
    """
    Create a video from the list of .xyz files using OpenCV and ffmpeg.

    Args:
        xyz_files (list): List of paths to .xyz files.
        color_map (dict): Dictionary mapping atom types to unique colors.
        video_path (str): Path where the video will be saved.
    """
    if not xyz_files:
        print("No .xyz files found.")
        return
    
    first_img = visualize_xyz(xyz_files[0], color_map)
    first_img = cv2.cvtColor(first_img, cv2.COLOR_BGR2RGB)
    height, width, _ = first_img.shape
    
    frameCnt = len(xyz_files)
    frameClipCnt = max(frameCnt/1000, 1)
    maxFps = 30
    seconds = max(60 * frameClipCnt, 60)# make every animation a minute
    # fps = 1 # min(1 / min_frame_duration, frameCnt / seconds)
    fps = round(min(frameCnt/seconds, maxFps), 2)
    
    print(f"Total frames: {frameCnt}, FPS: {fps}, Seconds: {seconds}, Frame Clip Count: {frameClipCnt}")


    process = (
        ffmpeg
        .input('pipe:', format='rawvideo', pix_fmt='rgb24', s='{}x{}'.format(width, height), framerate=fps)
        .output(video_path, pix_fmt='yuv420p', vcodec='libx264', r=fps)
        .overwrite_output()
        .run_async(pipe_stdin=True)
    )

    for file in xyz_files:
        img = visualize_xyz(file, color_map)
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
        process.stdin.write(img.astype(np.uint8).tobytes())
            
    process.stdin.close()
    process.wait()

    
    print(f"MPEG video saved at {video_path}")

def get_video_duration(video_path):
    """Get video duration in seconds using ffprobe"""
    result = subprocess.run([
        "ffprobe", "-v", "error", "-show_entries",
        "format=duration", "-of", "json", video_path
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    duration = json.loads(result.stdout)["format"]["duration"]
    return float(duration)

def compress_video_to_target_size(input_path, output_path, target_size_mb):
    duration = get_video_duration(input_path)

    # Convert target size to kilobits (MB → bytes → bits → kilobits)
    target_kbits = target_size_mb * 1024 * 8
    bitrate_kbps = target_kbits / duration

    print(f"Compressing to ~{bitrate_kbps:.2f} kbps for {duration:.2f} seconds to fit {target_size_mb} MB")

    # Run ffmpeg compression
    subprocess.run([
        "ffmpeg", "-i", input_path,
        "-b:v", f"{int(bitrate_kbps)}k",
        "-bufsize", f"{int(bitrate_kbps*2)}k",  # smoother rate control
        "-maxrate", f"{int(bitrate_kbps)}k",
        "-vcodec", "libx264", "-preset", "slow",  # slower = better compression
        "-y", output_path
    ])

    

if __name__ == "__main__":
    
    imgSize = 500  # Size of the output image
    
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
            # print(lines[1].strip().split()[-1][:-2])
            cellSize = float(lines[1].strip().split()[-1][:-2])  # Assuming the second line contains the cell size
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
                x, y, z = correctCoordinate(x, y, z, cellSize)
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

    # first_img = visualize_xyz(xyz_files[0], color_map)
    # first_img = cv2.cvtColor(first_img, cv2.COLOR_BGR2RGB)
    # height, width, _ = first_img.shape
    
    video_path = os.path.join(directory_to_search, "visualization.mp4")
    
    makeVideo(xyz_files, color_map, video_path)
    
    # frameCnt = len(xyz_files)
    # frameClipCnt = max(frameCnt/1000, 1)
    # maxFps = 30
    # seconds = max(60 * frameClipCnt, 60)# make every animation a minute
    # # fps = 1 # min(1 / min_frame_duration, frameCnt / seconds)
    # fps = round(min(frameCnt/seconds, maxFps), 2)
    
    # print(f"Total frames: {frameCnt}, FPS: {fps}, Seconds: {seconds}, Frame Clip Count: {frameClipCnt}")

    # process = (
    #     ffmpeg
    #     .input('pipe:', format='rawvideo', pix_fmt='rgb24', s='{}x{}'.format(width, height), framerate=fps)
    #     .output(video_path, pix_fmt='yuv420p', vcodec='libx264', r=fps)
    #     .overwrite_output()
    #     .run_async(pipe_stdin=True)
    # )

    # for file in xyz_files:
    #     img = visualize_xyz(file, color_map)
    #     img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    #     process.stdin.write(img.astype(np.uint8).tobytes())

    # process.stdin.close()
    # process.wait()

    
    # print(f"MPEG video saved at {video_path}")
    
    videoSize = os.path.getsize(video_path)/(1024 * 1024)  # Size in MB
    
    
    # if it is larger than 25 MB for email purposes make a condensed mp4
    if videoSize > 25:
        
        newVideoPath = os.path.join(directory_to_search, "visualization_compressed.mp4")
    
        targetSize = 25 # Size in MB
        
        compress_video_to_target_size(video_path, newVideoPath, targetSize)

    # Save the frames as a GIF
    # rename the gif with the directory info if there is directory info
    # gif_path = os.path.join(directory_to_search, "visualization.gif")
    # Image.fromarray(frames[0]).save(
    # gif_path, save_all=True, append_images=[Image.fromarray(frame) for frame in frames[1:]], duration=500, loop=0
    # )
    # print(f"GIF saved at {gif_path}")