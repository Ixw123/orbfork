from ase.io import read, write
from ase.visualize import view

from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

from io import StringIO
from PIL import Image, ImageChops
import os
import glob
import argparse

def get_most_recent_directory(base_path):
    # Get all entries in base_path that are directories
    dirs = [d for d in glob.glob(os.path.join(base_path, '*')) if os.path.isdir(d) and "output" in d]
    if not dirs:
        return None
    most_recent = max(dirs, key=os.path.getmtime)
    return most_recent

def save_atoms_png(atoms, filename, radii=0.5):
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    plot_atoms(atoms, ax, radii=radii, rotation=("90x,90y,90z"))
    ax.set_axis_off()
    fig.tight_layout(pad=0) 
    fig.savefig(filename, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close(fig)

def crop_image(png_file):
    img = Image.open(png_file)
    
    if img.mode != "RGB":
        print(f"Converting {png_file} from {img.mode} to RGB")
        img = img.convert("RGB")
    
    bg = Image.new(img.mode, img.size, img.getpixel((0, 0)))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox:
        img = img.crop(bbox)
        # Force safe save
    try:
        img.save(png_file, format="PNG")
    except Exception as e:
        print(f"Error saving image {png_file}: {e}")
        raise

# Function to extract a specific frame from the .xyz file
def extract_frame(lines, frame_number):
    # Parse the file
    num_atoms = int(lines[0].strip())  # Number of atoms from the first line
    # print(f"Number of frames: {len(lines) // (num_atoms + 2)}")
    frame_start = frame_number * (num_atoms + 2)  # Start of the desired frame
    frame_lines = lines[frame_start : frame_start + num_atoms + 2]  # Extract frame
    return "".join(frame_lines)


# Function to visualize the trajectory
def visualize_trajectory(xyz_file, reaction_name):
    # Extract the specific frame from the .xyz file

    with open(xyz_file, "r") as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip())
    num_frames = len(lines) // (num_atoms + 2)

    # Read all files in and use largest bounding box or fit to unit cell

    png_files = []
    for i in range(num_frames):
        frame_data = extract_frame(lines, i)
        file_wrapper = StringIO(frame_data)

        atoms = read(file_wrapper, format="xyz")

        # write("test.png", atoms, format="png")
        # Image.open("test.png").show()
        # print("Mode:", Image.open("test.png").mode)

        # del atoms[
        #     [atom.index for atom in atoms if atom.index not in [327, 1154, 546, 547]]
        # ]


        if not os.path.exists("reaction_images"):
            os.makedirs("reaction_images")

        # Create an ASE view
        f_name = f"reaction_images/{reaction_name}_{i:03}.png"

        # atoms.center(vacuum=5.0)  # Add vacuum padding
        # fix output based on arbitrary size of system
        # scale = 1.0 / max(atoms.get_cell().lengths())
        # atoms.set_positions(atoms.get_positions() * scale)
        # atoms.set_cell(atoms.get_cell() * scale)
        # atoms.wrap()
        # full jpeg
        # first_frame = frames[0]
        write(f_name, atoms, format="png")
        # save_atoms_png(atoms, f_name)
        # crop_image(f_name)
        png_files.append(f_name)

        # write('frame0.png', atoms, format='png')
        # Image.open(f_name).show()

    # Print the first atoms info 
    file_wrapper.seek(0)
    atoms = read(file_wrapper, format="xyz", index=0)
    write('first_frame.xyz', atoms, format='xyz')
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    for i, pos in enumerate(positions):
        print(f"Atom {i}: x={pos[0]:.6f}, y={pos[1]:.6f}, z={pos[2]:.6f}")
        print(f"Symbol {symbols[i]}")

    # print(atoms.get_chemical_symbols())
    print(atoms.get_cell())

    return png_files


def create_gif(png_files, gif_name="output.gif", duration=200):
    # print(images)
    print(gif_name)
    """
    Create a GIF from a list of PNG files.

    Args:
        png_files (list): List of PNG file paths.
        gif_name (str): Name of the output GIF file.
        duration (int, optional): Duration of each frame in milliseconds. Defaults to 200.
    """
    images = []
    for png_file in png_files:
        img = Image.open(png_file).convert("RGB")
        images.append(img)

    if images:
        images[0].save(
            gif_name,
            save_all=True,
            append_images=images[1:],
            duration=duration,
            loop=0,
            disposal=2,
        )
        print(f"GIF created: {gif_name}")
    else:
        print("No PNG images found in the folder.")


# examples/output/.xyz
def main(reaction_name):
    # Path to the .xyz file
    basePath = os.getcwd()
    print(basePath)

    # Get most recently modified dir hack
    latestRun = get_most_recent_directory(basePath)
    print(latestRun)

    xyz_file = os.path.join(latestRun, reaction_name + ".xyz")

    print(xyz_file)
    # xyz_file = f"output/{reaction_name}.xyz"

    # Visualize the trajectory
    png_files = visualize_trajectory(xyz_file, reaction_name)
    # Create a GIF from the PNG files
    create_gif(png_files, gif_name=latestRun + "/" + reaction_name + ".gif", duration=200)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize trajectory and create GIF.")
    parser.add_argument(
        "reaction_name", type=str, help="Name of the reaction to visualize."
    )
    args = parser.parse_args()
    reaction_name = args.reaction_name

    # assert os.path.exists(
    #     f"output/{reaction_name}.xyz"
    # ), f"File output/{reaction_name}.xyz does not exist."
    main(reaction_name)
