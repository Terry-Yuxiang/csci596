import numpy as np
from pymol import cmd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Read temperature data
def read_all_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)
    times = data[:, 0]
    layers = data[:, 1:]
    return times, layers

# Map temperatures to colors
def temperature_to_color(temp_values, vmin, vmax, cmap="coolwarm"):
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    colormap = plt.get_cmap(cmap)
    return [colormap(norm(temp))[:3] for temp in temp_values]

# Randomly move water molecules
def apply_random_displacement_whole_molecule(matrix_size=20, max_displacement=0.1):
    num_molecules = matrix_size * matrix_size
    for i in range(num_molecules):
        atom_indices = [3 * i + 1, 3 * i + 2, 3 * i + 3]
        displacement = np.random.uniform(-max_displacement, max_displacement, size=3)
        cmd.translate(list(displacement), selection=f"id {atom_indices[0]}-{atom_indices[2]}")

# Assign colors based on temperature
def assign_colors_to_pymol(layer_temps, vmin, vmax, matrix_size=20):
    rows_per_layer = matrix_size // len(layer_temps)
    index = 1
    colors = temperature_to_color(layer_temps, vmin, vmax)
    for layer_idx, color in enumerate(colors):
        start_row = layer_idx * rows_per_layer
        end_row = start_row + rows_per_layer
        for row in range(start_row, end_row):
            for col in range(matrix_size):
                cmd.set_color(f"color_{index}", color)
                cmd.color(f"color_{index}", f"id {index}")
                cmd.color(f"color_{index}", f"id {index + 1}")
                cmd.color(f"color_{index}", f"id {index + 2}")
                index += 3

# Generate images for each frame
def generate_images(file_path, output_dir="frames", matrix_size=20, max_frames=500, vmin=0, vmax=8, cmap="coolwarm"):
    times, layers = read_all_data(file_path)
    cmd.load("water_matrix_20x20.pdb", "water_matrix")

    cmd.show("spheres", "water_matrix")
    cmd.set("sphere_scale", 0.6, "water_matrix")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for frame_idx, time in enumerate(times[:max_frames]):
        print(f"Generating frame {frame_idx + 1}/{max_frames}, time: {time:.3f}")
        apply_random_displacement_whole_molecule(matrix_size, max_displacement=0.1)
        assign_colors_to_pymol(layers[frame_idx], vmin, vmax, matrix_size)
        image_path = os.path.join(output_dir, f"frame_{frame_idx:04d}.png")
        cmd.png(image_path, width=1920, height=1080, dpi=300)
        cmd.refresh()

file_path = "layers_res.txt"
generate_images(file_path, max_frames=500, vmin=0, vmax=8, cmap="coolwarm")
