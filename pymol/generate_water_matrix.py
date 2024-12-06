def generate_water_matrix(output_file, rows, cols, spacing=3.0):
    """
    Generate a water molecule matrix and save to a PDB file.
    """
    with open(output_file, 'w') as f:
        atom_id = 1
        res_id = 1
        for i in range(rows):
            for j in range(cols):
                x = i * spacing
                y = j * spacing
                z = 0.0

                f.write(f"ATOM  {atom_id:5d}  O   HOH A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                atom_id += 1

                f.write(f"ATOM  {atom_id:5d}  H1  HOH A{res_id:4d}    {x + 0.957:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                atom_id += 1

                f.write(f"ATOM  {atom_id:5d}  H2  HOH A{res_id:4d}    {x:8.3f}{y + 0.957:8.3f}{z:8.3f}  1.00  0.00           C\n")
                atom_id += 1

                res_id += 1

generate_water_matrix("water_matrix_20x20.pdb", rows=20, cols=20, spacing=3.0)
print("PDB file generated: water_matrix_20x20.pdb")
