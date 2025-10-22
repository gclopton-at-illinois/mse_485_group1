import numpy as np
import sys, json, argparse, os, itertools
import ase.io

seed = 44592

def main():

    # Parsing Arguments
    parser = argparse.ArgumentParser(description = "Dimensions of the CeO_2 supercell and the percentage of Oxygens (d/2 in O_{2-d}) to remove")
    parser.add_argument("x_count", type = int, help="Number of unit cells in x-direction")
    parser.add_argument("y_count", type = int, help="Number of unit cells in y-direction")
    parser.add_argument("z_count", type = int, help="Number of unit cells in z-direction")
    parser.add_argument("p_O2_remove", type = float, help="Percentage of Oxygen atoms to remove (ratio of oxygen vacancies)")

    args = parser.parse_args()

    # Checking that arguments are within our bounds
    assert args.x_count > 0, "Number of unit cells in x-direction must be greater than 0" 
    assert args.y_count > 0, "Number of unit cells in y-direction must be greater than 0"
    assert args.z_count > 0, "Number of unit cells in z-direction must be greater than 0"
    assert (args.p_O2_remove >= 0 and args.p_O2_remove <= 1), "Percentage of Oxygen atoms removed must be between 0 and 1 (inclusive)"

    script_dir = os.path.dirname(__file__)
    relative_json_path = os.path.join('..', 'inputs', 'structure', 
                                      f"data.ceo2_{args.x_count}x{args.y_count}x{args.z_count}nm.json")

    full_json_path = os.path.join(script_dir, relative_json_path)

    # Check the associated JSON file for its parameter values
    with open(full_json_path, "r") as f:
        supercell_params = json.load(f)
    
    O_type_num = supercell_params["type_mapping"]["O"]
    init_num_O2 = supercell_params["composition_counts"]["O"]
    num_Ce = supercell_params["composition_counts"]["Ce"]

    # generate atom ids to remove
    rng = np.random.default_rng(seed = seed)
    O2_remove = rng.choice(range(num_Ce + 1, num_Ce + init_num_O2 + 1), 
                           size = int(init_num_O2 * args.p_O2_remove), replace = False)

    # remove those atoms from the lmp file
    relative_lmp_path = os.path.join('..', 'inputs', 'structure', f"data.ceo2_{args.x_count}x{args.y_count}x{args.z_count}nm.lmp")
    relative_lmp_temp = os.path.join('..', 'inputs', 'structure', 
                                     f"data.ceo2_{args.x_count}x{args.y_count}x{args.z_count}nm_pO2_{args.p_O2_remove}.lmp")
    relative_lmp_json = os.path.join('..', 'inputs', 'structure', 
                                     f"data.ceo2_{args.x_count}x{args.y_count}x{args.z_count}nm_pO2_{args.p_O2_remove}.json")


    full_lmp_path = os.path.join(script_dir, relative_lmp_path)
    full_lmp_temp = os.path.join(script_dir, relative_lmp_temp)
    full_json_temp = os.path.join(script_dir, relative_lmp_json)
    
    ##### FILE UPDATING -- REMOVING O ATOMS #####

    lammps_atom_style = "charge"

    # 1. Read the .lmp file into an ASE 'Atoms' object. The 'style' argument is essential for LAMMPS data files.
    try:
        atoms = ase.io.read(full_lmp_path, format = 'lammps-data', atom_style = lammps_atom_style)
    
    except Exception as e:
        print(f"\n--- ERROR ---")
        print(f"ASE failed to read the file. This is often due to a "
            f"wrong 'lammps_atom_style'.")
        print(f"You set the style to '{lammps_atom_style}'. Is that correct?")
        print(f"Details: {e}")
        sys.exit()

    print(f"Original atom count: {len(atoms)}")

    # 2. Convert 1-based LAMMPS IDs to 0-based ASE indices (LAMMPS ID 1 is ASE index 0, LAMMPS ID 2 is ASE index 1, etc.)
    indices_to_delete = [i - 1 for i in O2_remove]

    # 3. Sort the indices in REVERSE order for ASE to use properly
    indices_to_delete.sort(reverse=True)

    print(f"Removing {len(indices_to_delete)} atoms...")

    # 4. Delete the atoms one by one from the 'atoms' object
    for index in indices_to_delete:
        try:
            del atoms[index]
        except IndexError:
            print(f"Warning: Could not remove atom with ID {index + 1} "
                f"(index {index}). It might have already been "
                f"removed or was out of bounds.")

    print(f"New atom count: {len(atoms)}")

    # 5. Write the modified 'atoms' object back to a new .lmp file
    #    ASE automatically re-numbers all Atom IDs from 1 to N,
    #    updates the header counts, and re-maps all Bond/Angle/etc. lists.
    ase.io.write(full_lmp_temp, atoms, format='lammps-data', atom_style = lammps_atom_style)

    print(f"\nSuccessfully wrote new file to {full_lmp_temp}")
    
    supercell_params["n_atoms"] = len(atoms)
    supercell_params["composition_counts"]["O"] -= len(O2_remove)

    with open(full_json_temp, "w") as f:
        json.dump(supercell_params, f, indent = 2)

if __name__ == "__main__" :
    main()