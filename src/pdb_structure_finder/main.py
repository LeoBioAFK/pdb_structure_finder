from pathlib import Path

from read_json import json_data
from utility import parse_commandLine_args


def manage_pdb_file(store_pdb, pdb_dir):
    if store_pdb:
        return assign_pdb_folder(pdb_dir=pdb_dir)
    else:
        print("-[NO SAVE MODE]- pdb/ciff files won't be saved\n")
        print("Run with --save_pdb for saving them instead!")
        return None


def assign_pdb_folder(pdb_dir):
    return pdb_dir


def main():
    args = parse_commandLine_args()
    PDB_DIR = manage_pdb_file(store_pdb=args.store_pdb, pdb_dir=args.pdb_dir)
    if PDB_DIR:
        Path.mkdir(PDB_DIR)
        print(f"-> {PDB_DIR} created")

    # json_file
    json_instruction = json_data(json_file_path=args.json_file)
    print(json_instruction)


if "__main__" == __name__:
    main()
