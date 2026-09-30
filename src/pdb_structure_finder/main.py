import json
import sys
from pathlib import Path

import httpx
from db_requests import read_proteinDataBank, read_uniprot
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

    # ================================================
    # PARSE COMMAND LINE ARGUMENT
    # ================================================
    args = parse_commandLine_args()
    PDB_DIR = manage_pdb_file(store_pdb=args.store_pdb, pdb_dir=args.pdb_dir)
    if PDB_DIR:
        Path.mkdir(PDB_DIR)
        print(f"-> {PDB_DIR} created")

    # =================================================
    # READ JSON
    # =================================================
    try:
        json_content = json_data(json_file_path=args.json_file)
    except FileNotFoundError:
        print(f"ERROR: JSON file not found: {args.json_file}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"ERROR: Invalid JSON format: {args.json_file}")
        sys.exit(1)

    # =================================================
    # UNIPROT JSON INFO
    # =================================================
    try:
        structures = read_uniprot(json_content)
    except httpx.HTTPStatusError as error:
        print(f"HTTP error: {error}")
    except httpx.RequestError as error:
        print(f"Request errro: {error}")

    # =================================================
    # PROTEIN DATA BANK
    # =================================================
    try:
        pippo = read_proteinDataBank(
            df_structure=structures, store_pdb=args.store_pdb, pdb_dir=PDB_DIR
        )
    except httpx.HTTPStatusError as error:
        print(f"HTTP error: {error}")
    except httpx.RequestError as error:
        print(f"Request errro: {error}")


if "__main__" == __name__:
    main()
