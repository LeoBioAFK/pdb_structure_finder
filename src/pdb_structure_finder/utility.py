import argparse


def parse_commandLine_args():

    parser = argparse.ArgumentParser(
        prog="PDB_Hunter",
        description="""Search for protein structure that satisfy your requests |
        Database -> Protein Data Bank (PDB) |
        You can even request some specific ligand/metal encoded by CCD code""",
        epilog="""Input protein list + parameter must be a json file |
        For help run 'interactive_json.py' for an interactive path for create |
        your own json of input""",
    )
    parser.add_argument(
        "json_file",
        help="Path to json file of input containing the list of uniprot ids "
        "and optional parameters",
    )
    parser.add_argument(
        "--store_pdb",
        default=False,
        action="store_true",
        help="Bool value (default is false) for saving pdbs when found",
    )
    parser.add_argument(
        "--pdb_dir",
        default="pdb_dir",
        help="explicit the folder where you want to save the pdbs (when found)",
    )
    return parser.parse_args()
