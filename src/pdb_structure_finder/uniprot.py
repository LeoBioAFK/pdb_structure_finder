from pathlib import Path

import gemmi
import httpx
import pandas as pd

url = "https://rest.uniprot.org/uniprotkb/P02768.json"
response = httpx.get(url)
print(response)
if response.status_code == 200:
    content = response.json()

# condensed requests
PDB_infos = [
    elem for elem in content["uniProtKBCrossReferences"] if elem["database"] == "PDB"
]
print("one row database scraped")
print(PDB_infos[0]["properties"])

report = []

for info in PDB_infos:
    uniprot_id = "P00450"
    database_id = info["id"]
    method = next(x["value"] for x in info["properties"] if x["key"] == "Method")
    resolution = next(
        x["value"] for x in info["properties"] if x["key"] == "Resolution"
    )
    chain = next(x["value"] for x in info["properties"] if x["key"] == "Chains")

    report.append(
        {
            "uniprot_id": uniprot_id,
            "database": "PDB",
            "database_id": database_id,
            "method": method,
            "resolution": resolution,
            "chain": chain,
        }
    )

# NOTE: dataframe containing all entries from Structure field (Uniprot)
df = pd.DataFrame(report)
print(df)

# Try to use PDB API for getting the sctructures
# Still try with P00450
pdb_code = df["database_id"].iloc[0]
url_pdb = f"https://files.rcsb.org/download/{pdb_code}.cif"
response_pdb = httpx.get(url_pdb)
print(response_pdb)
if response_pdb.status_code == 200:
    content_pdb = response_pdb.text

# save cif file
DEST_FOLDER = Path("../../tests")
with open(f"{DEST_FOLDER}/{pdb_code}.cif", "w") as k:
    k.write(content_pdb)

# check non-polymer entity
# Using GraphQL
import requests

url_graph = "https://data.rcsb.org/graphql"
query = """
    query GetEntry($pdb_id: String!) {
        entry(entry_id: $pdb_id) {
            rcsb_id
            nonpolymer_entities {
                rcsb_id
                rcsb_nonpolymer_entity_container_identifiers {
                    entry_id
                    entity_id
                    auth_asym_ids
                    asym_ids
                    nonpolymer_comp_id
                }
                nonpolymer_comp {
                    chem_comp {
                        id
                        formula_weight
                        name
                        formula
                    }
                }
            }
        } 
    }
"""
variables = {"pdb_id": pdb_code}
response_graph = requests.post(url_graph, json={"query": query, "variables": variables})

data = response_graph.json()
print(data)
entry = data["data"]["entry"]
nonpolymer_entities = entry["nonpolymer_entities"]
metal_entity = None

for entity in nonpolymer_entities:
    comp_id = entity["rcsb_nonpolymer_entity_container_identifiers"][
        "nonpolymer_comp_id"
    ]

    if comp_id == "CU":
        metal_entity = entity
        break

if metal_entity is None:
    raise ValueError(f"Nessuna entita' CU trovata per {pdb_code}")

container = metal_entity["rcsb_nonpolymer_entity_container_identifiers"]

asym_ids = container["asym_ids"]

print("Metal:", container["nonpolymer_comp_id"])
print("Asym IDs:", asym_ids)

url_atoms = f"https://models.rcsb.org/v1/{pdb_code}/atoms"
params = {
    "label_asym_id": asym_ids[0],
    "encoding": "cif",
    "copy_all_categories": "false",
    "download": "false",
}
response_atoms = httpx.get(
    url_atoms,
    params=params,
)
response_atoms.raise_for_status()
cif_text = response_atoms.text
# print(cif_text)

# read residues surrounding directly from cif

doc = gemmi.cif.read_string(cif_text)
block = doc.sole_block()

data = block.get_mmcif_category("_struct_conn")
idx = data["id"].index("metalc21")
connection = {key: values[idx] for key, values in data.items()}
print(connection)


def search_nonpolymer_f(df_structure: pd.DataFrame):
    """
    Check if structure from PDB bind the
    ligand/s of interest
    """
    for entry in df_structure.itertuples:
        pdb_code = entry.database_id
        url_pdb = f"https://files.rcsb.org/download/{pdb_code}.cif"
        response_pdb = httpx.get(url_pdb)

        if response_pdb.status_code == 200:
            content_pdb = (
                response_pdb.text
            )  # FIXME: this will be used in the savefile section below

        # TODO: insert a function here with a flag (default=false) for
        # TODO saving the files, variable read from argument parser input

        nonpolymer_entities = extract_graphsql(pdb_code=pdb_code)

        # transform nonpolymer_entities into a big list/dict and return everything, handled by an external function (to be defined)
        ligand_ids = read_nonpolymer(non_pol=nonpolymer_entities)

        # HACK search atoms around
        a_around = search_surround(pdb_code=pdb_code, ligands_ids=ligand_ids)


def extract_graphsql(pdb_code: str):

    url_graph = "https://data.rcsb.org/graphql"
    query = """
        query GetEntry($pdb_id: String!) {
            entry(entry_id: $pdb_id) {
                rcsb_id
                nonpolymer_entities {
                    rcsb_id
                    rcsb_nonpolymer_entity_container_identifiers {
                        entry_id
                        entity_id
                        auth_asym_ids
                        asym_ids
                        nonpolymer_comp_id
                    }
                    nonpolymer_comp {
                        chem_comp {
                            id
                            formula_weight
                            name
                            formula
                        }
                    }
                }
            } 
        }
    """
    variables = {"pdb_id": pdb_code}
    response_graph = requests.post(
        url_graph, json={"query": query, "variables": variables}
    )

    data = response_graph.json()
    return data["data"]["entry"]["nonpolymer_entities"]


def read_nonpolymer(non_pol):

    for f in non_pol:
        comp_id = f["rcsb_nonpolymer_entity_container_identifiers"][
            "nonpolymer_comp_id"
        ]

        if comp_id == "CU":
            metal_entity = f
            print("Metal:", comp_id)
            print("Asym IDs:", comp_id["asym_ids"])

    if metal_entity is None:
        raise ValueError(f"Nessuna entita' CU trovata per {pdb_code}")

    return comp_id["asym_ids"]


def search_surround(pdb_code, ligands_ids):

    url_atoms = f"https://models.rcsb.org/v1/{pdb_code}/atoms"

    for l in ligands_ids:
        params = {
            "label_asym_id": asym_ids[0],
            "encoding": "cif",
            "copy_all_categories": "false",
            "download": "false",
        }
        response_atoms = httpx.get(
            url_atoms,
            params=params,
        )
        response_atoms.raise_for_status()
        cif_text = response_atoms.text
        doc = gemmi.cif.read_string(cif_text)
        block = doc.sole_block()

        data = block.get_mmcif_category("_struct_conn")
        idx = data["id"].index("metalc21")
        connection = {key: values[idx] for key, values in data.items()}

        return connection
