import gemmi
import httpx
import pandas as pd
import requests


def read_uniprot(json_content):
    for entry in json_content:
        uniprot_json(id_uniprot=entry["id_uniprot"])


def pdb_info(id_uniprot, pdb_codes):
    report = []

    for info in pdb_codes:
        database_id = info["id"]
        method = next(x["value"] for x in info["properties"] if x["key"] == "Method")
        resolution = next(
            x["value"] for x in info["properties"] if x["key"] == "Resolution"
        )
        chain = next(x["value"] for x in info["properties"] if x["key"] == "Chains")

        report.append(
            {
                "id_uniprot": id_uniprot,
                "database": "PDB",
                "database_id": database_id,
                "method": method,
                "resolution": resolution,
                "chain": chain,
            }
        )
    return pd.DataFrame(report)


def uniprot_json(id_uniprot):

    url = f"https://rest.uniprot.org/uniprotkb/{id_uniprot.upper()}.json"

    try:
        response = httpx.get(url, timeout=10)
        response.raise_for_status()
        content = response.json()

    except httpx.HTTPStatusError:
        raise

    except httpx.RequestError:
        raise

    pdb_codes = [
        elem
        for elem in content["uniProtKBCrossReferences"]
        if elem["database"] == "PDB"
    ]

    return pdb_info(id_uniprot=id_uniprot, pdb_codes=pdb_codes)


def read_proteinDataBank(df_structure: pd.DataFrame, store_pdb, pdb_dir):
    for entry in df_structure.itertuples:
        search_nonpolymer_f(
            pdb_code=entry.database_id, store_pdb=store_pdb, pdb_dir=pdb_dir
        )


def search_nonpolymer_f(pdb_code, store_pdb, pdb_dir):
    """
    Check if structure from PDB bind the
    ligand/s of interest
    """
    if store_pdb:
        url_pdb = f"https://files.rcsb.org/download/{pdb_code.upper()}.cif"

        try:
            response_pdb = httpx.get(url_pdb)
            response_pdb.raise_for_status()
            content_pdb = response_pdb.text

        except httpx.HTTPStatusError:
            raise
        except httpx.RequestError:
            raise

        content_pdb = response_pdb.text

        with open(f"{pdb_dir}/{pdb_code.upper()}.cif", "w") as k:
            k.write(content_pdb)

    nonpolymer_entities = extract_graphsql(pdb_code=pdb_code)

    # transform nonpolymer_entities into a big list/dict and return everything, handled by an external function (to be defined)
    ligand_ids = read_nonpolymer(non_pol=nonpolymer_entities, pdb_code=pdb_code)

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


def read_nonpolymer(non_pol, pdb_code):

    for f in non_pol:
        comp_id = f["rcsb_nonpolymer_entity_container_identifiers"][
            "nonpolymer_comp_id"
        ]

        if comp_id == "CU":  # FIXME: sostituire CU con il metallo letto dal json
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
