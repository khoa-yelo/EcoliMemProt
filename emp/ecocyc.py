"""Module to query data from Ecocyc"""
import sys

directory_to_add = r"C:\Users\khoah\CovertLabProjects\EcoliMemProt"
sys.path.append(directory_to_add)

import requests
import xmltodict
import pandas as pd
from emp.utils import find_values_by_key, join_list

LINKED_DATABASE = "dblink"


class EcocycQuerier:
    query_str = "https://websvc.biocyc.org/getxml?id=ECOLI:{0}&detail=full"
    failed_queries = []

    def __init__(self, email="", password=""):
        self.session = self.initiate_sesson(email, password)

    def initiate_sesson(self, email="", password=""):
        if not email:
            email = input("Enter Email: ")
        if not password:
            password = input("Enter Password: ")
        session = requests.Session()
        session.post(
            "https://websvc.biocyc.org/credentials/login/",
            data={"email": email, "password": password},
        )

        return session

    def query(self, entity_id, return_dict=False, query_str=""):
        if not query_str:
            query_str = self.query_str.format(entity_id)
        query_res = self.session.get(query_str)
        if query_res.status_code != 200:
            self.failed_queries.append(entity_id)
            print("Error ", entity_id, query_res.status_code)
        else:
            content = xmltodict.parse(query_res.content)
        self.content = content
        if return_dict:
            return content

    def look_for(self, key):
        query_content = find_values_by_key(self.content, key)
        query_content = join_list(query_content)
        df_res = pd.DataFrame(query_content)
        protein_data_result = {}
        if key == LINKED_DATABASE:
            for db in ["PDB", "ALPHAFOLD"]:
                ids = (
                    df_res.groupby("dblink-db")
                    .agg(list)
                    .to_dict()["dblink-oid"]
                    .get(db, [])
                )
                protein_data_result[db] = ids
        if protein_data_result:
            return protein_data_result
        else:
            return df_res


if __name__ == "__main__":
    # Input email and password to avoid manual input:
    # eq = EcocycQuerier("email@stanford.edu", "password")
    # Add your email and password below before running
    eq = EcocycQuerier("@stanford.edu", "")
    eq.query("NUOA-MONOMER")
    res = eq.look_for(LINKED_DATABASE)
    print(res)
