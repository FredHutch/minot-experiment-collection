"""Functions for getting data from the KEGG API."""

import requests
from functools import lru_cache

# KEGG names
@lru_cache(maxsize=None)
def get_kegg_name(ko):
    r = requests.get("http://rest.kegg.jp/list/{}".format(ko))
    return r.text.split("\t")[-1].rstrip("\n")
