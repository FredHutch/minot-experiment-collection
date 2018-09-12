"""Functions for getting information from the eggNOG API."""

import requests
from functools import lru_cache


@lru_cache(maxsize=None)
def get_species_for_eggnog_cluster(eggnog_cluster):
    # Strip out everything before the '.'
    eggnog_cluster = eggnog_cluster.split(".", 1)[-1].upper()
    
    # Get the NOGNAME
    data = {
        "desc": "",
        "seqid": "*@{}".format(eggnog_cluster),
        "target_species": "",
        "level": "",
        "nognames": "",
        "page": 0
    }
    r = requests.post("http://eggnogapi.embl.de/meta_search", data=data)
    assert r.status_code == 200, "Could not fetch data for {}".format(eggnog_cluster)

    # Now get the cluster members
    r = requests.get(
        "http://eggnogapi.embl.de/nog_data/json/extended_members/{}".format(r.json()["matches"][0]["nogname"]))

    dat = r.json()
    
    assert isinstance(dat, dict), "Returned data is formatted unexpectedly ({})".format(eggnog_cluster)
    assert "members" in dat, "Returned data is formatted unexpectedly ({})".format(eggnog_cluster)

    species_list = set()

    for m in dat["members"].values():
        assert isinstance(m, list)
        species_list.add(m[0])

    return list(species_list)                
