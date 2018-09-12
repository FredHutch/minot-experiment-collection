"""Functions for getting information from the eggNOG API."""

import requests
from functools import lru_cache


@lru_cache(maxsize=None)
def get_species_for_eggnog_cluster(eggnog_cluster):
    # Strip out everything before the '.'
    eggnog_cluster = eggnog_cluster.split(".", 1)[-1].upper()
    
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

    dat = r.json()

    assert isinstance(dat, dict), "Returned data is formatted unexpectedly ({})".format(eggnog_cluster)
    assert "matches" in dat, "Returned data is formatted unexpectedly ({})".format(eggnog_cluster)

    species_list = set()

    for match in dat["matches"]:
        assert isinstance(match, dict)
        assert "quick_members" in match
        for member in match["quick_members"]:
            assert isinstance(member, list)
            assert len(member) == 4
            species_list.add(member[2])

    return list(species_list)                
