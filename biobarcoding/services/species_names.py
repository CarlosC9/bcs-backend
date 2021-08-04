import urllib.parse
from typing import List
import requests

from biobarcoding.db_models.metadata import SpeciesNameToCanonical

species_names_map = {}


def get_canonical_species_names(sess, in_: List[str]) -> List[str]:
    """
    Canonicalize species names

    :param in_: List of species names to canonicalize
    :return: List of canonicalized species names. None if it was not possible to do
    """

    _ = []
    for sn in in_:
        lsn = sn.lower().strip()
        found = False
        if lsn not in species_names_map:
            # Check existence in database
            species_name = sess.query(SpeciesNameToCanonical).filter(SpeciesNameToCanonical.name == lsn).first()
            if species_name:
                # The entry may exist but it can be wrong
                if species_name.canonical_name == "":
                    pass  # Activate the following line to avoid repeating the lookup into the database
                    # species_names_map[lsn] = None
                else:
                    found = True
                    species_names_map[lsn] = species_name.canonical_name
            else:
                # Match using GBIF
                url = f"https://api.gbif.org/v1/species/match?name={urllib.parse.urlencode(sn)}"
                response = requests.request("GET", url)
                if response.status_code == 200:
                    species_name = SpeciesNameToCanonical()
                    species_name.name = lsn
                    r = response.json()
                    if r["matchType"] is "NONE" or (r["matchType"] is "FUZZY" and r["confidence"] <= 70):
                        species_name.canonical_name = ""
                        species_name.scientific_name = ""
                    elif r["matchType"] is "EXACT" or (r["matchType"] is "FUZZY" and r["confidence"] > 70):
                        species_name.rank = r["rank"]
                        species_name.synonym = r["synonym"]
                        species_name.canonical_name = r["canonicalName"]
                        species_name.scientific_name = r["scientificName"]
                    sess.add(species_name)
        else:
            found = True
        if found:
            v = species_names_map[sn.lower()]
        else:
            v = None
        _.append(v)
    return v
