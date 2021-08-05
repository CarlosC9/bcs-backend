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
    from biobarcoding import engine
    if sess is None:
        if engine is None:
            from biobarcoding.common.pg_helpers import create_pg_database_engine
            from biobarcoding.db_models import DBSession
            from biobarcoding import get_global_configuration_variable
            db_connection_string = get_global_configuration_variable('DB_CONNECTION_STRING')
            engine = create_pg_database_engine(db_connection_string, "bcs", recreate_db=False)
            DBSession.configure(bind=engine)  # reconfigure the sessionmaker used by this scoped_session
        sess = DBSession()

    _ = []
    any_gbif_request = False
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
                url = f"https://api.gbif.org/v1/species/match?name={urllib.parse.quote_plus(sn)}"
                response = requests.request("GET", url)
                if response.status_code == 200:
                    species_name = SpeciesNameToCanonical()
                    species_name.name = lsn
                    r = response.json()
                    if r["matchType"] == "NONE" or (r["matchType"] == "FUZZY" and r["confidence"] <= 70):
                        species_name.canonical_name = ""
                        species_name.scientific_name = ""
                    elif r["matchType"] == "EXACT" or (r["matchType"] == "FUZZY" and r["confidence"] > 70):
                        species_name.rank = r["rank"]
                        species_name.synonym = r["synonym"]
                        species_name.canonical_name = r["canonicalName"]
                        species_name.scientific_name = r["scientificName"]
                        species_names_map[lsn] = species_name.canonical_name
                        found = True
                    any_gbif_request = True
                    sess.add(species_name)
        else:
            found = True
        if found:
            v = species_names_map[sn.lower()]
        else:
            v = None
        _.append(v)
    if any_gbif_request:
        sess.commit()

    return _


# lst = ["Euphorbia lathyris L.",
#        'Lobularia libyca (Viv.) Meisn.',
#        'Schizogyne sericea (L. f.) DC.',
#        'Wahlenbergia lobelioides (L. f.) Link subsp. lobelioides',
#        'Cistus monspeliensis L. subsp. canariensis Rivas-Mart., Martín-Osorio & Wildpret']
# s = get_canonical_species_names(None, lst)
# print(s)
# s = get_canonical_species_names(None, lst)
# print(s)
