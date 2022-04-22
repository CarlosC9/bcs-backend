import re
import urllib.parse
from typing import List

import requests

from ..db_models.metadata import SpeciesNameToCanonical

GBIF_URL = 'https://api.gbif.org/v1/species/'
GBIF_MATCH_URL = GBIF_URL + 'match?name='

species_names_map = {}


def get_taxon_from_gbif(sn: str):
    # Match using GBIF
    url = f"{GBIF_MATCH_URL}{urllib.parse.quote_plus(sn)}"
    r = requests.request("GET", url)
    if r.status_code == 200:
        return r.json()
    return None


def get_lineage_from_gbif(key: int):
    # Match using GBIF
    url = f"{GBIF_URL}{key}/parents"
    r = requests.request("GET", url)
    if r.status_code == 200:
        return r.json()
    return None


def get_canonical_species_names(sess, in_: List[str], underscores=False) -> List[str]:
    """
    Canonicalize species names

    :param sess: Database session to access the cache table
    :param in_: List of species names to canonicalize
    :param underscores: if True, replace whitespace and "-" by "_"
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
                url = f"{GBIF_MATCH_URL}{urllib.parse.quote_plus(sn)}"
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
            v = species_names_map[lsn]
            if underscores and v:
                v = re.sub("[, -]", "_", v)
        else:
            v = None
        _.append(v)
    if any_gbif_request:
        sess.commit()

    return _


def get_canonical_species_info(sess, in_: List[str]) -> List[dict]:
    """
    Canonicalize species names and provide gbif data

    :param sess: Database session to access the cache table
    :param in_: List of species names to canonicalize
    :return: List of gbif data. None if it was not possible to do
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
        r = get_taxon_from_gbif(sn)
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
            elif r:
                species_name = SpeciesNameToCanonical()
                species_name.name = lsn
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
            v = r
        else:
            v = None
        _.append(v)
    if any_gbif_request:
        sess.commit()

    return _


def get_canonical_species_lineages(sess, in_: List[str]) -> List[List[dict]]:
    """
    Canonicalize species names and provide gbif data

    :param sess: Database session to access the cache table
    :param in_: List of species names to canonicalize
    :return: List of gbif data. None if it was not possible to do
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
        r = get_taxon_from_gbif(sn)
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
            elif r:
                species_name = SpeciesNameToCanonical()
                species_name.name = lsn
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
        if found and r.get('usageKey'):
            v = get_lineage_from_gbif(r["usageKey"])
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
