from sqlalchemy import create_engine, text, Table, Column, BigInteger, ForeignKey, Index, MetaData, UniqueConstraint
import os, sys, getopt, subprocess

CHADO_XML_FOLDER = "chado_xml"
DB_RELATIONSHIP_COMMENT = """Specifies relationships between databases.  This is 
particularly useful for ontologies that use multiple prefix IDs for it\'s vocabularies. For example,
the EDAM ontology uses the prefixes \"data\", \"format\", \"operation\" and others. Each of these would
have a record in the db table.  An \"EDAM\" record could be added for the entire ontology to the
db table and the previous records could be linked as \"part_of\" EDAM.  As another example
databases housing cross-references may have sub databases such as NCBI (e.g. Taxonomy, SRA, etc).
This table can use a \'part_of\' record to link all of them to NCBI."""

def main(host):
    engine = create_engine('postgresql://postgres:postgres@' + host + ':5432/')
    with engine.connect().execution_options(autocommit=True) as conn:
        edam_id = conn.execute(text("select db_id from db where name=\'EDAM\'")).fetchone()
        if edam_id is None:
            meta = MetaData()
            meta.reflect(engine, only=["db", "cvterm"])#this get information about the db and cvterm tables in the existing database
            #CREATE db_relationship table to relate EDAM submodules with EDAM itself
            db_relationship = Table('db_relationship', meta,
                                    Column('db_relationship_id', BigInteger, primary_key=True),
                                    Column('type_id', 
                                            BigInteger,
                                            ForeignKey('cvterm.cvterm_id', ondelete="CASCADE", initially="DEFERRED"),
                                            nullable=False),
                                    Column('subject_id', 
                                            BigInteger,
                                            ForeignKey('db.db_id', ondelete="CASCADE", initially="DEFERRED"),
                                            nullable=False),
                                    Column('object_id',
                                            BigInteger,
                                            ForeignKey('db.db_id', ondelete="CASCADE", initially="DEFERRED"),
                                            nullable=False),
                                    UniqueConstraint('type_id', 'subject_id', 'object_id', name='db_relationship_c1'),
                                    Index('db_relationship_idx1', 'type_id'),
                                    Index('db_relationship_idx2', 'subject_id'),
                                    Index('db_relationship_idx3', 'object_id'),
                                    comment=DB_RELATIONSHIP_COMMENT
                                  )
            db_relationship.create(bind=conn)
            
            #INSERT THE CHADO XML (EDAM) IN THE DATABASE
            command = ["stag-storenode.pl", "-d", "dbi:Pg:dbname=postgres;host=" + host + ";port=5432",
                        "--user", "postgres", "--password", "postgres", os.path.join(CHADO_XML_FOLDER, "EDAM.chado")]

            subprocess.run(command)

            #FILL THE db_relationship TABLE
            edam_id = conn.execute(text("select db_id from db where name=\'EDAM\'")).fetchone()[0]
            db_ids = conn.execute(text("select db_id from db where db_id > " + str(edam_id))).fetchall()
            type_id = conn.execute(text("select cvterm_id from cvterm join cv on cvterm.cv_id = cv.cv_id where cvterm.name = \'part_of\' and cv.name = \'relationship\'")).fetchone()[0]

            for db_id in db_ids:
              conn.execute(text("INSERT INTO db_relationship (type_id,subject_id,object_id)" +
                                  "VALUES(" + str(type_id) + "," + str(db_id[0]) + "," + str(edam_id) + ")"))

            
        else:
            print("EDAM was already inserted")

if __name__ == "__main__":
    try:
      argv = sys.argv[1:]
      opts, args = getopt.getopt(argv,"h:")
    except getopt.GetoptError:
      print('edam_insertion.py -h <hostname>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         main(arg)
         sys.exit()