import zlib

from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import mapper, sessionmaker


class MatrixDigest(object):
    pass


def inflate(data):
    decompress = zlib.decompressobj()
    inflated = decompress.decompress(data)
    inflated += decompress.flush()
    return inflated  #.decode("utf-8")


def dump_matrices(conn_string: str, md_path: str):
    engine = create_engine(conn_string, echo=True)
    metadata = MetaData(engine)
    ms_matrix_digests = Table("matrix_digests", metadata, autoload=True)
    mapper(MatrixDigest, ms_matrix_digests)
    Session = sessionmaker(bind=engine)
    session = Session()
    matrices = session.query(MatrixDigest).all()
    for matrix in matrices:
        print(f"ID: {matrix.id}, Individuals: {matrix.individuals}, Loci: {matrix.number_of_loci}")
        dec = inflate(matrix.data)
        with open(f"{md_path}/{matrix.id}.xml", "wb") as f:
            f.write(dec)
    session.close()


if __name__ == '__main__':
    # Assuming a temporary MySQL in Docker, with the Demiurge database:
    #   docker run -d -p 33060:3306 --name mysql-db -e MYSQL_ROOT_PASSWORD=secret mysql
    dump_matrices("mysql+pymysql://root:secret@localhost:33060/demiurgene_production",
                  "/home/rnebot/Downloads/borrame/demiurgene_matrix_digests/")
