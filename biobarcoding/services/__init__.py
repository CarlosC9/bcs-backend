import os


##
# CONNECTIONS
##

def conn_chado():
    from flask import current_app
    from chado import ChadoInstance
    conn = ChadoInstance(dbhost=current_app.config["CHADO_HOST"],
        dbname=current_app.config["CHADO_DATABASE"],
        dbuser=current_app.config["CHADO_USER"],
        dbpass=current_app.config["CHADO_PASSWORD"],
        dbschema=current_app.config["CHADO_SCHEMA"],
        dbport=current_app.config["CHADO_PORT"])
    return conn


##
# EXECUTIONS
##

def exec_cmds(*args):
    import subprocess
    out = err = []
    for cmd in args:
        print(cmd)
        process = subprocess.Popen(cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True)
        o, e = process.communicate()
        out.append(o)
        print(f'OUT: \n{o.decode("utf-8")}\n')
        if e:
            print(f'ERROR: \n{e.decode("utf-8")}')
            err.append(e)
    return out, err


##
# VISUALIZATIONS
##

def unfolded_print(obj, level=2):
    try:
        d = obj if isinstance(obj, dict) else obj.__dict__
        for key in d.keys():
            val = d.get(key)
            print(f"{level * ' '}{key} ({type(val)}):")
            unfolded_print(val, level=level + 2)
    except Exception as e:
        if isinstance(obj, (tuple, list, set)):
            print(f"{level * ' '}[")
            for i in obj:
                unfolded_print(i, level=level + 2)
                print(f"{level * ' '} ,")
            print(f"{level * ' '}]")
        else:
            print(f"{level*' '} {obj}")


##
# CONVERTERS
##

def orm2json(row):
    # TODO: missing inherited fields
    d = {}
    try:
        for column in row.__table__.columns:
            d[column.name] = str(getattr(row, column.name))
    except:
        pass
    return d


##
# SQLALCHEMY OPS
##

def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if not instance:
        params = dict((k, v) for k, v in kwargs.items())
        instance = model(**params)
        session.add(instance)
        session.flush()
    return instance


def get_simple_query(session, model, **kwargs):
    kwargs = {k:v for k,v in kwargs.items() if v is not None}
    return session.query(model).filter_by(**kwargs)


##
# LOGS
##

def log_exception(e):
    import os, sys
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    print(e)


##
# FILES MANAGEMENT
##

def get_encoding(file):
    import chardet
    with open(file, 'rb') as fo:
        enc = chardet.detect(fo.read())['encoding']
    return enc


def get_bioformat(input_file, format):
    ext = format if format else os.path.splitext(input_file)[1][1:]
    return {'frn': 'fasta', 'fna': 'fasta', 'faa': 'fasta', 'fas': 'fasta', 'fasta': 'fasta',
            'gb': 'genbank', 'gbf': 'genbank', 'gbk': 'genbank', 'genbank': 'genbank',
            'gff': 'gff3', 'gff3': 'gff3',
            'nex': 'nexus', 'nxs': 'nexus', 'nexus': 'nexus',
            'aln': 'clustal', 'clustal': 'clustal',
            'phy': 'phylip', 'phylip': 'phylip',
            'nwx': 'newick', 'tree': 'newick', 'newick': 'newick'}.get(ext)


def tsv_csv_parser(file):   # with csv lib
    with open(file, encoding=get_encoding(file)) as fo:
        import csv
        tab = csv.reader(fo, delimiter='\t')
        keys = next(tab)
        for rec in tab:
            yield dict(zip(keys, rec))


def tsv_pd_parser(file):   # with pandas lib
    import pandas as pd
    df = pd.read_csv(file, sep='\t', header=0, error_bad_lines=False,
                     encoding=get_encoding(file))
    for i, r in df.iterrows():
        yield r.to_dict()


def gff_parser(file):
    from BCBio import GFF
    with open(file, encoding=get_encoding(file)) as f:
        for rec in GFF.parse(f):
            yield rec


def seqs_parser(file, format='fasta'):
    if format == 'gff':
        return gff_parser(file)
    if format == 'tsv':
        return tsv_pd_parser(file)
    from Bio import SeqIO
    with open(file, encoding=get_encoding(file)) as fo:
        for rec in SeqIO.parse(fo, format):
            yield rec
