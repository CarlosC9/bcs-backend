import os


##
# CONNECTIONS
##
from biobarcoding.rest import filter_parse, order_parse


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
# CONVERTIONS
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


def get_bioformat(file, format):
    ext = format if format else os.path.splitext(file)[1][1:]
    return {'frn': 'fasta', 'fna': 'fasta', 'faa': 'fasta', 'fas': 'fasta', 'fasta': 'fasta',
            'gb': 'genbank', 'gbf': 'genbank', 'gbk': 'genbank', 'genbank': 'genbank',
            'gff': 'gff3', 'gff3': 'gff3',
            'nex': 'nexus', 'nxs': 'nexus', 'nexus': 'nexus',
            'aln': 'clustal', 'clustal': 'clustal',
            'phy': 'phylip', 'phylip': 'phylip',
            'nhx': 'newick', 'nwx': 'newick', 'tree': 'newick', 'newick': 'newick'}.get(ext)


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


##
# SQLALCHEMY OPS
##

def get_orm_params(orm, **params):
    return dict([(k, v) for k, v in params.items() if k in orm.__table__.columns])


def get_or_create(session, model, **params):
    # params = get_orm_params(model, **kwargs)
    instance = session.query(model).filter_by(**params).first()
    if not instance:
        instance = model(**params)
        session.add(instance)
        session.flush()
    return instance


def get_simple_query(session, model, **kwargs):
    # kwargs = {k: v for k, v in kwargs.items() if v is not None}
    params = get_orm_params(model, **kwargs)
    return session.query(model).filter_by(**params)


def paginator(query, pagination):
    if 'pageIndex' in pagination and 'pageSize' in pagination:
        page = pagination.get('pageIndex')
        page_size = pagination.get('pageSize')
        return query \
            .offset((page - 1) * page_size) \
            .limit(page_size)
    return query


def get_query(session, orm, query=None, id=None, aux_filter=None, aux_order=None, **kwargs):
    """
     reserved keywords in kwargs:
       'value': specific values of orm fields to filter
       'filter': advanced filtering clauses (see also filter_parse)
       'order': advanced ordering clauses (see also order_parse)
       'pagination': pageIndex and pageSize to paginate
       'searchValue': full-text search value (hopefully)
     otherwise it will be treated as 'value'
    """
    query = query or session.query(orm)
    count = 0
    if id:
        query = query.filter(orm.id == id)
    else:
        if not kwargs.get('value'):
            kwargs['value'] = {}
        for k, v in kwargs.items():
            if not k in ['value', 'filter', 'order', 'pagination', 'searchValue'] and v:
                kwargs['value'][k] = v
        if kwargs.get('value'):
            query = query.filter_by(**get_orm_params(orm, **kwargs.get('value')))
        if kwargs.get('filter'):
            query = query.filter(filter_parse(orm, kwargs.get('filter'), aux_filter))
        if kwargs.get('order'):
            query = query.order_by(order_parse(orm, kwargs.get('order'), aux_order))
        count = query.count()
        if kwargs.get('pagination'):
            query = paginator(query, kwargs.get('pagination'))
    return query, count