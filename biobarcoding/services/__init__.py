

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
from ..rest import filter_parse, order_parse


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
# CONVERSIONS
##

def force_underscored(text: str, replace: str = " ,.-") -> str:
    if not replace:
        replace = "".join(set(c for c in text if not c.isalnum()))
    return text.translate({ord(i): "_" for i in replace})


def orm2json(row):
    # TODO: replace for generate_json ?
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
    import os
    ext = format if format else os.path.splitext(file)[1][1:]
    return {'frn': 'fasta', 'fna': 'fasta', 'faa': 'fasta', 'fas': 'fasta', 'fasta': 'fasta',
            'gb': 'genbank', 'gbf': 'genbank', 'gbk': 'genbank', 'genbank': 'genbank',
            'gff': 'gff3', 'gff3': 'gff3',
            'nex': 'nexus', 'nxs': 'nexus', 'nexus': 'nexus',
            'aln': 'clustal', 'clustal': 'clustal',
            'ph': 'phylip', 'phy': 'phylip', 'phylip': 'phylip',
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
    schema = orm.Schema._declared_fields or orm.__table__.columns
    return dict([(k, v) for k, v in params.items() if k in schema])


def get_or_create(session, model, **params):
    # params = get_orm_params(model, **kwargs)
    instance = session.query(model).filter_by(**params).first()
    if not instance:
        instance = model(**params)
        session.add(instance)
        session.flush()
        # TODO: create ACLs if possible too ?
        try:
            from ..db_models import DBSession
            # acl
            from ..db_models.sysadmin import ACL
            acl = get_or_create(DBSession, ACL,
                                    object_uuid=instance.uuid,
                                    object_type=instance.obj_type_id)
            # acldetail
            from ..db_models.sysadmin import ACLDetail
            from flask import g
            from ..db_models.sysadmin import PermissionType
            acldetail = get_or_create(DBSession, ACLDetail,
                                          acl_id=acl.id,
                                          authorizable_id=g.n_session.identity.id,
                                          permission_id=DBSession.query(PermissionType.id).filter(PermissionType.name=='delete').one())
        except Exception as e:
            pass
    return instance


def get_simple_query(session, model, null_sensitive=True, **kwargs):
    if not null_sensitive:
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
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


def get_query(session, orm, query=None, id=None, aux_filter=None, default_filter=None, aux_order=None, **kwargs):
    """
     reserved keywords in kwargs:
       'values': specific values of orm fields to filter
       'filter': advanced filtering clauses (see also filter_parse)
       'order': advanced ordering clauses (see also order_parse)
       'pagination': pageIndex and pageSize to paginate
       'searchValue': full-text search value (hopefully)
     otherwise it will be treated as 'values'
    """
    query = query or session.query(orm)
    count = 0
    if id:
        query = query.filter(orm.id == id)
    else:
        if not kwargs.get('values'):
            kwargs['values'] = {}
        if default_filter:
            if "filter" in kwargs:
                kwargs["filter"].append(default_filter)
            else:
                kwargs["filter"] = [default_filter]
        for k, v in kwargs.items():
            if not k in ['values', 'filter', 'order', 'pagination', 'searchValue'] and v:
                kwargs['values'][k] = v
        if kwargs.get('values'):
            query = query.filter_by(**get_orm_params(orm, **kwargs.get('values')))
        if kwargs.get('filter'):
            query = query.filter(filter_parse(orm, kwargs.get('filter', {}), aux_filter, session))
        if kwargs.get('searchValue', '') != '':
            if hasattr(orm, '_ts_vector'):
                query = query.filter(orm._ts_vector.match(kwargs.get('searchValue')))
        if kwargs.get('order'):
            for o in order_parse(orm, kwargs.get('order'), aux_order):
                query = query.order_by(o)
        count = query.count()
        if kwargs.get('pagination'):
            query = paginator(query, kwargs.get('pagination'))
    return query, count
