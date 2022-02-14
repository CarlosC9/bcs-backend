"""From empty database to current

Revision ID: 15183552bca7
Revises: 
Create Date: 2022-01-11 14:50:20.479711

"""
from alembic import op
import sqlalchemy as sa
import sqlalchemy_utils
import biobarcoding
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = '15183552bca7'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('bo_bar_coding_regions',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.UnicodeText(), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('bo_specimens',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.UnicodeText(), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('fs_folders',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('fs_objects',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('fso_type', sa.String(length=6), nullable=False),
    sa.Column('name', sa.String(length=1024), nullable=True),
    sa.Column('full_name', sa.String(length=2048), nullable=True),
    sa.Column('parent_id', sa.BigInteger(), nullable=True),
    sa.ForeignKeyConstraint(['parent_id'], ['fs_folders.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_foreign_key(None, 'fs_folders', 'fs_objects', ['id'], ['id'])
    op.create_table('geo_regions',
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('geo_id', sa.Integer(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('identity_id', sa.Integer(), nullable=True),
    sa.Column('attributes', sa.JSON(), nullable=True),
    sa.Column('style', sa.JSON(), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('geo_id')
    )
    op.create_table('hie_hierarchy_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_algorithm_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_job_mgmt_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_job_statuses',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_processes',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.Text(), nullable=True),
    sa.Column('schema_inputs', sa.JSON(), nullable=True),
    sa.Column('schema_outputs', sa.JSON(), nullable=True),
    sa.Column('execution', sa.JSON(), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('md_ontologies',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('md_publications',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('md_species_names',
    sa.Column('name', sa.String(length=500), nullable=False),
    sa.Column('rank', sa.String(length=20), nullable=True),
    sa.Column('synonym', sa.Boolean(), nullable=True),
    sa.Column('canonical_name', sa.String(length=200), nullable=True),
    sa.Column('scientific_name', sa.String(length=200), nullable=True),
    sa.Column('is_canonical', sa.Boolean(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.PrimaryKeyConstraint('name')
    )
    op.create_table('md_taxa',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('md_taxonomies',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('ngd_port_types',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('name', sa.String(length=128), nullable=True),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('object_types',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=False),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('sa_annotation_form_item',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('name', sa.String(length=80), nullable=False),
    sa.Column('description', sa.String(length=512), nullable=True),
    sa.Column('cvterm_id', sa.Integer(), nullable=True),
    sa.Column('dbxref_id', sa.Integer(), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('cvterm_id'),
    sa.UniqueConstraint('name')
    )
    op.create_table('sa_auth_authenticators',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('validation_endpoint', sa.String(length=1024), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('sa_auth_authorizables',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('authorizable_type_id', sa.Integer(), nullable=False),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('sa_auth_collections',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('name', sa.String(length=512), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_functions',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=512), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('sa_auth_permission_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('rank', sa.Integer(), nullable=False),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_task_statuses',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('fs_files',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('content_type', sa.String(length=256), nullable=True),
    sa.Column('embedded_content', sa.LargeBinary(), nullable=True),
    sa.Column('content_size', sa.Integer(), nullable=True),
    sa.Column('content_location', sa.JSON(), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['fs_objects.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('hie_hierarchies',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('h_type_id', sa.Integer(), nullable=True),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['h_type_id'], ['hie_hierarchy_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_algorithms',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.Text(), nullable=True),
    sa.Column('algorithm_type_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['algorithm_type_id'], ['jobs_algorithm_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_compute_resources',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('jm_type_id', sa.Integer(), nullable=False),
    sa.Column('jm_location', sa.JSON(), nullable=False),
    sa.Column('jm_credentials', sa.JSON(), nullable=False),
    sa.Column('jm_params', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('agreement', sa.JSON(), nullable=True),
    sa.Column('consumption_counters', sa.JSON(), nullable=True),
    sa.ForeignKeyConstraint(['jm_type_id'], ['jobs_job_mgmt_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_process_configurations',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.Text(), nullable=True),
    sa.Column('process_id', sa.Integer(), nullable=False),
    sa.Column('params', sa.JSON(), nullable=True),
    sa.ForeignKeyConstraint(['process_id'], ['jobs_processes.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_wfs_job_mgmt_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('process_id', sa.Integer(), nullable=False),
    sa.Column('jm_type_id', sa.Integer(), nullable=False),
    sa.Column('input_mapper', sa.JSON(), nullable=True),
    sa.Column('output_mapper', sa.JSON(), nullable=True),
    sa.ForeignKeyConstraint(['jm_type_id'], ['jobs_job_mgmt_types.id'], ),
    sa.ForeignKeyConstraint(['process_id'], ['jobs_processes.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('ngd_fos',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('do_type_id', sa.Integer(), nullable=True),
    sa.Column('native_id', sa.BigInteger(), nullable=True),
    sa.Column('native_table', sa.String(length=80), nullable=True),
    sa.Column('name', sa.String(length=300), nullable=True),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['do_type_id'], ['object_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('native_table', 'native_id', name='ngd_fos_c1'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('sa_annotation_form_field',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('type', sa.String(length=32), nullable=False),
    sa.Column('range', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('view_type', sa.String(length=32), nullable=True),
    sa.Column('multiple', sa.Boolean(), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_form_item.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_form_item_object_type',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('form_item_id', sa.Integer(), nullable=False),
    sa.Column('object_type_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['form_item_id'], ['sa_annotation_form_item.id'], ),
    sa.ForeignKeyConstraint(['object_type_id'], ['object_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('form_item_id', 'object_type_id', name='sa_annotation_form_item_object_type_c1')
    )
    op.create_table('sa_annotation_form_template',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_form_item.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_collections_detail',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('collection_id', sa.Integer(), nullable=False),
    sa.Column('object_type', sa.Integer(), nullable=False),
    sa.Column('object_uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.ForeignKeyConstraint(['collection_id'], ['sa_auth_collections.id'], ),
    sa.ForeignKeyConstraint(['object_type'], ['object_types.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_groups',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_auth_authorizables.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_identities',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=255), nullable=True),
    sa.Column('email', sa.String(length=255), nullable=True),
    sa.Column('can_login', sa.Boolean(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.Column('deactivation_time', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_auth_authorizables.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_obj_types_perm_types',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('object_type_id', sa.Integer(), nullable=False),
    sa.Column('permission_type_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['object_type_id'], ['object_types.id'], ),
    sa.ForeignKeyConstraint(['permission_type_id'], ['sa_auth_permission_types.id'], ),
    sa.PrimaryKeyConstraint('id', 'object_type_id', 'permission_type_id')
    )
    op.create_table('sa_auth_organizations',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_auth_authorizables.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_permissions',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('object_type', sa.Integer(), nullable=False),
    sa.Column('object_uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.ForeignKeyConstraint(['object_type'], ['object_types.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('object_uuid', 'object_type', name='sa_auth_permissions_c1')
    )
    op.create_table('sa_auth_roles',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_auth_authorizables.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_task_instances',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('description', sa.String(length=80), nullable=True),
    sa.Column('params', sa.JSON(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.Column('finalization_time', sa.DateTime(), nullable=True),
    sa.Column('status', sa.Integer(), nullable=True),
    sa.Column('log', sa.Text(), nullable=True),
    sa.ForeignKeyConstraint(['status'], ['sa_task_statuses.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('bo_msas',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('bo_phylo_trees',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('bo_seq_sims',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('bo_sequences',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('specimen_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.ForeignKeyConstraint(['specimen_id'], ['bo_specimens.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('fs_fos_in_files',
    sa.Column('file_id', sa.BigInteger(), nullable=False),
    sa.Column('fos_id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['file_id'], ['fs_files.id'], ),
    sa.ForeignKeyConstraint(['fos_id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('file_id', 'fos_id')
    )
    op.create_table('hie_hierarchy_levels',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('name', sa.String(length=128), nullable=True),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('hierarchy_id', sa.Integer(), nullable=False),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['hierarchy_id'], ['hie_hierarchies.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_algorithm_configurations',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('description', sa.Text(), nullable=True),
    sa.Column('params', sa.JSON(), nullable=True),
    sa.Column('algorithm_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['algorithm_id'], ['jobs_algorithms.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_jobs',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('process_id', sa.Integer(), nullable=False),
    sa.Column('resource_id', sa.Integer(), nullable=True),
    sa.Column('identity_id', sa.Integer(), nullable=False),
    sa.Column('inputs', sa.JSON(), nullable=True),
    sa.Column('status', sa.String(length=80), nullable=True),
    sa.Column('log', sa.Text(), nullable=True),
    sa.Column('outputs', sa.JSON(), nullable=True),
    sa.Column('execution_state', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('deleted', sa.Boolean(), nullable=True),
    sa.ForeignKeyConstraint(['identity_id'], ['sa_auth_identities.id'], ),
    sa.ForeignKeyConstraint(['process_id'], ['jobs_processes.id'], ),
    sa.ForeignKeyConstraint(['resource_id'], ['jobs_compute_resources.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_process_algorithms',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('process_id', sa.Integer(), nullable=False),
    sa.Column('algorithm_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['algorithm_id'], ['jobs_algorithms.id'], ),
    sa.ForeignKeyConstraint(['process_id'], ['jobs_processes.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('jobs_wfs_compute_resources',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('process_id', sa.Integer(), nullable=False),
    sa.Column('resource_id', sa.Integer(), nullable=False),
    sa.Column('native_process_id', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['process_id'], ['jobs_processes.id'], ),
    sa.ForeignKeyConstraint(['resource_id'], ['jobs_compute_resources.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('ngd_case_studies',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('ngd_datasets',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('identity_id', sa.Integer(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.Column('structure', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('provenance', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('ngd_processes',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('process_type', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_form_attribute',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_form_field.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_form_relationship',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_form_field.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_form_tag',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_form_field.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_form_template_field',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('form_field_id', sa.Integer(), nullable=False),
    sa.Column('form_template_id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=True),
    sa.Column('rank', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['form_field_id'], ['sa_annotation_form_field.id'], ),
    sa.ForeignKeyConstraint(['form_template_id'], ['sa_annotation_form_template.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('form_template_id', 'rank', name='sa_annotation_form_template_field_c1')
    )
    op.create_table('sa_annotation_item',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('object_uuid', biobarcoding.db_models.GUID(), nullable=False),
    sa.Column('type', sa.String(length=80), nullable=False),
    sa.Column('name', sa.String(length=80), nullable=False),
    sa.Column('rank', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['object_uuid'], ['ngd_fos.uuid'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('object_uuid', 'name', name='sa_annotation_item_c2'),
    sa.UniqueConstraint('object_uuid', 'rank', name='sa_annotation_item_c1')
    )
    op.create_table('sa_auth_groups_identities',
    sa.Column('group_id', sa.Integer(), nullable=False),
    sa.Column('identity_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['group_id'], ['sa_auth_groups.id'], ),
    sa.ForeignKeyConstraint(['identity_id'], ['sa_auth_identities.id'], ),
    sa.PrimaryKeyConstraint('group_id', 'identity_id')
    )
    op.create_table('sa_auth_identities_authenticators',
    sa.Column('identity_id', sa.Integer(), nullable=False),
    sa.Column('authenticator_id', sa.Integer(), nullable=False),
    sa.Column('email', sa.String(length=255), nullable=True),
    sa.Column('name', sa.String(length=255), nullable=True),
    sa.Column('authenticator_info', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('last_login_time', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['authenticator_id'], ['sa_auth_authenticators.id'], ),
    sa.ForeignKeyConstraint(['identity_id'], ['sa_auth_identities.id'], ),
    sa.PrimaryKeyConstraint('identity_id', 'authenticator_id')
    )
    op.create_index('index_IdentityAuthenticator_on_AuthenticatorInfo_gin', 'sa_auth_identities_authenticators', ['authenticator_info'], unique=False, postgresql_using='gin')
    op.create_table('sa_auth_organizations_identities',
    sa.Column('organization_id', sa.Integer(), nullable=False),
    sa.Column('identity_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['identity_id'], ['sa_auth_identities.id'], ),
    sa.ForeignKeyConstraint(['organization_id'], ['sa_auth_organizations.id'], ),
    sa.PrimaryKeyConstraint('organization_id', 'identity_id')
    )
    op.create_table('sa_auth_permissions_expression',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('acl_id', sa.Integer(), nullable=True),
    sa.Column('expression', sa.String(length=500), nullable=True),
    sa.Column('validity_start', sa.DateTime(), nullable=True),
    sa.Column('validity_end', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['acl_id'], ['sa_auth_permissions.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_roles_identities',
    sa.Column('role_id', sa.Integer(), nullable=False),
    sa.Column('identity_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['identity_id'], ['sa_auth_identities.id'], ),
    sa.ForeignKeyConstraint(['role_id'], ['sa_auth_roles.id'], ),
    sa.PrimaryKeyConstraint('role_id', 'identity_id')
    )
    op.create_table('sa_browser_filters',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('name', sa.String(length=80), nullable=False),
    sa.Column('type', sa.String(length=80), nullable=False),
    sa.Column('values', sa.JSON(), nullable=True),
    sa.Column('user_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['user_id'], ['sa_auth_identities.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('name', 'type', 'user_id', name='sa_browser_filters_c1'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('geo_layers',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('wks', sa.String(length=80), nullable=True),
    sa.Column('geoserver_name', sa.String(length=80), nullable=True),
    sa.Column('wms_url', sa.String(length=255), nullable=True),
    sa.Column('properties', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.Column('published', sa.Boolean(), nullable=True),
    sa.Column('in_postgis', sa.Boolean(), nullable=True),
    sa.Column('is_deleted', sa.Boolean(), nullable=True),
    sa.Column('layer_type', sa.String(length=80), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['ngd_datasets.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_index('index_GeographicLayer_on_Properties_gin', 'geo_layers', ['properties'], unique=False, postgresql_using='gin')
    op.create_table('hie_hierarchy_nodes',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('name', sa.String(length=128), nullable=True),
    sa.Column('uuid', biobarcoding.db_models.GUID(), nullable=True),
    sa.Column('hierarchy_id', sa.Integer(), nullable=False),
    sa.Column('parent_node_id', sa.Integer(), nullable=True),
    sa.Column('level_id', sa.Integer(), nullable=True),
    sa.Column('reference_node_id', sa.Integer(), nullable=True),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['hierarchy_id'], ['hie_hierarchies.id'], ),
    sa.ForeignKeyConstraint(['level_id'], ['hie_hierarchy_levels.id'], ),
    sa.ForeignKeyConstraint(['parent_node_id'], ['hie_hierarchy_nodes.id'], ),
    sa.ForeignKeyConstraint(['reference_node_id'], ['hie_hierarchy_nodes.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('uuid')
    )
    op.create_table('ngd_case_studies_functional_objects',
    sa.Column('case_study_id', sa.BigInteger(), nullable=False),
    sa.Column('functional_object_id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['case_study_id'], ['ngd_case_studies.id'], ),
    sa.ForeignKeyConstraint(['functional_object_id'], ['ngd_fos.id'], ),
    sa.PrimaryKeyConstraint('case_study_id', 'functional_object_id')
    )
    op.create_table('ngd_dataset_port_types',
    sa.Column('dataset_id', sa.BigInteger(), nullable=False),
    sa.Column('port_type_id', sa.BigInteger(), nullable=False),
    sa.ForeignKeyConstraint(['dataset_id'], ['ngd_datasets.id'], ),
    sa.ForeignKeyConstraint(['port_type_id'], ['ngd_port_types.id'], ),
    sa.PrimaryKeyConstraint('dataset_id', 'port_type_id')
    )
    op.create_table('ngd_process_instances',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('instantiated_process_id', sa.BigInteger(), nullable=True),
    sa.Column('status', sa.String(length=10), nullable=True),
    sa.Column('job_id', sa.BigInteger(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.Column('hash_of_canonical', sa.String(length=64), nullable=True),
    sa.Column('params', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['ngd_fos.id'], ),
    sa.ForeignKeyConstraint(['instantiated_process_id'], ['ngd_processes.id'], ),
    sa.ForeignKeyConstraint(['job_id'], ['jobs_jobs.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('hash_of_canonical')
    )
    op.create_table('ngd_processes_ports',
    sa.Column('id', sa.BigInteger(), autoincrement=True, nullable=False),
    sa.Column('name', sa.String(length=128), nullable=True),
    sa.Column('process_id', sa.BigInteger(), nullable=True),
    sa.Column('port_type_id', sa.BigInteger(), nullable=True),
    sa.Column('input', sa.Boolean(), nullable=False),
    sa.Column('attributes', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['port_type_id'], ['ngd_port_types.id'], ),
    sa.ForeignKeyConstraint(['process_id'], ['ngd_processes.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_field',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('form_field_id', sa.Integer(), nullable=False),
    sa.Column('value', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['form_field_id'], ['sa_annotation_form_field.id'], ),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_item.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_template',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('form_template_id', sa.Integer(), nullable=False),
    sa.Column('value', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
    sa.ForeignKeyConstraint(['form_template_id'], ['sa_annotation_form_template.id'], ),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_item.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_annotation_text',
    sa.Column('id', sa.BigInteger(), nullable=False),
    sa.Column('value', sa.String(length=512), nullable=True),
    sa.ForeignKeyConstraint(['id'], ['sa_annotation_item.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('sa_auth_permissions_detail',
    sa.Column('id', sa.Integer(), autoincrement=True, nullable=False),
    sa.Column('acl_id', sa.Integer(), nullable=False),
    sa.Column('acl_expression_id', sa.Integer(), nullable=True),
    sa.Column('authorizable_id', sa.Integer(), nullable=False),
    sa.Column('permission_id', sa.Integer(), nullable=False),
    sa.Column('validity_start', sa.DateTime(), nullable=True),
    sa.Column('validity_end', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['acl_expression_id'], ['sa_auth_permissions_expression.id'], ),
    sa.ForeignKeyConstraint(['acl_id'], ['sa_auth_permissions.id'], ),
    sa.ForeignKeyConstraint(['authorizable_id'], ['sa_auth_authorizables.id'], ondelete='CASCADE'),
    sa.ForeignKeyConstraint(['permission_id'], ['sa_auth_permission_types.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('ngd_process_instances_ports',
    sa.Column('process_instance_id', sa.BigInteger(), nullable=False),
    sa.Column('port_id', sa.BigInteger(), nullable=False),
    sa.Column('dataset_id', sa.BigInteger(), nullable=True),
    sa.ForeignKeyConstraint(['dataset_id'], ['ngd_datasets.id'], ),
    sa.ForeignKeyConstraint(['port_id'], ['ngd_processes_ports.id'], ),
    sa.ForeignKeyConstraint(['process_instance_id'], ['ngd_process_instances.id'], ),
    sa.PrimaryKeyConstraint('process_instance_id', 'port_id')
    )
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('ngd_process_instances_ports')
    op.drop_table('sa_auth_permissions_detail')
    op.drop_table('sa_annotation_text')
    op.drop_table('sa_annotation_template')
    op.drop_table('sa_annotation_field')
    op.drop_table('ngd_processes_ports')
    op.drop_table('ngd_process_instances')
    op.drop_table('ngd_dataset_port_types')
    op.drop_table('ngd_case_studies_functional_objects')
    op.drop_table('hie_hierarchy_nodes')
    op.drop_index('index_GeographicLayer_on_Properties_gin', table_name='geo_layers', postgresql_using='gin')
    op.drop_table('geo_layers')
    op.drop_table('sa_browser_filters')
    op.drop_table('sa_auth_roles_identities')
    op.drop_table('sa_auth_permissions_expression')
    op.drop_table('sa_auth_organizations_identities')
    op.drop_index('index_IdentityAuthenticator_on_AuthenticatorInfo_gin', table_name='sa_auth_identities_authenticators', postgresql_using='gin')
    op.drop_table('sa_auth_identities_authenticators')
    op.drop_table('sa_auth_groups_identities')
    op.drop_table('sa_annotation_item')
    op.drop_table('sa_annotation_form_template_field')
    op.drop_table('sa_annotation_form_tag')
    op.drop_table('sa_annotation_form_relationship')
    op.drop_table('sa_annotation_form_attribute')
    op.drop_table('ngd_processes')
    op.drop_table('ngd_datasets')
    op.drop_table('ngd_case_studies')
    op.drop_table('jobs_wfs_compute_resources')
    op.drop_table('jobs_process_algorithms')
    op.drop_table('jobs_jobs')
    op.drop_table('jobs_algorithm_configurations')
    op.drop_table('hie_hierarchy_levels')
    op.drop_table('fs_fos_in_files')
    op.drop_table('bo_sequences')
    op.drop_table('bo_seq_sims')
    op.drop_table('bo_phylo_trees')
    op.drop_table('bo_msas')
    op.drop_table('sa_task_instances')
    op.drop_table('sa_auth_roles')
    op.drop_table('sa_auth_permissions')
    op.drop_table('sa_auth_organizations')
    op.drop_table('sa_auth_obj_types_perm_types')
    op.drop_table('sa_auth_identities')
    op.drop_table('sa_auth_groups')
    op.drop_table('sa_auth_collections_detail')
    op.drop_table('sa_annotation_form_template')
    op.drop_table('sa_annotation_form_item_object_type')
    op.drop_table('sa_annotation_form_field')
    op.drop_table('ngd_fos')
    op.drop_table('jobs_wfs_job_mgmt_types')
    op.drop_table('jobs_process_configurations')
    op.drop_table('jobs_compute_resources')
    op.drop_table('jobs_algorithms')
    op.drop_table('hie_hierarchies')
    op.drop_table('fs_files')
    op.drop_table('sa_task_statuses')
    op.drop_table('sa_auth_permission_types')
    op.drop_table('sa_auth_functions')
    op.drop_table('sa_auth_collections')
    op.drop_table('sa_auth_authorizables')
    op.drop_table('sa_auth_authenticators')
    op.drop_table('sa_annotation_form_item')
    op.drop_table('object_types')
    op.drop_table('ngd_port_types')
    op.drop_table('md_taxonomies')
    op.drop_table('md_taxa')
    op.drop_table('md_species_names')
    op.drop_table('md_publications')
    op.drop_table('md_ontologies')
    op.drop_table('jobs_processes')
    op.drop_table('jobs_job_statuses')
    op.drop_table('jobs_job_mgmt_types')
    op.drop_table('jobs_algorithm_types')
    op.drop_table('hie_hierarchy_types')
    op.drop_table('geo_regions')
    op.drop_table('fs_objects')
    op.drop_table('fs_folders')
    op.drop_table('bo_specimens')
    op.drop_table('bo_bar_coding_regions')
    # ### end Alembic commands ###
