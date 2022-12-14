#!/usr/bin/perl

use strict;
use warnings;
#use lib '/Users/cain/cvs_stuff/schema/trunk/chado/lib';
use Bio::FeatureIO;
use Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use URI::Escape;
use Carp;
use Bio::GMOD::DB::Adapter;
use ExtUtils::MakeMaker;  #to get prompt
use Module::Load;

=head1 NAME

$0 - Bulk loads gff3 files into a chado database.

=head1 SYNOPSIS

  % $0 [options]
  % cat <gff-file> | $0 [options]

=head1 OPTIONS

 --gfffile         The file containing GFF3 (optional, can read 
                     from stdin)
 --fastafile       Fasta file to load sequence from
 --organism        The organism for the data 
                    (use the value 'fromdata' to read from GFF organism=xxx)
 --dbprofile       Database config profile name
 --dbname          Database name
 --dbuser          Database user name
 --dbpass          Database password
 --dbhost          Database host
 --dbport          Database port
 --analysis        The GFF data is from computational analysis
 --noload          Create bulk load files, but don't actually load them.
 --nosequence      Don't load sequence even if it is in the file
 --notransact      Don't use a single transaction to load the database
 --drop_indexes    Drop indexes of affected tables before starting load
                     and recreate after load is finished; generally
                     does not help performance.
 --validate        Validate SOFA terms before attempting insert (can
                     cause script startup to be slow, off by default)
 --ontology        Give directions for handling misc Ontology_terms
 --skip_vacuum     Skip vacuuming the tables after the inserts (default)
 --no_skip_vaccum  Don't skip vacuuming the tables
 --inserts         Print INSERT statements instead of COPY FROM STDIN
 --noexon          Don't convert CDS features to exons (but still create
                     polypeptide features) 
 --recreate_cache  Causes the uniquename cache to be recreated
 --remove_lock     Remove the lock to allow a new process to run
 --save_tmpfiles   Save the temp files used for loading the database
 --random_tmp_dir  Use a randomly generated tmp dir (the default is
                     to use the current directory)
 --no_target_syn   By default, the loader adds the targetId in 
                     the synonyms list of the feature. This flag 
                     desactivate this.
 --unique_target   Trust the unicity of the target IDs. IDs are case 
                     sensitive. By default, the uniquename of a new target 
                     will be 'TargetId_PrimaryKey'. With this flag, 
                     it will be 'TargetId'. Furthermore, the Name of the 
                     created target will be its TargetId, instead of the
                     feature's Name.
 --dbxref          Use either the first Dbxref annotation as the
                     primary dbxref (that goes into feature.dbxref_id),
                     or if an optional argument is supplied, the first
                     dbxref that has a database part (ie, before the ':')
                     that matches the supplied pattern is used. 
 --delete          Instead of inserting features into the database,
                     use the GFF lines to delete features as though
                     the CRUD=delete-all option were set on all lines
                     (see 'Deletes and updates via GFF below'). The
                     loader will ask for confirmation before continuing.
 --delete_i_really_mean_it
                   Works like --delete except that it does not ask
                     for confirmation.
 --fp_cv           Name of the feature property controlled vocabulary
                     (defaults to 'feature_property').
 --noaddfpcv       By default, the loader adds GFF attribute types as
                     new feature_property cv terms when missing.  This flag
                     deactivates it.
   ** dgg note: should rename this flag: --[no]autoupdate 
            for Chado tables cvterm, cv, db, organism, analysis ...
   
 --manual          Detailed manual pages
 --custom_adapter  Use a custom subclass adaptor for Bio::GMOD::DB::Adapter
                     Provide the path to the adapter as an argument
 --private_schema  Load the data into a non-public schema.
 --use_public_cv   When loading into a non-public schema, load any cv and
                     cvterm data into the public schema
 --end_sql         SQL code to execute after the data load is complete
 --allow_external_parent 
                   Allow Parent tags to refer to IDs outside the current
                   GFF file

Note that all of the arguments that begin 'db' as well as organism can
be provided by default by Bio::GMOD::Config, which was installed when
'make install' was run.  Also note the the option dbprofile and all other
db* options are mutually exclusive--if you supply dbprofile, do not
supply any other db* options, as they will not be used.

=head1 DESCRIPTION

The GFF in the datafile must be version 3 due to its tighter control of
the specification and use of controlled vocabulary.  Accordingly, the names
of feature types must be exactly those in the Sequence Ontology Feature
Annotation (SOFA), not the synonyms and not the accession numbers (SO
accession numbers may be supported in future versions of this script).

Note that the ##sequence-region directive is not supported as a way of
declaring a reference sequence for a GFF3 file.  The ##sequence-region
directive is not expressive enough to define what type of thing the
sequence is (ie, is it a chromosome, a contig, an arm, etc?).  If
your GFF file uses a ##sequence-region directive in this way, you
must convert it to a full GFF3 line.  For example, if you have 
this line:

  ##sequence-region chrI 1 9999999

Then is should be converted to a GFF3 line like this:

  chrI	.	chromosome	1	9999999	.	.	.	ID=chrI

=head2 How GFF3 is stored in chado

Here is summary of how GFF3 data is stored in chado:

=over

=item Column 1 (reference sequence)

The reference sequence for the feature becomes the srcfeature_id
of the feature in the featureloc table for that feature.  That featureloc 
generally assigned a rank of zero if there are other locations associated
with this feature (for instance, for a match feature), the other locations
will be assigned featureloc.rank values greater than zero.

=item Column 2 (source)

The source is stored as a dbxref.  The chado instance must of an entry
in the db table named 'GFF_source'.  The script will then create a dbxref
entry for the feature's source and associate it to the feature via
the feature_dbxref table.

=item Column 3 (type)

The cvterm.cvterm_id of the SOFA type is stored in feature.type_id.

=item Column 4 (start)

The value of start minus 1 is stored in featureloc.fmin (one is subtracted
because chado uses interbase coordinates, whereas GFF uses base coordinates).

=item Column 5 (end)

The value of end is stored in featureloc.fmax.

=item Column 6 (score)

The score is stored in one of the score columns in the analysisfeature 
table.  The default is analysisfeature.significance.  See the
section below on analysis results for more information.

=item Column 7 (strand)

The strand is stored in featureloc.strand.

=item Column 8 (phase)

The phase is stored in featureloc.phase.  Note that there is currently
a problem with the chado schema for the case of single exons having 
different phases in different transcripts.  If your data has just such
a case, complain to gmod-schema@lists.sourceforge.net to find ways
to address this problem.

=item Column 9 (group)

Here is where the magic happens.

=over

=item Assigning feature.name, feature.uniquename

The values of feature.name and feature.uniquename are assigned 
according to these simple rules:

=over 

=item If there is an ID tag, that is used as feature.uniquename

otherwise, it is assigned a uniquename that is equal to
'auto' concatenated with the feature_id.

=item If there is a Name tag, it's value is set to feature.name;

otherwise it is null.

Note that these rules are much more simple than that those that
Bio::DB::GFF uses, and may need to be revisited.

=back

=item Assigning feature_relationship entries

All Parent tagged features are assigned feature_relationship
entries of 'part_of' to their parent features.  Derived_from
tags are assigned 'derived_from' relationships.  Note that
parent features must appear in the file before any features
use a Parent or Derived_from tags referring to that feature.

=item Alias tags

Alias values are stored in the synonym table, under
both synonym.name and synonym.synonym_sgml, and are
linked to the feature via the feature_synonym table.

=item Dbxref tags

Dbxref values must be of the form 'db_name:accession', where 
db_name must have an entry in the db table, with a value of 
db.name equal to 'DB:db_name'; several database names were preinstalled
with the database when 'make prepdb' was run.  Execute 'SELECT name
FROM db' to find out what databases are already availble.  New dbxref
entries are created in the dbxref table, and dbxrefs are linked to
features via the feature_dbxref table.

=item Gap tags

Currently is mostly ignored--the value is stored as a featureprop,
but otherwise is not used yet.

=item Note tags

The values are stored as featureprop entries for the feature.

=item Any custom (ie, lowercase-first) tags

Custom tags are supported.  If the tag doesn't already exist in 
the cvterm table, it will be created.  The value will stored with its 
associated cvterm in the featureprop table.

=item Ontology_term

When the Ontology_term tags are used, items from the Gene Ontology
and Sequence Ontology will be processed automatically when the standard
DB:accession format is used (e.g. GO:0001234).  To use other ontology
terms, you must specify that mapping of the DB indentifiers in the GFF
file and the name of the ontologies in the cv table as a comma separated
tag=value pairs.  For example, to use plant and cell ontology terms,
you would supply on the command line:

  --ontology 'PO=plant ontology,CL=cell ontology'

where 'plant ontology' and 'cell ontology' are the names in the cv table
exactly as they appear.

=item Target tags

Proper processing of Target tags requires that there be two source features
already available in the database, the 'primary' source feature (the
chromosome or contig) and the 'subject' from the similarity analysis,
like an EST, cDNA or syntenic chromosome.  If the subject feature is not
present, the loader will attempt to create a placeholder feature object
in its place.  If you have a fasta file the contains the subject, you can
use the perl script, L<gmod_fasta2gff3.pl>, that comes with this distribution
to make a GFF3 file suitable for loading into chado before loading your
analysis results.

=item CDS and UTR features

The way CDS features are represented in Chado is as an intersection of
a transcript's exons and the transcripts polypeptide feature.  To allow
proper translation of a GFF3 file's CDS features, this loader will 
convert CDS and UTR feature lines to corresponding exon features (and add
a featureprop note that the exon was inferred from a GFF3 CDS and/or UTR line),
and create a polypeptide feature that spans the genomic region from
the start of translation to the stop.

If your GFF3 file contains both exon and CDS/UTR features, then you will
want to suppress the creation of the exon features and instead will
only want a polypeptide feature to be created.  To do this, use the
--noexon option.  In this case, the CDS and UTR features will still 
be converted to exon features as described above.

Note that in the case where your GFF file contains CDS and/or UTR features
that do not belong to 'central dogma' genes (that is, that have a
gene, transcript and CDS/exon features), none of the above will happen
and the features will be stored as is.

=back

=back

=head2 NOTES

=over

=item Loading fasta file

When the --fastafile is provided with an argument that is the path to
a file containing fasta sequence, the loader will attempt to update the
feature table with the sequence provided.  Note that the ID provided in the
fasta description line must exactly match what is in the feature table
uniquename field.  Be careful if it is possible that the uniquename of the
feature was changed to ensure uniqueness when it was loaded from the
original GFF.  Also note that when loading sequence from a fasta file, 
loading GFF from standard in is disabled.  Sorry for any inconvenience.


=item ##sequence-region

This script does not use sequence-region directives for anything.
If it represents a feature that needs to be inserted into the database,
it should be represented with a full GFF line.  This includes the
reference sequence for the features if it is not already in the database,
like a chromosome.  For example, this:

  ##sequence-region chr1 1	213456789

should change to this:

  chr1	UCSC	chromosome	1	213456789	.	.	.	ID=chr1

=item Transactions

This application will, by default, try to load all of the data at
once as a single transcation.  This is safer from the database's
point of view, since if anything bad happens during the load, the 
transaction will be rolled back and the database will be untouched.  
The problem occurs if there are many (say, greater than a 2-300,000)
rows in the GFF file.  When that is the case, doing the load as 
a single transcation can result in the machine running out of memory
and killing processes.  If --notranscat is provided on the commandline,
each table will be loaded as a separate transaction.

=item SQL INSERTs versus COPY FROM

This bulk loader was originally designed to use the PostgreSQL
COPY FROM syntax for bulk loading of data.  However, as mentioned
in the 'Transactions' section, memory issues can sometimes interfere
with such bulk loads.  In another effort to circumvent this issue,
the bulk loader has been modified to optionally create INSERT statements
instead of the COPY FROM statements.  INSERT statements will load
much more slowly than COPY FROM statements, but as they load and
commit individually, they are more more likely to complete successfully.
As an indication of the speed differences involved, loading 
yeast GFF3 annotations (about 16K rows), it takes about 5 times
longer using INSERTs versus COPY on my laptop.

=item Deletes and updates via GFF

There is rudimentary support for modifying the features in an
existing database via GFF.  Currently, there is only support for
deleting.  In order to delete, the GFF line must have a custom
tag in the ninth column, 'CRUD' (for Create, Replace, Update and
Delete) and have a recognized value.  Currently the two recognized
values are CRUD=delete and CRUD=delete-all.  

IMPORTANT NOTE: Using the delete operations have the potential of creating
orphan features (eg, exons whose gene has been deleted).  You should be
careful to make sure that doesn't happen. Included in this distribution is a
PostgreSQL trigger (written in plpgsql) that will delete all orphan children
recursively, so if a gene is deleted, all transcripts, exons and 
polypeptides that belong to that gene will be deleted too.  See
the file modules/sequence/functions/delete-trigger.plpgsql for
more information.

=over

=item delete

The delete option will delete one and only one feature for which the 
name, type and organism match what is in the GFF line with what is 
in the database.  Note that feature.uniquename are not considered, nor
are the coordinates presented in the GFF file.  This is so that 
updates via GFF can be done on the coordinants.  If there is more than 
one feature for which the name, type and organism match, the loader will
print an error message and stop.  If there are no features that match
the name, type and organism, the loader will print a warning message
and continue.

=item delete-all

The delete-all option works similarly to the delete option, except that 
it will delete all features that match the name, type and organism in the
GFF line (as opposed to allowing only one feature to be deleted).  If there
are no features that match, the loader will print a warning message and
continue.

=back

=item The run lock

The bulk loader is not a multiuser application.  If two separate
bulk load processes try to load data into the database at the same
time, at least one and possibly all loads will fail.  To keep this from
happening, the bulk loader places a lock in the database to prevent
other gmod_bulk_load_gff3.pl processes from running at the same time.
When the application exits normally, this lock will be removed, but if
it crashes for some reason, the lock will not be removed.  To remove the
lock from the command line, provide the flag --remove_lock.  Note that
if the loader crashed necessitating the removal of the lock, you also
may need to rebuild the uniquename cache (see the next section).

=item The uniquename cache

The loader uses the chado database to create a table that caches
feature_ids, uniquenames, type_ids, and organism_ids of the features
that exist in the database at the time the load starts and the
features that will be added when the load is complete.  If it is possilbe
that new features have been added via some method that is not this
loader (eg, Apollo edits or loads with XORT) or if a previous load using
this loader was aborted, then you should supply
the --recreate_cache option to make sure the cache is fresh.

=item Sequence

By default, if there is sequence in the GFF file, it will be loaded
into the residues column in the feature table row that corresponds
to that feature.  By supplying the --nosequence option, the sequence
will be skipped.  You might want to do this if you have very large
sequences, which can be difficult to load.  In this context, "very large"
means more than 200MB.

Also note that for sequences to load properly, the GFF file must have
the ##FASTA directive (it is required for proper parsing by Bio::FeatureIO),
and the ID of the feature must be exactly the same as the name of the
sequence following the > in the fasta section.

=item The ORGANISM table

This script assumes that the organism table is populated with information
about your organism.  If you are unsure if that is the case, you can
execute this command from the psql command-line:

  select * from organism;

If you do not see your organism listed, execute this command to insert it:

  insert into organism (abbreviation, genus, species, common_name)
                values ('H.sapiens', 'Homo','sapiens','Human');

substituting in the appropriate values for your organism.

=item Parents/children order

Parents must come before children in the GFF file.

=item Analysis

If you are loading analysis results (ie, blat results, gene predictions), 
you should specify the -a flag.  If no arguments are supplied with the
-a, then the loader will assume that the results belong to an analysis
set with a name that is the concatenation of the source (column 2) and
the method (column 3) with an underscore in between.  Otherwise, the
argument provided with -a will be taken as the name of the analysis
set.  Either way, the analysis set must already be in the analysis
table.  The easist way to do this is to insert it directly in the
psql shell:

  INSERT INTO analysis (name, program, programversion)
               VALUES  ('genscan 2005-2-28','genscan','5.4');

There are other columns in the analysis table that are optional; see
the schema documentation and '\d analysis' in psql for more information.

Chado has four possible columns for storing the score in the GFF score 
column; please use whichever is most appropriate and identifiy it 
with --score_col flag (significance is the default). Note that the name
of the column can be shortened to one letter.  If you have more than
one score associated with each feature, you can put the other scores in
the ninth column as a tag=value pair, like 'identity=99', and the
bulk loader will put it in the featureprop table (provided there
is a cvterm for identity; see the section above concerning custom
tags).  Available options are:

=over

=item significance (default)

=item identity

=item normscore

=item rawscore

=back

A planned addtion to the functionality of handling analysis results
is to allow "mixed" GFF files, where some lines are analysis results
and some are not.  Additionally, one will be able to supply lists
of types (optionally with sources) and their associated entry in the
analysis table.  The format will probably be tag value pairs:

  --analysis match:Rice_est=rice_est_blast, \
             match:Maize_cDNA=maize_cdna_blast, \
             mRNA=genscan_prediction,exon=genscan_prediction

=item Grouping features by ID

The GFF3 specification allows features like CDSes and match_parts to
be grouped together by sharing the same ID.  This loader does not support
this method of grouping.  Instead the parent feature must be explicitly
created before the parts and the parts must refer to the parent with the 
Parent tag.

=item External Parent IDs

The GFF3 specification states that IDs are only valid within a single
GFF file, so you can't have Parent tags that refer to IDs in another 
file.  By specificifying the "allow_external_parent" flag, you can
relax this restriction.  A word of warning however: if the parent feature's
uniquename/ID was modified during loading (to make it unique), this
functionality won't work, as it won't beable to find the original
feature correctly.  Actually, it may be worse than not working,
it may attach child features to the wrong parent.  This is why it is
a bad idea to use this functionality!  Please use with caution.

=back

=head1 AUTHORS

Allen Day E<lt>allenday@ucla.eduE<gt>, Scott Cain E<lt>scain@cpan.orgE<gt>

Copyright (c) 2011

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

my ($ORGANISM, $GFFFILE,$FASTAFILE,$DBPROFILE, $DBNAME, $DBUSER, $DBPASS,$DBHOST, $DBPORT, 
    $ANALYSIS, $ANALYSIS_GROUP, $GLOBAL_ANALYSIS, $NOLOAD, $VALIDATE, $INSERTS,
    $NOTRANSACT, $NOSEQUENCE, $SCORE_COL, $ONTOLOGY, $SKIP_VACUUM,
    $DROP_INDEX, $NOEXON, $NOUNIQUECACHE, $RECREATE_CACHE, $SAVE_TMPFILES,$RANDOM_TMP_DIR,
    $NO_TARGET_SYN, $UNIQUE_TARGET, $DBXREF, $FP_CV, $NO_ADDFP_CV, 
    $MANPAGE, $DEBUG, $DELETE, $DELETE_CONFIRM, $CUSTOM_ADAPTER,
    $PRIVATE_SCHEMA, $USE_PUBLIC_CV, $NO_SKIP_VACUUM, $END_SQL, $REMOVE_LOCK,
    $ALLOW_EXTERNAL_PARENT, );

GetOptions(
    'organism=s' => \$ORGANISM,
    'gfffile=s'  => \$GFFFILE,
    'fastafile=s'=> \$FASTAFILE,
    'dbprofile=s'=> \$DBPROFILE,
    'dbname=s'   => \$DBNAME,
    'dbuser=s'   => \$DBUSER,
    'dbpass=s'   => \$DBPASS,
    'dbhost=s'   => \$DBHOST,
    'dbport=s'   => \$DBPORT,
    'analysis:s' => \$ANALYSIS, # = means it is required, : means optional
    'noload'     => \$NOLOAD,
    'validate'   => \$VALIDATE,
    'notransact' => \$NOTRANSACT,
    'nosequence' => \$NOSEQUENCE,
    'score_col=s'=> \$SCORE_COL,
    'ontology=s' => \$ONTOLOGY,
    'skip_vacuum'=> \$SKIP_VACUUM,
    'no_skip_vacuum' => \$NO_SKIP_VACUUM,
    'drop_indexes'=>\$DROP_INDEX,
    'inserts'    => \$INSERTS,
    'noexon'     => \$NOEXON,
    'nouniquecache'=> \$NOUNIQUECACHE,
    'recreate_cache'=> \$RECREATE_CACHE,
    'remove_lock'   => \$REMOVE_LOCK,
    'save_tmpfiles'=>\$SAVE_TMPFILES,
    'random_tmp_dir'=>\$RANDOM_TMP_DIR,
    'no_target_syn'=> \$NO_TARGET_SYN,
    'unique_target' => \$UNIQUE_TARGET,
    'dbxref:s'      => \$DBXREF,
    'fp_cv:s'       => \$FP_CV,
    'noaddfpcv'     => \$NO_ADDFP_CV,
    'manual'   => \$MANPAGE,
    'debug'   => \$DEBUG,
    'delete'  => \$DELETE,
    'delete_i_really_mean_it' => \$DELETE_CONFIRM,
    'custom_adapter=s' => \$CUSTOM_ADAPTER,
    'private_schema=s' => \$PRIVATE_SCHEMA,
    'use_public_cv'    => \$USE_PUBLIC_CV,
    'end_sql:s'        => \$END_SQL,
    'allow_external_parent' => \$ALLOW_EXTERNAL_PARENT,
) 
# or ( system( 'pod2text', $0 ), exit -1 );
or pod2usage(-verbose => 1, -exitval => 1);
pod2usage(-verbose => 2, -exitval => 1) if $MANPAGE;

$SIG{__DIE__} = $SIG{INT} = 'cleanup_handler';

my $ORGANISM_FROM_CMDLINE = $ORGANISM;

unless ($DBNAME) {
    if (eval {require Bio::GMOD::Config;
          Bio::GMOD::Config->import();
          require Bio::GMOD::DB::Config;
          Bio::GMOD::DB::Config->import();
          1;  } ) {
        my $gmod_conf = $ENV{'GMOD_ROOT'} || "/var/lib/gmod" ?
                  Bio::GMOD::Config->new($ENV{'GMOD_ROOT'} || "/var/lib/gmod") :
                  Bio::GMOD::Config->new();

        my $profile = $DBPROFILE || 'default';
        my $db_conf = Bio::GMOD::DB::Config->new($gmod_conf,$profile);
        $DBNAME    = $db_conf->name();
        $DBUSER    = $db_conf->user();
        $DBPASS    = $db_conf->password();
        $DBHOST    = $db_conf->host();
        $DBPORT    = $db_conf->port();
        $ORGANISM ||= $db_conf->organism();
    }
}

$GFFFILE  ||='stdin';  #nobody better name their file 'stdin'
#$DBNAME   ||='chado';
$DBPASS   ||='';
$DBHOST   ||='localhost';
$DBPORT   ||='5432';
$VALIDATE ||=0;
$NOTRANSACT ||=0;
$NOSEQUENCE ||=0;
$INSERTS    ||=0;
$SCORE_COL  ||='significance';
$FP_CV      ||='feature_property';
$NO_ADDFP_CV ||=0;

die "You must supply a database name" unless $DBNAME;

#die "You must supply an organism" unless $ORGANISM;

$GLOBAL_ANALYSIS=0;
if ((defined $ANALYSIS) and ($ANALYSIS eq '')) { 
  $ANALYSIS = 1; #ie, it was specified on the command line with no arg
} elsif ($ANALYSIS) {
  $GLOBAL_ANALYSIS = 1;
  $ANALYSIS_GROUP = $ANALYSIS; # analysis group specified on the command line
  $ANALYSIS = 1;
} else {
  $ANALYSIS = 0;
}

if ((defined $DBXREF) and ($DBXREF eq '')) {
  $DBXREF=1; #ie, it was on the command line with no arg
}

$SKIP_VACUUM = $NO_SKIP_VACUUM ? 0 : 1;

my %argv;
  $argv{organism}     = $ORGANISM;
  $argv{gfffile}      = $GFFFILE;
  $argv{fastafile}    = $FASTAFILE;
  $argv{dbprofile}    = $DBPROFILE;
  $argv{dbname}       = $DBNAME;
  $argv{dbuser}       = $DBUSER;
  $argv{dbpass}       = $DBPASS;
  $argv{dbhost}       = $DBHOST;
  $argv{dbport}       = $DBPORT;
  $argv{is_analysis}  = $ANALYSIS;
  $argv{noload}       = $NOLOAD;
  $argv{validate}     = $VALIDATE;
  $argv{notransact}   = $NOTRANSACT;
  $argv{nosequence}   = $NOSEQUENCE;
  $argv{score_col}    = $SCORE_COL;
  $argv{ontology}     = $ONTOLOGY;
  $argv{skip_vacuum}  = $SKIP_VACUUM;
  $argv{drop_indexes} = $DROP_INDEX;
  $argv{inserts}      = $INSERTS;
  $argv{global_analysis}=$GLOBAL_ANALYSIS;
  $argv{analysis_group}=$ANALYSIS_GROUP;
  $argv{noexon}       = $NOEXON;
  $argv{nouniquecache}= $NOUNIQUECACHE;
  $argv{recreate_cache}=$RECREATE_CACHE;
  $argv{save_tmpfiles}= $SAVE_TMPFILES;
  $argv{random_tmp_dir}=$RANDOM_TMP_DIR;
  $argv{no_target_syn}= $NO_TARGET_SYN;
  $argv{unique_target}= $UNIQUE_TARGET;
  $argv{dbxref}       = $DBXREF;
  $argv{fp_cv}        = $FP_CV;
  $argv{addpropertycv}= ! $NO_ADDFP_CV; #dgg
  $argv{private_schema} = $PRIVATE_SCHEMA;
  $argv{use_public_cv}= $USE_PUBLIC_CV;
  $argv{allow_external_parent} = $ALLOW_EXTERNAL_PARENT;

#allow a feature_id to be referenced by multiple featureloc.feature_id's
my %locgroup = ();

########################

my $adapter = 'Bio::GMOD::DB::Adapter';
if ($CUSTOM_ADAPTER) {
  $adapter = "$adapter\:\:$CUSTOM_ADAPTER";
  unless (eval {load $adapter;
            1; } ) {
    warn "\n\nWhile attempting to load Bio::GMOD::DB::Adapter::$CUSTOM_ADAPTER\n";
    warn "an error was encountered.  Please check that this is a valid perl module\n";
    warn "and that it is in perl's search path.\n";

    warn "\nError was: $@\n";

    warn "Exiting...\n\n";
    exit(-1);
  }
}

my $chado = $adapter->new(%argv);

$chado->remove_lock(force => 1) if $REMOVE_LOCK;
$chado->place_lock();
my $lock = 1;

#if we need custom ontology mapping, cache them here
if ($ONTOLOGY) {
  my @pairs = split /\,/, $ONTOLOGY;
  foreach (@pairs) {
    my ($tag, $value) = split/\=/;
    $chado->cache('ontology',$tag, $value);
  }
}

$chado->file_handles();

my $gffio;
if ($GFFFILE eq 'stdin' and !$FASTAFILE) {
    $gffio = Bio::FeatureIO->new(-fh   => \*STDIN , 
                                 -format => 'gff', 
                                 -ignore_seq_region => 1,
                                 -validate_terms => $VALIDATE);
} elsif ($GFFFILE and $GFFFILE ne 'stdin') {
    $gffio = Bio::FeatureIO->new(-file => $GFFFILE, 
                                 -format => 'gff', 
                                 -ignore_seq_region => 1,
                                 -validate_terms => $VALIDATE);
}

warn "Preparing data for inserting into the $DBNAME database\n";
warn "(This may take a while ...)\n";

my $seen_cds = my $seen_exon = my $seen_bad_cds = 0;

my $ORGANISM_FROMDATA= ($ORGANISM =~ /fromdata/);
if($ORGANISM_FROMDATA) {
    $ORGANISM= "null"; # is this useful?
} 
elsif ($ORGANISM_FROM_CMDLINE) {
    $ORGANISM = $ORGANISM_FROM_CMDLINE;
    $chado->organism($ORGANISM);
    $chado->organism_id($ORGANISM)
       or die "$ORGANISM organism not found in the database";
}
elsif (defined $gffio && $gffio->organism) {
    $ORGANISM = $gffio->organism;
    $chado->organism($ORGANISM);
    $chado->organism_id($ORGANISM)
       or die "$ORGANISM organism not found in the database";

}
else {
    $chado->organism($ORGANISM);
    $chado->organism_id($ORGANISM)
       or die "$ORGANISM organism not found in the database";
}

if ($DELETE && !$DELETE_CONFIRM) {
    my $really_delete 
   = prompt("Do you really want me to delete features using this GFF file",'N');
    if (lc $really_delete eq 'y') {
        $DELETE_CONFIRM = $DELETE;
    }
    else {
        warn "OK, I'm stopping instead.\n\n";
        $chado->remove_lock();
        exit(0);
    }
}

my ($analysis_err_str,$cds_err_str);
my $processing_CDS = 0;
my $feature_iterator; my $itern=0; my %seen_organism;
FEATURE:
while(my $feature = 
  (defined $feature_iterator && $feature_iterator->next_feature) 
    || (defined $gffio && $gffio->next_feature)){
  $chado->primary_dbxref('');
  
  my $featuretype = $feature->type->name;

  # dgg: pull organism from 1st feature??
  # * may be many per gff-file; e.g. uniprot input
  if($ORGANISM_FROMDATA) {
    my($gff_organism) = 
               defined($feature->annotation->get_Annotations('organism'))
              ? ($feature->annotation->get_Annotations('organism'))[0]
              : '';
    if($gff_organism && $gff_organism ne $chado->organism()) {
      #next 4 lines patch by Alexie P.
      $ORGANISM = "$gff_organism";# is it pesky bperl object? # sadly yes
      if (ref($gff_organism)){
         $ORGANISM = $gff_organism->value();
      }
      warn "Organism $ORGANISM from data\n" unless($seen_organism{$ORGANISM}++);
      $chado->organism($ORGANISM);
      $chado->organism_id($ORGANISM) #? dont die if many orgs? auto-add *
        or die "$ORGANISM organism not found in the database";
      }
  }

  if($feature->annotation->get_Annotations('CRUD')||$DELETE_CONFIRM ) {
    my $continue = $chado->handle_crud($feature, $DELETE_CONFIRM);
    next if ($continue == 1 || $DELETE_CONFIRM);  
                             #that is, this was a delete or update operation
                             #and it is done
  }
  
  $itern++;
  warn(join(",","f".$itern,$featuretype,$feature->seq_id),"\n") if $DEBUG;

  $seen_exon= 1 if $featuretype =~ /exon$/ and !$processing_CDS;
  if ($featuretype =~ /(CDS|UTR)/) {
    my $continue_on = $chado->handle_CDS($feature);
    $seen_cds = 1 if !$seen_cds && $featuretype =~ /CDS/;
    next FEATURE unless ($continue_on == 0); 

    if (!$seen_bad_cds ) {
      warn <<END;

This GFF file has CDS and/or UTR features that do not belong to a 
'central dogma' gene (ie, gene/transcript/CDS).  The features of 
this type are being stored in the database as is.

END
;
      $seen_bad_cds = 1;
    }
  }

  if ( !$cds_err_str && $seen_cds && $seen_exon && !$NOEXON && !$processing_CDS) {
    $cds_err_str = 
        "\n\nThere are both CDS and exon features in this file, but\n"
       ."you did not set the --noexon option, which you probably want.\n"
       ."Please see `perldoc gmod_bulk_load_gff3.pl` for more information.\n\n";
    warn $cds_err_str;
  }
  my $nextfeature    = $chado->nextfeature();
  my $nextfeatureloc = $chado->nextfeatureloc();

  my $type           = $chado->get_type($featuretype);
  my ($src, $seqlen) = $chado->get_src_seqlen($feature);

  if(!$src){
    $src = $chado->src_second_chance($feature);
  }
  die "no feature for ".$feature->seq_id unless $src;

  if($feature->annotation->get_Annotations('Parent')){
    $chado->handle_parent($feature);
  }

  if($feature->annotation->get_Annotations('Derives_from')){
    $chado->handle_derives_from($feature);
  }

  my $source      = defined ($feature->source) ?
                       $feature->source->value
                       : '.';

  my($uniquename) = defined(($feature->annotation->get_Annotations('ID'))[0]) ?
                ($feature->annotation->get_Annotations('ID'))[0]->value
                : "auto".$nextfeature;
  $uniquename     = $uniquename->value if ref($uniquename);

  $uniquename     = $chado->uniquename_validation( $uniquename,
                                       $type,
                                       $chado->organism_id,
                                       $nextfeature);

  if (defined(($feature->annotation->get_Annotations('ID'))[0]) &&
      ($feature->annotation->get_Annotations('ID'))[0]->value ne $uniquename) {
    #need to keep a temporary map of modified uniquenames
    $chado->modified_uniquename(
          orig_id => ($feature->annotation->get_Annotations('ID'))[0]->value,
          modified_id => $uniquename,
          organism_id=>$chado->organism_id);
  }

  my($name)    = defined(($feature->annotation->get_Annotations('Name'))[0]) ?
                   ($feature->annotation->get_Annotations('Name'))[0]->value :
                 defined(($feature->annotation->get_Annotations('ID'))[0])   ?
                   ($feature->annotation->get_Annotations('ID'))[0]->value   :
                 "$featuretype-$uniquename";

  $name           = $name->value if ref($name);

  if ($uniquename eq $feature->seq_id or 
      (defined( $chado->modified_uniquename(modified_id => $uniquename, organism_id => $chado->organism_id)) and
        $chado->modified_uniquename(modified_id => $uniquename, organism_id => $chado->organism_id) eq $feature->seq_id)) {

    $chado->reftype_property($featuretype,$type); #  a reference sequence?? yes 
    } 
  
  my $current_feature_id=0;
  if($chado->cache('feature',$uniquename)){
    #seen this feature before
    $locgroup{$uniquename}++;
  }
  else {
    $chado->cache('feature',$uniquename,$nextfeature);
    $locgroup{$uniquename}       = 0;
    $current_feature_id=$nextfeature;
  }

  #if there are Targets, match types or scores and this is not ANALYSIS,
  #there is a problem.
  #
  if (!$analysis_err_str && !$ANALYSIS && (
      ((scalar($feature->annotation->get_Annotations('Target'))>0) 
         and
       ((($feature->annotation->get_Annotations('Target'))[0]->can('value') 
        && ($feature->annotation->get_Annotations('Target'))[0]->value)
           or
         (($feature->annotation->get_Annotations('Target'))[0]->can('display_text')
        && ($feature->annotation->get_Annotations('Target'))[0]->display_text)))
      or
      (defined($feature->score) and $feature->score ne '.') 
      or
      $featuretype =~ /match/ 
      ) ) {
      my @errs;
      push @errs, '* Has Target attributes'
           if (scalar($feature->annotation->get_Annotations('Target'))>0 and
              ($feature->annotation->get_Annotations('Target'))[0]->as_text);
      push @errs, '* Has scores'
           if (defined($feature->score) and $feature->score ne '.');
      push @errs, '* Has a match feature type'
           if ($featuretype =~ /match/);

      $analysis_err_str = join("\n", @errs);
      warn "\nThis file was not declared as analysis results (with the --analysis flag,\nbut this file contains attributes that imply that it is:\n$analysis_err_str\nYou can safely ignore this message if you don't need to access scores\nassociated with these features.\n\n";
  }


  if ($ANALYSIS 
      && $featuretype =~ /match/  
      && !defined($feature->annotation->get_Annotations('Target'))) {
    if (($feature->annotation->get_Annotations('ID'))[0]->can('value')) {
        $chado->cache('feature',
                     ($feature->annotation->get_Annotations('ID'))[0]->value,
                     $nextfeature);
    }
    elsif (($feature->annotation->get_Annotations('ID'))[0]->can('display_text')) {
        $chado->cache('feature',
                     ($feature->annotation->get_Annotations('ID'))[0]->display_text,
                     $nextfeature); 
    }

  }

#don't write a featureloc entry for srcfeatures
  unless ($src eq '\N' or $src == $nextfeature) {
#need to convert from base to interbase coords
    my $start = $feature->start eq '.' ? '\N' : ($feature->start - 1);
    my $end   = $feature->end   eq '.' ? '\N' : defined($feature->end) ? $feature->end : '\N';
    my $phase = ($feature->phase eq '.' or $feature->phase eq '') ? '\N' : $feature->phase;

    $chado->print_floc($nextfeatureloc, $chado->cache('feature',$uniquename), $src, $start, $end, $feature->strand, $phase,'0',$locgroup{$uniquename});
  }

  if ($feature->annotation->get_Annotations('Gap')) {
    $chado->handle_gap($feature,$uniquename);
  }


  if ($feature->annotation->get_Annotations('Note')) {
    $chado->handle_note($feature,$uniquename);
  }

#try to put unreserved tags in featureprop
#this requires that the terms exist in cvterm (and therefore that they
#have a dbxref)
  my @unreserved_tags = grep {/^[a-z]/} $feature->annotation->get_all_annotation_keys();
  if ( @unreserved_tags > 0 ) {
    $chado->handle_unreserved_tags($feature,$uniquename,@unreserved_tags);
  }

  if ( $chado->{const}{source_success} && $source && $source ne '.') {
    $chado->handle_source($feature,$uniquename,$source);
  }

  if ($feature->annotation->get_Annotations('Ontology_term')) {
    $chado->handle_ontology_term($feature,$uniquename);
  }

  if ($feature->annotation->get_Annotations('Dbxref')
       or $feature->annotation->get_Annotations('dbxref')) {
    $chado->handle_dbxref($feature,$uniquename);
  }

  my @aliases;
  if ($feature->annotation->get_Annotations('Alias')) {
    @aliases = map {$_->value} $feature->annotation->get_Annotations('Alias');
  }

  #if the uniquename was modified, put the orig ID in the alias list
  push @aliases, $chado->modified_uniquename(modified_id=>$uniquename,organism_id=>$chado->organism_id) 
      if $chado->modified_uniquename(modified_id=>$uniquename,organism_id=>$chado->organism_id);

  #Un-denormalizing the synonym table
  #if ($name ne '\N') {
  #  push @aliases, $name;
  #}
  #push @aliases, $uniquename;

  #need to unique-ify the list
  my %count;
  my @ualiases = grep {++$count{$_} < 2} @aliases;

  foreach my $alias (@ualiases) {
    $chado->synonyms($alias,$chado->cache('feature',$uniquename));
  }

  if($current_feature_id or $chado->cache('srcfeature',$nextfeature)) {
    $chado->print_f($nextfeature,$chado->organism_id,$name,$uniquename,$type,$seqlen,$chado->primary_dbxref);
  }

  if ($ANALYSIS && !defined(($feature->annotation->get_Annotations('Target'))[0])) {
    $chado->handle_nontarget_analysis($feature,$uniquename);
  }

  $chado->nextfeatureloc('++');
  #now deal with creating another feature for targets

  if (!$ANALYSIS && defined(($feature->annotation->get_Annotations('Target'))[0])) {
    die "Features in this GFF file have Target tags, but you did not indicate\n"
    ."--analysis on the command line";
  }
  elsif (defined(($feature->annotation->get_Annotations('Target'))[0])) {
      $chado->handle_target($feature, $uniquename,$name,$featuretype,$type);
  }
  $chado->nextfeature('++');
}

if ($feature_iterator = $chado->process_CDS() ) {
    $processing_CDS=1;
    goto FEATURE;
}

$chado->end_files();

#$search_uniquename->finish;
#$validate_uniquename->finish;

#deal with sequence 
unless ($NOSEQUENCE or !defined $gffio) { #ugh--reversed unless logic
  while (my $seq = $gffio->next_seq) {
    my $string = $seq->seq();
    my $name   = $seq->display_id();
    $chado->print_seq($name,$string);
  }
}

if ($FASTAFILE) {
    #use SeqIO to parse the fasta file
    my $in = Bio::SeqIO->new(-file   => $FASTAFILE,
                             -format => 'fasta');
    while (my $seq = $in->next_seq) {
        $chado->print_fasta($seq->display_id,$seq->seq);
    }
}

$chado->flush_caches();

$chado->load_data() unless $NOLOAD;

if ($END_SQL) {
  $chado->dbh->do($END_SQL);
}

$chado->remove_lock();
exit(0);


sub cleanup_handler {
    warn "@_\nAbnormal termination, trying to clean up...\n\n" if @_;  #gets the message that the die signal sent if there is one
    if ($chado && $chado->dbh->ping) {
        #$chado->dbh->{AutoCommit} = 1;
        $chado->cleanup_tmp_table;
        if ($lock) {
            warn "Trying to remove the run lock (so that --remove_lock won't be needed)...\n";
            $chado->remove_lock; #remove the lock only if we've set it
        }
        #$chado->dbh->rollback;
        print STDERR "Exiting...\n";
    }
    exit(1);
}

