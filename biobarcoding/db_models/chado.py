from biobarcoding.db_models import ORMBase


class Analysis(ORMBase):
    __tablename__ = "analysis"
    __table_args__ = {'autoload':True}

    # analysis_id = Column(BigInteger(), primary_key=True, nullable=False)
    # name = Column(String())
    # description = Column(Text())
    # program = Column(String())
    # programversion = Column(String())
    # algorithm = Column(String())
    # sourcename = Column(String())
    # sourceversion = Column(String())
    # sourceuri = Column(Text())
    # timeexecuted = Column(TIMESTAMP())


class Organism(ORMBase):
    __tablename__ = "organism"
    __table_args__ = {'autoload':True}

    # organism_id = Column(BigInteger(), primary_key=True, nullable=False)
    # abbreviation = Column(String())
    # genus = Column(Text())
    # species = Column(String())
    # common_name = Column(String())
    # infraspecific_name = Column(String())
    # type_id = Column(BigInteger())
    # comment = Column(Text())


class Feature(ORMBase):
    __tablename__ = "feature"
    __table_args__ = {'autoload':True}

    # feature_id = Column(BigInteger(), primary_key=True, nullable=False)


class AnalysisFeature(ORMBase):
    __tablename__ = "analysisfeature"
    __table_args__ = {'autoload':True}

    # analysisfeature_id = Column(BigInteger(), primary_key=True, nullable=False)
    # feature_id = Column(BigInteger(), ForeignKey('feature.feature_id'), nullable=False)
    # analysis_id = Column(BigInteger(), ForeignKey('analysis.analysis_id'), nullable=False)
    # rawscore = Column(Numeric())
    # normscore = Column(Numeric())
    # significance = Column(Numeric())
    # identity = Column(Numeric())


class Db(ORMBase):
    __tablename__ = "db"
    __table_args__ = {'autoload':True}

    # db_id = Column(BigInteger(), primary_key=True, nullable=False)
    # name = Column(String(), nullable=False)
    # description = Column(String())
    # urlprefix = Column(String())
    # url = Column(String())


# Remain

class Acquisition(ORMBase):
    __tablename__ = "acquisition"
    __table_args__ = {'autoload':True}


class AcquisitionRelationship(ORMBase):
    __tablename__ = "acquisition_relationship"
    __table_args__ = {'autoload':True}


class Acquisitionprop(ORMBase):
    __tablename__ = "acquisitionprop"
    __table_args__ = {'autoload':True}


class AnalysisCvterm(ORMBase):
    __tablename__ = "analysis_cvterm"
    __table_args__ = {'autoload':True}


class AnalysisDbxref(ORMBase):
    __tablename__ = "analysis_dbxref"
    __table_args__ = {'autoload':True}


class AnalysisPub(ORMBase):
    __tablename__ = "analysis_pub"
    __table_args__ = {'autoload':True}


class AnalysisRelationship(ORMBase):
    __tablename__ = "analysis_relationship"
    __table_args__ = {'autoload':True}


class Analysisfeatureprop(ORMBase):
    __tablename__ = "analysisfeatureprop"
    __table_args__ = {'autoload':True}


class Analysisprop(ORMBase):
    __tablename__ = "analysisprop"
    __table_args__ = {'autoload':True}


class Arraydesign(ORMBase):
    __tablename__ = "arraydesign"
    __table_args__ = {'autoload':True}


class Arraydesignprop(ORMBase):
    __tablename__ = "arraydesignprop"
    __table_args__ = {'autoload':True}


class Assay(ORMBase):
    __tablename__ = "assay"
    __table_args__ = {'autoload':True}


class AssayBiomaterial(ORMBase):
    __tablename__ = "assay_biomaterial"
    __table_args__ = {'autoload':True}


class AssayProject(ORMBase):
    __tablename__ = "assay_project"
    __table_args__ = {'autoload':True}


class Assayprop(ORMBase):
    __tablename__ = "assayprop"
    __table_args__ = {'autoload':True}


class Biomaterial(ORMBase):
    __tablename__ = "biomaterial"
    __table_args__ = {'autoload':True}


class BiomaterialDbxref(ORMBase):
    __tablename__ = "biomaterial_dbxref"
    __table_args__ = {'autoload':True}


class BiomaterialRelationship(ORMBase):
    __tablename__ = "biomaterial_relationship"
    __table_args__ = {'autoload':True}


class BiomaterialTreatment(ORMBase):
    __tablename__ = "biomaterial_treatment"
    __table_args__ = {'autoload':True}


class Biomaterialprop(ORMBase):
    __tablename__ = "biomaterialprop"
    __table_args__ = {'autoload':True}


class CellLine(ORMBase):
    __tablename__ = "cell_line"
    __table_args__ = {'autoload':True}


class CellLineCvterm(ORMBase):
    __tablename__ = "cell_line__cvterm"
    __table_args__ = {'autoload':True}


class CellLineCvtermprop(ORMBase):
    __tablename__ = "cell_line_cvtermprop"
    __table_args__ = {'autoload':True}


class CellLineDbxref(ORMBase):
    __tablename__ = "cell_line_dbxref"
    __table_args__ = {'autoload':True}


class CellLineFeature(ORMBase):
    __tablename__ = "cell_line_feature"
    __table_args__ = {'autoload':True}


class CellLineLibrary(ORMBase):
    __tablename__ = "cell_line_library"
    __table_args__ = {'autoload':True}


class CellLinePub(ORMBase):
    __tablename__ = "cell_line_pub"
    __table_args__ = {'autoload':True}


class CellLineRelationship(ORMBase):
    __tablename__ = "cell_line_relationship"
    __table_args__ = {'autoload':True}


class CellLineSynonym(ORMBase):
    __tablename__ = "cell_line_synonym"
    __table_args__ = {'autoload':True}


class CellLineprop(ORMBase):
    __tablename__ = "cell_lineprop"
    __table_args__ = {'autoload':True}


class CellLinepropPub(ORMBase):
    __tablename__ = "cell_lineprop_pub"
    __table_args__ = {'autoload':True}


class Chadoprop(ORMBase):
    __tablename__ = "chadoprop"
    __table_args__ = {'autoload':True}


class Channel(ORMBase):
    __tablename__ = "channel"
    __table_args__ = {'autoload':True}


class Contact(ORMBase):
    __tablename__ = "contact"
    __table_args__ = {'autoload':True}


class ContactRelationship(ORMBase):
    __tablename__ = "contact_relationship"
    __table_args__ = {'autoload':True}


class Contactprop(ORMBase):
    __tablename__ = "contactprop"
    __table_args__ = {'autoload':True}


class Control(ORMBase):
    __tablename__ = "control"
    __table_args__ = {'autoload':True}


class Cv(ORMBase):
    __tablename__ = "cv"
    __table_args__ = {'autoload':True}


class Cvprop(ORMBase):
    __tablename__ = "cvprop"
    __table_args__ = {'autoload':True}


class Cvterm(ORMBase):
    __tablename__ = "cvterm"
    __table_args__ = {'autoload':True}


class CvtermDbxref(ORMBase):
    __tablename__ = "cvterm_dbxref"
    __table_args__ = {'autoload':True}


class CvtermRelationship(ORMBase):
    __tablename__ = "cvterm_relationship"
    __table_args__ = {'autoload':True}


class Cvtermpath(ORMBase):
    __tablename__ = "cvtermpath"
    __table_args__ = {'autoload':True}


class Cvtermprop(ORMBase):
    __tablename__ = "cvtermprop"
    __table_args__ = {'autoload':True}


class Cvtermsynonym(ORMBase):
    __tablename__ = "cvtermsynonym"
    __table_args__ = {'autoload':True}


class Dbprop(ORMBase):
    __tablename__ = "dbprop"
    __table_args__ = {'autoload':True}


class Dbxref(ORMBase):
    __tablename__ = "dbxref"
    __table_args__ = {'autoload':True}


class Dbxrefprop(ORMBase):
    __tablename__ = "dbxrefprop"
    __table_args__ = {'autoload':True}


class Eimage(ORMBase):
    __tablename__ = "eimage"
    __table_args__ = {'autoload':True}


class Element(ORMBase):
    __tablename__ = "element"
    __table_args__ = {'autoload':True}


class ElementRelationship(ORMBase):
    __tablename__ = "element_relationship"
    __table_args__ = {'autoload':True}


class Elementresult(ORMBase):
    __tablename__ = "elementresult"
    __table_args__ = {'autoload':True}


class ElementresultRelationship(ORMBase):
    __tablename__ = "elementresult_relationship"
    __table_args__ = {'autoload':True}


class Environment(ORMBase):
    __tablename__ = "environment"
    __table_args__ = {'autoload':True}


class EnvironmentCvterm(ORMBase):
    __tablename__ = "environment_cvterm"
    __table_args__ = {'autoload':True}


class Expression(ORMBase):
    __tablename__ = "expression"
    __table_args__ = {'autoload':True}


class ExpressionCvterm(ORMBase):
    __tablename__ = "expression_cvterm"
    __table_args__ = {'autoload':True}


class ExpressionCvtermprop(ORMBase):
    __tablename__ = "expression_cvtermprop"
    __table_args__ = {'autoload':True}


class ExpressionImage(ORMBase):
    __tablename__ = "expression_image"
    __table_args__ = {'autoload':True}


class ExpressionPub(ORMBase):
    __tablename__ = "expression_pub"
    __table_args__ = {'autoload':True}


class Expressionprop(ORMBase):
    __tablename__ = "expressionprop"
    __table_args__ = {'autoload':True}


class FeatureContact(ORMBase):
    __tablename__ = "feature_contact"
    __table_args__ = {'autoload':True}


class FeatureCvterm(ORMBase):
    __tablename__ = "feature_cvterm"
    __table_args__ = {'autoload':True}


class FeatureCvtermDbxref(ORMBase):
    __tablename__ = "feature_cvterm_dbxref"
    __table_args__ = {'autoload':True}


class FeatureCvtermPub(ORMBase):
    __tablename__ = "feature_cvterm_pub"
    __table_args__ = {'autoload':True}


class FeatureCvtermprop(ORMBase):
    __tablename__ = "feature_cvtermprop"
    __table_args__ = {'autoload':True}


class FeatureDbxref(ORMBase):
    __tablename__ = "feature_dbxref"
    __table_args__ = {'autoload':True}


class FeatureExpression(ORMBase):
    __tablename__ = "feature_expression"
    __table_args__ = {'autoload':True}


class FeatureExpressionprop(ORMBase):
    __tablename__ = "feature_expressionprop"
    __table_args__ = {'autoload':True}


class FeatureGenotype(ORMBase):
    __tablename__ = "feature_genotype"
    __table_args__ = {'autoload':True}


class FeaturePhenotype(ORMBase):
    __tablename__ = "feature_phenotype"
    __table_args__ = {'autoload':True}


class FeaturePub(ORMBase):
    __tablename__ = "feature_pub"
    __table_args__ = {'autoload':True}


class FeaturePubprop(ORMBase):
    __tablename__ = "feature_pubprop"
    __table_args__ = {'autoload':True}


class FeatureRelationship(ORMBase):
    __tablename__ = "feature_relationship"
    __table_args__ = {'autoload':True}


class FeatureRelationshipPub(ORMBase):
    __tablename__ = "feature_relationship_pub"
    __table_args__ = {'autoload':True}


class FeatureRelationshipprop(ORMBase):
    __tablename__ = "feature_relationshipprop"
    __table_args__ = {'autoload':True}


class FeatureRelationshippropPub(ORMBase):
    __tablename__ = "feature_relationshipprop_pub"
    __table_args__ = {'autoload':True}


class FeatureSynonym(ORMBase):
    __tablename__ = "feature_synonym"
    __table_args__ = {'autoload':True}


class Featureloc(ORMBase):
    __tablename__ = "featureloc"
    __table_args__ = {'autoload':True}


class FeaturelocPub(ORMBase):
    __tablename__ = "featureloc_pub"
    __table_args__ = {'autoload':True}


class Featuremap(ORMBase):
    __tablename__ = "featuremap"
    __table_args__ = {'autoload':True}


class FeaturemapContact(ORMBase):
    __tablename__ = "featuremap_contact"
    __table_args__ = {'autoload':True}


class FeaturemapDbxref(ORMBase):
    __tablename__ = "featuremap_dbxref"
    __table_args__ = {'autoload':True}


class FeaturemapOrganism(ORMBase):
    __tablename__ = "featuremap_organism"
    __table_args__ = {'autoload':True}


class FeaturemapPub(ORMBase):
    __tablename__ = "featuremap_pub"
    __table_args__ = {'autoload':True}


class Featuremapprop(ORMBase):
    __tablename__ = "featuremapprop"
    __table_args__ = {'autoload':True}


class Featurepos(ORMBase):
    __tablename__ = "featurepos"
    __table_args__ = {'autoload':True}


class Featureposprop(ORMBase):
    __tablename__ = "featureposprop"
    __table_args__ = {'autoload':True}


class Featureprop(ORMBase):
    __tablename__ = "featureprop"
    __table_args__ = {'autoload':True}


class FeaturepropPub(ORMBase):
    __tablename__ = "featureprop_pub"
    __table_args__ = {'autoload':True}


class Featurerange(ORMBase):
    __tablename__ = "featurerange"
    __table_args__ = {'autoload':True}


class Genotype(ORMBase):
    __tablename__ = "genotype"
    __table_args__ = {'autoload':True}


class Genotypeprop(ORMBase):
    __tablename__ = "genotypeprop"
    __table_args__ = {'autoload':True}


class Library(ORMBase):
    __tablename__ = "library"
    __table_args__ = {'autoload':True}


class LibraryContact(ORMBase):
    __tablename__ = "library_contact"
    __table_args__ = {'autoload':True}


class LibraryCvterm(ORMBase):
    __tablename__ = "library_cvterm"
    __table_args__ = {'autoload':True}


class LibraryDbxref(ORMBase):
    __tablename__ = "library_dbxref"
    __table_args__ = {'autoload':True}


class LibraryExpression(ORMBase):
    __tablename__ = "library_expression"
    __table_args__ = {'autoload':True}


class LibraryExpressionprop(ORMBase):
    __tablename__ = "library_expressionprop"
    __table_args__ = {'autoload':True}


class LibraryFeature(ORMBase):
    __tablename__ = "library_feature"
    __table_args__ = {'autoload':True}


class LibraryFeatureprop(ORMBase):
    __tablename__ = "library_featureprop"
    __table_args__ = {'autoload':True}


class LibraryPub(ORMBase):
    __tablename__ = "library_pub"
    __table_args__ = {'autoload':True}


class LibraryRelationship(ORMBase):
    __tablename__ = "library_relationship"
    __table_args__ = {'autoload':True}


class LibraryRelationshipPub(ORMBase):
    __tablename__ = "library_relationship_pub"
    __table_args__ = {'autoload':True}


class LibrarySynonym(ORMBase):
    __tablename__ = "library_synonym"
    __table_args__ = {'autoload':True}


class Libraryprop(ORMBase):
    __tablename__ = "libraryprop"
    __table_args__ = {'autoload':True}


class LibrarypropPub(ORMBase):
    __tablename__ = "libraryprop_pub"
    __table_args__ = {'autoload':True}


class Magedocumentation(ORMBase):
    __tablename__ = "magedocumentation"
    __table_args__ = {'autoload':True}


class Mageml(ORMBase):
    __tablename__ = "mageml"
    __table_args__ = {'autoload':True}


class NdExperiment(ORMBase):
    __tablename__ = "nd_experiment"
    __table_args__ = {'autoload':True}


class NdExperimentAnalysis(ORMBase):
    __tablename__ = "nd_experiment_analysis"
    __table_args__ = {'autoload':True}


class NdExperimentContact(ORMBase):
    __tablename__ = "nd_experiment_contact"
    __table_args__ = {'autoload':True}


class NdExperimentDbxref(ORMBase):
    __tablename__ = "nd_experiment_dbxref"
    __table_args__ = {'autoload':True}


class NdExperimentGenotype(ORMBase):
    __tablename__ = "nd_experiment_genotype"
    __table_args__ = {'autoload':True}


class NdExperimentPhenotype(ORMBase):
    __tablename__ = "nd_experiment_phenotype"
    __table_args__ = {'autoload':True}


class NdExperimentProject(ORMBase):
    __tablename__ = "nd_experiment_project"
    __table_args__ = {'autoload':True}


class NdExperimentProtocol(ORMBase):
    __tablename__ = "nd_experiment_protocol"
    __table_args__ = {'autoload':True}


class NdExperimentPub(ORMBase):
    __tablename__ = "nd_experiment_pub"
    __table_args__ = {'autoload':True}


class NdExperimentStock(ORMBase):
    __tablename__ = "nd_experiment_stock"
    __table_args__ = {'autoload':True}


class NdExperimentStockDbxref(ORMBase):
    __tablename__ = "nd_experiment_stock_dbxref"
    __table_args__ = {'autoload':True}


class NdExperimentStockprop(ORMBase):
    __tablename__ = "nd_experiment_stockprop"
    __table_args__ = {'autoload':True}


class NdExperimentprop(ORMBase):
    __tablename__ = "nd_experimentprop"
    __table_args__ = {'autoload':True}


class NdGeolocation(ORMBase):
    __tablename__ = "nd_geolocation"
    __table_args__ = {'autoload':True}


class NdGeolocationprop(ORMBase):
    __tablename__ = "nd_geolocationprop"
    __table_args__ = {'autoload':True}


class NdProtocol(ORMBase):
    __tablename__ = "nd_protocol"
    __table_args__ = {'autoload':True}


class NdProtocolReagent(ORMBase):
    __tablename__ = "nd_protocol_reagent"
    __table_args__ = {'autoload':True}


class NdProtocolprop(ORMBase):
    __tablename__ = "nd_protocolprop"
    __table_args__ = {'autoload':True}


class NdReagent(ORMBase):
    __tablename__ = "nd_reagent"
    __table_args__ = {'autoload':True}


class NdReagentRelationship(ORMBase):
    __tablename__ = "nd_reagent_relationship"
    __table_args__ = {'autoload':True}


class NdReagentprop(ORMBase):
    __tablename__ = "nd_reagentprop"
    __table_args__ = {'autoload':True}


class OrganismCvterm(ORMBase):
    __tablename__ = "organism_cvterm"
    __table_args__ = {'autoload':True}


class OrganismCvtermprop(ORMBase):
    __tablename__ = "organism_cvtermprop"
    __table_args__ = {'autoload':True}


class OrganismDbxref(ORMBase):
    __tablename__ = "organism_dbxref"
    __table_args__ = {'autoload':True}


class OrganismPub(ORMBase):
    __tablename__ = "organism_pub"
    __table_args__ = {'autoload':True}


class OrganismRelationship(ORMBase):
    __tablename__ = "organism_relationship"
    __table_args__ = {'autoload':True}


class Organismprop(ORMBase):
    __tablename__ = "organismprop"
    __table_args__ = {'autoload':True}


class OrganismpropPub(ORMBase):
    __tablename__ = "organismprop_pub"
    __table_args__ = {'autoload':True}


class Phendesc(ORMBase):
    __tablename__ = "phendesc"
    __table_args__ = {'autoload':True}


class Phenotype(ORMBase):
    __tablename__ = "phenotype"
    __table_args__ = {'autoload':True}


class PhenotypeComparison(ORMBase):
    __tablename__ = "phenotype_comparison"
    __table_args__ = {'autoload':True}


class PhenotypeComparisonCvterm(ORMBase):
    __tablename__ = "phenotype_comparison_cvterm"
    __table_args__ = {'autoload':True}


class PhenotypeCvterm(ORMBase):
    __tablename__ = "phenotype_cvterm"
    __table_args__ = {'autoload':True}


class Phenotypeprop(ORMBase):
    __tablename__ = "phenotypeprop"
    __table_args__ = {'autoload':True}


class Phenstatement(ORMBase):
    __tablename__ = "phenstatement"
    __table_args__ = {'autoload':True}


class Phylonode(ORMBase):
    __tablename__ = "phylonode"
    __table_args__ = {'autoload':True}


class PhylonodeDbxref(ORMBase):
    __tablename__ = "phylonode_dbxref"
    __table_args__ = {'autoload':True}


class PhylonodeOrganism(ORMBase):
    __tablename__ = "phylonode_organism"
    __table_args__ = {'autoload':True}


class PhylonodePub(ORMBase):
    __tablename__ = "phylonode_pub"
    __table_args__ = {'autoload':True}


class PhylonodeRelationship(ORMBase):
    __tablename__ = "phylonode_relationship"
    __table_args__ = {'autoload':True}


class Phylonodeprop(ORMBase):
    __tablename__ = "phylonodeprop"
    __table_args__ = {'autoload':True}


class Phylotree(ORMBase):
    __tablename__ = "phylotree"
    __table_args__ = {'autoload':True}


class PhylotreePub(ORMBase):
    __tablename__ = "phylotree_pub"
    __table_args__ = {'autoload':True}


class Phylotreeprop(ORMBase):
    __tablename__ = "phylotreeprop"
    __table_args__ = {'autoload':True}


class Project(ORMBase):
    __tablename__ = "project"
    __table_args__ = {'autoload':True}


class ProjectAnalysis(ORMBase):
    __tablename__ = "projectanalysis"
    __table_args__ = {'autoload':True}


class ProjectContact(ORMBase):
    __tablename__ = "project_contact"
    __table_args__ = {'autoload':True}


class ProjectDbxref(ORMBase):
    __tablename__ = "project_dbxref"
    __table_args__ = {'autoload':True}


class ProjectFeature(ORMBase):
    __tablename__ = "project_feature"
    __table_args__ = {'autoload':True}


class ProjectPub(ORMBase):
    __tablename__ = "project_pub"
    __table_args__ = {'autoload':True}


class ProjectRelationship(ORMBase):
    __tablename__ = "project_relationship"
    __table_args__ = {'autoload':True}


class ProjectStock(ORMBase):
    __tablename__ = "project_stock"
    __table_args__ = {'autoload':True}


class Projectprop(ORMBase):
    __tablename__ = "projectprop"
    __table_args__ = {'autoload':True}


class Protocol(ORMBase):
    __tablename__ = "protocol"
    __table_args__ = {'autoload':True}


class Protocolparam(ORMBase):
    __tablename__ = "protocolparam"
    __table_args__ = {'autoload':True}


class Pub(ORMBase):
    __tablename__ = "pub"
    __table_args__ = {'autoload':True}


class PubDbxref(ORMBase):
    __tablename__ = "pub_dbxref"
    __table_args__ = {'autoload':True}


class PubRelationship(ORMBase):
    __tablename__ = "pub_relationship"
    __table_args__ = {'autoload':True}


class Pubauthor(ORMBase):
    __tablename__ = "pubauthor"
    __table_args__ = {'autoload':True}


class PubauthorContact(ORMBase):
    __tablename__ = "pubauthor_contact"
    __table_args__ = {'autoload':True}


class Pubprop(ORMBase):
    __tablename__ = "pubprop"
    __table_args__ = {'autoload':True}


class Quantification(ORMBase):
    __tablename__ = "quantification"
    __table_args__ = {'autoload':True}


class QuantificationRelationship(ORMBase):
    __tablename__ = "quantification_relationship"
    __table_args__ = {'autoload':True}


class Quantificationprop(ORMBase):
    __tablename__ = "quantificationprop"
    __table_args__ = {'autoload':True}


class Stock(ORMBase):
    __tablename__ = "stock"
    __table_args__ = {'autoload':True}


class StockCvterm(ORMBase):
    __tablename__ = "stock_cvterm"
    __table_args__ = {'autoload':True}


class StockCvtermprop(ORMBase):
    __tablename__ = "stock_cvtermprop"
    __table_args__ = {'autoload':True}


class StockDbxref(ORMBase):
    __tablename__ = "stock_dbxref"
    __table_args__ = {'autoload':True}


class StockDbxrefprop(ORMBase):
    __tablename__ = "stock_dbxrefprop"
    __table_args__ = {'autoload':True}


class StockFeature(ORMBase):
    __tablename__ = "stock_feature"
    __table_args__ = {'autoload':True}


class StockFeaturemap(ORMBase):
    __tablename__ = "stock_featuremap"
    __table_args__ = {'autoload':True}


class StockGenotype(ORMBase):
    __tablename__ = "stock_genotype"
    __table_args__ = {'autoload':True}


class StockLibrary(ORMBase):
    __tablename__ = "stock_library"
    __table_args__ = {'autoload':True}


class StockPub(ORMBase):
    __tablename__ = "stock_pub"
    __table_args__ = {'autoload':True}


class StockRelationship(ORMBase):
    __tablename__ = "stock_relationship"
    __table_args__ = {'autoload':True}


class StockRelationshipCvterm(ORMBase):
    __tablename__ = "stock_relationship_cvterm"
    __table_args__ = {'autoload':True}


class StockRelationshipPub(ORMBase):
    __tablename__ = "stock_relationship_pub"
    __table_args__ = {'autoload':True}


class Stockcollection(ORMBase):
    __tablename__ = "stockcollection"
    __table_args__ = {'autoload':True}


class StockcollectionDb(ORMBase):
    __tablename__ = "stockcollection_db"
    __table_args__ = {'autoload':True}


class StockcollectionStock(ORMBase):
    __tablename__ = "stockcollection_stock"
    __table_args__ = {'autoload':True}


class Stockcollectionprop(ORMBase):
    __tablename__ = "stockcollectionprop"
    __table_args__ = {'autoload':True}


class Stockprop(ORMBase):
    __tablename__ = "stockprop"
    __table_args__ = {'autoload':True}


class StockpropPub(ORMBase):
    __tablename__ = "stockprop_pub"
    __table_args__ = {'autoload':True}


class Study(ORMBase):
    __tablename__ = "study"
    __table_args__ = {'autoload':True}


class StudyAssay(ORMBase):
    __tablename__ = "study_assay"
    __table_args__ = {'autoload':True}


class Studydesign(ORMBase):
    __tablename__ = "studydesign"
    __table_args__ = {'autoload':True}


class Studydesignprop(ORMBase):
    __tablename__ = "studydesignprop"
    __table_args__ = {'autoload':True}


class Studyfactor(ORMBase):
    __tablename__ = "studyfactor"
    __table_args__ = {'autoload':True}


class Studyfactorvalue(ORMBase):
    __tablename__ = "studyfactorvalue"
    __table_args__ = {'autoload':True}


class Studyprop(ORMBase):
    __tablename__ = "studyprop"
    __table_args__ = {'autoload':True}


class StudypropFeature(ORMBase):
    __tablename__ = "studyprop_feature"
    __table_args__ = {'autoload':True}


class Synonym(ORMBase):
    __tablename__ = "synonym"
    __table_args__ = {'autoload':True}


class Tableinfo(ORMBase):
    __tablename__ = "tableinfo"
    __table_args__ = {'autoload':True}


class Treatment(ORMBase):
    __tablename__ = "treatment"
    __table_args__ = {'autoload':True}
