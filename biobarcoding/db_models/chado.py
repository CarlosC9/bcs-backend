from sqlalchemy.ext.hybrid import hybrid_property
from . import ORMBaseChado


class Organism(ORMBaseChado):
    __tablename__ = "organism"
    __table_args__ = {'autoload': True}

    @hybrid_property
    def name(self):
        return Organism.genus + ' ' + Organism.species


class Acquisition(ORMBaseChado):
    __tablename__ = "acquisition"
    __table_args__ = {'autoload': True}


class AcquisitionRelationship(ORMBaseChado):
    __tablename__ = "acquisition_relationship"
    __table_args__ = {'autoload': True}


class Acquisitionprop(ORMBaseChado):
    __tablename__ = "acquisitionprop"
    __table_args__ = {'autoload': True}


class Analysis(ORMBaseChado):
    __tablename__ = "analysis"
    __table_args__ = {'autoload': True}


class AnalysisCvterm(ORMBaseChado):
    __tablename__ = "analysis_cvterm"
    __table_args__ = {'autoload': True}


class AnalysisDbxref(ORMBaseChado):
    __tablename__ = "analysis_dbxref"
    __table_args__ = {'autoload': True}


class AnalysisFeature(ORMBaseChado):
    __tablename__ = "analysisfeature"
    __table_args__ = {'autoload': True}


class AnalysisPub(ORMBaseChado):
    __tablename__ = "analysis_pub"
    __table_args__ = {'autoload': True}


class AnalysisRelationship(ORMBaseChado):
    __tablename__ = "analysis_relationship"
    __table_args__ = {'autoload': True}


class Analysisfeatureprop(ORMBaseChado):
    __tablename__ = "analysisfeatureprop"
    __table_args__ = {'autoload': True}


class Analysisprop(ORMBaseChado):
    __tablename__ = "analysisprop"
    __table_args__ = {'autoload': True}


class Arraydesign(ORMBaseChado):
    __tablename__ = "arraydesign"
    __table_args__ = {'autoload': True}


class Arraydesignprop(ORMBaseChado):
    __tablename__ = "arraydesignprop"
    __table_args__ = {'autoload': True}


class Assay(ORMBaseChado):
    __tablename__ = "assay"
    __table_args__ = {'autoload': True}


class AssayBiomaterial(ORMBaseChado):
    __tablename__ = "assay_biomaterial"
    __table_args__ = {'autoload': True}


class AssayProject(ORMBaseChado):
    __tablename__ = "assay_project"
    __table_args__ = {'autoload': True}


class Assayprop(ORMBaseChado):
    __tablename__ = "assayprop"
    __table_args__ = {'autoload': True}


class Biomaterial(ORMBaseChado):
    __tablename__ = "biomaterial"
    __table_args__ = {'autoload': True}


class BiomaterialDbxref(ORMBaseChado):
    __tablename__ = "biomaterial_dbxref"
    __table_args__ = {'autoload': True}


class BiomaterialRelationship(ORMBaseChado):
    __tablename__ = "biomaterial_relationship"
    __table_args__ = {'autoload': True}


class BiomaterialTreatment(ORMBaseChado):
    __tablename__ = "biomaterial_treatment"
    __table_args__ = {'autoload': True}


class Biomaterialprop(ORMBaseChado):
    __tablename__ = "biomaterialprop"
    __table_args__ = {'autoload': True}


class CellLine(ORMBaseChado):
    __tablename__ = "cell_line"
    __table_args__ = {'autoload': True}


class CellLineCvterm(ORMBaseChado):
    __tablename__ = "cell_line_cvterm"
    __table_args__ = {'autoload': True}


class CellLineCvtermprop(ORMBaseChado):
    __tablename__ = "cell_line_cvtermprop"
    __table_args__ = {'autoload': True}


class CellLineDbxref(ORMBaseChado):
    __tablename__ = "cell_line_dbxref"
    __table_args__ = {'autoload': True}


class CellLineFeature(ORMBaseChado):
    __tablename__ = "cell_line_feature"
    __table_args__ = {'autoload': True}


class CellLineLibrary(ORMBaseChado):
    __tablename__ = "cell_line_library"
    __table_args__ = {'autoload': True}


class CellLinePub(ORMBaseChado):
    __tablename__ = "cell_line_pub"
    __table_args__ = {'autoload': True}


class CellLineRelationship(ORMBaseChado):
    __tablename__ = "cell_line_relationship"
    __table_args__ = {'autoload': True}


class CellLineSynonym(ORMBaseChado):
    __tablename__ = "cell_line_synonym"
    __table_args__ = {'autoload': True}


class CellLineprop(ORMBaseChado):
    __tablename__ = "cell_lineprop"
    __table_args__ = {'autoload': True}


class CellLinepropPub(ORMBaseChado):
    __tablename__ = "cell_lineprop_pub"
    __table_args__ = {'autoload': True}


class Chadoprop(ORMBaseChado):
    __tablename__ = "chadoprop"
    __table_args__ = {'autoload': True}


class Channel(ORMBaseChado):
    __tablename__ = "channel"
    __table_args__ = {'autoload': True}


class Contact(ORMBaseChado):
    __tablename__ = "contact"
    __table_args__ = {'autoload': True}


class ContactRelationship(ORMBaseChado):
    __tablename__ = "contact_relationship"
    __table_args__ = {'autoload': True}


class Contactprop(ORMBaseChado):
    __tablename__ = "contactprop"
    __table_args__ = {'autoload': True}


class Control(ORMBaseChado):
    __tablename__ = "control"
    __table_args__ = {'autoload': True}


class Cv(ORMBaseChado):
    __tablename__ = "cv"
    __table_args__ = {'autoload': True}


class Cvprop(ORMBaseChado):
    __tablename__ = "cvprop"
    __table_args__ = {'autoload': True}


class Cvterm(ORMBaseChado):
    __tablename__ = "cvterm"
    __table_args__ = {'autoload': True}


class CvtermDbxref(ORMBaseChado):
    __tablename__ = "cvterm_dbxref"
    __table_args__ = {'autoload': True}


class CvtermRelationship(ORMBaseChado):
    __tablename__ = "cvterm_relationship"
    __table_args__ = {'autoload': True}


class Cvtermpath(ORMBaseChado):
    __tablename__ = "cvtermpath"
    __table_args__ = {'autoload': True}


class Cvtermprop(ORMBaseChado):
    __tablename__ = "cvtermprop"
    __table_args__ = {'autoload': True}


class Cvtermsynonym(ORMBaseChado):
    __tablename__ = "cvtermsynonym"
    __table_args__ = {'autoload': True}


class Db(ORMBaseChado):
    __tablename__ = "db"
    __table_args__ = {'autoload': True}


class Dbprop(ORMBaseChado):
    __tablename__ = "dbprop"
    __table_args__ = {'autoload': True}


class Dbxref(ORMBaseChado):
    __tablename__ = "dbxref"
    __table_args__ = {'autoload': True}


class Dbxrefprop(ORMBaseChado):
    __tablename__ = "dbxrefprop"
    __table_args__ = {'autoload': True}


class Eimage(ORMBaseChado):
    __tablename__ = "eimage"
    __table_args__ = {'autoload': True}


class Element(ORMBaseChado):
    __tablename__ = "element"
    __table_args__ = {'autoload': True}


class ElementRelationship(ORMBaseChado):
    __tablename__ = "element_relationship"
    __table_args__ = {'autoload': True}


class Elementresult(ORMBaseChado):
    __tablename__ = "elementresult"
    __table_args__ = {'autoload': True}


class ElementresultRelationship(ORMBaseChado):
    __tablename__ = "elementresult_relationship"
    __table_args__ = {'autoload': True}


class Environment(ORMBaseChado):
    __tablename__ = "environment"
    __table_args__ = {'autoload': True}


class EnvironmentCvterm(ORMBaseChado):
    __tablename__ = "environment_cvterm"
    __table_args__ = {'autoload': True}


class Expression(ORMBaseChado):
    __tablename__ = "expression"
    __table_args__ = {'autoload': True}


class ExpressionCvterm(ORMBaseChado):
    __tablename__ = "expression_cvterm"
    __table_args__ = {'autoload': True}


class ExpressionCvtermprop(ORMBaseChado):
    __tablename__ = "expression_cvtermprop"
    __table_args__ = {'autoload': True}


class ExpressionImage(ORMBaseChado):
    __tablename__ = "expression_image"
    __table_args__ = {'autoload': True}


class ExpressionPub(ORMBaseChado):
    __tablename__ = "expression_pub"
    __table_args__ = {'autoload': True}


class Expressionprop(ORMBaseChado):
    __tablename__ = "expressionprop"
    __table_args__ = {'autoload': True}


class Feature(ORMBaseChado):
    __tablename__ = "feature"
    __table_args__ = {'autoload': True}


class FeatureContact(ORMBaseChado):
    __tablename__ = "feature_contact"
    __table_args__ = {'autoload': True}


class FeatureCvterm(ORMBaseChado):
    __tablename__ = "feature_cvterm"
    __table_args__ = {'autoload': True}


class FeatureCvtermDbxref(ORMBaseChado):
    __tablename__ = "feature_cvterm_dbxref"
    __table_args__ = {'autoload': True}


class FeatureCvtermPub(ORMBaseChado):
    __tablename__ = "feature_cvterm_pub"
    __table_args__ = {'autoload': True}


class FeatureCvtermprop(ORMBaseChado):
    __tablename__ = "feature_cvtermprop"
    __table_args__ = {'autoload': True}


class FeatureDbxref(ORMBaseChado):
    __tablename__ = "feature_dbxref"
    __table_args__ = {'autoload': True}


class FeatureExpression(ORMBaseChado):
    __tablename__ = "feature_expression"
    __table_args__ = {'autoload': True}


class FeatureExpressionprop(ORMBaseChado):
    __tablename__ = "feature_expressionprop"
    __table_args__ = {'autoload': True}


class FeatureGenotype(ORMBaseChado):
    __tablename__ = "feature_genotype"
    __table_args__ = {'autoload': True}


class FeaturePhenotype(ORMBaseChado):
    __tablename__ = "feature_phenotype"
    __table_args__ = {'autoload': True}


class FeaturePub(ORMBaseChado):
    __tablename__ = "feature_pub"
    __table_args__ = {'autoload': True}


class FeaturePubprop(ORMBaseChado):
    __tablename__ = "feature_pubprop"
    __table_args__ = {'autoload': True}


class FeatureRelationship(ORMBaseChado):
    __tablename__ = "feature_relationship"
    __table_args__ = {'autoload': True}


class FeatureRelationshipPub(ORMBaseChado):
    __tablename__ = "feature_relationship_pub"
    __table_args__ = {'autoload': True}


class FeatureRelationshipprop(ORMBaseChado):
    __tablename__ = "feature_relationshipprop"
    __table_args__ = {'autoload': True}


class FeatureRelationshippropPub(ORMBaseChado):
    __tablename__ = "feature_relationshipprop_pub"
    __table_args__ = {'autoload': True}


class FeatureSynonym(ORMBaseChado):
    __tablename__ = "feature_synonym"
    __table_args__ = {'autoload': True}


class Featureloc(ORMBaseChado):
    __tablename__ = "featureloc"
    __table_args__ = {'autoload': True}


class FeaturelocPub(ORMBaseChado):
    __tablename__ = "featureloc_pub"
    __table_args__ = {'autoload': True}


class Featuremap(ORMBaseChado):
    __tablename__ = "featuremap"
    __table_args__ = {'autoload': True}


class FeaturemapContact(ORMBaseChado):
    __tablename__ = "featuremap_contact"
    __table_args__ = {'autoload': True}


class FeaturemapDbxref(ORMBaseChado):
    __tablename__ = "featuremap_dbxref"
    __table_args__ = {'autoload': True}


class FeaturemapOrganism(ORMBaseChado):
    __tablename__ = "featuremap_organism"
    __table_args__ = {'autoload': True}


class FeaturemapPub(ORMBaseChado):
    __tablename__ = "featuremap_pub"
    __table_args__ = {'autoload': True}


class Featuremapprop(ORMBaseChado):
    __tablename__ = "featuremapprop"
    __table_args__ = {'autoload': True}


class Featurepos(ORMBaseChado):
    __tablename__ = "featurepos"
    __table_args__ = {'autoload': True}


class Featureposprop(ORMBaseChado):
    __tablename__ = "featureposprop"
    __table_args__ = {'autoload': True}


class Featureprop(ORMBaseChado):
    __tablename__ = "featureprop"
    __table_args__ = {'autoload': True}


class FeaturepropPub(ORMBaseChado):
    __tablename__ = "featureprop_pub"
    __table_args__ = {'autoload': True}


class Featurerange(ORMBaseChado):
    __tablename__ = "featurerange"
    __table_args__ = {'autoload': True}


class Genotype(ORMBaseChado):
    __tablename__ = "genotype"
    __table_args__ = {'autoload': True}


class Genotypeprop(ORMBaseChado):
    __tablename__ = "genotypeprop"
    __table_args__ = {'autoload': True}


class Library(ORMBaseChado):
    __tablename__ = "library"
    __table_args__ = {'autoload': True}


class LibraryContact(ORMBaseChado):
    __tablename__ = "library_contact"
    __table_args__ = {'autoload': True}


class LibraryCvterm(ORMBaseChado):
    __tablename__ = "library_cvterm"
    __table_args__ = {'autoload': True}


class LibraryDbxref(ORMBaseChado):
    __tablename__ = "library_dbxref"
    __table_args__ = {'autoload': True}


class LibraryExpression(ORMBaseChado):
    __tablename__ = "library_expression"
    __table_args__ = {'autoload': True}


class LibraryExpressionprop(ORMBaseChado):
    __tablename__ = "library_expressionprop"
    __table_args__ = {'autoload': True}


class LibraryFeature(ORMBaseChado):
    __tablename__ = "library_feature"
    __table_args__ = {'autoload': True}


class LibraryFeatureprop(ORMBaseChado):
    __tablename__ = "library_featureprop"
    __table_args__ = {'autoload': True}


class LibraryPub(ORMBaseChado):
    __tablename__ = "library_pub"
    __table_args__ = {'autoload': True}


class LibraryRelationship(ORMBaseChado):
    __tablename__ = "library_relationship"
    __table_args__ = {'autoload': True}


class LibraryRelationshipPub(ORMBaseChado):
    __tablename__ = "library_relationship_pub"
    __table_args__ = {'autoload': True}


class LibrarySynonym(ORMBaseChado):
    __tablename__ = "library_synonym"
    __table_args__ = {'autoload': True}


class Libraryprop(ORMBaseChado):
    __tablename__ = "libraryprop"
    __table_args__ = {'autoload': True}


class LibrarypropPub(ORMBaseChado):
    __tablename__ = "libraryprop_pub"
    __table_args__ = {'autoload': True}


class Magedocumentation(ORMBaseChado):
    __tablename__ = "magedocumentation"
    __table_args__ = {'autoload': True}


class Mageml(ORMBaseChado):
    __tablename__ = "mageml"
    __table_args__ = {'autoload': True}


class NdExperiment(ORMBaseChado):
    __tablename__ = "nd_experiment"
    __table_args__ = {'autoload': True}


class NdExperimentAnalysis(ORMBaseChado):
    __tablename__ = "nd_experiment_analysis"
    __table_args__ = {'autoload': True}


class NdExperimentContact(ORMBaseChado):
    __tablename__ = "nd_experiment_contact"
    __table_args__ = {'autoload': True}


class NdExperimentDbxref(ORMBaseChado):
    __tablename__ = "nd_experiment_dbxref"
    __table_args__ = {'autoload': True}


class NdExperimentGenotype(ORMBaseChado):
    __tablename__ = "nd_experiment_genotype"
    __table_args__ = {'autoload': True}


class NdExperimentPhenotype(ORMBaseChado):
    __tablename__ = "nd_experiment_phenotype"
    __table_args__ = {'autoload': True}


class NdExperimentProject(ORMBaseChado):
    __tablename__ = "nd_experiment_project"
    __table_args__ = {'autoload': True}


class NdExperimentProtocol(ORMBaseChado):
    __tablename__ = "nd_experiment_protocol"
    __table_args__ = {'autoload': True}


class NdExperimentPub(ORMBaseChado):
    __tablename__ = "nd_experiment_pub"
    __table_args__ = {'autoload': True}


class NdExperimentStock(ORMBaseChado):
    __tablename__ = "nd_experiment_stock"
    __table_args__ = {'autoload': True}


class NdExperimentStockDbxref(ORMBaseChado):
    __tablename__ = "nd_experiment_stock_dbxref"
    __table_args__ = {'autoload': True}


class NdExperimentStockprop(ORMBaseChado):
    __tablename__ = "nd_experiment_stockprop"
    __table_args__ = {'autoload': True}


class NdExperimentprop(ORMBaseChado):
    __tablename__ = "nd_experimentprop"
    __table_args__ = {'autoload': True}


class NdGeolocation(ORMBaseChado):
    __tablename__ = "nd_geolocation"
    __table_args__ = {'autoload': True}


class NdGeolocationprop(ORMBaseChado):
    __tablename__ = "nd_geolocationprop"
    __table_args__ = {'autoload': True}


class NdProtocol(ORMBaseChado):
    __tablename__ = "nd_protocol"
    __table_args__ = {'autoload': True}


class NdProtocolReagent(ORMBaseChado):
    __tablename__ = "nd_protocol_reagent"
    __table_args__ = {'autoload': True}


class NdProtocolprop(ORMBaseChado):
    __tablename__ = "nd_protocolprop"
    __table_args__ = {'autoload': True}


class NdReagent(ORMBaseChado):
    __tablename__ = "nd_reagent"
    __table_args__ = {'autoload': True}


class NdReagentRelationship(ORMBaseChado):
    __tablename__ = "nd_reagent_relationship"
    __table_args__ = {'autoload': True}


class NdReagentprop(ORMBaseChado):
    __tablename__ = "nd_reagentprop"
    __table_args__ = {'autoload': True}


class OrganismCvterm(ORMBaseChado):
    __tablename__ = "organism_cvterm"
    __table_args__ = {'autoload': True}


class OrganismCvtermprop(ORMBaseChado):
    __tablename__ = "organism_cvtermprop"
    __table_args__ = {'autoload': True}


class OrganismDbxref(ORMBaseChado):
    __tablename__ = "organism_dbxref"
    __table_args__ = {'autoload': True}


class OrganismPub(ORMBaseChado):
    __tablename__ = "organism_pub"
    __table_args__ = {'autoload': True}


class OrganismRelationship(ORMBaseChado):
    __tablename__ = "organism_relationship"
    __table_args__ = {'autoload': True}


class Organismprop(ORMBaseChado):
    __tablename__ = "organismprop"
    __table_args__ = {'autoload': True}


class OrganismpropPub(ORMBaseChado):
    __tablename__ = "organismprop_pub"
    __table_args__ = {'autoload': True}


class Phendesc(ORMBaseChado):
    __tablename__ = "phendesc"
    __table_args__ = {'autoload': True}


class Phenotype(ORMBaseChado):
    __tablename__ = "phenotype"
    __table_args__ = {'autoload': True}


class PhenotypeComparison(ORMBaseChado):
    __tablename__ = "phenotype_comparison"
    __table_args__ = {'autoload': True}


class PhenotypeComparisonCvterm(ORMBaseChado):
    __tablename__ = "phenotype_comparison_cvterm"
    __table_args__ = {'autoload': True}


class PhenotypeCvterm(ORMBaseChado):
    __tablename__ = "phenotype_cvterm"
    __table_args__ = {'autoload': True}


class Phenotypeprop(ORMBaseChado):
    __tablename__ = "phenotypeprop"
    __table_args__ = {'autoload': True}


class Phenstatement(ORMBaseChado):
    __tablename__ = "phenstatement"
    __table_args__ = {'autoload': True}


class Phylonode(ORMBaseChado):
    __tablename__ = "phylonode"
    __table_args__ = {'autoload': True}


class PhylonodeDbxref(ORMBaseChado):
    __tablename__ = "phylonode_dbxref"
    __table_args__ = {'autoload': True}


class PhylonodeOrganism(ORMBaseChado):
    __tablename__ = "phylonode_organism"
    __table_args__ = {'autoload': True}


class PhylonodePub(ORMBaseChado):
    __tablename__ = "phylonode_pub"
    __table_args__ = {'autoload': True}


class PhylonodeRelationship(ORMBaseChado):
    __tablename__ = "phylonode_relationship"
    __table_args__ = {'autoload': True}


class Phylonodeprop(ORMBaseChado):
    __tablename__ = "phylonodeprop"
    __table_args__ = {'autoload': True}


class Phylotree(ORMBaseChado):
    __tablename__ = "phylotree"
    __table_args__ = {'autoload': True}


class PhylotreePub(ORMBaseChado):
    __tablename__ = "phylotree_pub"
    __table_args__ = {'autoload': True}


class Phylotreeprop(ORMBaseChado):
    __tablename__ = "phylotreeprop"
    __table_args__ = {'autoload': True}


class Project(ORMBaseChado):
    __tablename__ = "project"
    __table_args__ = {'autoload': True}


class ProjectAnalysis(ORMBaseChado):
    __tablename__ = "project_analysis"
    __table_args__ = {'autoload': True}


class ProjectContact(ORMBaseChado):
    __tablename__ = "project_contact"
    __table_args__ = {'autoload': True}


class ProjectDbxref(ORMBaseChado):
    __tablename__ = "project_dbxref"
    __table_args__ = {'autoload': True}


class ProjectFeature(ORMBaseChado):
    __tablename__ = "project_feature"
    __table_args__ = {'autoload': True}


class ProjectPub(ORMBaseChado):
    __tablename__ = "project_pub"
    __table_args__ = {'autoload': True}


class ProjectRelationship(ORMBaseChado):
    __tablename__ = "project_relationship"
    __table_args__ = {'autoload': True}


class ProjectStock(ORMBaseChado):
    __tablename__ = "project_stock"
    __table_args__ = {'autoload': True}


class Projectprop(ORMBaseChado):
    __tablename__ = "projectprop"
    __table_args__ = {'autoload': True}


class Protocol(ORMBaseChado):
    __tablename__ = "protocol"
    __table_args__ = {'autoload': True}


class Protocolparam(ORMBaseChado):
    __tablename__ = "protocolparam"
    __table_args__ = {'autoload': True}


class Pub(ORMBaseChado):
    __tablename__ = "pub"
    __table_args__ = {'autoload': True}


class PubDbxref(ORMBaseChado):
    __tablename__ = "pub_dbxref"
    __table_args__ = {'autoload': True}


class PubRelationship(ORMBaseChado):
    __tablename__ = "pub_relationship"
    __table_args__ = {'autoload': True}


class Pubauthor(ORMBaseChado):
    __tablename__ = "pubauthor"
    __table_args__ = {'autoload': True}


class PubauthorContact(ORMBaseChado):
    __tablename__ = "pubauthor_contact"
    __table_args__ = {'autoload': True}


class Pubprop(ORMBaseChado):
    __tablename__ = "pubprop"
    __table_args__ = {'autoload': True}


class Quantification(ORMBaseChado):
    __tablename__ = "quantification"
    __table_args__ = {'autoload': True}


class QuantificationRelationship(ORMBaseChado):
    __tablename__ = "quantification_relationship"
    __table_args__ = {'autoload': True}


class Quantificationprop(ORMBaseChado):
    __tablename__ = "quantificationprop"
    __table_args__ = {'autoload': True}


class Stock(ORMBaseChado):
    __tablename__ = "stock"
    __table_args__ = {'autoload': True}


class StockCvterm(ORMBaseChado):
    __tablename__ = "stock_cvterm"
    __table_args__ = {'autoload': True}


class StockCvtermprop(ORMBaseChado):
    __tablename__ = "stock_cvtermprop"
    __table_args__ = {'autoload': True}


class StockDbxref(ORMBaseChado):
    __tablename__ = "stock_dbxref"
    __table_args__ = {'autoload': True}


class StockDbxrefprop(ORMBaseChado):
    __tablename__ = "stock_dbxrefprop"
    __table_args__ = {'autoload': True}


class StockFeature(ORMBaseChado):
    __tablename__ = "stock_feature"
    __table_args__ = {'autoload': True}


class StockFeaturemap(ORMBaseChado):
    __tablename__ = "stock_featuremap"
    __table_args__ = {'autoload': True}


class StockGenotype(ORMBaseChado):
    __tablename__ = "stock_genotype"
    __table_args__ = {'autoload': True}


class StockLibrary(ORMBaseChado):
    __tablename__ = "stock_library"
    __table_args__ = {'autoload': True}


class StockPub(ORMBaseChado):
    __tablename__ = "stock_pub"
    __table_args__ = {'autoload': True}


class StockRelationship(ORMBaseChado):
    __tablename__ = "stock_relationship"
    __table_args__ = {'autoload': True}


class StockRelationshipCvterm(ORMBaseChado):
    __tablename__ = "stock_relationship_cvterm"
    __table_args__ = {'autoload': True}


class StockRelationshipPub(ORMBaseChado):
    __tablename__ = "stock_relationship_pub"
    __table_args__ = {'autoload': True}


class Stockcollection(ORMBaseChado):
    __tablename__ = "stockcollection"
    __table_args__ = {'autoload': True}


class StockcollectionDb(ORMBaseChado):
    __tablename__ = "stockcollection_db"
    __table_args__ = {'autoload': True}


class StockcollectionStock(ORMBaseChado):
    __tablename__ = "stockcollection_stock"
    __table_args__ = {'autoload': True}


class Stockcollectionprop(ORMBaseChado):
    __tablename__ = "stockcollectionprop"
    __table_args__ = {'autoload': True}


class Stockprop(ORMBaseChado):
    __tablename__ = "stockprop"
    __table_args__ = {'autoload': True}


class StockpropPub(ORMBaseChado):
    __tablename__ = "stockprop_pub"
    __table_args__ = {'autoload': True}


class Study(ORMBaseChado):
    __tablename__ = "study"
    __table_args__ = {'autoload': True}


class StudyAssay(ORMBaseChado):
    __tablename__ = "study_assay"
    __table_args__ = {'autoload': True}


class Studydesign(ORMBaseChado):
    __tablename__ = "studydesign"
    __table_args__ = {'autoload': True}


class Studydesignprop(ORMBaseChado):
    __tablename__ = "studydesignprop"
    __table_args__ = {'autoload': True}


class Studyfactor(ORMBaseChado):
    __tablename__ = "studyfactor"
    __table_args__ = {'autoload': True}


class Studyfactorvalue(ORMBaseChado):
    __tablename__ = "studyfactorvalue"
    __table_args__ = {'autoload': True}


class Studyprop(ORMBaseChado):
    __tablename__ = "studyprop"
    __table_args__ = {'autoload': True}


class StudypropFeature(ORMBaseChado):
    __tablename__ = "studyprop_feature"
    __table_args__ = {'autoload': True}


class Synonym(ORMBaseChado):
    __tablename__ = "synonym"
    __table_args__ = {'autoload': True}


class Tableinfo(ORMBaseChado):
    __tablename__ = "tableinfo"
    __table_args__ = {'autoload': True}


class Treatment(ORMBaseChado):
    __tablename__ = "treatment"
    __table_args__ = {'autoload': True}
