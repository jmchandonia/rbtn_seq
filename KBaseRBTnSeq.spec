#include <KBaseFeatureValues.spec>

/*
This module is for storing microbial growth phenotype data, e.g., from
FEBA project or other ENIGMA experiments
*/

module KBaseRBTnSeq {

typedef int bool;
    
/*
@id ws KBaseGenomes.Genome
*/
typedef string genome_ref;

/*
@id ws KBaseGenomes.Contig
*/
typedef string contig_ref;

/*
@id ws KBaseAssembly.SingleEndLibrary
*/
typedef string reads_ref;

/*
TnSeq barcode model, showing how the barcodes / transposon
were designed
*/
typedef tuple <string read_template, string flanking_sequence> tnseq_model;

/*
A single TnSeq mapped read
*/
typedef tuple <string read_name, string barcode, contig_ref scaffold, int insert_pos, string strand, bool is_unique, int hit_start, int hit_end, float bit_score, float pct_identity> mapped_read;

/*
A MappedReads object stores the mapping of reads to a genome
*/
typedef structure {
    genome_ref genome;
    reads_ref reads;
    tnseq_model model;
    list<mapped_read> mapped_reads;
} MappedReads;

/*
@id ws KBaseCommunities.Sample
*/
typedef string sample_ref;

/*
@id ws KBaseBiochem.Media
*/
typedef string media_ref;

/*
	Reference to a Feature object of a genome in the workspace
	@id subws KBaseGenomes.Genome.features.[*].id
*/
typedef string feature_ref;

/*
enum: insertion, deletion, substitution
The latter is not strictly necessary, but convenient to avoid storing
two separate events.
*/
typedef string change_type;

/*
  A Delta is a description of a single change to a strain.  A series
  of Deltas defines the transition from one Strain to another.  For
  sequenced insertions or substitutions, give the 0-indexed position
  on the contig where the insertion/substitution begins, in the +
  direction.  For sequenced deletions and substitutions, give the
  position and length.  The position of all Deltas should be
  calculated relative to the parent strain (derived_from_strain), so
  that the Deltas could be applied in any order.
@optional change feature contig sequence position length
*/
typedef structure {
    string description;
    change_type change;
    feature_ref feature;
    contig_ref contig;
    string sequence;
    int position;
    int length;
} Delta;

/*
@id ws KBaseRBTnSeq.Strain
*/
typedef string strain_ref;

/*
  A Strain is a particular genetic variant of an organism.  Optionally,
  it may be:
    * derived from another Strain (e.g., as an engineered mutant)
    * sequenced
    * isolated from a community
    * a wild-type example of a Genome
  If a strain is "wild type" it should have a non-null genome_ref and a
  null derived_from_strain.  If not wild type, genome_ref should be
  set to the "original" parent strain in KBase, if it exists, or null
  if it does not exist or is unknown.
@optional description genome derived_from_strain deltas isolated_from
*/
typedef structure {
    string name;
    string description;
    genome_ref genome;
    strain_ref derived_from_strain;
    list<Delta> deltas;
    sample_ref isolated_from;
} Strain;

/*
@id ws KBaseRBTnSeq.Pool
*/
typedef string pool_ref;

/*
  A Pool is a collection of barcoded strains.  Barcodes, tags, etc should
  be stored as Deltas in each strain.
*/
typedef structure {
    list<tuple<Strain strain, int count>> strains;
} Pool;

/*
 Reference to a compound object in a biochemistry
 @id subws KBaseBiochem.Biochemistry.compounds.[*].id
*/
typedef string compound_ref;
    
/*
  A Condition is something that is added to particular aliquots in
  a growth experiment, in addition to the media.  e.g., it may be a stress
  condition, or a nutrient.  Compound is needed if the condition is
  addition of a chemical in the KBase Biochemistry database.
@optional concentration units compound
*/
typedef structure {
    string name;
    float concentration;
    string units;
    compound_ref compound;
} Condition;

/*
@id ws KBaseRBTnSeq.Condition
*/
typedef string condition_ref;

/*
  GrowthParameters describes all the conditions a particular aliquot
  was subjected to in an experiment
@optional description gDNA_plate gDNA_well index media growth_method group temperature pH isLiquid isAerobic shaking growth_plate_id growth_plate_wells startOD endOD total_generations
*/
typedef structure {
    string description;
    string gDNA_plate;
    string gDNA_well;
    string index;
    media_ref media;
    string growth_method;
    string group;
    float temperature;
    float pH;
    bool isLiquid;
    bool isAerobic;
    string shaking;
    string growth_plate_id;
    string growth_plate_wells;
    float startOD;
    float endOD;
    float total_generations;
} GrowthParameters;

/*
@id ws KBaseRBTnSeq.GrowthParameters
*/
typedef string growth_parameters_ref;

/*
  A BarSeqExperiment is an experiment in which a pool is grown in
  several parallel aliquots (e.g., wells or tubes), each potentially
  treated with a different set of conditions (but usually it is just one
  Condition)
  BarSeqExperiment object represents only one such condition
@optional person mutant_lib_name tnseq_pool
*/
typedef structure {
    string name;
    string person;
    string mutant_lib_name;
    string start_date;
    string sequenced_at;
    GrowthParameters growth_parameters;
    list<Condition> conditions; /* usually one or two conditions */
    pool_ref tnseq_pool;
} BarSeqExperiment;

/*
@id ws KBaseRBTnSeq.BarSeqExperiment
*/
typedef string barseq_experiment_ref;

/*
  Number of times a barcode (i.e. a strain) was detected by sequencing a pool at beginning (refernce state)
  and at the end of GrowthParameters, and a calculated log ratio of strain abundance relative to a starting
  condition.
  feature_index - Genome.features[index] optional, if strain_index is NA and log ratio corresponds to a gene, not a strain
  strain_index - index of a strain in Pool.strains list
  count_begin - number of instances of the strain identified from sequencing at the beginning of experiment
  count_end - at the end of experiment
  norm_log_ratio - normalized log ratio between count_end and count_begin
*/
typedef tuple<int feature_index,int strain_index,int count_begin,int count_end,float norm_log_ratio> bar_seq_result;

/*
 A single (i.e. one condition) BarSeq experiment
 experiment - describes the BarSeq experiment
 bar_seq_result - list of counts. Can be per gene or per strain.
*/
typedef tuple<BarSeqExperiment experiment, list< bar_seq_result > results> bar_seq_exp;

/*
  BarSeqExperimentResults stores the log ratios calculated from
  a BarSeqExperiment.  There is one log ratio per strain or per gene per
  GrowthParameters.
  May contain a number of BarSeqExperiment done on same species (i.e. 
  usually it would be dozens and even more than a hundred barseq experiments 
  performed under various conditions starting from a the same library).

  The raw data, including start/end mutant counts, is stored in:
     list<bar_seq_exp> experiments
  It does not require to have values for all genes/loci and the number of features does not have 
  to be the same for all experiments.

  features_by_experiments - is a 2D matrix that contains only normalized log ratios from experiments[*]['results' = 1][*]['norm_log_ratio' = 5 ]
  Therefore, there is a redundancy, we store log ratios twice. 'experiments' is used to store all raw data and allows quickly retrieve a single condition. 
  'features_by_experiments' is used to support visualization widgets and other methods that work with KBaseFeatureValues.FloatMatrix2D. 
  It also allows to quickly retrieve log ratios for all conditions per gene.
*/
typedef structure {
    genome_ref genome;
    mapping<int feature_index, string feature_id> feature_index_to_id;
    list<bar_seq_exp> experiments;

    KBaseFeatureValues.FloatMatrix2D features_by_experiments;    
    mapping<string, int> col_to_index; /* column name of 'features_by_experiments' to index within 'experiments' mapping */
    mapping<string, int> row_to_index; /* row names (usually genes) to feature_index, i.e. reverse of 'feature_index_to_id' */
} BarSeqExperimentResults;

/*
@id ws KBaseRBTnSeq.BarSeqExperimentResults
*/
typedef string barseq_experiment_results_ref;

/*
  Computes essential genes from a TnSeq pool
  Input: pool_ref - reference to a Pool
  Output: list of genes 
*/
funcdef essential_genes(pool_ref) returns (list<feature_ref>) authentication required;

/*
  Computes essential genes from a BarSeq experiment
  Input: barseq_experiment_results_ref - reference to a BarSeqExperimentResults object
  Output: list of genes 
*/
funcdef essential_genes(barseq_experiment_results_ref) returns (list<feature_ref>) authentication required;

/*
  Computes gene fitness within a TnSeq pool
  Input: tnseq_library_ref - reference to a TnSeqLibrary
  Output: list of genes with their fitness
*/
funcdef gene_fitness(pool_ref) returns (list<tuple<feature_ref gene, float fitness>>) authentication required;

/*
  Computes gene fitness from a BarSeq experiment
  Input: barseq_experiment_results_ref - reference to a BarSeqExperimentResults object
  Output: list of genes with their fitness
*/
funcdef gene_fitness(barseq_experiment_results_ref) returns (list<tuple<feature_ref gene, float fitness>>) authentication required;

};
