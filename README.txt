# rbtn_seq

###################################################
#register new types with workspace
###################################################

ws-typespec-register -u jmc -t KBaseRBTnSeq.spec --add "Strain;Delta;Pool;MappedReads" --commit
ws-typespec-register -u jmc -release KBaseRBTnSeq

###################################################
#upload data using a script
###################################################

#from 'scripts' dir:
#should be logged in to kbase to execute the script
#params:
#file_with_meta_data file_with_logratios <-w workspace> <-n object_name> <-g genome_name>
#log ratios: RCH2_fit_logratios_good.tab.txt
#meta data: expsUsed.txt 
#-w WS id
#-g genome name, should be present in a given WS (e.g. ps_rch2_v3)
#-n name of the new object (e.g. rbRCH2)
#working example:
./run_in_kb_env ./upload-data.pl ../test/data/expsUsed.txt ../test/data/RCH2_fit_logratios_good.tab.txt -w mshatsky:1439854430561 -g ps_rch2_v3 -n rbRCH2



Widgets:

widgets/kbaseBarSeqExperimentResults.js
should be placed into:
narrative/src/notebook/ipython_profiles/profile_narrative/kbase_templates/static/kbase/js/widgets/function_output/

add:
<script src="{{ static_url("kbase/js/widgets/function_output/kbaseBarSeqExperimentResults.js") }}" type="text/javascript" charset="utf-8"></script>
to:
narrative/src/notebook/ipython_profiles/profile_narrative/kbase_templates/notebook.html


Additional dependencies:

widgets/extra/types/KBaseRBTnSeq.BarSeqExperimentResults/
should be placed into 
narrative_method_specs_ci/types/

widgets/extra/methods/view_barseq_results/
should be placed into
narrative_method_specs_ci/methods

