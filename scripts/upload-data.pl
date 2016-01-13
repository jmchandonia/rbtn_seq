#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::KBase::workspace::ScriptHelpers qw(get_ws_client workspace printObjectInfo);

my $usage = 
"file_with_meta_data file_with_logratios <-w workspace> <-n object_name> <-g genome_name>
";

my $workspace = "";
my $object_name = "new";
my $genome_name = "ps_rch2_v2";

(GetOptions('w=s' => \$workspace, 
	    'n=s' => \$object_name,
            'g=s' => \$genome_name) and scalar(@ARGV) >= 2)    || die $usage;

my ($FNM_META, $FNM_LOGRT) = @ARGV;

die "Can't open file $FNM_META\n" if ! -e $FNM_META;
die "Can't open file $FNM_LOGRT\n" if ! -e $FNM_LOGRT;

my $serv = get_ws_client();

sub createObject($$);
sub createObjects($$@);
sub createConditionObject($$$$);
sub createMediaObject($);
sub createGrowthParamsObj($$$$$$$$$$$$$$$$$$);
sub createBarSeqExperimentObject($$$$$$$);
sub createObjectsForMissingRefs($$$);

sub getIndexOfElemExactMatch($$);
sub formKBname(@);



#####################################################
#get genome to access genes (Features)
#####################################################
my $genome_ref = "";
my $output = $serv->list_objects( {"workspaces" => [($workspace)], "type" => "KBaseGenomes.Genome"} );

if ( scalar(@$output)>0 ) {
    foreach my $object_info (@$output) {
	#---print  "Existing objects in ws: ".$object_info->[1]."\n";
	$genome_ref = $object_info->[6]."/".$object_info->[0]."/".$object_info->[4]
	    if ($object_info->[1] eq  $genome_name);
    }
}
#print "Ref: $genome_ref\n";


#get genome object
my $genome = $serv->get_object({
    id => $genome_name, 
    type => "KBaseGenomes.Genome",
    workspace => $workspace
    #instance => $jobdata->{TranscriptSet_inst},
    #auth => $job->{auth}
});
#foreach(@{$genome->{metadata}}){
#    print "$_\n";
#}



my %Aliases2FeatID = ();
my %FeatID2index = ();
my $FeatIndex2id = {};

print "Genome: ",$genome->{data}->{scientific_name}, "\n";
print "Number of features : ",scalar(@{$genome->{data}->{features}}), "\n";
for(my $i=0; $i< scalar(@{$genome->{data}->{features}}); ++$i){
#foreach my $f (@{$genome->{data}->{features}}){
    my $f = $genome->{data}->{features}->[$i];

    if($f->{type} eq "CDS"){
	#print "F ".$f->{id}." type : ".$f->{type}."\n"; 
	my @aliases = @{$f->{aliases}};
	foreach(@aliases){
	    #print "alias: $_ -> ".$f->{id}."\n";
	    $Aliases2FeatID{ $_ } = $f->{id};
	}
	$FeatIndex2id->{ $i } = $f->{id};
	$FeatID2index{ $f->{id} } = $i;
    }
}

#####################################################
#get existing Media objects in WS to save creating 
#new versions for the same thing.
#as currently the media object is just the name.
#####################################################
my %Media2objref = ();

$output = $serv->list_objects( {"workspaces" => [($workspace)], "type" => "KBaseBiochem.Media"} );

if ( scalar(@$output)>0 ) {
    foreach my $object_info (@$output) {
	#---print  "Existing Media objects in ws: ".$object_info->[1]."\n";
	$Media2objref{ $object_info->[1] } = $object_info->[6]."/".$object_info->[0]."/".$object_info->[4];
    }
}


#####################################################
#read in file with metadata
#####################################################
open FILE, $FNM_META or die $!;

#get header
my @header = split /\t/, <FILE>;
die "Wrong number of columns ".scalar(@header)." expected 35\n" if scalar(@header)!=35;

#find experiment name index
my $ExpNameIndex = getIndexOfElemExactMatch(\@header, 'name');

#find Media index 
my $MediaIndex = getIndexOfElemExactMatch(\@header, 'Media');

my %Exprname2cond1Name = ();
my %Exprname2cond2Name = ();
my %Exprname2mediaName = ();
my %Cond2objref = ();

my @meta = (); #store all lines from the file for the following passes
while(<FILE>){
    chomp;
    push @meta, $_;

    my @l = split /\t/, $_;
    die "Wrong number of columns in meta file, line $_\n",scalar(@l)," expected 35\n" if scalar(@l)!=35;
    
    if(length($l[$MediaIndex]) > 0){ #some experiments don't have Media
	my $medname = formKBname( $l[$MediaIndex] );

	#ws_client, workspace, string_media_name 
	$Media2objref{ $medname } = createMediaObject($medname)           
            if ! exists $Media2objref{ $medname };
	
	die "Two experiments with the same name ".$l[$ExpNameIndex]."\n" if exists $Exprname2mediaName{ $l[$ExpNameIndex] };
	$Exprname2mediaName{ $l[$ExpNameIndex] } = $medname; 
    }

    #get conditions 1 and 2
    my $c1 = $l[ getIndexOfElemExactMatch(\@header, 'Condition_1') ];
    my $conc1 = $l[ getIndexOfElemExactMatch(\@header, 'Concentration_1') ];
    my $u1 = $l[ getIndexOfElemExactMatch(\@header, 'Units_1') ];
    my $cname1 = formKBname($c1, $conc1, $u1);
	
    if($conc1 !~ /NA/){
    	$Cond2objref{ $cname1 } = createConditionObject($cname1, $c1, $conc1, $u1 ) 
	    if ! exists $Cond2objref{ $cname1 };

	$Exprname2cond1Name{ $l[$ExpNameIndex] } = $cname1;
    }

    my $c2 = $l[ getIndexOfElemExactMatch(\@header, 'Condition_2') ];
    my $conc2 = $l[ getIndexOfElemExactMatch(\@header, 'Concentration_2') ];
    my $u2 = $l[ getIndexOfElemExactMatch(\@header, 'Units_2') ];
    my $cname2 = formKBname($c2, $conc2, $u2);

    if($conc2 !~ /NA/){
	$Cond2objref{ $cname2 } = createConditionObject($cname2, $c2, $conc2, $u2 ) 
		if ! exists $Cond2objref{ $cname2 };

	$Exprname2cond2Name{ $l[$ExpNameIndex] } = $cname2;
    }
}
close FILE;


#####################################################
#create Condition objs in ws
#####################################################
#we don't create Condition objects in ws anymore
#createObjectsForMissingRefs($serv, $workspace, \%Cond2objref);


#####################################################
#create Media objs in ws
#####################################################
createObjectsForMissingRefs($serv, $workspace, \%Media2objref);

#####################################################
#make a second pass and build GrowthParams objs
#####################################################
my %Grwth2objref = ();

foreach (@meta){
    my @l = split /\t/, $_;
  
    my $grwthname = formKBname( 
	$l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ], 
	$l[$ExpNameIndex],
	$l[ getIndexOfElemExactMatch(\@header, 'Description') ]
    );

    my $growthobj = createGrowthParamsObj(
	$grwthname,
	$l[ getIndexOfElemExactMatch(\@header, 'Description') ],
	$l[ getIndexOfElemExactMatch(\@header, 'gDNA.plate') ],
	$l[ getIndexOfElemExactMatch(\@header, 'gDNA.well') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Index') ],
	exists $Exprname2mediaName{ $l[$ExpNameIndex] } ? $Media2objref{ $Exprname2mediaName{ $l[$ExpNameIndex] } } : 'NA',
	$l[ getIndexOfElemExactMatch(\@header, 'Growth.Method') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Group') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Temperature') ],
	$l[ getIndexOfElemExactMatch(\@header, 'pH') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Liquid.v..solid') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Aerobic_v_Anaerobic') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Shaking') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Growth.Plate.ID') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Growth.Plate.wells') ],
	$l[ getIndexOfElemExactMatch(\@header, 'StartOD') ],
	$l[ getIndexOfElemExactMatch(\@header, 'EndOD') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Total.Generations') ]
	);
    $Grwth2objref{ $l[$ExpNameIndex] } = $growthobj;
}

#we don't create GrowthParameters objects in ws anymore
#createObjectsForMissingRefs($serv, $workspace, \%Grwth2objref);

#####################################################
#make a third pass and build BarSeqExpr objs
#####################################################
my %Barseq2objref = ();

foreach (@meta){
    my @l = split /\t/, $_;

    my @conds = ();

    push @conds, $Cond2objref{ $Exprname2cond1Name{ $l[$ExpNameIndex] } }->{data} if exists $Exprname2cond1Name{ $l[$ExpNameIndex] };
    push @conds, $Cond2objref{ $Exprname2cond2Name{ $l[$ExpNameIndex] } }->{data} if exists $Exprname2cond2Name{ $l[$ExpNameIndex] };
   
    my $barseqname = formKBname( 
	$l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ],
	$l[$ExpNameIndex],
	$l[ getIndexOfElemExactMatch(\@header, 'Description') ]
    );

    #create BarSeqExperiment obj
    my $barseqobj = createBarSeqExperimentObject(
	$barseqname, 
	$l[ getIndexOfElemExactMatch(\@header, 'Person') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Date_pool_expt_started') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Sequenced.At') ],
	$Grwth2objref{ $l[$ExpNameIndex] }->{data},
	[ @conds ]
	);
    
    die "Two BarSeq experiments with the same name : ".$l[$ExpNameIndex]."\n" if
        exists $Barseq2objref{ $l[$ExpNameIndex] };

    $Barseq2objref{ $l[$ExpNameIndex] } = $barseqobj;
    #for(my $i=0; $i<= $#l; ++$i){
#	$meta{ $l[$ExpNameIndex] }{ $header[$i] } = $l[ $i ];	
 #   }
}

#we don't create these objects in ws anymore, too many
#createObjectsForMissingRefs($serv, $workspace, \%Barseq2objref);

####################################################
#create BarSeqExperimentResults
#
#typedef tuple<feature_ref,int strain_index,int count_begin,int count_end,float norm_log_ratio> bar_seq_result;
#
#/*
#  BarSeqExperimentResults stores the log ratios calculated from
#  a BarSeqExperiment.  There is one log ratio per strain per
#  GrowthParameters.
#*/
#typedef structure {
#    barseq_experiment_ref experiment;
#    list<bar_seq_result> results;
#} BarSeqExperimentResults;
##################################################### 

#print "Test: ",$Aliases2FeatID{ "Psest_4147" }, " : ", $FeatID2index{ $Aliases2FeatID{ "Psest_4147" } } , "\n";

# my $elem = 
#     [( #$genome_ref."/features/id/".$Aliases2FeatID{ "Psest_4147" } ,
#       int($FeatID2index{ $Aliases2FeatID{ "Psest_1061" } }),
#       -1,
#       -1,
#       -1,
#       -2
#     )];

# my $name = "test1";
# my $params = {
#     "name" => $name,
#     "type" => "KBaseRBTnSeq.BarSeqExperimentResults",
# };

# $params->{data}->{genome} = $genome_ref;
# $params->{data}->{experiments} = [ ( [( $Barseq2objref{ (keys %Barseq2objref)[0] }, [ ($elem) ] )]   )]; 
# $params->{data}->{feature_index_to_id} = $FeatIndex2id; 

# print "Test: ",$params->{name}, " : ", $params->{data}->{genome}, "\n";
# my %BarseqRes2objref = ();
# $BarseqRes2objref{ $name } = $params;
# createObjectsForMissingRefs($serv, $workspace, \%BarseqRes2objref);



#####################################################
#readin file with log ratios
##################################################### 

#FloatMatrix2D object
my $matrix2D = { 
    "row_ids" => [],
    "col_ids" => [],
    "values"  => []
};
my %BarseqLongName2short = ();
my %BarseqShortName2orig = ();
my $row_to_index = {};

my %data = ();
open FILE, $FNM_LOGRT or die $!;

#get header
my @headerData = split /\t/, <FILE>;
die "Too few columns in data file ".scalar(@headerData)."\n" if scalar(@headerData)<5;

#trim experiment annotations, keep only experiment names
foreach (@headerData){
    ($_) = split;
}

for(my $i=4; $i<= $#headerData; ++$i){
    die "Experiment ".$headerData[$i]." is not in meta file\n" if !exists $Barseq2objref{ $headerData[$i] };
    my @t = split /\|/,  $Barseq2objref{ $headerData[$i] }{name};
    my $shortname = $t[ $#t ];
    $shortname =~ s/_/ /g;

    #short names can be non unique, add index
    if(exists $BarseqShortName2orig{ $shortname }){
	my $counter = 2;
	while( exists $BarseqShortName2orig{ $shortname."$counter" } ){
	    ++$counter;
	}
	$shortname = $shortname."$counter";
    }
    #print "Test: name : ", $Barseq2objref{ $headerData[$i] }{name}, " : ", $shortname, "\n";
    $BarseqLongName2short{ $Barseq2objref{ $headerData[$i] }{name} } = $shortname;
    $BarseqShortName2orig{ $shortname } = $headerData[$i];
    push @{$matrix2D->{ "col_ids" }}, $shortname;#$headerData[$i];
}
#foreach (@headerData){
#    print "test: $_\n";
#    print "ok\n" if exists $meta{$_};
#}


my %barseqdata=(); #hash of barseq_experiment_ref -> ref to list< bar_seq_result >
my $gcounter = 0;


while(<FILE>){
    chomp;
    my @l = split /\t/, $_;
    die "Wrong number of columns in data file, line $_\n",scalar(@l)," expected ", scalar(@headerData), "as in the first line.\n" 
	if scalar(@l) != scalar(@headerData);
    
    my ($locusId, $sysName, $desc, $comb, @lratios) = @l;

    if( !exists $FeatID2index{ $sysName } ){
	if( exists $Aliases2FeatID{ $sysName } ){
		$sysName = $Aliases2FeatID{ $sysName };
	}else{
	    print STDERR "Error: $sysName not found in genome $genome_name. Skipping this gene.\n";
	    die "Error: $sysName not found in genome $genome_name\n";
	}
    }
    die "Error: with index: for $sysName not found in genome $genome_name\n"
        if !exists $FeatID2index{ $sysName };

    my $feat_index=int($FeatID2index{ $sysName });

    ++$gcounter;
    for(my $i=0; $i<= $#lratios; ++$i){
	#fill in data
	#typedef tuple<int feature_index,int strain_index,int count_begin,int count_end,float norm_log_ratio> bar_seq_result;
	#typedef tuple<barseq_experiment_ref experiment, bar_seq_result results> bar_seq_exp;

	my @res = ($feat_index, 0,0,0,0+$lratios[ $i ]);
	push @{$barseqdata{ $headerData[ $i + 4 ] }}, [ @res ];	
    }

    push @{$matrix2D->{ "values" }}, [ @lratios ];
    push @{$matrix2D->{ "row_ids" }}, $sysName ;
    $row_to_index->{ $sysName } = $feat_index;
}
close FILE;
print "Saving data for $gcounter genes\n";

#prepare BarSeqExperimentResults object;
my $params = {
    "name" => $object_name,
    "type" => "KBaseRBTnSeq.BarSeqExperimentResults",
};

$params->{data}->{genome} = $genome_ref;
$params->{data}->{feature_index_to_id} = $FeatIndex2id; 
$params->{data}->{experiments}  = [];
$params->{data}->{features_by_experiments}  = $matrix2D;
$params->{data}->{col_to_index}  = {};
$params->{data}->{row_to_index}  = $row_to_index;

#fill in experiment field
#foreach (keys %barseqdata){
#use the same order as in matrix2D
foreach my $shortname( @{$matrix2D->{ "col_ids" }} ){
    my $name = $BarseqShortName2orig{ $shortname };

    push @{ $params->{data}->{experiments} }, [(  $Barseq2objref{ $name }->{data} , $barseqdata{ $name } )];

    $params->{data}->{col_to_index}->{  $shortname } = $#{ $params->{data}->{experiments} };
}

#print "Test: ",$params->{name}, " : ", $params->{data}->{genome}, "\n";
my %BarseqRes2objref = ();
$BarseqRes2objref{ $object_name } = $params;
createObjectsForMissingRefs($serv, $workspace, \%BarseqRes2objref);

exit(0);



###########################################################
# functions
###########################################################

#join array of strings into one string and subs 
#illegal chars with something else
sub formKBname(@){
    my $out = join('|', @_);

    #subs white spaces with '_'
    $out =~ s/\s+/_/g;

    #subst , with '_'
    $out =~ s/,/_/g;

    #subst : with '_'
    $out =~ s/:/_/g;

    #subst % with 'pct'
    $out =~ s/%/pct/g;

    #subst / with '_over_'
    $out =~ s/\//_over_/g;

    #subst () with '_'
    $out =~ s/\(/_/g;
    $out =~ s/\)/_/g;
    
    #subst ; with '_'
    $out =~ s/\;/_/g;

    #subst ? with '_'
    $out =~ s/\?/_/g;

    #subst non-ascii with '_'
    $out =~ s/[^[:ascii:]]/_/g;

    #collapse multiple '_' into one
    $out =~ s/_+/_/g;

    return $out;
}

#creates a new KBaseRBTnSeq.Condition object
#input:  obj name, condition name, concentration, units
sub createConditionObject($$$$){
	my $params = {
		"name" => $_[0],
		"type" => "KBaseRBTnSeq.Condition",
	};
	$params->{data}->{name} = $_[1];
	$params->{data}->{concentration} =   $_[2]+0;
	$params->{data}->{units} = $_[3];

	return $params;
}


#creates a new KBaseBiochem.Media object
#input:  obj_name 
sub createMediaObject($){
	my $params = {
		"name" => $_[0],
		"type" => "KBaseBiochem.Media",
	};
	$params->{data}->{name} = $_[0];
	$params->{data}->{id} =   $_[0];#"kb|type.0";#$_[2];
	$params->{data}->{isDefined} = 0;
	$params->{data}->{type} = "custom";
	$params->{data}->{isMinimal} = 0;
	$params->{data}->{mediacompounds} = [];
	
	return $params;
}

#creates a new KBaseRBTnSeq.GrowthParameters object
#input: 
#0    obj name 
#1    string description;
#2    string gDNA_plate;
#3    string gDNA_well;
#4    string index;
#5    media_ref media;
#6    string growth_method;
#7    string group;
#8    float temperature;
#9    float pH;
#10    bool isLiquid;
#11    bool isAerobic;
#12    string shaking;
#13    string growth_plate_id;
#14    string growth_plate_wells;
#15    float startOD;
#16    float endOD;
#17    float total_generations;
sub createGrowthParamsObj($$$$$$$$$$$$$$$$$$){
    #foreach(@_){
    #	print "param:$_:\n";
    #}
    my $params = {
	"name" => $_[0],
	"type" => "KBaseRBTnSeq.GrowthParameters",
    };
    $params->{data}->{description} = $_[1];
    $params->{data}->{gDNA_plate}  = $_[2];
    $params->{data}->{gDNA_well} = $_[3];
    $params->{data}->{index} = $_[4];
    $params->{data}->{media} = $_[5] if $_[5] ne 'NA';
    $params->{data}->{growth_method} = $_[6];
    $params->{data}->{group} = $_[7];
    $params->{data}->{temperature} = $_[8]+0 if $_[8] !~ 'NA' and $_[8]>0;
    $params->{data}->{pH} = $_[9]+0 if $_[9] !~ 'NA' and $_[9]>0;
    if($_[10] =~ /liquid/i){
	$params->{data}->{isLiquid} = 1;
    }elsif($_[10] =~ /solid/i){
	$params->{data}->{isLiquid} = 0;
    }
    if($_[11] =~ /aerobic/i and $_[11] !~ /anaerobic/i){
	$params->{data}->{isAerobic} = 1;
    }elsif($_[11] =~ /anaerobic/i){ 
	$params->{data}->{isAerobic} = 0;
    }
    $params->{data}->{shaking} = $_[12];
    $params->{data}->{growth_plate_id} = $_[13] if length($_[13])>0;
    $params->{data}->{growth_plate_wells} = $_[14] if length($_[14])>0;;
    $params->{data}->{startOD} = $_[15]+0 if $_[15] !~ 'NA' and $_[15]>0;
    $params->{data}->{endOD} = $_[16]+0   if $_[16] !~ 'NA' and $_[16]>0;
    $params->{data}->{total_generations} = $_[17]+0 if $_[17] !~ 'NA' and $_[17]>0;
    
    #print "Ref to media obj:".$params->{data}->{media}.":\n";
    return $params;
}


#create BarSeqExperiment obj
#input:  
#0 obj_name 
#1 person
#2 mutant_lib_name
#3 start_date
#4 sequenced_at
#5 growth_parameters_ref growth_parameters
#6 list<condition_ref> conditions;
sub createBarSeqExperimentObject($$$$$$$){
	my $params = {
		"name" => $_[0],
		"type" => "KBaseRBTnSeq.BarSeqExperiment",
	};
	$params->{data}->{name} = $_[0];
	$params->{data}->{person} = $_[1];
	$params->{data}->{mutant_lib_name} =   $_[2];
	$params->{data}->{start_date} = $_[3];
	$params->{data}->{sequenced_at} = $_[4];
	$params->{data}->{growth_parameters} = $_[5];
	$params->{data}->{conditions} = $_[6];

	return $params;
}


#creates new objects
#input: ws_client, ws_name, array of $params
sub createObjects($$@){
	my ($serv, $workspace, @params) = @_;

	#print "Saving object named : ". $params->{name}. " : Typed : ".$params->{type}. " :\n";

	my $ws_ver = '';
	eval { $ws_ver = $serv->ver(); };
	if($@) {
	    print "Object could not be saved! Error connecting to the WS server.\n";
	    print STDERR $@->{message}."\n";
	    if(defined($@->{status_line})) {print STDERR $@->{status_line}."\n" };
	    print STDERR "\n";
	    exit 1;
	}
	
	my $saveObjectsParams = {
	    "workspace" => $workspace,
	    "objects" => [ ]  #add objects next
	};
	
	foreach my $p (@params){
	    push @{$saveObjectsParams->{objects}},
	    {
		"data"  => $p->{data},
		"name"  => $p->{name},
		"type"  => $p->{type},
		"meta"  => $p->{metadata},
		"provenance" => 
		    [ {
			"service"=>"Workspace",
			"service_ver"=>$ws_ver,
			"script"=>"upload-data.pl",
			"script_command_line"=> "@ARGV"
		      } ]           
	    };
	}
	
	my $output;
	eval { $output = $serv->save_objects( $saveObjectsParams ); };
	if($@) {
	    print "Object could not be saved!\n";
	    print STDERR $@->{message}."\n";
	    if(defined($@->{status_line})) {print STDERR $@->{status_line}."\n" };
	    print STDERR "\n";
	    exit 1;
	}

	my @refs = ();
	#Report the results
	print "Object saved.  Details:\n";
	if (scalar(@$output)>0) {
	    foreach my $object_info (@$output) {
		#return reference to the created object
		push @refs, $object_info->[6]."/".$object_info->[0]."/".$object_info->[4];
		#return $object_info->[6]."/".$object_info->[0]."/".$object_info->[4];
		#return 
		#printObjectInfo($object_info);
		#print $object_info,"\n";
	    }
	} else {
	    print STDERR "No details returned: ws save_objects\n";
	}
	
	return @refs;
}



#creates a new object
#input: $params, ws_client
sub createObject($$){
	my ($params, $serv) = @_;

	print "Saving object named : ". $params->{name}. " : Typed : ".$params->{type}. " :\n";

	my $ws_ver = '';
	eval { $ws_ver = $serv->ver(); };
	if($@) {
    		print "Object could not be saved! Error connecting to the WS server.\n";
    		print STDERR $@->{message}."\n";
    		if(defined($@->{status_line})) {print STDERR $@->{status_line}."\n" };
    		print STDERR "\n";
    		exit 1;
	}

	# set provenance info
	my $PA = {
                "service"=>"Workspace",
                "service_ver"=>$ws_ver,
                "script"=>"upload-data.pl",
                "script_command_line"=> "@ARGV"
          };
	$params->{provenance} = [ $PA ];

	my $saveObjectsParams = {
		"workspace" => $params->{workspace},
                "objects" => [
                           {
	                        "data"  => $params->{data},
                                "name"  => $params->{name},
                                "type"  => $params->{type},
                                "meta"  => $params->{metadata},
                                "provenance" => $params->{provenance}
                           }
                        ]
        };

	my $output;
	eval { $output = $serv->save_objects( $saveObjectsParams ); };
	if($@) {
    		print "Object could not be saved!\n";
    		print STDERR $@->{message}."\n";
    		if(defined($@->{status_line})) {print STDERR $@->{status_line}."\n" };
    		print STDERR "\n";
    		exit 1;
	}

	#Report the results
	print "Object saved.  Details:\n";
	if (scalar(@$output)>0) {
        	foreach my $object_info (@$output) {
			#return reference to the created object
			return $object_info->[6]."/".$object_info->[0]."/".$object_info->[4];
			#return 
                	#printObjectInfo($object_info);
			#print $object_info,"\n";
        	}
	} else {
        	print STDERR "No details returned: ws save_objects\n";
	}
}

#input: array of strings, string to search for
#returns index of the query string
#dies if not found or multiple copies exist
sub getIndexOfElemExactMatch($$){
    my @array = @{$_[0]};
    my $str = $_[1];

    my @index = grep { $array[$_] =~ /^$str$/ } 0..$#array;
    die "Multiple or non existent field '$str' in array\n" if scalar(@index)!=1;
    
    return $index[0];
}


#input: 
#$serv, $workspace, ref to hash
#check whether hash data is a ref to a HASH (i.e. $params)
#then creates ws objects for all these
sub createObjectsForMissingRefs($$$){
    my ($serv, $workspace) = ($_[0], $_[1]);
    my $h =  $_[2];

    my @params = ();
    my @nms = ();
    
    foreach ( keys %{$h} ){
	if(  ref( $h->{$_} ) eq "HASH" ){ #i.e. obj reference is not yet created
	    push @params, $h->{$_};
	    push @nms, $_;

	    print "Name: ".$h->{$_}->{name}." : Type : ". $h->{$_}->{type}." :\n";
	}
    }

    return if scalar(@params) <= 0;

    print "Saving ".scalar(@params)." objects in ws\n";
    my @refs = createObjects($serv, $workspace, @params);
    die "Wrong number of refs for input params\n" if scalar(@params) != scalar(@refs);
    
    for(my $i=0; $i<=$#params; ++$i){
	$h->{ $nms[$i] } = $refs[ $i ];
    }
}
