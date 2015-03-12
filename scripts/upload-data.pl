#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = 
"file_with_meta_data file_with_logratios <-w workspace> <-n object_name>
";

my $workspace = "";
my $object_name = "new";

(GetOptions('w=s' => \$workspace, 
	    'n=s' => \$object_name) and scalar(@ARGV) >= 2)    || die $usage;

my ($FNM_META, $FNM_LOGRT) = @ARGV;

die "Can't open file $FNM_META\n" if ! -e $FNM_META;
die "Can't open file $FNM_LOGRT\n" if ! -e $FNM_LOGRT;

my $serv = get_ws_client();

sub createObject($$);
sub createConditionObject($$$);
sub createMediaObject($$$);
sub createGrowthParamsObj($$$$$$$$$$$$$$$$$$$$);
sub createBarSeqExperimentObject($$$$$$$$$);

sub getIndexOfElemExactMatch($$);

#####################################################
#get existing Media objects in WS to save creating 
#new versions for the same thing.
#as currently the media object is just the name.
#####################################################
my %Media2objref = ();

my $output = $serv->list_objects( {"workspaces" => [($workspace)], "type" => "KBaseBiochem.Media"} );

if ( scalar(@$output)>0 ) {
    foreach my $object_info (@$output) {
	#---
	print  "Existing objects in ws: ".$object_info->[1]."\n";
	$Media2objref{ $object_info->[1] } = $object_info->[6]."/".$object_info->[0]."/".$object_info->[4]."\n";
    }
}

#####################################################
#get existing Conditions objects in WS to save creating 
#new versions for the same thing.
#####################################################
my %Cond2objref = ();

$output = $serv->list_objects( {"workspaces" => [($workspace)], "type" => "KBaseRBTnSeq.Condition"} );

if ( scalar(@$output)>0 ) {
    foreach my $object_info (@$output) {
	#---
	print  "Existing objects in ws: ".$object_info->[1]."\n";
	$Cond2objref{ $object_info->[1] } = $object_info->[6]."/".$object_info->[0]."/".$object_info->[4]."\n";
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

my %Brseq2objref = ()

while(<FILE>){
    chomp;
    my @l = split /\t/, $_;
    die "Wrong number of columns in meta file, line $_\n",scalar(@l)," expected 35\n" if scalar(@l)!=35;
    
    #ws_client, workspace, string_media_name 
    $Media2objref{ $l[$MediaIndex] } = createMediaObject($serv, $workspace, $l[$MediaIndex])
           if ! exists $Media2objref{ $l[$MediaIndex] } and ! exists $Media2objref{ $l[$MediaIndex] };

    #get conditions 1 and 2
    my $c1 = $l[ getIndexOfElemExactMatch(\@header, 'Condition_1') ];
    my $con1 = $l[ getIndexOfElemExactMatch(\@header, 'Concentration_1') ];
    my $u1 = $l[ getIndexOfElemExactMatch(\@header, 'Units_1') ];

    $Cond2objref{ "$c1 : $con1 : $u1" } = createConditionObject($serv, $workspace, "$c1 : $con1 : $u1", "$c1 : $con1 : $u1", $c1, $con1, $u1 ) 
	    if ! exists $Cond2objref{ "$c1 : $con1 : $u1" };

    my $c2 = $l[ getIndexOfElemExactMatch(\@header, 'Condition_2') ];
    my $con2 = $l[ getIndexOfElemExactMatch(\@header, 'Concentration_2') ];
    my $u2 = $l[ getIndexOfElemExactMatch(\@header, 'Units_2') ];

    $Cond2objref{ "$c2 : $con2 : $u2" } = createConditionObject($serv, $workspace, "$c2 : $con2 : $u2",  "$c2 : $con2 : $u2", $c2, $con2, $u2 ) 
	    if ! exists $Cond2objref{ "$c2 : $con2 : $u2" };

    #create Condition obj
    my $growthobj = createGrowthParamsObj(
	 $serv, $workspace, $l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ]." : ".$l[$ExpNameIndex]." : ".$l[ getIndexOfElemExactMatch(\@header, 'Description') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Description') ],
	$l[ getIndexOfElemExactMatch(\@header, 'gDNA.plate') ],
	$l[ getIndexOfElemExactMatch(\@header, 'gDNA.well') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Index') ],
	$Media2objref{ $l[$MediaIndex] },
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


    #create BarSeqExperiment obj
    my $barseqobj = createBarSeqExperimentObj(
	$serv, $workspace, $l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ]." : ".$l[$ExpNameIndex]." : ".$l[ getIndexOfElemExactMatch(\@header, 'Description') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Person') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Mutant.Library') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Date_pool_expt_started') ],
	$l[ getIndexOfElemExactMatch(\@header, 'Sequenced.At') ],
	$growthobj,
	"c2" eq "NA" ? [( $Cond2objref{ "$c1 : $con1 : $u1" } )] : 
	               [( $Cond2objref{ "$c1 : $con1 : $u1" }, $Cond2objref{ "$c2 : $con2 : $u2" } )],
	);
    
    die "Two BarSeq experiments with the same name : ".$l[$ExpNameIndex]."\n" if
        exists $Brseq2objref{ $l[$ExpNameIndex] };

    $Brseq2objref{ $l[$ExpNameIndex] } = $barseqobj;
    #for(my $i=0; $i<= $#l; ++$i){
#	$meta{ $l[$ExpNameIndex] }{ $header[$i] } = $l[ $i ];	
 #   }
}
close FILE;

exit(0);


#####################################################
#readin file with log ratios
##################################################### 
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
    die "Experiment ".$headerData[$i]." is not in meta file\n" if !exists $Brseq2objref{ $headerData[$i] };
}
#foreach (@headerData){
#    print "test: $_\n";
#    print "ok\n" if exists $meta{$_};
#}


#typedef structure {
#    barseq_experiment_ref experiment;
#    list<bar_seq_result> results;
#} BarSeqExperimentResults;
my %brseqdata=(); #hash of many BarSeqExperimentResults

while(<FILE>){
    chomp;
    my @l = split /\t/, $_;
    die "Wrong number of columns in data file, line $_\n",scalar(@l)," expected ", scalar(@headerData), "as in the first line.\n" 
	if scalar(@l) != scalar(@headerData);
    
    my ($locusId, $sysName, $desc, $comb, @lratios) = @l;

    for(my $i=0; $i<= $#lratios; ++$i){
	#fill in data
	#typedef tuple<int strain_index,int count_begin,int count_end,float norm_log_ratio> bar_seq_result;
	my @res = (0,0,0,0+$lratios[ $i ]);
	push @{$brseqdata{ $headerData[ $i + 4 ] }->{ results }}, [ @res ];	
    }
}
close FILE;

#fill in experiment field
foreach (keys %brseqdata){
    $brseqdata{ $_ }->{ experiment } = "$_";
}



###################################################
my $params = {
       "id" => $object_name,
       "type" => "KBaseRBTnSeq.BarSeqExperimentResults-0.1",
       workspace => $workspace,
       metadata => ""
};

#try the first one
$params->{data} = $brseqdata{  $headerData[4] };

print "T: ", $params->{data}->{experiment}, "\n";
#print "T: ", $params->{data}->{results}->[0]->[0], " ", $params->{data}->{results}->[0]->[3], "\n";
#print "T: ", $params->{data}->{results}->[9]->[0], " ", $params->{data}->{results}->[9]->[3], "\n";

###################################################
use Bio::KBase::workspace::ScriptHelpers qw(get_ws_client workspace printObjectInfo);

#lookup version number of WS Service that will be loading the data
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


# setup the new save_objects parameters
my $saveObjectsParams = {
                "objects" => [
                           {
                                "data"  => $params->{data},
                                "name"  => $params->{id},
                                "type"  => $params->{type},
                                #"meta"  => $params->{metadata},
                                "provenance" => $params->{provenance}
                           }
                        ]
        };
if ($params->{workspace} =~ /^\d+$/ ) { #is ID
        $saveObjectsParams->{id}=$workspace+0;
} else { #is name
        $saveObjectsParams->{workspace}=$workspace;
}

#Calling the server
my $output;
eval { $output = $serv->save_objects($saveObjectsParams); };
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
                printObjectInfo($object_info);
        }
} else {
        print "No details returned!\n";
}
print "\n";
exit(0);


#creates a new KBaseRBTnSeq.Condition object
#input:  ws_client, workspace, obj name 
#condition name, concentration, units ) 
sub createConditionObject($$$$$$){
	my $params = {
		"name" => $_[2],
		"type" => "KBaseRBTnSeq.Condition",
		"workspace" => $_[1],
	};
	$params->{data}->{name} = $_[3];
	$params->{data}->{concentration} =   $_[4];
	$params->{data}->{units} = $_[5];

	return createObject($params, $_[0]);
}


#creates a new KBaseBiochem.Media object
#input:  ws_client, workspace, obj_name 
sub createMediaObject($$$){
	my $params = {
		"name" => $_[2],
		"type" => "KBaseBiochem.Media",
		"workspace" => $_[1],
	};
	$params->{data}->{name} = $_[2];
	$params->{data}->{id} =   $_[2];
	$params->{data}->{isDefined} = 0;
	$params->{data}->{type} = "custom";
	$params->{data}->{isMinimal} = 0;
	$params->{data}->{mediacompounds} = [];
	
	return createObject($params, $_[0]);
}

#creates a new KBaseRBTnSeq.GrowthParameters object
#input:  ws_client, workspace, obj name 
#3    string description;
#4    string gDNA_plate;
#5    string gDNA_well;
#6    string index;
#7    media_ref media;
#8    string growth_method;
#9    string group;
#10    float temperature;
#11    float pH;
#12    bool isLiquid;
#13    bool isAerobic;
#14    string shaking;
#15    string growth_plate_id;
#16    string growth_plate_wells;
#17    float startOD;
#18    float endOD;
#19    float total_generations;
sub createGrowthParamsObj($$$$$$$$$$$$$$$$$$$$){
	my $params = {
		"name" => $_[2],
		"type" => "KBaseRBTnSeq.GrowthParameters",
		"workspace" => $_[1],
	};
	$params->{data}->{description} = $_[3];
	$params->{data}->{gDNA_plate}  = $_[4];
	$params->{data}->{gDNA_well} = $_[5];
	$params->{data}->{index} = $_[6];
	$params->{data}->{media} = $_[7];
	$params->{data}->{growth_method} = $_[8];
	$params->{data}->{group} = $_[9];
	$params->{data}->{temperature} = $_[10];
	$params->{data}->{pH} = $_[11];
	$params->{data}->{isLiquid} = $_[12];
	$params->{data}->{isAerobic} = $_[13];
	$params->{data}->{shaking} = $_[14];
	$params->{data}->{growth_plate_id} = $_[15];
	$params->{data}->{growth_plate_wells} = $_[16];
	$params->{data}->{startOD} = $_[17];
	$params->{data}->{endOD} = $_[18];
	$params->{data}->{total_generations} = $_[19];

	return createObject($params, $_[0]);
}


#create BarSeqExperiment obj
#input:  ws_client, workspace, obj_name 
#3 person
#4 mutant_lib_name
#5 start_date
#6 sequenced_at
#7 growth_parameters_ref growth_parameters
#8 list<condition_ref> conditions;
sub createBarSeqExperimentObject($$$$$$$$$){
	my $params = {
		"name" => $_[2],
		"type" => "KBaseRBTnSeq.BarSeqExperiment",
		"workspace" => $_[1],
	};
	$params->{data}->{person} = $_[3];
	$params->{data}->{mutant_lib_name} =   $_[4];
	$params->{data}->{start_date} = $_[5];
	$params->{data}->{sequenced_at} = $_[6];
	$params->{data}->{growth_parameters} = $_[7];
	$params->{data}->{conditions} = $_[8];

	return createObject($params, $_[0]);
}


#creates a new object
#input: $params, ws_client
sub createObject($$){
	my ($params, $serv) = @_;

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
			return $object_info->[6]."/".$object_info->[0]."/".$object_info->[4]."\n";
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



 
