#!/usr/bin/perl -w


use strict;
use warnings;

no warnings qw/uninitialized/;

use Getopt::Long qw(GetOptionsFromArray :config pass_through no_auto_abbrev bundling);
use File::Basename;


#######################################################################
##			ENV variables				     ##
##   Modify CYTOFPIPE_HOME to point to cytofpipe master directory    ##
#######################################################################

my $cytofpipe_home = '/shared/ucl/depts/cancer/apps/cytofpipe/v1.3';
my $cytofpipe_array_home = '/shared/ucl/depts/cancer/apps/cytofpipe_array/v0.1';

my $rand_id=`od -N 4 -t uL -An /dev/urandom | tr -d " " | tr -d "\n"`;

$ENV{'CYTOFPIPE_HOME'} = $cytofpipe_home;
$ENV{'CYTOFPIPE_ARRAY_HOME'} = $cytofpipe_array_home;

$ENV{'R_MAX_NUM_DLLS'} = 153;
$ENV{'RAND_ID'} = $rand_id;
#######################################################################




my @args=@ARGV;

my $inputdir=''; my $outputdir=''; my $markersfile='';
my $configfile='';
my $groupsfile='';
my $flow='0';
my $cytof='0';
my $displayall='';
my $all='';
my $downsample='';
my $randomtsneSeed='';
my $randomsampleSeed='';
my $randomflowSeed='';
my $iter='';	
GetOptionsFromArray (
    \@args,
    "i=s" => \$inputdir,
    "o=s"   => \$outputdir,
    "m=s"   => \$markersfile,
    "config=s"   => \$configfile,
    "groups=s"   => \$groupsfile,
    "flow"   => \$flow,
    "cytof"   => \$cytof,
    "displayAll"   => \$displayall,
    "all"   => \$all,
    "downsample=i"   => \$downsample,
    "randomSampleSeed"   => \$randomsampleSeed,
    "randomTsneSeed"   => \$randomtsneSeed,
    "randomFlowSeed"   => \$randomflowSeed,
    "iter=i"   => \$iter,
    "<>"   => \&print_usage
) or die "\n";


if ($inputdir eq '' || $outputdir eq '' || $markersfile eq '' || $iter eq ''){         
	usage_clustering("Please check that you are providing a inputdir (-i), outputdir (-o), markersfile (-m) and the desired number of cytofpipe runs (--iter)\n");
	return;
}
if (!-e "$inputdir") {
	usage_clustering("Can't find directory with fcs files <$inputdir>\n");
	return;
}		
if (-e "$outputdir") {
	usage_clustering("The outputdir <$outputdir> already exists, please choose a different outputdir\n");
	return;
}
if (!-e "$markersfile") {
	usage_clustering("Can't find markers file <$markersfile>\n");
	return;
}
if(!&isnum($iter) || ($iter < 1 || $iter > 1000)){
	usage_clustering("<$iter> is not a valid iterations number. Please insert a number between 1 and 1,000\n");
	return;
}	

if($flow eq '1' && $cytof eq '1') {
	usage_clustering("These two parameters [--flow, --cytof] can not be used jointly, please choose one of them (or none for default options)\n");
	return;
}
if($downsample ne ''){
	if($all eq '1' && $downsample ne ''){
		usage_clustering("These two parameters [--all, --downsample NUM] can not be used jointly, please choose one of them (or none for default options)\n");
       	      	return;
	}
	if(!&isnum($downsample) || ($downsample < 100 || $downsample > 100000)){
		usage_clustering("<$downsample> is not a valid downsample numebr. Please insert a number between 100 and 100,000\n");
		return;
	}	
}
if ($configfile ne ''){
	if (!-e "$configfile") {
		usage_clustering("Can't find config file <$configfile>\n");
		return;
	}else{
		check_config_clustering($configfile)
	}
}
if ($groupsfile ne ''){
	if (!-e "$groupsfile") {
		usage_clustering("Can't find groups file <$groupsfile>\n");
		return;
	}
}


#-----------------------------------
#- Remve the --iter param from ARGV
#-----------------------------------
my $index_iter;
my $index_out;
for (my $i=0; $i<scalar @ARGV; $i++){
	if($ARGV[$i] eq '--iter'){
		$index_iter=$i;
	}elsif($ARGV[$i] eq '-o'){
                $index_out=$i;
        }
}
splice @ARGV, $index_iter, 2;


#---------------------------------------------------------------
#- Run first iteration in normal mode and rest in array mode
#---------------------------------------------------------------

my $outoputdir=$ARGV[($index_out+1)];

open(OUTF, ">$ENV{PWD}/run_j${rand_id}.txt");
print OUTF "$ENV{CYTOFPIPE_HOME}/cytofpipe.qsub --clustering "."@ARGV"."\n";
for (my $iter_i=2;$iter_i<=$iter;$iter_i++){
	$ARGV[($index_out+1)]=$outputdir."_".$iter_i;
	print OUTF "$ENV{CYTOFPIPE_HOME}/cytofpipe.qsub --clustering "."@ARGV"." --array\n";
}
print OUTF "$ENV{CYTOFPIPE_ARRAY_HOME}/post_array.pl $ENV{CYTOFPIPE_ARRAY_HOME}\n";
close(OUTF);

my $post_iter=$iter+1;

`qsub -t 1-${iter} -N j${rand_id} -v RAND_ID=j${rand_id} $ENV{CYTOFPIPE_ARRAY_HOME}/submit_cytofpipe_jobs.sh`;
`qsub -t ${post_iter} -hold_jid j${rand_id} -N post_j${rand_id} -cwd -v RAND_ID=j${rand_id} -v OUT=$outputdir -v ITER=$iter $ENV{CYTOFPIPE_ARRAY_HOME}/submit_cytofpipe_jobs.sh`;





#---------------------------------
#- USAGE
#---------------------------------

sub print_usage {

	my @a=@_;

	my $usage0="";
        my $usage1="Program: Cytofpipe_array";
        my $usage2= "Version: 0.1";
        my $usage3= "Contact: Lucia Conde <l.conde\@ucl.ac.uk>";
        my $usage4="";
        my $usage5="Usage:   cytofpipe_array -i DIR -o DIR -m FILE --iter NUM [options]";
        my $usage6="";
        my $usage7="Required: -i DIR    Input directory with the FCS files";
        my $usage8="	      -o DIR    Output directory where results will be generated";
        my $usage9="	      -m FILE   File with markers that will be selected for clustering";
        my $usage10="	      --iter NUM   Number of cytofpipe --clustering runs";
        my $usage11="Options: --config FILE                     Configuration file to customize the analysis";
        my $usage12="	      --flow | --cyto           Flow cytometry data (transformation = autoLgcl) or Cytof data (transformation = cytofAsinh) [--cytof]";
        my $usage13="	      --all | --downsample NUM  Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]";
        my $usage14="	      --displayAll                      Display all markers in output files [NULL]";
        my $usage15="	      --randomSampleSeed                    Use a random sampling seed instead of default seed used for reproducible expression matrix merging [NULL]";
        my $usage16="	      --randomTsneSeed                      Use a random tSNE seed instead of default seed used for reproducible tSNE results [NULL]";
        my $usage17="	      --randomFlowSeed                      Use a random flowSOM seed instead of default seed used for reproducible flowSOM results [NULL]";
        my $usage18="	      --groups FILE                      Get marker level plots for groups of samples [NULL]";
        my $usage19="";
    
	print "$usage0\n$usage1\n$usage2\n$usage3\n";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n$usage11\n";
	print "$usage12\n$usage13\n$usage14\n$usage15\n";
	print "$usage16\n$usage17\n$usage18\n$usage19\n";

	die "ERROR: Invalid argument '@a'. Please check above the usage of cytofpipe_array\n";
}


sub usage_clustering {

	my $error=shift;

	die qq(
	Program: Cytofpipe_array
	Version: 0.1
	Contact: Lucia Conde <l.conde\@ucl.ac.uk>

	Usage:   cytofpipe_array -i DIR -o DIR -m FILE --iter NUM [options]	
	Required: -i DIR    Input directory with the FCS files
		      -o DIR    Output directory where results will be generated
		      -m FILE   File with markers that will be selected for clustering
		      --iter NUM   Number of cytofpipe --clustering runs
	Options: --config FILE                     Configuration file to customize the analysis
		      --flow | --cyto           Flow cytometry data (transformation = autoLgcl) or Cytof data (transformation = cytofAsinh) [--cytof]
		      --all | --downsample NUM  Use all events in the analysis or downsample each FCS file to the specified number of events (with no replacement for sample with events < NUM) [--downsample 10000]
		      --displayAll                      Display all markers in output files [NULL]
		      --randomSampleSeed                    Use a random sampling seed instead of default seed used for reproducible expression matrix merging [NULL]
		      --randomTsneSeed                      Use a random tSNE seed instead of default seed used for reproducible tSNE results [NULL]
		      --randomFlowSeed                      Use a random flowSOM seed instead of default seed used for reproducible flowSOM results [NULL]
		      --groups FILE                      Get marker level plots for groups of samples [NULL]

	ERROR: $error

	);

}

sub usage_clustering_config {

	my $error=shift;
	die qq(

        Program: Cytofpipe_array
	Version: 0.1
	Contact: Lucia Conde <l.conde\@ucl.ac.uk>

	------------------
	CONFIG file format
	------------------
	[ clustering ]  		#-- MANDATORY FIELD, IT SHOULD BE THE FIRST LINE OF THE CONFIG FILE

	GATING = yes|no  		#-- MANDATORY FIELD

	TRANSFORM = autoLgcl, cytofAsinh, logicle, arcsinh or none  #-- MANDATORY FIELD
	MERGE = ceil, all, min, or fixed			  	#-- MANDATORY FIELD
	DOWNSAMPLE = number between 100 and 100000  		#-- MANDATORY FIELD if MERGE = fixed or ceil

	#- DimRed method (tSNE) parameters:
	PERPLEXITY = 30
	THETA = 0.5
	MAX_ITER = 1000

	#- Clustering methods:
	PHENOGRAPH = yes|no
	CLUSTERX = yes|no
	DENSVM = yes|no
	FLOWSOM = yes|no
	FLOWSOM_K = number between 2 and 50  		#-- MANDATORY FIELD if FLOWSOM = YES

	#- Additional visualization methods:
	PCA = yes|no
	ISOMAP = yes|no

	#- Other:
	DISPLAY_ALL = yes|no
	RANDOM_SAMPLE_SEED = yes|no
	RANDOM_TSNE_SEED = yes|no
	RANDOM_FLOW_SEED = yes|no
	---------------------------------

	ERROR: $error

	);

}


sub isnum ($) {
    return 0 if $_[0] eq '';
    $_[0] ^ $_[0] ? 0 : 1
}

sub check_config_clustering {
	my $config = shift(@_);

	my $gating;my $transform;my $merge;my $downsample;my $flowsom;my $flowsom_k;
	my $perplexity; my $theta; my $max_iter;

	local $/ = undef;
	open(INF, "$config");
	my $content = <INF>;
	my @lines = split /\r\n|\n|\r/, $content;

	my $first_line=1;
	foreach my $line(@lines){
		chomp $_;
		if($first_line == 1 && $line !~ /^\s*\[\s*clustering\s*\]\s*$/i){
			usage_clustering_config("Invalid config file. Please make sure that the first line of the config file is \"[ clustering ]\"");
			return;
		}
		if($line=~/^GATING\s*\=\s*(.*)\s*$/){
			$gating=$1;
			if($gating !~/^yes$/i && $gating !~/^no$/i){
				usage_clustering_config("\"$gating\" **************** is not a valid GATING option in <$config>. Please correct the config file and choose if you want to include atomatic gating on your analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^TRANSFORM\s*\=\s*(.*)\s*$/){
			$transform=$1;
			if($transform ne "autoLgcl" && $transform ne "cytofAsinh" && $transform ne "logicle" && $transform ne "arcsinh" && $transform ne "none"){
				usage_clustering_config("Can't recognize \"$transform\" as a valid transformation method in <$config>. Please correct the config file and choose one of the available methods (\"autoLgcl\", \"cytofAsinh\", \"logicle\", \"arcsinh\", \"none\"\) or omit the config file to run with default parameters");
				return;
			}
		}
		if($line=~/^MERGE\s*\=\s*(.*)\s*$/){
			$merge=$1;
			if($merge ne "ceil" && $merge ne "all" && $merge ne "min" && $merge ne "fixed"){
				usage_clustering_config("Can't recognize \"$merge\" as a valid merge method in <$config>. Please correct the config file and choose one of the available methods (\"ceil\", \"all\", \"min\", \"fixed\"\) or omit the config file to run with default parameters");
				return;
			}
		}
		if($line=~/^PCA\s*\=\s*(.*)\s*$/){
			my $pca=$1;
			if($pca !~/^yes$/i && $pca !~/^no$/i){
				usage_clustering_config("\"$pca\" is not a valid PCA option in <$config>. Please correct the config file and choose if you want to include PCA for cluster visualization in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^ISOMAP\s*\=\s*(.*)\s*$/){
			my $isomap=$1;
			if($isomap !~/^yes$/i && $isomap !~/^no$/i){
				usage_clustering_config("\"$isomap\" is not a valid ISOMAP option in <$config>. Please correct the config file and choose if you wnat to include ISOMAP for cluster visualization in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^PHENOGRAPH\s*\=\s*(.*)\s*$/){
			my $phenograph=$1;
			if($phenograph !~/^yes$/i && $phenograph !~/^no$/i){
				usage_clustering_config("\"$phenograph\" is not a valid PHENOGRAPH option in <$config>. Please correct the config file and choose if you want to include PHENOGRAPH as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^CLUSTERX\s*\=\s*(.*)\s*$/){
			my $clusterx=$1;
			if($clusterx !~/^yes$/i && $clusterx !~/^no$/i){
				clustering_config("\"$clusterx\" is not a valid CLUSTERX option in <$config>. Please correct the config file and choose if you want to include CLUSTERX as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^DENSVM\s*\=\s*(.*)\s*$/){
			my $densvm=$1;
			if($densvm !~/^yes$/i && $densvm !~/^no$/i){
				usage_clustering_config("\"$densvm\" is not a valid DENSVM option in <$config>. Please correct the config file and choose if you want to include DENSVM as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^FLOWSOM\s*\=\s*(.*)\s*$/){
			$flowsom=$1;
			if($flowsom !~/^yes$/i && $flowsom !~/^no$/i){
				usage_clustering_config("\"$flowsom\" is not a valid FLOWSOM option in <$config>. Please correct the config file and choose if you want to include FLOWSOM as clustering method in the analysis \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^DOWNSAMPLE\s*\=\s*(.*)\s*$/){
			$downsample=$1;
		}
		if($line=~/^FLOWSOM_K\s*\=\s*(.*)\s*$/){
			$flowsom_k=$1;
		}
		if($line=~/^PERPLEXITY\s*\=\s*(.*)\s*$/){
			$perplexity=$1;
			if($perplexity !~/^\d+$/ || $perplexity < 5 || $perplexity > 50){
				usage_clustering_config("\"$perplexity\" is not a valid PERPLEXITY option in <$config>. Please correct the config file and choose a value between 5 and 50");
				return;
			}
		}
		if($line=~/^THETA\s*\=\s*(.*)\s*$/){
			$theta=$1;
			if($theta !~/([0-9]*[.])?[0-9]+$/ || $theta < 0 || $theta > 1){
				usage_clustering_config("\"$theta\" is not a valid THETA option in <$config>. Please correct the config file and choose a value between 0 and 1");
				return;
			}
		}
		if($line=~/^MAX_ITER\s*\=\s*(.*)\s*$/){
			$max_iter=$1;
			if($max_iter !~/^\d+$/ || $max_iter < 100 || $max_iter > 5000){
				usage_clustering_config("\"$max_iter\" is not a valid MAX_ITER option in <$config>. Please correct the config file and choose a value between 100 and 5000");
				return;
			}
		}
		if($line=~/^DISPLAY_ALL\s*\=\s*(.*)\s*$/){
			my $display=$1;
			if($display !~/^yes$/i && $display !~/^no$/i){
				usage_clustering_config("\"$display\" is not a valid DISPLAY_ALL option in <$config>. Please correct the config file and choose if you want to display all the markers in the output files \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^RANDOM_SAMPLE_SEED\s*\=\s*(.*)\s*$/){
			my $sampleSeed=$1;
			if($sampleSeed !~/^yes$/i && $sampleSeed !~/^no$/i){
				usage_clustering_config("\"$sampleSeed\" is not a valid SAMPLE_SEED option in <$config>. Please correct the config file and choose if you want to use the default sampling seed for reproducible expression matrix merging \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^RANDOM_TSNE_SEED\s*\=\s*(.*)\s*$/){
			my $tsneSeed=$1;
			if($tsneSeed !~/^yes$/i && $tsneSeed !~/^no$/i){
				usage_clustering_config("\"$tsneSeed\" is not a valid TSNE_SEED option in <$config>. Please correct the config file and choose if you want to use the default tSNE seed for repreducible tSNE results \(\"YES or NO\"\)");
				return;
			}
		}
		if($line=~/^RANDOM_FLOW_SEED\s*\=\s*(.*)\s*$/){
			my $flowSeed=$1;
			if($flowSeed !~/^yes$/i && $flowSeed !~/^no$/i){
				usage_clustering_config("\"$flowSeed\" is not a valid FLOW_SEED option in <$config>. Please correct the config file and choose if you want to use the default flowSOM seed for repreducible flowSOM results \(\"YES or NO\"\)");
				return;
			}
		}
		$first_line++;
	}
	close(INF);
	if(!$gating || $gating eq ''){
		usage_clustering_config("Gating parameter not found in <$config>. Please correct the config file and enter a valid Gating option \(\"YES\" or \"NO\"\) or omit the config file to run with default parameters \(GATING = NO\)");
		return;
	}
	if(!$transform || $transform eq ''){
		usage_clustering_config("Transformation parameter not found in <$config>. Please correct the config file and enter a valid transformation method \(\"autoLgcl\", \"cytofAsinh\", \"logicle\", \"arcsinh\", \"none\"\) or omit the config file to run with default parameters \(TRANSFORMATION = arcsinh\)");
		return;
	}
	if(!$merge || $merge eq ''){
		usage_clustering_config("Merge parameter not found in <$config>. Please correct the config file and enter a valid merge method \(\"ceil\", \"all\", \"min\", \"fixed\"\) or omit the config file to run with default parameters \(MERGE = fixed, DOWNSAMPLE = 10000\)");
		return;
	}
	if($merge =~ /^fixed$/i || $merge =~ /^ceil$/i){
		if(!$downsample || $downsample eq '' || $downsample !~/^\d+$/ || $downsample < 100 || $downsample > 100000){
			if(!$downsample || $downsample eq ''){
				usage_clustering_config("Downsample parameter not found in <$config>. Please correct the config file and enter a valid size between 50 and 100000");
				return;
			}elsif($downsample !~/^\d+$/ || $downsample < 100 || $downsample > 100000){
				usage_clustering_config("Can't recognize \"$downsample\" as a valid downsample number in <$config>. Please correct the config file and choose a downsample size between 100 and 100000");
				return;
			}
		}
	}
	if($flowsom =~ /^yes$/i){
		if(!$flowsom_k || $flowsom_k eq '' || $flowsom_k !~/^\d+$/ || $flowsom_k < 2 || $flowsom_k > 50){
			if(!$flowsom_k || $flowsom_k eq ''){
				usage_clustering_config("FlowSOM_k parameter not found in <$config>. Please correct the config file and enter a valid cluster number between 2 and 50");
				return;
			}elsif($flowsom_k !~/^\d+$/ || $flowsom_k < 2 || $flowsom_k > 50){
				usage_clustering_config("Can't recognize \"$flowsom_k\" as a valid number of FlowSOM clusters in <$config>. Please correct the config file and choose a cluster number between 2 and 50");
				return;
			}
		}
	}

}

