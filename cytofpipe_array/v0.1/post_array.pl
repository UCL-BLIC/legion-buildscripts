#!/usr/bin/perl -w


use strict;


my $rand_id=$ENV{RAND_ID};
my $out=$ENV{OUT};
my $iter=$ENV{ITER};

#my $cytofpipe_array_home=$ENV{CYTOFPIPE_ARRAY_HOME};
my $cytofpipe_array_home=$ARGV[0];

system("mkdir $ENV{PWD}/${out}_array_files");

my %clusteringmethod=();
for (my $i=2; $i<=$iter; $i++){
	opendir(DIR, "$ENV{PWD}/${out}_$i");
	my @files= grep { /cytofpipe_.*_cluster_.*.csv/ } readdir(DIR);
	foreach my $file(@files){
		if($file =~ /(cytofpipe_(.*)_cluster_(.*)).csv/){
			my $cl=$2; my $type=$3; my $id=$1;
			$clusteringmethod{$cl}=1;
			if(!-e "$ENV{PWD}/${out}_array_files/$cl/$type"){
				system("mkdir -p $ENV{PWD}/${out}_array_files/$cl/$type");
			}
			system("perl -i -pe \'s/\"\"/\"clusterID\"/ if 1\' $ENV{PWD}/${out}_$i/$file");	# header of first column was empty, add header
			system("perl -i -pe \'s/\",/_${i}\",/ if \$\. > 1\'  $ENV{PWD}/${out}_$i/$file"); # add iteration number to cluster id
			system("mv $ENV{PWD}/${out}_$i/$file $ENV{PWD}/${out}_array_files/$cl/$type/${id}_${i}.csv");
		}
	}
	closedir(DIR);
	system("rm -rf $ENV{PWD}/${out}_$i");
}

my @clust=keys(%clusteringmethod);
my $clust=join(",",@clust);

system("OUT=$ENV{PWD}/$out ITER=$iter CLUST=$clust R CMD BATCH ${cytofpipe_array_home}/merge_files.R");
system("rm -rf merge_files.Rout");
system("rm -rf $ENV{PWD}/*${rand_id}*");

