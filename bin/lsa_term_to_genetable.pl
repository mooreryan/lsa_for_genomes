#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Provide term idx's of interest (one per line) and get back
#a tab-delimited map of which genomes contain genes within that term.
#v.0.1.0 by Sean McAllister

# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("t:s:m:g:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-t = term idx of interest (one per line)\n";
        print "-s = seqs.clu.tsv file from MMseqs2.0 clustering\n";
	print "-m = term map file \"td_matrix.idx_to_term_map.txt\"\n";
	print "-g = genome order for output (one per line)\n";
	print "-h = This help message\n\n";
	die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %MMseqsClusters;
my %TermIdx;
my %Genomes;

open(IN, "<$options{s}") or die "\n\nFile $options{s} does not exist or was not given. Try -h for the help file.\n\n";
my @clu_data = <IN>; close(IN);
foreach my $line (@clu_data)
	{   	chomp($line);
		my @data = split('\t', $line);
		my $clu_rep = $data[0];
		my $gene_id = $data[1];
		chomp($clu_rep); chomp($gene_id);
		$MMseqsClusters{$clu_rep}{'cluster_rep'} = $clu_rep;
		$MMseqsClusters{$clu_rep}{'genes_in_cluster'} .= $gene_id.";";
		$MMseqsClusters{$clu_rep}{'count'} += 1;
	}

open(IN2, "<$options{m}") or die "\n\nFile $options{m} does not exist or was not given. Try -h for the help file.\n\n";
my @term_idx = <IN2>; close(IN2);
foreach my $line2 (@term_idx)
	{	chomp($line2);
		my @data2 = split('\t', $line2);
		my $term_idx = $data2[0];
		my $term_rep = $data2[1];
		chomp($term_idx); chomp($term_rep);
		$MMseqsClusters{$term_rep}{'term_idx'} = $term_idx;
		$TermIdx{$term_idx}{'cluster_rep'} = $term_rep;
	}

open(GENOME, "<$options{g}") or die "\n\nFile $options{g} does not exist or was not given. Try -h for the help file.\n\n";
my @genomes = <GENOME>; close(GENOME);
my $unid = 1000001;
foreach my $genome (@genomes)
	{	chomp($genome);
		$Genomes{$unid}{'original'} = $genome;
		$genome =~ s/\W/_/g;
		$genome =~ s/\s/_/g;
		$Genomes{$unid}{'cleaned'} = $genome;
		$unid += 1;
	}

#Print line 1 (headers)
print "Cluster representative annotation\tTerm ID\t";
foreach my $i (sort keys %Genomes)
	{	print "$Genomes{$i}{'cleaned'}\t";		
	}
print "\n";

open(TERM, "<$options{t}") or die "\n\nFile $options{t} does not exist or was not given. Try -h for the help file.";
my @terms_of_interest = <TERM>; close(TERM);
foreach my $term (@terms_of_interest)
	{	my %TempGeneDirectory;
		chomp($term);
		my $prep_cluster_rep = $TermIdx{$term}{'cluster_rep'};
		my @clust2 = split('~', $prep_cluster_rep);
		print "$clust2[1]\t$term\t";
		my $genes_of_interest = $MMseqsClusters{$prep_cluster_rep}{'genes_in_cluster'};
		my @gene_list = split(';', $genes_of_interest);
		foreach my $j (@gene_list)
			{	my @split_gene = split('~', $j);
				$TempGeneDirectory{$split_gene[0]}{'count'} += 1;
			}
		foreach my $genome (sort keys %Genomes)
			{	if (exists $TempGeneDirectory{$Genomes{$genome}{'cleaned'}})
					{	print "$TempGeneDirectory{$Genomes{$genome}{'cleaned'}}{'count'}X\t"	
					}
				else
					{	print "\t";
					}
			}
		print "\n";
	}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -