#!/usr/bin/perl -w
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;

my ($input_file,$output_file);

GetOptions(
	'i=s'    => \$input_file,
	'o=s'    => \$output_file,
    'n=i'    => \$length_threshold,
);

defined($input_file) or die("use -i to specify input file");
defined($output_file) or die("use -o to specify output file");
defined($length_threshold) or die("use -n to specify length threshold");

print "The sequences shorter than $length_threshold will be removed!\n";

my $catchseq_seqio_obj = Bio::SeqIO->new(-file=>"$input_file", -format=>'fasta');

open (O,">$output_file");

while(my $seq_obj = $catchseq_seqio_obj->next_seq)
{
	my $display_name = $seq_obj->display_name;
	my $desc = $seq_obj->desc;
	my $seq = $seq_obj->seq;
	my $seq_type = $seq_obj->alphabet;
	my $seq_length = $seq_obj->length;

    my $seqContent = "";
	if($seq_length >= $length_threshold) {
        $seqContent = ">".$display_name." ".$desc."\n".$seq."\n";
        print O $seqContent;
   	} else {
        $seqContent = ">".$display_name." ".$desc."\n".$seq."\n";
        print "The length of the following sequence is shorter than $length_threshold :\n";
        print $seqContent;
        print "\n";
        print "This sequence has been removed!\n";
        print "\n";
    }
}

close O;

print "Success!\n";
