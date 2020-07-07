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
);

defined($input_file) or die("use -i to specify input file");
defined($output_file) or die("use -o to specify output file");

#print "The sequences containing illegal characters, such as 'B', 'J', 'O', 'U', 'X' and 'Z', will be removed!\n";

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

    if((index($seq, "B")==-1) and (index($seq, "J")==-1) and (index($seq, "O")==-1) and (index($seq, "U")==-1) and (index($seq, "X")==-1) and (index($seq, "Z")==-1))
    {
        $seqContent = ">".$display_name." ".$desc."\n".$seq."\n";
        print O $seqContent;
    }
    else
    {
		$seqContent = ">".$display_name." ".$desc."\n".$seq."\n";

		print "The following sequence contains illegal characters, such as 'B', 'J', 'O', 'U', 'X' and 'Z' :\n";
        print $seqContent;
        print "\n";
        print "This sequence has been removed!\n";
        print "\n";
    }
}

close O;

print "Success!\n";
