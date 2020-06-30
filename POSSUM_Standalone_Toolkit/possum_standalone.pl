#!/usr/bin/perl
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

#use strict;
#use warnings;

use Getopt::Long;
use File::Path;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;


my ($input_file,$output_file,$type,$pssmdir,$firstArgument,$secondArgument,$header);

GetOptions(
	'i=s'    => \$input_file,
	'o=s'    => \$output_file,
	't=s'    => \$type,
	'p=s'    => \$pssmdir,
	'a=i'    => \$firstArgument,
	'b=i'    => \$secondArgument,
	'h=s'	 => \$header,
);

defined($input_file) or die("use -i to specify input file");
defined($output_file) or die("use -o to specify output file");
defined($type) or die("use -t to specify encoding algorithm.\nAvailable algorithms:\naac_pssm, d_fpssm, smoothed_pssm, ab_pssm, pssm_composition, rpm_pssm, s_fpssm, dpc_pssm, k_separated_bigrams_pssm, tri_gram_pssm, eedp, tpc, edp, rpssm, pse_pssm, dp_pssm, pssm_ac, pssm_cc, aadp_pssm, aatp or medp");
defined($pssmdir) or die("use -p to specify the directory of pssm files");


my $catchseq_seqio_obj = Bio::SeqIO->new(-file=>"$input_file", -format=>'fasta');

print "Begin to check input sequences \n";
while(my $seq_obj = $catchseq_seqio_obj->next_seq)
{
	my $display_name = $seq_obj->display_name;
	my $desc = $seq_obj->desc;
	my $seq = $seq_obj->seq;
	my $seq_type = $seq_obj->alphabet;
	my $seq_length = $seq_obj->length;

	if($seq_length < 50)
	{
        die("This input file contains sequence(s) shorter than 50, please remove it/them using utils/removeShortSequences.pl and try again!\n\n Usage example: \n perl removeShortSequences.pl -i example.fasta -o example_corrected.fasta -n 50\n");
   	}

	if((index($seq,"B")!=-1) or (index($seq,"J")!=-1) or (index($seq,"O")!=-1) or (index($seq,"U")!=-1) or (index($seq,"X")!=-1) or (index($seq,"Z")!=-1))
	{
        die("Sequence(s) should not contain illegal characters, such as 'B', 'J', 'O', 'U', 'X' and 'Z'. Please remove it/them using utils/removeIllegalSequences.pl and try again!\n\n Usage example: \n perl removeIllegalSequences.pl -i example.fasta -o example_corrected.fasta\n");
   	}
}
print "End to check input sequences \n";

my $my_env_osname = "$^O";
if($my_env_osname eq 'MSWin32'){
	print "\n";
	print "Current platform is Windows OS \n";
	print "\n";
} else {
	print "\n";
	print "Current platform is Unix/Linux/Mac OS X \n";
	print "\n";
}


my $output_file_prefix = "";
my $output_file_suffix = "";

@output_file_vec = split(/\./,$output_file);
$output_file_prefix = $output_file_vec[0];
$output_file_suffix = ".".$output_file_vec[1];

my @suffixes = $output_file_suffix;
my ($ouputFileName,$dirName,$suffixName) = fileparse($output_file, $output_file_suffix);

print "The specified output file name is $ouputFileName\n";
print "The specified output file directory is $dirName\n";
#print "The specified output file suffix is $suffixName\n";
print "\n";

print "Start to create $dirName directory (if needed)\n";
mkpath($dirName);
print "Success to create $dirName directory (if needed)\n";

if($header eq ''){
	$header="T";
}
if($type eq "aac_pssm") {
	if($header eq "T"){
		print "Start to extract aac_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aac_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";

			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract aac_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aac_pssm feature\n";
	}
} elsif($type eq "d_fpssm") {
	if($header eq "T"){
		print "Start to extract d_fpssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract d_fpssm feature\n";

		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract d_fpssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Sucess to extract d_fpssm feature\n";
	}
} elsif($type eq "smoothed_pssm"){
	if($firstArgument==''){
		$firstArgument=7;
	}
	if($secondArgument==''){
		$secondArgument=50;
	}
	my $dimension2 = $secondArgument*20;
	if($header eq "T"){
		print "Start to extract smoothed_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b $secondArgument\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b $secondArgument`;
		print "Success to extract smoothed_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension2 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension2 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif($header eq "F"){
		print "Start to extract smoothed_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b $secondArgument\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b $secondArgument`;
		print "Success to extract smoothed_pssm feature\n";
	}
} elsif($type eq "ab_pssm") {

	if($header eq "T"){
		print "Start to extract ab_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract ab_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract ab_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract ab_pssm feature\n";
	}
}

elsif($type eq "pssm_composition")
{
	if($header eq "T"){
		print "Start to extract pssm_composition feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract pssm_composition feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type  -n  400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type   -n  400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract pssm_composition feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract pssm_composition feature\n";
	}
}

elsif($type eq "rpm_pssm")
{
	if($header eq "T"){
		print "Start to extract rpm_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract rpm_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract rpm_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract rpm_pssm feature\n";
	}
}

elsif($type eq "s_fpssm")
{
	if($header eq "T"){
		print "Start to extract s_fpssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract s_fpssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract s_fpssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract s_fpssm feature\n";
	}
}

elsif($type eq "dpc_pssm")
{
	if($header eq "T"){
		print "Start to extract dpc_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract dpc_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract dpc_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract dpc_pssm feature\n";
	}
}

elsif($type eq "k_separated_bigrams_pssm")
{
	if($firstArgument==''){
		$firstArgument=1;
	}
	if($secondArgument==''){
		$secondArgument=0;
	}

	if($header eq "T"){
		print "Start to extract k_separated_bigrams_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract k_separated_bigrams_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;

			}
		}
	} elsif ($header eq "F"){
		print "Start to extract k_separated_bigrams_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract k_separated_bigrams_pssm feature\n";
	}

}

elsif($type eq "tri_gram_pssm")
{
	if($header eq "T"){
		print "Start to extract tri_gram_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract tri_gram_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 8000 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 8000 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract tri_gram_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract tri_gram_pssm feature\n";
	}
}

elsif($type eq "eedp")
{
	if($header eq "T"){
		print "Start to extract eedp feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract eedp feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract eedp feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract eedp feature\n";
	}
}

elsif($type eq "tpc")
{
	if($header eq "T"){
		print "Start to extract tpc feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract tpc feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 400 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract tpc feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract tpc feature\n";
	}
}

elsif($type eq "edp")
{
	if($header eq "T"){
		print "Start to extract edp feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract edp feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 20 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract edp feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract edp feature\n";
	}
}

elsif($type eq "rpssm")
{
	if($header eq "T"){
		print "Start to extract rpssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract rpssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 110 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 110 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract rpssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract rpssm feature\n";
	}
}

elsif($type eq "pse_pssm")
{
	if($firstArgument==''){
		$firstArgument=1;
	}
	if($secondArgument==''){
		$secondArgument=0;
	}

	if($header eq "T"){
		print "Start to extract pse_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pse_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 40 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 40 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}

	} elsif ($header eq "F"){
		print "Start to extract pse_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pse_pssm feature\n";
	}
}

elsif($type eq "dp_pssm")
{
	if($firstArgument==''){
		$firstArgument=5;
	}
	if($secondArgument==''){
		$secondArgument=0;
	}

	my $dimension15 = ($firstArgument+1)*40;
	if($header eq "T"){
		print "Start to extract dp_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract dp_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension15 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension15 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract dp_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract dp_pssm feature\n";
	}
}

elsif($type eq "pssm_ac")
{
	if($firstArgument==''){
		$firstArgument=10;
	}
	if($secondArgument==''){
		$secondArgument=0;
	}

	my $dimension16 = $firstArgument*20;
	if($header eq "T"){
		print "Start to extract pssm_ac feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pssm_ac feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension16 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension16 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract pssm_ac feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pssm_ac feature\n";
	}
}

elsif($type eq "pssm_cc")
{
	if($firstArgument==''){
		$firstArgument=10;
	}
	if($secondArgument==''){
		$secondArgument=0;
	}

	my $dimension17 = $firstArgument*380;
	if($header eq "T"){
		print "Start to extract pssm_cc feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pssm_cc feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension17 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n $dimension17 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract pssm_cc feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a $firstArgument -b 0`;
		print "Success to extract pssm_cc feature\n";
	}
}

elsif($type eq "aadp_pssm")
{
	if($header eq "T"){
		print "Start to extract aadp_pssm feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aadp_pssm feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract aadp_pssm feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aadp_pssm feature\n";
	}
}

elsif($type eq "aatp")
{
	if($header eq "T"){
		print "Start to extract aatp feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aatp feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract aatp feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract aatp feature\n";
	}
}

elsif($type eq "medp")
{
	if($header eq "T"){
		print "Start to extract medp feature\n";
		print "python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o tmp/$ouputFileName\_no_header -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract medp feature\n";
		if (-e "tmp/$ouputFileName\_no_header"){
			print "\n";
			print "Adding header\n";
			print "python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file\n";
			print `python src/headerHandler.py -i tmp/$ouputFileName\_no_header -p $type -n 420 -o $output_file`;
			if($my_env_osname eq 'MSWin32'){
				print "\n";
				print "Deleting tmp\\$ouputFileName\_no_header\n";
				print "del tmp\\$ouputFileName\_no_header\n";
				print `del tmp\\$ouputFileName\_no_header`;
			} else {
				print "\n";
				print "Deleting tmp/$ouputFileName\_no_header\n";
				print "rm tmp/$ouputFileName\_no_header\n";
				print `rm tmp/$ouputFileName\_no_header`;
			}
		}
	} elsif ($header eq "F"){
		print "Start to extract medp feature\n";
		print "python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0\n";
		print `python src/possum.py -i $input_file -o $output_file -t $type -p $pssmdir -a 0 -b 0`;
		print "Success to extract medp feature\n";
	}
}

else {
	die("Please specify an available algorithm:\naac_pssm, d_fpssm, smoothed_pssm, ab_pssm, pssm_composition, rpm_pssm, s_fpssm, dpc_pssm, k_separated_bigrams_pssm, tri_gram_pssm, eedp, tpc, edp, rpssm, pse_pssm, dp_pssm, pssm_ac, pssm_cc, aadp_pssm, aatp or medp");
}
