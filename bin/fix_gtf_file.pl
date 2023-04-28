#! /usr/bin/perl -w

#  Copyright (C) 2014-02-20 Stefan Lang

#  This program is free software; you can redistribute it 
#  and/or modify it under the terms of the GNU General Public License 
#  as published by the Free Software Foundation; 
#  either version 3 of the License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful, 
#  but WITHOUT ANY WARRANTY; without even the implied warranty of 
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License 
#  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1 get_sequences_4_bed_file.pl

This script gets you all sequences for a given bed file. Make sure you use the right genome for that.

To get further help use 'get_sequences_4_bed_file.pl -help' at the comman line.

=cut

use Getopt::Long;
use strict;
use warnings;

use stefans_libs::fastaDB;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my $includes = "-I " . join( " -I ", @INC );

my ( $help, $debug, $file, $outfile );

Getopt::Long::GetOptions(
	 "-file=s"    => \$file,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $file) {
	$error .= "the cmd line switch -file is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	print helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
 	return "
 $errorMessage
 command line switches for get_fa_sequences_4_bed_file.pl

   -bed_file   :the bed file containing the positions for the fasta table
   -fa         :the fasta genome database
   -outfile    :the fasta db file
      
   -help       :print this help
   -debug      :verbose output
   

"; 
}


my ( $task_description);

$task_description .= 'perl '.$includes.' '.$plugin_path .'/fix_get_file.pl';
$task_description .= " -file $bed_file" if (defined $bed_file);
$task_description .= " -outfile $outfile" if (defined $outfile);

open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );

my ($genomeDB, $bed, @tmp, $acc, $start, $end, $name );

if ( $bed_file =~ m/\.gz$/ ) {
	open( IN, "zcat $bed_file |" )
	  or die "I could not gunzip the zipped file on the fly\n";
}
else {
	open( IN, "<$bed_file" )
	  or die "I could not open file '$bed_file'\n$!\n";
}
my ( $seqname, $source, $feature, $start, $end, $score,	$strand, $frame, $attribute);

my $zero;

open( OUT, ">$outfile" ) or die "I could not create the outfile '$outfile'\n$!\n";

foreach (<IN>) {
	chomp($_);
	if ( substr( $_, 0, 1 ) eq "#" ) {
		next;
	}
	( $seqname, $source, $feature, $start, $end, $score,	$strand, $frame, $attribute) = split ("\t", $_);

	if ($feature == "gene"){
		$zero = $start -1;
	}
	if ( $attribute =~m/gene_name "([\w\d\-\._]*)";/ ){
		$seqname = $1;
	}
	$start -= $zero;
	$end -=  $zero;
	print OUT join( "\t", $seqname, $source, $feature, $start, $end, $score,	$strand, $frame, $attribute)."\n";

}

close ( IN );
close ( OUT );
