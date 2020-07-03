#!/export/projects/bioinfo2/to16r/bin/conda/bin/perl -w
#
# File: little.renameContigs.pl
# Time-stamp: <17-Nov-2015 10:18:29 tdo>
# $Id: $
#
# (c) Copyright: Genome Research Ltd. 2017
#
# Author: Thomas Dan Otto
# Licencse: GNU General Public License v3.0
#
#
# Description: Renames contigs / chromosomes for 
#

use strict;


my $prefix=shift;
my $num=shift;
if (!defined($num)){
	$num=1;
}

while (<STDIN>) {
  if ((/>NODE/) || /^>unitig/ || />scaf/ || />scf/ || />newseq/ || /_00_/) {
	print ">$prefix\_00_$num\n";
	$num++
  }  elsif ((/>.*(_MT)/ || />.*(_API)/)) {

        print ">$prefix$1\n";
  }
  elsif ((/>.*_(\d+)_/)) {

        print ">$prefix\_$1\n";
  }
elsif ((/>.*_(\d+)/)) {
	
	print ">$prefix\_$1\n";
  }
 elsif ((/>.*chr(\d+)/)) {

        print ">$prefix\_$1\n";
  }
  elsif ((/>/)) {
	print ">$prefix\_$num\n";
	$num++
  }
  else {
	print
  }
}
