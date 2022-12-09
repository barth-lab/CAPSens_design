#!/usr/bin/perl

###################################################################
# This script determines what Rosetta score represents the
# cutoff necessary to keep a specified fraction of the top scores
# from a run.
#
# It writes formatted output to STDOUT which contains the cutoff
# determined and the number of scores evaluated so far.
# First 10 characters = cutoff; next 8 = number of scores evaluated.
#
# SEM 3/02; rev. JJG 4/02
#
# 6/02: revised to also filter out bumpiest fraction of scores 
#
#  usage: smart_scorfilter.pl score_fraction bkrep_fraction [scorefile_path]
#
###################################################################



my $fraction_to_keep = shift @ARGV; # real number; should be >0 & <1
my $bump_fraction_to_keep = shift @ARGV; # real number; should be >0 & <1
my $scpath = shift @ARGV | '.';     

$bump_limit = 99999.9;
$score_limit = 99999.9;
$Nstruct = 0;

# Put all the scorefiles in one big array.
@scorefiles=`ls $scpath/*sc`;

# make sure scorefiles exist
if ( @scorefiles != 0 ) {

    # filtering on bumps, find the column
    chomp(@scorefiles);
    $bump_column=`findCol.pl bk_rep $scorefiles[0] | cut -c16-`;
    #printf STDERR "col $bump_column\n";
    
    # get the scores and bumps out of the scorefiles
    foreach $scfile (@scorefiles) {
    	 open (SCOREFILE,"$scfile");
    	 while (<SCOREFILE>)
    	 {
    	     # only count decoys and those rejected due to scores or bumps
    	     next unless /output_decoy|score_failure|bkrep_failure/;
    
    	     $line_score = (split)[1];
    	     push(@scores, $line_score);
    
    	     $line_bumps = (split)[$bump_column];
    	     push(@bumps, $line_bumps);
    	 }
    }
    $Nstruct = $#scores+1;
    
    if ($bump_fraction_to_keep<1) { 
    # eliminate bumpy scores
    	 @sorted_bumps = sort { $a <=> $b } @bumps;
    	 $index = $bump_fraction_to_keep * scalar(@sorted_bumps);
    	 $bump_limit = $sorted_bumps[$index];
    #    printf STDERR "$bump_limit at $index\n";
    
    	 for ($i=0; $i<@scores; $i++) {
    	     push(@filtered_scores,$scores[$i]) if ($bumps[$i] < $bump_limit);
    	 }
    	 @scores = @filtered_scores;	
    }	
    
    # Turn 'em around and cut 'em off.
    @scores = sort {$a <=> $b} @scores;
    $score_limit = $scores[$fraction_to_keep * $#scores];

}


#output
if ($output_file) { 
    open(OUT,">$output_file");
    printf OUT "%15s%7.5f: %10.2f\n","Score limit of ",
           $fraction_to_keep,$score_limit;
    printf OUT "%15s%7.5f: %10.2f\n","BKrep limit of ",
           $bump_fraction_to_keep,$bump_limit;
    printf OUT "%22s: %10d\n","Structures",$Nstruct;
}
else {
    printf "%15s%7.5f: %10.2f\n","Score limit of ",
           $fraction_to_keep,$score_limit;
    printf "%15s%7.5f: %10.2f\n","BKrep limit of ",
           $bump_fraction_to_keep,$bump_limit;
    printf "%22s: %10d\n","Structures",$Nstruct;
}

