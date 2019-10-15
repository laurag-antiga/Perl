#!/usr/bin/perl

use strict;
use warnings;

  ##################################################
  # PERL PROJECT                                   #
  # Laura González Antiga                          #
  # MSc Bioinformatics for Health Science, 2018-19 #
  ##################################################

#Declaring variables:
my (@values, @header, %data, %secondary_matrix, @order, %matrix, $i, $j, %indexes, $key, $position);


### STORING THE INPUT DATA ###
#Opening files: 
open (INPUT, $ARGV[0]) || die "File not opened; please pass the input as a text file in $ARGV[0]";

while (<INPUT>) {
    chomp;
    if ($_ =~ /^(\d)/) {            #All rows starting with a number will be part of the matrix
        @values= split(/\s+/, $_);
      pop @values;                  #We delete the letter of the last column
      $key = shift @values;         #We separate the first element (number of row) as the key of the has
      $data{$key} = [ @values ];    #And save an array of the scores  as a value of each key
        }
    elsif ($_ =~ /^P0/) {
        @header= split ('\s+', $_); #We keep the order of letters by column in a different variable
    }
  }

close INPUT;




shift @header;      #We just need the order of letters, not the row ID
print "Original matrix:\n@header\n\n";

# A loop just to check that our data-matrix was correctly filled:
# (It won't be sorted, but the correct row order is kept as key)
for my $row (keys %data ) {
    print "@{ $data{$row} }\n"; 
}

#print "\n\nMatrix with sorted scores per row:\n";

my $reference= \%data; 

### SORTING LOOP ###
#The next loop will sort each row by decreasing order of scores
#The first commented print shows each row's values ordered
#The second commented print shows the indexes from 0 to 3 in which the letters were reordered
#Uncomment both to show relation between each row's scores and indexes

foreach my $key (keys %data) {
    $secondary_matrix{$key} = [ sort {$b <=> $a } @{$reference->{$key} }];
    # We save in a secondary matrix the same key for each row, but he scores sorted decreasingly
    #print "@{ $secondary_matrix{$key} }\n";
    
    #And then we save the same order, but indexed from 0-1 in a different matrix:
    #Uncomment the print to show the index rows
    $indexes{$key} = [ sort {$reference->{$key}->[$b] <=>  $reference->{$key}->[$a] } 0 .. $#header ];
    #print "Index: @{ $indexes{$key} }\n"
   }

 
print "\n\nMatrix with sorted rows and scores:\n";

### ORDER VECTOR ###
# Next, we order rows by decreasing order of the first column of each, 
# and save the key (number of row) ordered in @order variable, so we can have a vector
# with the order in which we have to process the data.

for my $value ( @order= sort { $secondary_matrix{$b}->[0] <=> $secondary_matrix{$a}->[0] }  keys %secondary_matrix ) {
    print "@{ $secondary_matrix{$value} }\n";
        
}
    
 # All the script works with the @order vector, that corresponds to the
# keys of the hash and to the number of the row were the scores were 
# stored in the initial matrix. But it will come handy to have the same
# vector storing the same positions in the same order, but decreased by one unit,
# because positions in arrays start at 0 and ends in n-1
# That's what this loop makes. Uncomment to visualize differences between @order and @positions
my @positions;

    for (my $i = 0; $i < scalar @order;  $i++) {
        $positions[$i] = $order[$i]-1;
    }   
    #print "Matrix order vector: @order\nArray positions vector:@positions\n\n";
print "\n The order is: @order\n\n";


print "Nucleotide matrix (function input):\n";

### SCORE-NUCLEOTIDE EQUIVALENCE ###
#Finally, this loop transforms each score in its corresponding nucleotide
#and it will be the function's input
#Uncomment the first print to see the corresponding score-nucleotide of each row
foreach my $key (keys %secondary_matrix) {
    #We save the sorted rows in a new matrix, to keep both the scores ones and 
    #the nucleotides one
    $matrix{$key}  = [@{ $secondary_matrix{$key} }] ;
     for ($j=0; $j < scalar @header ; $j++ ) {
        #And then we substitute each position of each row by the nucleotide
        #corresponding to the %index matrix's number of that position in the @header vector
        $matrix{$key}->[$j]  =~ s{(.*)}{$header[$indexes{$key}->[$j]]} ;
          }
     #print "first: @{ $secondary_matrix{$key} }\nsecond: ";
     print "@{ $matrix{$key} }\n";
 } 

print "\n\n\n";


##### Now we are setting some variables we will eventually use:

# Empty counter and set threshold
my $count = 0;
my $threshold = 10000;

### PERMUTATOR EXECUTION ###
my ($result, $counter) = permutator(\@order, \%matrix, \$count, \$threshold);


### AFTER THE FUNCTION ###
# Delete surplus
# This is needed because the function works with rows, so it addes 4 nucleotides each time
# To make an exact thereshold of 10000, some sequences will be incomplete
# With this cleaning, we get only the 10000 complete sequences (listed in the @total array)
my $residue = $threshold-scalar @$result ;
my @total = splice @$result ,0, $residue;

#print $residue."\n@total\n";


# We want to save the sequences and the scores in a new file:
open (RESULT, '>GONZALEZ_LAURA_PER1819_FINALSCRIPT.txt') || die "Result file not creatd";
print RESULT "  ##################################################
  # PERL PROJECT   Results                         #
  # Laura González Antiga                          #
  # MSc Bioinformatics for Health Science, 2018-19 #
  ##################################################\n\n\n";


### CALCULATING THE SCORE ###
# With our sequences listed, we now have to reorder the nucleotides so that they will go back 
#to it's position (the key of the matrix)
foreach my $var (@total){

   # We separate each nucleotide as an element
    my @line = split //, $var;
    print RESULT "Sorted by score line @line"."\n";

    my $score;
    

    # We calculate the score of each nucleotide:
    for (my $i = 0; $i < scalar @line;  $i++) {

        #For each position (nucleotide)
        my $letter = $line[$i];
            

        #And following the @header index:       
        if ($letter eq $header[0]) {$position=0;}
        elsif ($letter eq $header[1]) {$position=1;}
        elsif ($letter eq $header[2]) {$position=2;}
        else {$position=3;}
        
        # The value is the score stored in the $order[$i] row and the corresponding 
        #nucleotide ($position) column, and we add one score per nucleotide
        # of the seuqences each round of the loop
        my $value = ${$data{$order[$i]}}[$position];
        
        $score += $value;
        
    }

    # To end this loop, we print both the correct (ordered sequence) and the score:
   my @line_sorted = @line[@positions]; 
   print RESULT "Correct sequence     @line_sorted"."\n";
   print RESULT "Score = $score\n\n";

   
}


#A final print, to know how many sequences have been printed and 
# which was the thereshold
print "Threshold = ".$threshold."\n\n";
print "Amount of sequences = ". scalar @total. "\n";







####### PERMUTATOR SUBROUTINE #######

sub permutator {
    
    # Obtain variables
	my ($aref, $href, $counter, $threshold) = @_;
	my @order = @$aref;
	my %input = %$href;
	if ($$counter  < $$threshold){
    #	print $$counter;
    
  	#print "\n\n ++++ INITIALIZING permutator... ++++\n\n";
    
    # Check @order array size
    my $arrsize = scalar @order;
    #print " *** Order size = $arrsize *** \n  ";
    
    # The current line is stored, corresponding to the first element in @order
    my @current = @{$input{$order[0]}};

    # Depending on the size, re-call permutator or reurn baseline
    if ($arrsize > 1){

        # The first element in @order is discarded to be the input in the next permutator runs
        shift @order;
        
        # The new strings for all the nucleotides will be stored
        my @allstrings; 
        
        # For each nucleotide in my current line, I concatenate the result of a permutator run with lower-level elements
        foreach my $nucleotide (@current){
     #       print " \n --- We are working with $nucleotide in $arrsize --- \n";
            
            # The new strings for this nucleotide will be stored 
            my @newstrings ;
            
            # Next permutator run
            my ($run, $counter) = permutator(\@order, \%input, \$$counter, \$$threshold);
            
            # For each string in the permutator run
            foreach my $string (@$run){
                my $concat = $nucleotide.$string;
    #            print "Concatenated $concat\n";
                push @newstrings, $concat;
            }
            
            push @allstrings, @newstrings;
   #         print "\n --- $nucleotide in $arrsize is ended --- \n"t'
            
        }
        
        
        # Return current concatenated
        return (\@allstrings, $$counter);
  #      print "*** Ended current $arrsize = @current*** \n\n";
        

    }
    else{
        # The current line is returned
        $$counter = $$counter + 4;
        return (\@current, $$counter);
 #       print "*** Ended last current = @current***\n";
        #return @current;
        
    }
	}else{
	    my @dummy = ("");
        return (\@dummy, $$counter);
	}

}

