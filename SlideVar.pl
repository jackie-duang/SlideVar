#! /usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    perl $0 -in <input.fasta> -l <species.list> -w <window size> -con <conserved number> -div <changed number> -bn <background can not conserved species>

    -in input.fasta : sequence aligned in fasta format file
    -l species.list : species list with species marked 
    -w window size : default 20 , the window size for calculating the conserved region
    -con conserved number : default 18 , threshold for conserved region, identity >= conserved number is conserved 
    -div diverged number : default 12 , threshold for divergent region, identity <= diverged number is divergent
    -bn  : default 0 , how many species could be divergent in background species (in case of some species are not conserved in background species because of assembly error or other reasons)
	
    ---
    species list file format, one species per line and foreground species marked with '*':
    human
    mouse
    snake *
    frog 
    caecilian *
    ...

    fasta file format: # species name should not contain '.' '-' '\@' etc. ; '_' is allowed ; the first species should be the reference species
    >species_name
    AAGCTTGGG
    or 
    >species_name.seqId
    AAGCTTGGG
USAGE


my ($in,$window,$speciesList, $conservedNumber, $changedNumber, $canChanged);

GetOptions(
    "in=s" => \$in,
    "w=i" => \$window,
    "l=s" => \$speciesList,
    "con=i" => \$conservedNumber,
    "div=i" => \$changedNumber,
    "bn=i" => \$canChanged,
);

# set default value
$window = 20 if (!defined $window);
$conservedNumber = 18 if (!defined $conservedNumber);
$changedNumber = 12 if (!defined $changedNumber);
$canChanged = 0 if (!defined $canChanged);

# check options
die $usage if (!defined $in or !defined $speciesList);

my %info ;
my %seq ;
my $length = 0;
my %addInfo ;
my $line = 0;
my %mapSites ;
my $ref = '' ;

my $file_prefix = "$in.w${window}c${conservedNumber}d${changedNumber}b${canChanged}";
my $out = "$file_prefix.info";

$/ = ">";
my %original_pos = ();
my %changedInfo = ();
my $last_pos = 0 ;
open I , "< $in";
while (<I>){
    chomp;
    my @a = split /\n+/,$_ ;
    next if @a==0 ;

    # get id and sequence
    my $sp = shift @a ;
    $sp =~ /^(\w+)/;
    $sp = $1 ;
    my $seq = join '',@a ;
    $seq = uc($seq);

    # record the length of the longest id
    if (length($sp) > $length){
        $length = length($sp);
    }

    # record the number of nucleotides in each species
    my @nucls = split //,$seq;
    my $start = 1 ;
    if ($ref eq ''){
        $ref = $sp;
        my $ref_site = 0 ;

        # record the nucleotide at each site in the reference species
        for (my $i=0;$i<@nucls;$i++){
            $mapSites{$i} = $ref_site ;
            next if ($nucls[$i] eq '-');
            $ref_site ++ ;
            $seq{'ref'}{$i} = $nucls[$i];
            $original_pos{$i} = $start + $ref_site - 1;
            $last_pos = $i ;
        }
    }
    else{
        my $window_id = 0 ;
        my $step = 0 ;
        my %map = ();
        my ($I, $D, $S, $G) = (0,0,0,0);
        # record the nucleotide at each site in each species and calculate the number of conserved nucleotides in each window
        for (my $i=0;$i<@nucls;$i++){
            if ((exists $seq{'ref'}{$i}) || ($nucls[$i] ne '-')){
                $step ++ ;
                $map{$step} = $i ;
                if ((exists $seq{'ref'}{$i}) and ($nucls[$i] eq $seq{'ref'}{$i})){
                    $window_id ++ ;
                }
                else{
                    if (!exists $seq{'ref'}{$i} and $nucls[$i] ne '-'){
                        # insert
                        $I ++ ;
                    }
                    elsif (exists $seq{'ref'}{$i} and $nucls[$i] eq '-'){
                        # delete
                        $D ++ ;
                    }
                    else{
                        # substitution
                        $S ++ ;
                    }
                }
                if ($step >= $window){
                    $addInfo{$sp}{$map{$step-$window+1}} = $window_id;
                    $changedInfo{$sp}{$map{$step-$window+1}}{I} = $I ;
                    $changedInfo{$sp}{$map{$step-$window+1}}{D} = $D ;
                    $changedInfo{$sp}{$map{$step-$window+1}}{S} = $S ;
                    $changedInfo{$sp}{$map{$step-$window+1}}{M} = $window_id ;
                    if ((exists $seq{'ref'}{$map{$step-$window+1}}) and ($nucls[$map{$step-$window+1}] eq $seq{'ref'}{$map{$step-$window+1}})){
                        $window_id -- ;
                    }
                    elsif (!exists $seq{'ref'}{$map{$step-$window+1}} and $nucls[$map{$step-$window+1}] ne '-'){
                        # insert
                        $I -- ;
                    }
                    elsif (exists $seq{'ref'}{$map{$step-$window+1}} and $nucls[$map{$step-$window+1}] eq '-'){
                        # delete
                        $D -- ;
                    }
                    else{
                        # substitution
                        $S -- ;
                    }
                }
            }
        }
        # the final window is not complete, add it to the last position
        for (my $i=$step-$window+2;$i<=$step;$i++){
            if (exists $addInfo{$sp}{$map{$i}}){
                # print "$i is already added\n";
                next ;
            }
            $addInfo{$sp}{$map{$i}} = $window_id;
            $changedInfo{$sp}{$map{$i}}{I} = $I ;
            $changedInfo{$sp}{$map{$i}}{D} = $D ;
            $changedInfo{$sp}{$map{$i}}{S} = $S ;
            $changedInfo{$sp}{$map{$i}}{M} = $window_id ;
            if ((exists $seq{'ref'}{$map{$i}}) and ($nucls[$map{$i}] eq $seq{'ref'}{$map{$i}})){
                $window_id -- ;
            }
            elsif (!exists $seq{'ref'}{$map{$i}} and $nucls[$map{$i}] ne '-'){
                # insert
                $I -- ;
            }
            elsif (exists $seq{'ref'}{$i} and $nucls[$i] eq '-'){
                # delete
                $D -- ;
            }
            else{
                # substitution
                $S -- ;
            }
        }
    }
}
close I ;
$/ = "\n";

my %backgroundSpecies = ();
my %foregroundSpecies = ();

# read species list and record the foreground species and background species
open I, "< $speciesList" or die $!;
while (<I>){
    chomp;
    next if $_ =~ /^$/ ;
    my @a = split /\s+/, $_ ;
    next if @a>=2 and $a[1] =~ /pheno/i ;
    $line ++ ;
    $info{$line}{id} = $a[0] ;
    if (@a>1 and ($a[1] =~ /\*|#/ or $a[1] eq '0')){
        $foregroundSpecies{$a[0]} = 1 ;
    }
    else{
        next if ($a[0] eq $ref);
        $backgroundSpecies{$a[0]} = 1 ;
    }
}
close I ;

open O , "> $out";
printf O "%-${length}s", "RefPos";
foreach my $k (sort {$a <=> $b} keys %{$seq{'ref'}}){
    printf O " %-4d",$mapSites{$k} + 1;
}
print O "\n";

printf O "%-${length}s", "WithGapPos";
foreach my $k (sort {$a <=> $b} keys %{$seq{'ref'}}){
    printf O " %-4d",$k + 1;
}
print O "\n";

foreach my $line (sort {$a <=> $b} keys %info){
    my $id = $info{$line}{id};
    printf O "%-${length}s", $id;
    foreach my $k (sort {$a <=> $b} keys %{$seq{'ref'}}){
        if ($id eq $ref){
            printf O " %-4d", $window;
            next ;
        }
        my $title = sprintf "%-4d",$k + 1;
        if (exists $addInfo{$id}{$k}){
            printf O " %-4d", $addInfo{$id}{$k};
        }else{
            printf O (" %-4d", 0);
        }
    }
    print O "\n";
}
close O ;

my %hash = ();
my %tempMap = () ;
my @positions ;
open I , "< $out";
while (<I>){
    chomp;
    my @a = split /\s+/,$_;
    next if @a==0 ;
    if ($a[0] eq 'RefPos'){
        @positions = @a ;
        next ;
    }
    elsif ($a[0] eq 'WithGapPos'){
        my @temps = @a ;
        for (my $i=1;$i<@temps;$i++){
            $tempMap{$positions[$i]} = $temps[$i] - 1 ;
        }
        next ;
    }
    my $species = $a[0] ;
    next if $species eq $ref ;
    for (my $i=1;$i<@a;$i++){
        $hash{$tempMap{$positions[$i]}}{$species} = $a[$i] ;
    }
}
close I ;

open O , "> $out.list";
print O "file\toriginalPositionBasedRef\tTargetDivergedSpecies\tBackgroundNotConservedSpecies\n";
foreach my $pos (sort {$a <=> $b} keys %{$seq{'ref'}}){
    my $notConservedNumber = 0 ;
    my %notConservedSpecies = ();
    my $npos = $original_pos{$pos};

    foreach my $tpSp (sort keys %backgroundSpecies){
        if (exists $hash{$pos}{$tpSp}){
            if ($hash{$pos}{$tpSp} < $conservedNumber){
                $notConservedNumber ++ ;
                $notConservedSpecies{$tpSp} = $hash{$pos}{$tpSp} ;
            }
        }
        else{
            $notConservedNumber ++ ;
            $notConservedSpecies{$tpSp} = 0 ;
        }
    }

    if ($notConservedNumber > $canChanged){
        print O "$in\t$npos\t.\t$notConservedNumber\n";
        next ;
    }

    my %changedSpecies = ();
    foreach my $sp (sort keys %foregroundSpecies){
        if (exists $hash{$pos}{$sp}){
            if ($hash{$pos}{$sp} <= $changedNumber and $hash{$pos}{$sp} >= 0){
                $changedSpecies{$sp} = $hash{$pos}{$sp} ;
            }
        }
        else{
            $changedSpecies{$sp} = 0 ;
        }
    }

    my $light = 1;
    my $changedSpLine = "$in\t$npos\t";
    foreach my $sp (sort keys %changedSpecies){
        # get Insert Delete Substitution number
        my $I = '-' ;
        my $D = $window ;
        my $S = '-' ;
        my $M = '-' ;
        # my $window_count = $tempMap{$pos} ;
        my $window_count = $pos ;
        if (exists $changedInfo{$sp}{$window_count}){
            $I = $changedInfo{$sp}{$window_count}{I} ;
            $D = $changedInfo{$sp}{$window_count}{D} ;
            $S = $changedInfo{$sp}{$window_count}{S} ;
            $M = $changedInfo{$sp}{$window_count}{M} ;
        }
        else{
            $light = 0 ;
        }
        $changedSpLine .= "$sp:M${M}I${I}D${D}S${S};" ;
    }
    $changedSpLine =~ s/;$// ;
    # if ($changedSpLine ){
    if ($changedSpLine !~ /:M/){
        $changedSpLine .= '.';
    }
    my $notConservedSpLine = '#';
    foreach my $sp (sort keys %notConservedSpecies){
        # get Insert Delete Substitution number
        my $I = '-' ;
        my $D = $window ;
        my $S = '-' ;
        my $M = '-' ;
        my $window_count = $pos ;
        if (exists $changedInfo{$sp}{$window_count}){
            $I = $changedInfo{$sp}{$window_count}{I} ;
            $D = $changedInfo{$sp}{$window_count}{D} ;
            $S = $changedInfo{$sp}{$window_count}{S} ;
            $M = $changedInfo{$sp}{$window_count}{M} ;
        }
		else{
            $light = 0 ;
		}
        $notConservedSpLine .= "$sp:M${M}I${I}D${D}S${S};" ;
    }
    $notConservedSpLine =~ s/;$// ;
    if ($notConservedSpLine ne '#'){
        $notConservedSpLine =~ s/^#//;
    }else{
        $notConservedSpLine = '.';
    }
    last if $light == 0 ;
    print O "$changedSpLine\t$notConservedSpLine\n";
}
close O ;

