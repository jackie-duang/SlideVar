use warnings;
use strict;

if (@ARGV < 4){
    die "perl $0 <input.maf> <species.list> <min_length> <output_dir_prefix>\n";
}

my ($in , $speciesFile, $min_length, $outdir) = @ARGV;

my %sortSp ;
my $count = 0 ;
open I , "< $speciesFile";
while(<I>){
    chomp;
    my @a = split /\s+/;
    next if @a==0 ;
    $sortSp{$a[0]} = $count ;
    $count ++ ;
}
close I;

my %seq ;
my %length ;
my $block_num = 0 ;
my $length = 0 ;
$/ = "\na";
open I , "< $in";
while(<I>){
    chomp;
    my @lines = split /\n+/;
    next if @lines==0 ;
    my $light = 0 ;
    my $length = 0 ;
    my $ref_length = 0 ;
    my $pos = '';
    foreach my $line (@lines){
        my @a = split /\s+/,$line ;
        next if @a==0 ;
        next unless $a[0] eq 's';
        my ($sp, $start, $len, $strand, $srcLen, $seq) = @a[1,2,3,4,5,6];
        $sp =~ /^(\w+)\.(\S+)/;
        $sp = $1 ;
        my $chr = $2 ;
        my $newseq = $seq ;
        $newseq =~ s/-//g;
        if ($light == 0 and $length == 0){
            $length = length $newseq ;
            $ref_length = length $seq ;
            my $end = $start + $len ;
            $start = $start + 1 ;
            $pos = "$chr\t$start\t$end\t$sp";
        }
        last if ($length < $min_length);
        my $seq_len = length $seq ;
        if ($seq_len == $ref_length){
            $seq{$block_num}{$sp} = $seq ;
        }else{
            print STDERR "Warning: $line\n";
        }
        $light ++ ;
    }
    if ($light > 0){
        $length{$block_num}{len} = $ref_length ;
        $length{$block_num}{pos} = $pos ;
        $block_num ++ 
    }
}
close I;
$/ = "\n";

print STDERR "\033[31mTotal $block_num blocks\033[0m\n";

foreach my $b_n (sort {$a <=> $b} keys %seq){
    my $new_b_n = sprintf "%06d",$b_n;
    $new_b_n =~ /^(\d{3})(\d{3})$/;
    my $dir = "$outdir/$1/$2";
    `mkdir -p $dir` if (!-d $dir);
    open F , "> $dir/Block$new_b_n.info";
    print F "$length{$b_n}{pos}\tBlock$new_b_n\n";
    close F ;
    open O , "> $dir/Block$new_b_n.fasta";
    foreach my $sp (sort {$sortSp{$a} <=> $sortSp{$b}} keys %sortSp){
        if (exists $seq{$b_n}{$sp}){
            print O ">$sp\n$seq{$b_n}{$sp}\n";
        }
        else{
            my $length = $length{$b_n}{len} ;
            my $seq = '' ;
            $seq .= '-' x $length;
            print O ">$sp\n$seq\n";
        }
    }
    close O;
}