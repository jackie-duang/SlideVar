use warnings;
use strict;

if (@ARGV != 3) {
    die "perl $0 <output_dir_from_maf2fasta> <species.list> <output_file>\n";
}
my ($in, $list, $out) = @ARGV ;

my @sps = ();
open I , "< $list";
while (<I>){
    chomp;
    my @a = split /\s+/, $_ ;
    next if @a==0 ;
    next unless (@a>1 and ($a[1] eq '*' or $a[1] eq '#'));
    push @sps, $a[0] ;
}
close I ;
my %sps = map {$_ => 1} @sps ;

my @files = <"$in/*/*/*info.list">;

open O , "> $out";
print O "Block_id\tBackground_conserved_sites";
foreach my $sp (@sps){
    print O "\t$sp";
}
print O "\n";
foreach my $file (@files){
    open I , "< $file";
    <I>;
    my $all = 0 ;
    my %hash ;
    my %info ; 
    while (<I>){
        chomp;
        my @a = split ;
        next if @a==0 ;
        next if $a[3] =~ /^\d+$/;
        $all ++ ;
        if ($a[2] ne '.'){
            my @species = split /;/,$a[2] ;
            foreach my $species (@species){
                $species =~ /^(\w+):(\w+)$/;
                $hash{$1}{$a[1]} = $2 ;
                $info{$a[1]}{$1} = $2 ;
            }
        }
    }
    close I ;
    if (keys %hash == 0) {
        next ;
    }
    $file =~ /(Block\d+)\.fasta/;
    my $id = $1 ;
    # print O "$id\t$all";
    # foreach my $sp (@sps){
    #     if (exists $hash{$sp}){
    #         my $n = keys %{$hash{$sp}};
    #         print O "\t$n";
    #     }
    #     else{
    #         print O "\t0";
    #     }
    # }
    # print O "\n";

    foreach my $site (sort {$a <=> $b} keys %info){
        print O "$id\t$site";
        foreach my $sp (@sps){
            if (exists $info{$site}{$sp}){
                print O "\t$info{$site}{$sp}";
            }
            else{
                print O "\t-";
            }
        }
    print O "\n";
    }
}
close O ;