#!/usr/bin/perl -w
#count_venn.vcf

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'fasta|f=s','help|h','prefix|p=s');
my @vcffiles = @ARGV;

unless($opt{fasta}) {
    $opt{fasta} = 'hs38DH.fa';
}
unless($opt{prefix}) {
    $opt{prefix} = 'merge';
}
my $total = 0;
open FASTA, "<$opt{fasta}\.fai" or die "unable to find fai index for $opt{fasta}, please index with samtools";
while (my $line = <FASTA>) {
    chomp($line);
    my ($chr,$length,$offset,$linelen,$linebytes) = split(/\t/,$line);
    $total += $length;
}

foreach $vcf (@vcffiles) {
  $shuff = $vcf;
  $shuff =~ s/vcf+/shuff.vcf/;
  my $gfile = $shuff;
  if ($vcf ne $vcffiles[0] && ! -e $shuff) {
    if ($shuff =~m/gz$/) {
      system(qq{vcf-shuffle-cols -t $vcffiles[0] $vcf |bgzip > $shuff})
    }else {
      system(qq{vcf-shuffle-cols -t $vcffiles[0] $vcf > $shuff})
    }
  }elsif (! -e $shuff) {
     system(qq{ln -s $vcf $shuff});
  }
  if ($shuff =~m/gz$/) {
    system("tabix $shuff") if ($vcf =~m/gz$/);
  }else {
    system("bgzip $shuff") unless (-e "$vcf\.gz") ;
    system("tabix $shuff\.gz") unless ($vcf =~m/gz$/);
    $gfile = "$shuff\.gz";
  }
  push @zipvcf, $gfile;
  $fname = (split(/\//,$vcf))[-1];
  $fname =~ s/\.vcf//;
}
my $command = "vcf-compare ".join(" ",@zipvcf)." > $opt{prefix}\.vcf_compare.out";
system($command);
$command = "vcf-isec -f --prefix $opt{prefix}\.integ ".join(" ",@zipvcf);
system($command);

my %ct = ();
my %snpbins = ();
open CTS, ">$opt{prefix}\.cts" or die $!;
print CTS "Estimated sum: ".$total,"\n";
open VC, "<$opt{prefix}\.vcf_compare.out" or die $!;
while (my $line = <VC>) {
	my @key = ();
    next if($line =~ m/#/);
	if ($line =~ m/^VN\s+(.+)/) {
		my ($ct,@flist) = split(/\s+/,$1);
		foreach my $dset (@zipvcf) {
			if (grep(/$dset/,@flist)) {
				push @key, 1;
			}else {
		push @key, 0;
			}
		}
	$key = join("",@key);
	$ct{$key} = $ct;
	$total -= $ct;
    }
}
my @g = bits(scalar(@vcffiles));
$ct{$g[0]} = $total;
foreach (@g) {
    $ct{$_} = 0 unless ($ct{$_});
    print CTS join("\t",$_,$ct{$_}),"\n";
}

system("Rscript $lca -c $opt{prefix}.cts -s $opt{prefix}.stats");

my @key1 = split(//,$g[-1]);
my @key2;
foreach $i (0..$#key1) {
    push @key2, $i if ($key1[$i] > 0);
}
open STATS, "<$opt{prefix}.stats" or die $!;
my @keeppos;
while (my $line = <STATS>) {
    chomp($line);
    next unless ($line =~ m/postprobs\[(\d+)\]\s+(\S+)/);
    my ($index,$prob) = ($1,$2);
    next unless ($prob >= 0.8);
    my $key = $g[$index-1];
    @key1 = split(//,$key);
    @key2 = ();
    foreach $i (0..$#key1) {
	push @key2, $i if ($key1[$i] > 0);
    }
    my $subset =  $opt{prefix}.'.integ'.join('_',@key2).'.vcf.gz';
    push @keeppos, $subset if (-e $subset);
}

system("vcf-concat ".join(" ",@keeppos)." |vcf-sort |bgzip >  $opt{prefix}.baysic.vcf.gz");


sub bits { glob join "", map "{0,1}", 1..$_[0] }
