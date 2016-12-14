#!/usr/bin/perl -w
# Auther: Zechen Chong
$| = 1;
use strict;
use warnings;

my $kmerfile = shift or die "Usage: $0 <kmer.stat> <germline_novo_kmer_read.fq> \n";

my $rdfile = shift or die "Usage: $0 <kmer.stat> <germline_novo_kmer_read.fq> \n";

my ($start, $end);

my $id = 0;
my %kmer2id = ();
my %id2seq = ();
$start = time;
printf STDERR "begin kmer2id ...\n";
open IN, $kmerfile or die $!;
my $kmer;
while (<IN>) {
	next unless /SOMATIC/;
	my @e = split /\s+/,$_;
	next if defined $kmer2id{$e[0]};
	$kmer2id{$e[0]} = $id;
	$id++;
	$kmer = $e[0];
}
my $len = length $kmer;
my $kmer_num = $id;
$end = time;
printf STDERR "kmer2id takes %d seconds\n", $end - $start;
close IN;

my @ids = (); 
my @sz = ();
for (my $i = 0; $i < $id; $i++) {
	push @ids, $i;
	push @sz, 1;
}
printf STDERR "begin id2seq...\n";
$start = time;
open IN, $rdfile or die $!;

my @relations = ();
my $seq;
while (<IN>) {
	if ($. % 4 == 2) {
		chomp;
		$seq = $_;
		{
			$id2seq{$id} = $seq;
			push @ids, $id;
			push @sz, 1;
			for (my $i = 0; $i <= length($seq) - $len; $i++) {
				my $kmer = substr($seq, $i, $len);
				if (exists $kmer2id{$kmer}) {
					#print STDERR ($kmer, " ", $kmer2id{$kmer}, " ", $id, "\n");
					&union($kmer2id{$kmer}, $id);
				} elsif (exists $kmer2id{&rev_com($kmer)}) {
					#print STDERR ($kmer, " ", $kmer2id{&rev_com($kmer)}, " ", $id, "\n");
					&union($kmer2id{&rev_com($kmer)}, $id);
				}
			}
			$id++;
		}
	}
}

close IN;
$end = time;
printf STDERR "id2seq takes %d seconds\n", $end - $start;

#printf STDERR "begin clustering ...\n";

#printf STDERR "clustering takes %d seconds\n", $end - $start;

printf STDERR "begin sorting ids...\n";
$start = time;
my $line = 0;
my @idmap = ();
foreach my $i (@ids) {
	push @idmap, { newid => $i, oldid => $line };
	$line ++;
}

@idmap = sort {$a->{newid} <=> $b->{newid}} @idmap;
$end = time;
printf STDERR "sorting ids takes %d seconds\n", $end - $start;
printf STDERR "begin output results...\n";

$start = time;
my $group = -1;
my $pre = -1;

my @holder = ();

for my $i (@idmap) {
	if ($i->{newid} != $pre) {
		if (scalar @holder >= 3) {
			$group ++;
			&print_seq($group, \@holder);
		}
		@holder = ();
		$pre = $i->{newid};
		if ($i->{oldid} >= $kmer_num) {
			my $seq = $id2seq{$i->{oldid}};
			push @holder, $seq;
		} 
	} else {
		if ($i->{oldid} >= $kmer_num) {
			my $seq = $id2seq{$i->{oldid}};
			push @holder, $seq;
		} 
	}
}
		if (scalar @holder >= 3) {
			$group ++;
			&print_seq($group, \@holder);
		}
$end = time;
printf STDERR "Outputting results takes %d seconds\n", $end - $start;
printf STDERR "Finished\n";

1;

sub print_seq {
	my $gid = shift;
	my $group = shift;
	
	foreach my $seq (@$group) {
		print join("\t", ($gid, $seq)), "\n";
	}
}

sub rev_com {
	my $seq = shift;

	$seq =~ tr/ACGT/TGCA/;
	return reverse $seq;
}
sub union {
	my $p = shift;
	my $q = shift;
	my ($i, $j);
	for ($i = $p; $i != $ids[$i]; $i = $ids[$i]) {
		$ids[$i] = $ids[$ids[$i]];
	}
	for ($j = $q; $j != $ids[$j]; $j = $ids[$j]) {
		$ids[$j] = $ids[$ids[$j]];
	}

	return if ($i == $j);

	if ($sz[$j] + $sz[$i] < 1000) {
		if ($sz[$i] < $sz[$j]) {
			$ids[$i] = $j;
			$sz[$j] += $sz[$i];
		} else {
			$ids[$j] = $i;
			$sz[$i] += $sz[$j];
		}
	}
}
