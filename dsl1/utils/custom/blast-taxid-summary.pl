#!/usr/bin/env perl

use 5.016;
use strict;
use warnings;
use Bio::DB::Taxonomy;
use Data::Dumper;
use Getopt::Long;
use DBI;

# my quick and dirty script to download the latest livelist and gi2taxid data
# for local database lookups.  Can, and should, be cleaned up!

my ($file, $taxdb, $acc2taxid);
my $diff = 3; # 3% difference allowed for top hits

GetOptions(
    'file=s'         => \$file,
    'tax=s'          => \$taxdb,
    'acc2taxid=s'    => \$acc2taxid,
    'diff=i'         => \$diff
    # TODO: column for accession 
);

my @valid_ranks = qw(superkingdom kingdom phylum class order family genus species);

my $acc_dbh = DBI->connect("DBI:SQLite:dbname=$acc2taxid" ,"","", {RaiseError =>1})
    or die "Couldn't connect to database: " . DBI->errstr;

my $tax_dbh = Bio::DB::Taxonomy->new(-source    => 'sqlite', 
                                     -db        => $taxdb);

open(my $blastin, '<', $file) or die "Can't open file:$!";

my $ct = 0;

say join("\t", qw(QueryID QueryLen SubjectID SubjectLen SubjectTitle PercIdent 
    AlnLength Mismatch GapOpen QueryStart QueryEnd SubjectStart SubjectEnd QueryCov 
    EValue BitScore NCBI_TaxID SciName CommonName), @valid_ranks);

my @best;

my ($last_asv, $best_score);

HSP:
while (my $line = <$blastin>) {
    chomp $line;
    my @atts = split("\t", $line);
    
    my ($asv, $score) = @atts[0,5];
    if (!$last_asv or $last_asv ne $asv) {
        ($last_asv, $best_score) = ($asv, $score);
    } else {
        if ( ($best_score - $score) > $diff ) {
            next HSP;
        }
        if ($score > $best_score) {
            $best_score = $score # this does sometimes happen (we're not using bitscore)
        }
    }
    
    push @best, \@atts;
}

my %cache;

foreach my $atts (@best) {
    # strip out any ref-like stuff, just need the identifier
    (my $accver = $atts->[2]) =~ s/^[^\|]+\|([^\|]+)\|$/$1/;
    
    # need only accession
    my ($acc, $ver) = split(/\./, $accver);
    my $taxid = retrieve_taxid($acc, $acc_dbh);
        
    my $taxoninfo = exists $cache{$taxid} ? $cache{$taxid} : retrieve_taxinfo($taxid, $tax_dbh);
    
    if ($taxoninfo) {
        my ($taxon, $info) = @{$taxoninfo};
        say join("\t", @{$atts}, $taxid, $taxon->scientific_name(), join(',',$taxon->common_name()), @{$info}{@valid_ranks});
        $cache{$taxid} //= $taxoninfo;
    } else {
        warn "Can't find taxon for $taxid\n";
    }
}

sub retrieve_taxid {
    my ($id, $acc_dbh) = @_;
    my $sth = $acc_dbh->prepare_cached(<<SQL);
    SELECT taxon_id
    FROM acc2taxid
    WHERE
        accession = ?
SQL
    my ($taxid);
    $sth->execute($id) or die($sth->errstr);
    return $sth->fetchall_arrayref()->[0][0];
}

sub retrieve_taxinfo {
    my ($taxid, $tax_dbh) = @_;
    my $origtaxon = $tax_dbh->get_taxon(-taxonid => $taxid);
    unless ($origtaxon) {
        warn "No taxon returned for $taxid";
        return;
    }
    my %taxonranks = map {$_ => ''} @valid_ranks;
    if (exists $taxonranks{$origtaxon->rank}) {
        $taxonranks{$origtaxon->rank} = $origtaxon->node_name;
    }
    my $taxon = $origtaxon;
    while (my $anc = $tax_dbh->ancestor($taxon)) {
        if (exists $taxonranks{$anc->rank}) {
            $taxonranks{$anc->rank} = $anc->node_name;
        }
        $taxon = $anc;
    }
    [$origtaxon, \%taxonranks];
}

