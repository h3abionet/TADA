#!/usr/bin/env perl
use 5.016;
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Cwd;

# set up and load the latest NCBI taxonomy

my $db = Bio::DB::Taxonomy->new(-verbose    => 1, -force => 1,
                                -source    => 'sqlite' ,
                                -directory => Cwd::cwd,
                                -nodesfile => 'nodes.dmp',
                                -namesfile => 'names.dmp');

