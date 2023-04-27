#!perl
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'Stefanl_Libs::fastaDB' ) || print "Bail out!\n";
}

diag( "Testing Stefanl_Libs::fastDB $Stefanl_Libs::fastDB::VERSION, Perl $], $^X" );


