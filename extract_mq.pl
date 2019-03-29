#!/usr/bin/perl -w
use strict;

my $sf = shift || 0;  # Sequence field
my $ef = 2;  # Read field
my $MCOV = 1; # Minimum coverage
my $LCOV = 10; # At or below LCOV, any change is significant
my $LIMIT = 0.1; # Above MINCOV, a change is significant with a proportion higher
                 # than LIMIT
while (<>){
  my @f = split(/\s+/);
  $sf = set_seq_field(@f) unless ($sf);
  die("Could not determine seq field\n") unless ($sf);
  my $exp = $f[$ef];
  my $cov = $f[$sf-1];
  my $seq = $f[$sf];
  if ($exp ne 'N' && $exp ne 'n'){
    if ($exp eq '*'){
      print;
    }
    elsif ($cov > $MCOV && $seq =~ /[acgt]/io){
      if ($cov <= $LCOV){
        print;
      }
      else{
        my $n = $seq =~ tr/AaCcGgTt/AaCcGgTt/;
        print if ($n/$cov > $LIMIT);
      }
    }
  }
}

sub clean_seq{
  my $reads = shift;
  $reads =~ s/\^.//g;  # A read is born
  $reads =~ s/\$//g;   # A read is dead
  while ($reads =~ /[\+\-](\d+)/){ # Not interested in indels
    my $pattern = '[\+\-]'."$1\.\{$1\}";
    $reads =~ s/$pattern//;
  }
  return $reads;
}

sub set_seq_field{
  my $cov = shift;
  my $seq   = shift;
  my $result = 0;
  my $done = 0;
  while (!$done && scalar @_){
    my $qual = shift;
    my $cseq = clean_seq($seq);
    $result++;
    if (is_integer($cov) && $result>2){
      if (length $cseq == $cov && length $qual == $cov){
        if (!($cseq =~ /\d/)){
          $done = 1;
        }
      }
    }
    $cov = $seq;
    $seq = $qual;
  }
  unless ($done){
    tolog("Could not detect sequence field. Defaulting to fifth column");
    $result = 4;
  }
  return $result;
}

sub tolog{
  my $t = shift;
  print "$t\n";
}


sub is_integer{
  my $n = shift;
  my $result = 1;
  $result = 0 if ($n =~ /\D/);
  return $result;
}