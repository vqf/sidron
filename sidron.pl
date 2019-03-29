#!/usr/bin/perl -w
use strict;
#use vBAM::GetPileupRegion;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Storable;

my $VERBOSE = 1;
my $pobj = [];   # Contains file descriptors
while (scalar @ARGV >= 2){
  my $file = shift;
  my $tabl = shift;
  open (my $fh, $file) or die ("Cannot open $file: $!\n");
  my $bt = -B $file;
  if ($bt){
    tell_user("$file seems to be gzipped", 1);
    $fh = new IO::Uncompress::Gunzip $fh;
  }
  my $hsh = Storable::retrieve($tabl) or die("Could not open frequencies file: $!\n");
  my $fprobobj = Prob->new($hsh);  # Prob object for forward strand
  my $rprobobj = Prob->new($hsh);  # Prob object for reverse strand
  my $tmp = {
    'fh' => $fh,
    'fprob' => $fprobobj,
    'rprob' => $rprobobj,
    'chr'   => '',
    'pos'   => 0
  };
  push @$pobj, $tmp;
}
#my $VERBOSE = shift || 0;
my @QUALITY_GROUPS = (10, 20, 30, 40);
my @MAPQ_GROUPS    = (25, 37, 60);
my $MAXCOV = 1e4;
my $ALIGQUAL = 0;   # 1 to take alignment quality into account
my $WRITE_EVERY = 1e3;
my $I = 0;
my @temp = split /[\/\\]/, $0;
my $pname = pop @temp;
#die("Use: perl $pname tumor_pileup_file tumor_freq_file".
#    " [normal_pileup_file] [normal_freq_file] [verbosity]\n")
#  unless ($tfile && -e $tfile &&
#          $ttfile && -e $ttfile);

my %rev = (
  'A' => 'T',
  'C' => 'G',
  'G' => 'C',
  'T' => 'A',
  'N' => 'N'
);


my $chr_field = 0;
my $pos_field = 1;
my $exp_field = 2;
my $nseq_field = 0;

my $counter = 0;
while (get_next_line($pobj)){
  #my $nline = <$fhn>;
  #chomp $nline;
  $counter++;
  # Init objects
  my $clearly = 0;   # Store clearly hz separately (TODO or not TODO)
  my $nin = 0;       # Not in normal
  my $t_alleles = '';
  my $n_alleles = '';
  my @pileups = ();   # Pileup lines
  my @extra   = ();   # Sidron scores
  foreach my $f (@$pobj){
    # Set models
    my $cseq = $f->{'cseq'};
    my ($tref_base, $tsec_base) = get_ref($cseq);
    unless($tsec_base){
      $tsec_base = $rev{$tref_base};
    }
    $t_alleles = "$tref_base$tsec_base";
    my $hz_model_direct  = [[1,[$tref_base]]];
    my $het_model_direct = [[1,[$tref_base]],[1,[$tsec_base]]];
    my $tref_base_rev = $rev{uc($tref_base)};
    my $tsec_base_rev = $rev{uc($tsec_base)};
    my $hz_model_rev  = [[1,[$tref_base_rev]]];
    my $het_model_rev = [[1,[$tref_base_rev]],[1,[$tsec_base_rev]]];
    if ($tref_base){
      my $tprob = get_conf_prob($f,
                                $hz_model_direct, $hz_model_rev,
                                $het_model_direct, $het_model_rev);
      push @pileups, $f->{'line'};
      push @extra, $t_alleles, $tprob;
    }
  }
  my $outline = join("\t", @pileups)."\t".join("\t", @extra)."\t-";
  tell_user('.', 2);
  print "$outline\n";
}
#Close open file handles
foreach my $f (@$pobj){
  close $f->{'fh'};
}

sub get_next_line{
  my $fobj = shift;   # File descriptors: [file_handle, chr, pos, current_line, cache_pos]
  return 0 unless (ref($fobj) eq 'ARRAY' && scalar @$fobj);
  my $mfd = shift @$fobj;
  my $mfh = $mfd->{'fh'};
  my $oldchr = $mfd->{'chr'} || '';
  my $oldpos = $mfd->{'pos'} || 0;
  my $success = 0;
  while (!$success && ($mfd->{'line'} = <$mfh>)){
    chomp $mfd->{'line'};
    $success = 1;
    my ($newchr, $newpos) = split(/\s+/, $mfd->{'line'}, 3);
    $mfd->{'chr'} = $newchr;
    $mfd->{'pos'} = $newpos;
    foreach my $fd (@$fobj){
      my $fh = $fd->{'fh'};
      if ($newchr ne $oldchr){ # First file gets to new chr
        while ($fd->{'chr'} ne $newchr && ($fd->{'line'} = <$fh>)){
          chomp $fd->{'line'};
          ($fd->{'chr'}, $fd->{'pos'}) = split(/\s+/, $fd->{'line'});
        }
      }
      while ($fd->{'chr'} eq $newchr && $fd->{'pos'} < $newpos && ($fd->{'line'} = <$fh>)){
        ($fd->{'chr'}, $fd->{'pos'}) = split(/\s+/, $fd->{'line'});
      }
      unless ($fd->{'chr'} eq $newchr && $fd->{'pos'} == $newpos){
        $success = 0;
      }
    }
    $oldchr = $newchr if ($oldchr ne $newchr);
  }
  unshift @$fobj, $mfd;
  if ($success){
    # Fill useful fields
    foreach my $f (@$fobj){
      chomp $f->{'line'};
      my $line = $f->{'line'};
      my @f = split(/\s+/, $line);
      if (!exists($f->{'seq_field'})){
        $f->{'seq_field'} = set_seq_field(@f);
      }
      my $sf = $f->{'seq_field'};
      my $qf = $sf + 1;
      my $seq = $f[$sf];
      my $exp  = uc($f[$exp_field]);
      my $rexp = lc($exp);
      my $cseq = clean_seq($seq);
      $cseq =~ s/\./$exp/g;
      $cseq =~ s/\,/$rexp/g;
      $f->{'cseq'} = $cseq;
      $f->{'qual'} = $f[$qf];
    }
    return 1;
  }
  elsif(eof($mfh)){
    return 0;  # No more lines in $mfh
  }
  else{
    warn("Possible error: no success, but first file has not ended");
  }
}

sub div{
  my $num = shift;
  my $den = shift;
  my $result = 0;
  return 1;
  if ($num =~ /[\-\+]Inf/ && $den =~ /[\-\+]Inf/){
    $result = 'NaN';
  }
  elsif($num =~ /[\-\+]Inf/){
    $result = $num;
  }
  elsif ($den =~ /[\-\+]Inf/){
    $result = 0;
  }
  elsif ($num eq 'NaN' || $den eq 'NaN'){
    $result = 'NaN';
  }
  elsif ($den == 0){
    $result = '+Inf';
    $result = '-Inf' if (abs($num) != $num);
  }
  else{
    $result = $num - $den;
  }
  return $result;
}

sub get_conf_prob{
  my $fobj = shift;
  my $cseq = $fobj->{'cseq'};
  my $qual = $fobj->{'qual'};
  my $mapq = $fobj->{'mapq'} || '';
  my $f  = $fobj->{'fprob'};
  my $r  = $fobj->{'rprob'};
  my $hz_model_direct = shift;
  my $hz_model_rev = shift;
  my $het_model_direct = shift;
  my $het_model_rev = shift;
  #Init models
  my $tmodel = ['A', 'C', 'G', 'T'];
  $f->init($tmodel);
  $r->init($tmodel);
  #
  my @tbases = split(//, $cseq);
  my @tquals = split(//, $qual);
  my $n1 = scalar @tbases;
  my $n2 = scalar @tquals;
  my @mquals = split(//, $mapq) if ($ALIGQUAL);
  while (my $tbase = shift @tbases){
    my $q = shift @tquals;
    die($qual) unless (defined($q));
    #die ("$n1\t$n2\t$tbase\n") unless $q;
    my $m = shift @mquals if ($ALIGQUAL);
    #return 1 unless defined($m);
    my $bqual = ord($q) - 33;
    my $mqual = ord($m) - 33 if ($ALIGQUAL);
    my $qqual = set_groups($bqual, @QUALITY_GROUPS);
    my $tmqual = set_groups($mqual, @MAPQ_GROUPS) if ($ALIGQUAL);
    if (exists($rev{uc($tbase)})){
      if ($tbase eq uc($tbase)){  # Forward strand
        my $qstate = [[$tbase, 'qual', $qqual]];
        push @$qstate, [$tbase, 'mapq', $tmqual] if ($ALIGQUAL);
        $f->add_state($qstate);
      }
      else{ # Reverse strand
        my $qstate = [[$rev{uc($tbase)}, 'qual', $qqual]];
        push @$qstate, [$rev{uc($tbase)}, 'mapq', $tmqual] if ($ALIGQUAL);
        $r->add_state($qstate);
      }
    }
    else{
      tell_user("Line $fobj->{'line'}: found non standard base $tbase. Skipping.", 2);
    }
  }
  #my $prob_hz = $f->get_prob($hz_model_direct) +
  #              $r->get_prob($hz_model_rev);
  #my $prob_het = $f->get_prob($het_model_direct) +
  #               $r->get_prob($het_model_rev);
  #die($f->debug());
  my $prob_hz = $f->get_prob($hz_model_direct) +
                $r->get_prob($hz_model_rev);
  my $prob_het = $f->get_prob($het_model_direct) +
                 $r->get_prob($het_model_rev);
                 
  my $c = 0;#$f->get_c() + $r->get_c();
  #print $f->debug()."\n---\n".$r->debug();
  #exit(0);
  
  #print "--- $cseq\t$prob_hz\t$prob_het\n";
  #my $het_hz = $prob_het - $prob_hz;
  my $prob_hz_c = $prob_hz - $c;
  my $prob_het_c = $prob_het - $c;
  my $het_hz = $prob_het_c - $prob_hz_c;
  return $het_hz;
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
  open (APP, ">>sidrlog.txt");
  print APP "$t\n";
  close APP;
}

sub is_integer{
  my $n = shift;
  my $result = 1;
  $result = 0 if ($n =~ /\D/);
  return $result;
}

sub set_groups{
  my $val = shift;
  my @groups = sort{$a <=> $b} @_;
  my $result = 0;
  my $min = shift @groups;
  while ($val > $min && scalar @groups){
    $min = shift @groups;
    $result++;
  }
  return $result;
}


sub get_ref{
  my $cseq = shift;
  my $half = (length $cseq) >> 1;
  my @result = ();
  my $cutoff = 0;
  my @store = ();
  my $na = $cseq =~ tr/Aa/Aa/;
  push @store, ['A', $na];
  my $nc = $cseq =~ tr/Cc/Cc/;
  push @store, ['C', $nc];
  my $ng = $cseq =~ tr/Gg/Gg/;
  push @store, ['G', $ng];
  my $nt = $cseq =~ tr/Tt/Tt/;
  push @store, ['T', $nt];
  my @sorted = sort{$a->[1] <=> $b->[1]} @store;
  my $best = pop @sorted;
  my $second = pop @sorted;
  $second->[0] = '' unless ($second->[1]);
  @result = ($best->[0], $second->[0]);
  return @result;
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


sub tell_user{
  my $text = shift;
  my $v = shift;
  if ($VERBOSE >= $v){
    open (APN, ">>warnings.log");
    print APN "$text\n";
    close APN;
  }
  return 0;
}

package Prob;

  sub new{
    my $pkg = shift;
    my $table = shift;
    my $self = {
      'table' => $table,
      'configuration' => [],
      'setup' => {
        'level_sep' => '|',    # Reserved chars
        'element_sep' => '#'
      }
    };
    bless ($self, $pkg);
    # Init
    if (!exists($self->{'table'}{'processed'})){
      $self->get_frequencies($self->{'table'});
      $self->{'table'}{'processed'} = 1;
    }
    #die($self->debug());
    return $self;
  }
  
  sub init{
    my $self = shift;
    my $tmodel = shift;
    $self->{'configuration'} = [];
    $self->{'nrstates'} = [];
    if ($tmodel && ref($tmodel) eq 'ARRAY' && !exists($self->{'total_model'})){
      my $total_model = [];
      foreach my $base (@$tmodel){
        my $lfreq = $self->{'table'}{$base}{'total'};
        my $freq = 10**$lfreq;
        push @$total_model, [$freq, [$base]];
      }
      $self->{'total_model'}{'model'} = $total_model;
      $self->set_model($total_model);
      $self->{'total_model'}{'string'} = $self->{'current_model'};
    }
    return 0;
  }
  
  sub set_model{
    my $self = shift;
    my $model = shift;
    die("You need to provide a model to get a probability")
      unless ($model);
    # Normalize frequencies
    my $total = 0;
    foreach my $comp (@$model){
      my $p = $comp->[0];
      $total += $p;
    }
    if ($total != 1){
      foreach my $comp (@$model){
        $comp->[0] /= $total;
      }
    }
    # Check if cached
    my $string = $self->stringify($model);
    if (exists($self->{'cached_models'}{$string})){
      $self->{'current_model'} = $string;
    }
    else{
      # Check model for consistency and add to cache
      my $pointer = $self->{'table'};
      foreach my $comp (@$model){
        my $p = $comp->[0];
        my $c = $comp->[1];
        unless ($self->get_pointer($c)){
          die("Incorrect model: $string\n");
        }
      }
      $self->{'cached_models'}{$string} = $model;
      $self->{'current_model'} = $string;
    }
    return 0;
  }
  
  sub get_pointer{
    my $self = shift;
    my $path = shift;
    my $pointer = $self->{'table'};
    foreach my $key (@$path){
      if (ref($pointer) eq 'ARRAY'){
        if (exists($pointer->[$key])){
          $pointer = $pointer->[$key];
        }
        else{
          return 0;
        }
      }
      elsif (ref($pointer) eq 'HASH'){
        if (exists($pointer->{$key})){
          $pointer = $pointer->{$key};
        }
        else{
          return 0;
        }
      } 
    }
    return $pointer;
  }
  
  sub stringify{
    my $self = shift;
    my $struct = shift;
    my $result = '';
    foreach my $el (@$struct){
      if (ref($el) eq 'ARRAY'){
        $result .= $self->{'setup'}{'level_sep'}.$self->stringify($el).
                    $self->{'setup'}{'level_sep'};
      }
      else{
        $result .= $self->{'setup'}{'element_sep'}.
                    $el.$self->{'setup'}{'element_sep'};
      }
    }
    $result =~ s/\./\_/g;
    return $result;
  }
  
  sub add_state{
    my $self = shift;
    my $states = shift;  # Array of arrays. Each sub-array contains
                         # the path to a state
    push @{$self->{'configuration'}}, [$states, $self->stringify($states)];
    $self->{'nrstates'} = [] if (exists($self->{'nrstates'})
                                 && scalar @{$self->{'nrstates'}});
    return 0;
  }
  
  sub get_prob{
    # Gets log10 probability of current config given a model.
    my $self = shift;
    my $model = shift;
    #my $log   = shift;  # If true, prob is logarithmic
    $self->set_model($model);
    unless (! scalar @{$self->{'configuration'}} ||
            (exists($self->{'nrstates'}) && scalar @{$self->{'nrstates'}})
            ){
      my $nrstates = [];
      my @multinomial = ();
      my @states = sort {$a->[1] cmp $b->[1]} @{$self->{'configuration'}};
      my $ostate = shift @states;
      unshift @$nrstates, [1, $ostate];
      while (my $cstate = shift @states){
        if ($cstate->[1] eq $ostate->[1]){
          $nrstates->[0][0]++;
        }
        else{
          unshift @$nrstates, [1, $cstate];
          $ostate = $cstate;
        }
      }
      foreach my $el (@$nrstates){
        push @multinomial, $el->[0];
      }
      $self->{'nrstates'} = $nrstates;
      $self->{'multinomial'} = $self->multinomial(@multinomial);
    }
    #print $self->debug()."\n--\n";
    my $result = $self->calculate();
    #$calc .= "\n";
    #$pths .= "\n";
    #print $pths;
    #print $calc;
    #print "$result\n";
    return $result;
  }
  
  sub calculate{
    my $self = shift;
    my $nrstates = $self->{'nrstates'};
    my $factor = $self->{'multinomial'} || 0;
    my $result = $factor;
    foreach my $state (@$nrstates){
      my $n          = $state->[0];
      my $local_paths = $state->[1][0];
      my $current_model = $self->{'cached_models'}{$self->{'current_model'}};
      my $sp = 0;
      foreach my $m (@$current_model){
        my $pmodel = log10($m->[0]);  # Prob of the model
        my @path = @{$m->[1]};  # Model-independent path
        my $lp = 0;             # Conf probability
        foreach my $lpath (@$local_paths){
          my @tpath = ();
          push @tpath, @path, @$lpath;
          my $ptr = $self->get_pointer(\@tpath);
          if (ref($ptr) ne 'HASH'){
            #tell_user("Bad path: " . join(', ', @tpath)."\n", 1);
          }
          else{
            $lp += $ptr->{'total'};
          }
        }
        if ($sp){  #More than one factor
          # D'oh!!!!
          my $tsp = 10**$sp;
          my $tsm = 10**($pmodel + $lp);
          $sp = log10($tsp + $tsm);
        }
        else{  # Waste less time for the first factor
          $sp = $pmodel + $lp;
        }
      }
      $result += $sp * $n;
    }
    return $result;
  }
  
  sub get_c{
    # Gets total probability in the current configuration.
    my $self = shift;
    $self->{'current_model'} = $self->{'total_model'}{'string'};
    my $result = $self->calculate();
    return $result;
  }
  
  sub get_prob_c{
    # Gets log10 probability of current config given a model divided by the total
    # probability of the config.
    my $self = shift;
    my $model = shift;
    my $tmodel = shift;  # Total model, no frequencies given
    #my $log   = shift;  # If true, prob is logarithmic
    # Obtain frequencies from table to complete total model
    unless (exists($self->{'total_model'})){
      my $total_model = [];
      foreach my $base (@$tmodel){
        my $lfreq = $self->{'table'}{$base}{'total'};
        my $freq = 10**$lfreq;
        push @$total_model, [$freq, [$base]];
      }
      $self->{'total_model'}{'model'} = $total_model;
      $self->set_model($total_model);
      $self->{'total_model'}{'string'} = $self->{'current_model'};
    }
    my $result = $self->get_prob($model);
    $self->{'current_model'} = $self->{'total_model'}{'string'};
    my $t = $self->calculate();
    #my $t = $self->get_prob($self->{'total_model'}{'model'});
    $result -= $t;
    return $result;
  }
  
  
  sub get_frequencies{
    my $self = shift;
    my $hash = shift;
    my $total = log10($self->{'table'}{'total'});
    my @sorted = sort {
      if ($a eq 'total') { return 1; }
      elsif ($b eq 'total') { return -1; }
      else { return $a cmp $b; } } keys %$hash;
    foreach my $key (@sorted){
      if (exists($hash->{'total'}) && $key ne 'total'){
        $self->get_frequencies(\%{$hash->{$key}}, $hash->{'total'});
      }
      elsif ($total){
        $hash->{$key} = log10($hash->{$key}) - $total;
      }
    }
    return 0;
  }
  
  sub log10{
    my $n = shift;
    return 'NaN' unless ($n > 0);
    my $result = log($n)/log(10);
    return $result;
  }
  
  sub index_table{
    my $self = shift;
    my $hash = shift;
    foreach my $key (keys %$hash){
      if (exists($hash->{'total'}) && $key ne 'total'){
        $self->index_table(\%{$hash->{$key}});
      }
      if ($key eq 'total'){
        $hash->{'index'} = $self->{'i'};
        $self->{'index'}[$self->{'i'}] = $hash->{'total'};
        $self->{'i'}++;
      }
    }
  }
  
  sub binomial{
    # Logarithmic!
    my $self = shift;
    my $n    = shift;
    my $k    = shift;
    if ($k > $n ||
        $k < 0  ||
        $n < 0){
      return 0;
    }
    my $result = 0;
    for (my $i = 1; $i <= $k; $i++){
      $result += log10(1 + (($n - $k)/$i));
    }
    return $result;
  }
    
  sub multinomial{
    # Logarithmic!
    my $self = shift;
    my @coeffs = sort{$a <=> $b} @_;
    my $result = 0;
    my $cummul = 0;
    foreach my $coeff (@coeffs){
      my $upper = $coeff + $cummul;
      $result += $self->binomial($upper, $coeff);
      $cummul = $upper;
    }
    return $result;
  }

  
  sub debug{
    my $hash = shift;
    my $i = 0;
    my $level = 0;
    my $uid = "\&__$i\__";
    my $xml = '<root>'.$uid.'</root>';
    my $id_list = {$uid => $hash};
    while (scalar keys %$id_list){
      my $new_id_list = {};
      $level++;
      foreach my $id (keys %$id_list){
        my $temp_xml = '';
        my $href = $id_list->{$id};
        if (ref($href) eq 'ARRAY'){
          my $counter = 0;
          foreach my $val (@$href){
            $i++;
            $uid = "\&__$i\__";
            $new_id_list->{$uid} = $val;
            $temp_xml .= "\<c_$level\_$counter\>$uid\<\/c_$level\_$counter\>";
            $counter++;
          }
        }
        elsif (ref($href) eq 'HASH' || ref($href) eq __PACKAGE__){
          foreach my $key (keys %$href){
            $i++;
            $uid = "\&__$i\__";
            $new_id_list->{$uid} = $href->{$key};
            my $safe = '';
            if (substr($key,0,1) =~ /[^a-zA-Z]/){
              $safe = 'a';
            }
            $temp_xml .= "\<$safe$key\>$uid\<\/$safe$key\>";
          }
        }
        else{
          $href =~ s/[<>]//g;
          $href = '<![CDATA['.$href.']]>';
          $temp_xml .= $href;
        }
        $temp_xml = '_empty_' unless ($temp_xml);
        die ("$id\t$temp_xml\n") unless ($xml =~ /$id/);
        $xml =~ s/$id/$temp_xml/;
      }
      $id_list = $new_id_list;
    }
    return $xml;
  }

  
1;

=head1 NAME

Sidron, a program to predict genotypes from ultra-sequencing data

=head1 SYNOPSIS

perl sidron.pl pileup_file hash_table_file [pileup_file hash_table_file] ...

=head1 DESCRIPTION

Sidron sets a score for each pileup line in the input that can be used to
predict the genotype of each position. In its present form, only two genotypes
are taken into account: homozygous for the most represented base or heterozygous
for the most and second most represented base. The score represents the
probability to get the configuration in the pileup given the heterozygous
genotype divided by the probability to get the same configuration given the
homozygous genotype. To compute those probabilities, Sidron reads a table
from hash_table_file that estimates the probability to get a base read given
a genomic base. This table must be a hash stored with the module Storable.

=head1 DEPENDENCIES

Sidron uses strict, Storable, and, optionally, IO::Uncompress::Gunzip.

=head1 INPUT

The input is composed of one or several pairs of pileup file and its
corresponding table. The pileup file can be compressed if it can be
decompressed with IO::Uncompress::Gunzip. The hash table must be a hash with
the following structure:

 $hash{
   expected_base => {
     read_base => {
       'qual' => [
         0 => {'total' => t0},
         1 => {'total' => t1},
         ...
       ]
     }
     read_base => ...
   }
   expected_base => ...
 }
 

At this point, the vector under 'qual' has values from 0 to 3, encompassing
qualities 0-10, 11-20, 21-30, and 31-40, respectively. Values t0, t1, ...
mean the probability to get I<read_base with> a given quality if the genomic
base is I<expected_base>.

=head1 OUTPUT

The output is sent to stdout. Each line contains each of the pileup fields for
each input pileup_file plus each possible genotype and its corresponding
score. Genotypes are expressed as a two-letter field. The heterozygous
genotype is the field, and the homozygous genotype is the first letter repeated.

=head1 LICENSE

Sidron

(c) Universidad de Oviedo, 2011

Authors: Victor Quesada, Xose S. Puente, and Carlos Lopez-Otin

=cut
