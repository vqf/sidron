#!/usr/bin/perl -w

# Filter mutations found in polyN stretches
# A mutation is filtered away if the reads where it is
# present contain a polyN longer than or equal to $N

use strict;
use Storable;
use Bio::DB::HTS;

my $pfile  = shift;
my $bfile  = shift;
my $hfile  = shift || '';
my @temp = split /[\/\\]/, $0;
my $pname = pop @temp;
die("Use: perl $pname pileup_file bam_file\n")
  unless ($pfile && -e $pfile &&
          $bfile && -e $bfile);

my $DEBUG = 0;
my $TMP = 1;
my $CHECK_QDIST = 0;
my $OFFSET = 2;       # For blat filter
my $MINQUAL = 0;
my $DIST_LIMIT = 1e-3;
my $MAXP = 5e-5;      # Quality and strand filter
my $CLIMIT = 8;
my $QTRIM = 10;       # Max qual for trimming
my $MAXREADS = 1e4;
#my $samview_command = 'samtools view '.$bfile;
my $BLATPORT = '$BLATPORT';
if (-e 'config.txt'){
  open (IN, 'config.txt');
  while (<IN>){
    chomp;
    s/\s//g;
    my ($factor, $value) = split(/\=/);
    if ($factor eq 'BLATPORT'){
      $BLATPORT = $value;
      tolog("Setting port to $value\n");
    }
  }
  close IN;
}
my $blat_command = 'gfClient -out=blast8 localhost '. $BLATPORT . ' "" stdin stdout';
my $N = 10;
my $seq_field = 0;
my $isindel = 0;
my $CHUNK_SIZE = 1000;
my $qprobs = {};
#if ($hfile && -e $hfile){
#  $CHECK_QDIST = 1;
#  my $tmp = Storable::retrieve($hfile);  # If check_qdists, checks if the number of low quality mutant
#                                            # bases is significantly too high
#  foreach my $ebase (keys %$tmp){
#    if (ref($tmp->{$ebase}) eq 'HASH' && exists($tmp->{$ebase}{'total'})){
#      foreach my $fbase (keys %{$tmp->{$ebase}}){
#        if (ref($tmp->{$ebase}{$fbase}) eq 'HASH' && exists($tmp->{$ebase}{$fbase}{'total'})){
#          my $tt = $tmp->{$ebase}{$fbase}{'qual'}{'total'};
#          my $lq = $tmp->{$ebase}{$fbase}{'qual'}{0}{'total'} || 1;
#          my $hq = ($tt - $lq) || 1;
#          $qprobs->{$ebase}{$fbase}{'lq'} = $lq/$tt;
#          $qprobs->{$ebase}{$fbase}{'hq'} = $hq/$tt;
#        }
#      }
#    }
#  }
#}
my @QUALITY_GROUPS = (10, 20, 30);
if (exists($qprobs->{'quality_groups'}) && scalar(@{$qprobs->{'quality_groups'}})){
  @QUALITY_GROUPS = @{$qprobs->{'quality_groups'}};
}

open (IN, $pfile);
open (FLT, ">>filtered_out");
log("Reading file $pfile\n");
my $aggregate = 0;
while (! eof(IN)){
  my $info = {};
  my $transl = [];
  my $counter = 0;
  while (!eof(IN) && $counter < $CHUNK_SIZE){
    my $l = <IN>;
    $counter++;
    $aggregate++;
    my @fields = split(/\s+/, $l);
    my ($chr, $pos) = @fields[0..1];
    $isindel = 1 if ($fields[2] eq '*');
    if (!$seq_field && !$isindel){
      $seq_field = set_seq_qual_fields(@fields);
    }
    
    push @{$info->{$chr}}, $pos;
    push @$transl, \@fields;
  }
  my $eof = eof(IN);
  tolog("Finished reading file chunk. Read $aggregate lines. EOF is $eof\n");
  
  my $binfo = {};
  if ($DEBUG && -e 'dump.hsh'){
    tolog("Getting hash from dump.hsh\n");
    $binfo = Storable::retrieve('dump.hsh');
  }
  else{
    $binfo = get_bam_reads($info, $bfile);
    Storable::store($binfo, 'dump.hsh') if ($DEBUG);
  }
  
  my $median_filter = statFilter->new();
  
  foreach my $ps (@$transl){
    # Process new pileup for each pileup position
    my $chr = $ps->[0];
    my $pos = $ps->[1];
    my $ref = $ps->[2];
    my $seq = '';
    my $qua = '';
    my @tomedian = ();
    my $min = 0;
    my $max = 0;
    my $read_len = 0;
    my $lowqual = 0;
    my $n = 0;
    my $mut = '';
    my $nplus = 0;   #Number of reads in the plus strand
    my $nminus = 0;  #Number of reads in the minus strand
    my $names = {};  #Keep track of read names to avoid double callings
    my $QLIMIT = $MAXP;
    my $pcov = scalar(@{$binfo->{$chr}{$pos}});
    foreach my $bread (@{$binfo->{$chr}{$pos}}){
      my $rname = $bread->[0];
      if (exists($names->{$rname})){
        next;
      }
      $names->{$rname} = 1;
      my ($ns, $nq, $l, $v) = filter_breads($bread, $chr, $pos, $ref);
      if ($CHECK_QDIST && $ns ne '.' && $ns ne ','){
        my $qgroup = set_groups(ord($nq) - 33, @QUALITY_GROUPS);
        $lowqual++ unless ($qgroup);
        $n++;
        $mut = uc($ns) unless ($mut);
        if ($ns){
          if ($ns eq lc($ns)){ #Minus strand
            $nminus++;
          }
          elsif($ns eq uc($ns)){
            $nplus++;
          }
        }
      }
      if ((ord($nq) - 33) > $MINQUAL){
        $seq .= $ns;
        $qua .= $nq;
        if ($l && $v){
          $read_len = $l unless ($read_len);
          $min = $l unless ($min);
          my $corrected_pos = int($v*$read_len/$l);
          push @tomedian, $corrected_pos;
          $min = $corrected_pos if ($corrected_pos < $min);
          $max = $corrected_pos if ($corrected_pos > $max);
          #warn ("$v\t$read_len\t$l\t$corrected_pos\n");
        }
      }
    }
    my @nf = ();
    if ($isindel){
      @nf = @$ps;
    }
    else{
      @nf = (@$ps[0..$seq_field-2], length $seq, $seq, $qua);
                  #@$ps[$seq_field+2..(scalar @$ps)-1]);
    }
    if (scalar @tomedian){
      #my $med = median(@tomedian);
      #my ($low, $high) = $median_filter->get_limits($read_len, scalar @tomedian);
      #my $passed_median = 0;
      #if ($med >= $low && $med <= $high){
      #  $passed_median = 1;
      #}
      #else{
      #  warn join("\t", @nf)."\t$med, $low, $high\n";
      #}
      my $dlow = $median_filter->get_maxdist_limit($read_len, scalar @tomedian, $DIST_LIMIT);
      my $mdist = abs($max-$min);
      if ($mdist >= $dlow){
        print join("\t", @nf)."\n";
      }
      else{
        print FLT join("\t", @nf)."\tmaxdist: $mdist, mindist: $dlow\n";
      }
    }
    elsif($seq){
      print join("\t", @nf)."\n";
    }
  }
}
close FLT;
close IN;

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


sub median{
  my @v = sort {$a <=> $b} @_;
  my $total = scalar @v;
  my $mid = $total >> 1;
  return $v[$mid];
}

sub min{
  my $result = shift;
  foreach my $v (@_){
    $result = $v if ($v < $result);
  }
  return $result;
}

sub max{
  my $result = shift;
  foreach my $v (@_){
    $result = $v if ($v > $result);
  }
  return $result;
}



sub filter_breads{
  my $sam_fields = shift;
  my $chr        = shift;
  my $genome_pos = shift;
  my $ref_base   = shift;
  my $isindel = ($ref_base eq '*') ? 1 : 0;
  my @result = ('', '');
  my $poly = 0;
  my $diff  = $genome_pos - $sam_fields->[3];
  my $cigar = $sam_fields->[5];
  my @seq = split(//, $sam_fields->[9]);
  my @qua = split(//, $sam_fields->[10]);
  my $rtrim = pop @$sam_fields;
  my $ltrim = pop @$sam_fields;
  my $c = cigar_as_array($cigar);
  my $prev_mode = '';
  my $foundit = 0;
  my $cummul_length = 0;
  my $read_len = get_sam_len($cigar);
  while (!$foundit && defined(my $zone = shift @$c)){
    my ($mode, $length) = @$zone;
    if ($mode eq 'M'){
      $cummul_length += $length;
      if ($cummul_length > $diff){ # Bingo!
        
        $foundit = 1;
        # Find base
        # Trimmed?
        if ($diff <= $ltrim || $diff >= $rtrim){
          return @result;
        }
        my $found_base = $seq[$diff];
        return @result unless ($found_base);
        #print "$sam_fields->[0]\t$diff\n";
        my $found_qual = $qua[$diff];
        my $fbase = $found_base;
        #print "$found_base\t$ref_base\n";
        if ($found_base ne $ref_base){ # Base in SAM is mutated
          #BLAT filter
          #warn join("\t", @$sam_fields)."\n";
          #warn("$ltrim\t$rtrim\n");
          #warn "-----\n";
          my $comm = "echo \"\>tmp\n$sam_fields->[9]\" | $blat_command";
          my @brest = `$comm`;
          my @bres = sort {blatscore($b) <=> blatscore($a)} @brest;
          #print join("", @bres);
          #print "----\n";
          return @result unless ($bres[0]);
          my @f1 = split(/\s+/, $bres[0]);
          if (
                $f1[1] ne $chr ||
                ($f1[8]>$genome_pos + $OFFSET && $f1[9]>$genome_pos + $OFFSET) ||
                ($f1[8]<$genome_pos - $OFFSET && $f1[9]<$genome_pos - $OFFSET)
              ){
            return @result;
          }
          if (scalar @bres > 1){
            my @f2 = split(/\s+/, $bres[1]);
            my $bestqual = $f1[2];
            my $second   = $f2[2];
            if ($second == $bestqual){
              return @result;
            }
          }
          my $goon = 1;
          my $pos = $diff - 1;
          # Check backwards
          while ($goon && $pos){
            $goon = 0;
            my $this_base = substr $sam_fields->[9], $pos, 1;
            if ($this_base eq $found_base){
              $goon = 1;
              $poly++;
            }
            $pos--;
          }
          # Check forward
          $goon = 1;
          $pos = $diff + 1;
          $goon = 0 if $pos >= length($sam_fields->[9]);
          while ($goon && $pos){
            $goon = 0;
            my $this_base = substr $sam_fields->[9], $pos, 1;
            tolog(join("\n", @$sam_fields)."\n$pos\t".length($sam_fields->[9])) unless ($this_base);
            if ($this_base eq $found_base){
              $goon = 1;
              $poly++;
            }
            $pos++;
            $goon = 0 if $pos >= length($sam_fields->[9]);
          }
          if ($poly<$N){
            if (check_flags($sam_fields->[1], 16)){
              if ($found_base eq $ref_base){
                $fbase = ',';
              }
              else{
                $fbase = lc($fbase);
              }
            }
            elsif ($found_base eq $ref_base){
              $fbase = '.';
            }
            @result = ($fbase, $found_qual, $read_len, $diff + 1);
          }
          #else{
          #  print debug($sam_fields)."\n$genome_pos\t$poly\n";
          #}
        }
        else{
          if (check_flags($sam_fields->[1], 16)){
            if ($found_base eq $ref_base){
              $fbase = ',';
            }
            else{
              $fbase = lc($fbase);
            }
          }
          elsif ($found_base eq $ref_base){
            $fbase = '.';
          }
          @result = ($fbase, $found_qual);
        }
        
      }
      else{ # Tidy up and keep trying
        #$diff -= $length;
      }
      $prev_mode = 'M';
    }
    elsif ($mode eq 'N'){
      $diff -= $length;
      if ($diff <= 0){  # Position inside a gap
        return @result;
      }
      $prev_mode = 'N';
    }
    elsif ($mode eq 'I'){
      $diff += $length;
      #if ($cummul_length > $diff){ # Bingo!
      #  die('Insert');
      #}
    }
    elsif ($mode eq 'D'){
      $diff -= $length;
      if ($cummul_length > $diff){ # Bingo!
        return @result;
      }
    }
    elsif ($mode eq 'S'){
      #$diff += $length;
      #$cummul_length += $length;
      splice(@seq, $cummul_length, $length);
      splice(@qua, $cummul_length, $length);
      if ($cummul_length > $diff){ # Bingo!
        return @result;
      }
    }
  }
  return @result;
}

sub check_flags{
  my $flag = shift;
  my $mask = shift;
  return $flag & $mask;
}

sub blatscore{
  my $line = shift;
  my @t = split(/\s+/, $line);
  return $t[11];
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

sub is_integer{
  my $n = shift;
  my $result = 1;
  $result = 0 if ($n =~ /\D/);
  return $result;
}


sub set_seq_qual_fields{
  my $cov = shift;
  my $seq   = shift;
  my $result = 0;
  my $done = 0;
  while (!$done && scalar @_){
    my $qual = shift;
    my $cseq = clean_seq($seq);
    $result++;
    if (is_integer($cov) && $result>3){
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
    tolog("Could not detect sequence field. Defaulting to fifth column\n");
    $result = 4;
  }
  return $result;
}

sub tolog{
  my $t = shift;
  open (APP, ">>polylog.txt");
  print APP "$t\n";
  close APP;
}

sub get_sam_len{
  my $cigar = shift;
  my $len = 0;
  my $c = cigar_as_array($cigar);
  foreach my $cig (@$c){
    my ($mode, $length) = @$cig;
    if ($mode ne 'D' && $mode ne 'N'){
      $len += $length;
    }
  }
  return $len;
}

sub cigar_as_array{
  my $cigar = shift;
  my $result = [];
  my $oldcigar = 0;
  my @c = ();
  if ($cigar =~ /^\D/){
    @c = $cigar =~ /(\w\d+)/g;
    $oldcigar = 1;
  }
  else{
    @c = $cigar =~ /(\d+\w)/g;
  }
  foreach my $zone (@c){
    my ($length, $mode);
    ($length, $mode) = $zone =~ /(\d+)(\w)/ unless ($oldcigar);
    ($mode, $length) = $zone =~ /(\w)(\d+)/ if ($oldcigar);
    push @$result, [$mode, $length];
  }
  return $result;
}

sub get_bam_reads{
  my $info = shift;
  my $bfile = shift;
  my $bam = Bio::DB::HTS->new(-bam => $bfile, -autoindex => 1);
  my $result = {};
  foreach my $chr (keys %$info){
    tolog("Getting alignments from chromosome $chr\n");
    foreach my $pos (@{$info->{$chr}}){
      my $segment = $bam->segment($chr, $pos, $pos);
      my $iter = $segment->features(-iterator=>1);
      my $nseqs = 0;
      while (my $al = $iter->next_seq){
        # Build seq string
        my $flag = $al->flag();
        if (! ($flag & 0xD00)){  #Avoid duplicates, secondary and supplementary alignments
          my $qs = $al->query->dna;
          # Build score string
          my @sc = $al->qscore();
          # Set trim
          my $t1 = 0;  # Left trim
          my $t2 = scalar(@sc) - 1;  #Right trim
          # Left
          my $c = $t1;
          my $currentq = $sc[$c];
          while ($currentq < $QTRIM && $c < $t2){
            $c++;
            $currentq = $sc[$c];
          }
          $t1 = $c;
          #Right
          $c = $t2;
          $currentq = $sc[$c];
          while ($currentq < $QTRIM && $c > $t1){
            $c--;
            $currentq = $sc[$c];
          }
          $t2 = $c;
          my $qual = '';
          foreach my $sc (@sc){
            $qual .= chr($sc+33);
          }
          my $rarr = [$al->qname(), $al->flag(), $chr, $al->start(), $al->qual(), $al->cigar_str(), '=', 'not_det', 'not_det', $qs, $qual, $t1, $t2];
          push @{$result->{$chr}{$pos}}, $rarr;
          $nseqs++;
        }
        if ($nseqs > $MAXREADS){
          tolog("Cutting to $MAXREADS at $chr:$pos");
          last;
        }
      }
    }
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

package statFilter;

  sub new{
    my $pkg = shift;
    my $self = {};
    $self->{'cache'} = {};
    bless ($self, $pkg);
    return $self;
  }
  
  sub get_limits{
    my $self = shift;
    my $l    = shift;
    my $n    = shift;
    my $limit = shift || 1e-4;
    return 0 unless ($limit);
    my $result = 0;
    $n-- unless ($n%2);  # Simple model only valid for odd numbers
    my $L = $l + 1;
    my $str = "$l:$n:$limit";
    return 0 unless ($limit < 1);
    if (exists($self->{'cache'}{'median'}{$str})){
      $result = $self->{'cache'}{'median'}{$str};
    }
    else{
      $result = $self->get_median_limit($l, $n, $limit);
      $self->{'cache'}{'median'}{$str} = $result;
    }
    return ($result, $L - $result);
  }
  
  
  sub get_median_limit{
    my $self = shift;
    my $l = shift;
    my $n = shift;
    my $limit = shift || 1e-4;
    my $L = $l + 1;
    return 0 unless ($limit < 1);
    my $result = 0;
    my $cummul = 0;
    my $u = ($n-1)/2;
    my $v = 1;
    while($cummul<$limit){
      my $mult = $self->multinomial($u, $u, 1);
      my $p = $mult * (($v*($L-$v)/$L**2)**$u)/$L;
      $cummul += $p;
      $v++;
      $result++ if ($cummul<$limit);
    }
    return $result;
  }
  
  sub get_maxdist_limit{
    # Calculates the minimum maximal distance engulfing $n points in a read of length
    # $l with a confidence limit $limit
    my $self = shift;
    my $l = shift;
    my $n = shift;
    my $limit = shift || 1e-4;
    return 0 unless ($limit <= 1);
    my $result = 0;
    my $str = "$l:$n:$limit";
    if (exists($self->{'cache'}{'maxdist'}{$str})){
      $result = $self->{'cache'}{'maxdist'}{$str};
    }
    else{
      my $d = 0;
      my $p = 0;
      while($p<$limit){
        $p = ($l-$d) * (($d+1)/$l)**$n - ($l-$d-1) * ($d/$l)**$n;
        $d++;
        $result++ if ($p<$limit);
      }
      $self->{'cache'}{'maxdist'}{$str} = $result;
    }
    return $result;
  }
  
  sub binomial{
    my $self = shift;
    my $n    = shift;
    my $k    = shift;
    if ($k > $n ||
        $k < 0  ||
        $n < 0){
      return 0;
    }
    my $result = 1;
    for (my $i = 1; $i <= $k; $i++){
      $result *= 1 + (($n - $k)/$i);
    }
    return $result;
  }
  
  sub multinomial{
    my $self = shift;
    my @coeffs = sort{$a <=> $b} @_;
    my $result = 1;
    my $cummul = 0;
    foreach my $coeff (@coeffs){
      my $upper = $coeff + $cummul;
      $result *= $self->binomial($upper, $coeff);
      $cummul = $upper;
    }
    return $result;
  }
1;

