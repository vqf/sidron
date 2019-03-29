#!/usr/bin/perl
use strict;
use Cwd;


my @questions = qw/
  basename
  ref_genome
  read_folder
  file_pattern
  nproc
/;

my $defaults = {
  'read_folder' => getcwd(),
  'nproc' => 1,
  'PL' => 'ILLUMINA',
  'LB' => 'WGS'
};

my $tfile = shift;
my $info = ask(\@questions, $defaults, $tfile);

my $files = get_file_names($info->{'read_folder'}, $info->{'file_pattern'});
#if (!exists($files->{'paired'}) || !scalar(@{$files->{'paired'}})){
#  die("Could not find paired files in folder $info->{'read_folder'} with the pattern $info->{'file_pattern'}\n");
#}

#print join("\n", @{$files->{'paired'}[0]});
#print "\n----\n";
#print join("\n", @{$files->{'paired'}[1]});
#print "\n----\n";
#print join("\n", @{$files->{'unpaired'}});

my $header = <<END;
#!/bin/bash
#Generated with /home/vqf/software/alig.pl
END

foreach my $k (@questions){
  $header .= "#$k\t$info->{$k}\n";
}

$header .= "####\nbwa=bwa\n";
$header .= "samtools=samtools\n";
$header .= "nproc=$info->{'nproc'}\n####\n";

my $PL = $defaults->{'PL'};
my $LB = $defaults->{'LB'};

my $catcommand = 'zcat';
my $outfolder = $info->{'basename'};

if (exists($files->{'unpaired'}) && scalar(@{$files->{'unpaired'}})){
  my $uheader = $header;
  my $n = scalar(@{$files->{'unpaired'}});
  warn("Found $n files classified as unpaired\n");
  $uheader .= "#Create output folder if necessary\n";
  $uheader .= "if [[ ! -d $outfolder ]]; then\nmkdir $outfolder\nfi\n";
  $uheader .= "#Files classified as unpaired\n";
  my $rfolder = $info->{'read_folder'};
  $rfolder .= '/' if (!($rfolder =~ /\/$/));
  $uheader .= 'read_folder="' . $rfolder . '"' . "\n";
  $uheader .= 'fastq1="' . "$rfolder" . join(" $rfolder", @{$files->{'unpaired'}}) . "\"\n";
  $uheader .= "#RG data\n";
  $uheader .= "ID=$outfolder\n";
  $uheader .= "SM=$outfolder\n";
  $uheader .= "PL=$PL\n";
  $uheader .= "LB=$LB\n";
  $uheader .= "PU=$(zcat $files->{'unpaired'}[0] | head -1 | sed 's/[:].*//' | sed 's/@//')\n";

  
  $uheader .= "#Alignment\n";
  $uheader .= "\$bwa mem -t \$nproc -a $info->{'ref_genome'} " .
              "<($catcommand \$fastq1\) | " .
              #"samblaster --discordantFile $outfolder\/$outfolder.discordant.sam --splitterFile $outfolder\/$outfolder.splitreads.sam | " .
              "\$samtools view -Shu - | " .
              "\$samtools sort -@ \$nproc -m 6000000000 -T $outfolder\/tmp -o $outfolder\/$outfolder.sorted.bam\n";
  
  $uheader .= "# Index\n";
  $uheader .= "cd $outfolder\n";
  $uheader .= "\$samtools index $outfolder.sorted.bam\n";
  warn("Writing $outfolder\.sh\n");
  my $outfile = "unpaired_$outfolder\.sh";
  open (OUT, ">$outfile") or die("Error: Cannot open file, $!\n");
  print OUT $uheader;
  close OUT;
  my $ex_command = "chmod 755 $outfile";
  `$ex_command`;
}


if (exists($files->{'paired'}) && scalar(@{$files->{'paired'}})){
  my $n = scalar(@{$files->{'paired'}});
  warn("Found $n files classified as paired\n");
  my $pheader = $header;
  $pheader .= "#Create output folder if necessary\n";
  $pheader .= "if [[ ! -d $outfolder ]]; then\nmkdir $outfolder\nfi\n";
  $pheader .= "#Files classified as paired\n";
  my $rfolder = $info->{'read_folder'};
  $rfolder .= '/' if (!($rfolder =~ /\/$/));
  $pheader .= 'read_folder="' . $rfolder . '"' . "\n";
  $pheader .= 'fastq1="' . "$rfolder" . join(" $rfolder", @{$files->{'paired'}[0]}) . "\"\n";
  $pheader .= 'fastq2="' . "$rfolder" . join(" $rfolder", @{$files->{'paired'}[1]}) . "\"\n";

  $pheader .= "#RG data\n";
  $pheader .= "ID=$outfolder\n";
  $pheader .= "SM=$outfolder\n";
  $pheader .= "PL=$PL\n";
  $pheader .= "LB=$LB\n";
  $pheader .= "PU=\$(zcat $files->{'paired'}[0][0] | head -1 | sed 's/[:].*//' | sed 's/@//')\n\n";

  $pheader .= "#Alignment\n";
  $pheader .= "\$bwa mem -t \$nproc -a $info->{'ref_genome'} " .
              "<($catcommand \$fastq1\) <($catcommand \$fastq2) | " .
              #"samblaster --discordantFile $outfolder\/$outfolder.discordant.sam --splitterFile $outfolder\/$outfolder.splitreads.sam | " .
              "\$samtools view -Shu - | " .
              "\$samtools sort -@ \$nproc -m 600000000 -T $outfolder\/tmp -o $outfolder\/$outfolder.sorted.bam\n";
  $pheader .= "# Index\n";
  $pheader .= "cd $outfolder\n";
  $pheader .= "\$samtools index $outfolder.sorted.bam\n";
  $pheader .= "\$samtools flagstat $outfolder.sorted.bam >$outfolder.sorted.flagstat\n";
  #$pheader .= "\$samtools view -b $outfolder.splitreads.sam | samtools sort - -o $outfolder.splitreads.sorted.bam &\n";
  #$pheader .= "\$samtools view -b $outfolder.discordant.sam | samtools sort - -o $outfolder.discordant.sorted.bam &\n";

  $pheader .= "cd ..\n";
  
  warn("Writing $outfolder\.sh\n");
  my $outfile = "$outfolder\.sh";
  open (OUT, ">$outfile") or die("Error: Cannot open file, $!\n");
  print OUT $pheader;
  close OUT;
  my $ex_command = "chmod 755 $outfile";
  `$ex_command`;
}

sub get_file_names{
  my $folder = shift;
  my $expr   = shift;
  my $result = {};
  opendir(DIR, $folder) or die ("Could not open $folder: $!\n");
  my $files = {};
  while (my $f = readdir(DIR)){
    $files->{$f} = 1;
  }
  closedir DIR;
  foreach my $f (sort{$a cmp $b} keys %$files){
    if ($f =~ /^$expr.*gz$/){
      if ($f =~ /^$expr.*?\_\D*?(\d+).*gz$/){
        my $nread = $1;
        if ($nread == 1){  
          my $paired = $f;
          $paired =~ s/^($expr)(.*?)\_(\D*?)1/$1$2\_${3}2/;
          #print "$paired\n";
          if (exists($files->{$paired})){
            push @{$result->{'paired'}[0]}, $f;
            push @{$result->{'paired'}[1]}, $paired;
            warn("Found a paired one: $f, $paired\n");
            next;
          }
        }
        elsif ($nread ==2){
          my $paired = $f;
          $paired =~ s/^($expr)\_2/$1\_1/;
          next if (exists($files->{$paired}));
        }
        else{
          push @{$result->{'unpaired'}}, $f;
        }
      }
      else{
        push @{$result->{'unpaired'}}, $f;
      }
    }
  }
  return $result;
}

sub ask{
  my $questions = shift;
  my $def       = shift || {}; #Defaults
  my $file      = shift || '';
  my $info = {};
  my $result = {};
  my $isfile = 1;
  if ($file && -e $file){
    open (IN, $file) or ($isfile = 0);
    if ($isfile){
      while (<IN>){
        chomp;
        s/^#//;
        s/\"//g;
        my ($field, $val) = split(/[\s=,]+/);
        $info->{$field} = $val;
      }
    }
    close IN;
  }
  foreach my $q (@$questions){
    if (exists($info->{$q})){
      $result->{$q} = $info->{$q};
    }
    else{
      my $extra = '';
      if (exists($def->{$q})){
        $extra = " \($def->{$q}\)";
      }
      $result->{$q} = prompt("$q$extra");
    }
    if(!$result->{$q} && exists($def->{$q})){
      $result->{$q} = $def->{$q};
    }
  }
  return $result;
}

sub prompt{
  my $q = shift;
  print "$q\n";
  my $r = <STDIN>;
  chomp $r;
  return $r;
}

sub slurp{
  my $file = shift;
  my $h    = shift;
  local $/ = undef;
  open (IN, $file) or die ("Could not open file: $!");
  my $result = <IN>;
  close IN;
  return $result;
}