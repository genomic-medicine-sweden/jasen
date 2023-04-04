#!/usr/bin/env perl

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};
my $stats;

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use POSIX;
#use Math::Round qw/nearest/;

use Data::Dumper;
use threads;
use Thread::Queue;
use threads::shared;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDOUT "$0: ".(caller(1))[3].": @_\n";}

my %threadStatus:shared;
my @IOextensions=qw(.fastq .fastq.gz .fq);

exit(main());

sub main() {
  my $settings={
    #poly=>1, # number of reads per group (1=SE, 2=paired end)
    qualOffset=>33,
    # trimming
    min_quality=>30,
    bases_to_trim=>20, # max number of bases that will be trimmed on either side
    # cleaning
    min_avg_quality=>10,
    min_length=>62,# twice kmer length sounds good
    singletons=>1, # yes, output singletons
    'trim-on-average'=>0, # trimming until a base > min_quality is found
    trim=>1,              # yes perform trimming
    removeDots=>1,        # yes remove dots
  };
  
  GetOptions($settings,qw(poly=i infile=s@ outfile=s min_quality=i bases_to_trim=i min_avg_quality=i  min_length=i quieter trim! removeDots! debug qualOffset=i numcpus=i singletons! trim-on-average! auto)) or die "Error in command line arguments";
  $$settings{numcpus}||=getNumCPUs();
  $$settings{'zero-quality'}=chr($$settings{qualOffset}); # zero quality
  
  my $infile=$$settings{infile} or die "Error: need an infile\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: need an outfile\n".usage($settings);

  my $outfileDir;
  ($$settings{outBasename},$outfileDir,$$settings{outSuffix}) = fileparse($outfile,@IOextensions);
  my $singletonOutfile="$outfileDir/$$settings{outBasename}.singletons$$settings{outSuffix}";

  if($$settings{trim}){
    if($$settings{'trim-on-average'}){
      logmsg "Trimming on an expanding average";
    } else {
      logmsg "Trimming on a simple threshold";
    }
  } else {
    logmsg "I will not perform trimming because the user specified --notrim";
  }

  # trimming/filtering starts here
  my $printQueue=Thread::Queue->new;
  my $singletonPrintQueue=Thread::Queue->new;
  my $printThread=threads->new(\&printWorker,$outfile,$printQueue,$settings);
  my $singletonPrintThread=threads->new(\&singletonPrintWorker,$singletonOutfile,$singletonPrintQueue,$settings);
  
  my $entryCount=qualityTrimFastqPoly($infile,$printQueue,$singletonPrintQueue,$settings);
  
  $singletonPrintQueue->enqueue(undef); # term signal
  $singletonPrintThread->join;
  $printQueue->enqueue(undef); # term signal
  $printThread->join;
  
  my $numGood=sum(values(%threadStatus));
  my $freq_isClean=$numGood/$entryCount;
  $freq_isClean=nearest(0.01,$freq_isClean);
  logmsg "\nFinished! $freq_isClean pass rate.";

  return 0;
}

# returns a reference to an array of fastq entries.
# These entries are each a string of the actual entry
sub qualityTrimFastqPoly($;$){
  my($fastq,$printQueue,$singletonQueue,$settings)=@_;
  my $verbose=!$$settings{quieter};
  my $reportEvery=1000000;
    $reportEvery=1000 if($$settings{debug});

  # TODO check for poly with each fastq file
  my $poly=$$settings{poly}||checkPolyness($$fastq[0],$settings);
  $poly=1 if($poly<1); #sanity check
  $$settings{poly}=$poly; # this line is a band-aid until I pass the parameter directly the workers later
  
  # Choose automatic values before the threads are initialized.
  $settings=autoChooseParameters($$fastq[0],$settings) if($$settings{auto});
  # TODO insert these new settings into the threads
  # TODO choose new settings for each fastq

  # initialize the threads
  my (@t,$t);
  my $Q=Thread::Queue->new;
  for(0..$$settings{numcpus}-1){
    $t[$t++]=threads->create(\&trimCleanPolyWorker,$Q,$printQueue,$singletonQueue,$settings);
  }

  my $entryCount=0;
  my @fastq=@$fastq;
  for my $fastq (@fastq){

    logmsg "Trimming and cleaning a file $fastq with poly=$poly with a quality offset of $$settings{qualOffset}";
  
    # check the file before continuing
    my $readsToCheck=10;
    my $warnings=checkFirstXReads($fastq,$readsToCheck,$settings);
    warn "WARNING: There were $warnings warnings. Judging from the first $readsToCheck reads, it is possible that many or all sequences will be filtered out." if($warnings);

    # load all reads into the threads for analysis
    my $linesPerGroup=4*$poly;
    my $moreLinesPerGroup=$linesPerGroup-1; # calculate that outside the loop to save CPU
    if($$settings{inSuffix}=~/\.fastq$/){
      open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
    }
    elsif($$settings{inSuffix}=~/\.fastq\.gz/){
      open(FQ,"gunzip -c $fastq | ") or die "Could not open $fastq for reading: $!";
    }
    else{
      die "Could not determine the file type for reading based on your extension $$settings{inSuffix}";
    }

    my @entry; # buffer for the queue
    while(my $entry=<FQ>){
      $entry.=<FQ> for(2..$linesPerGroup); # e.g. 8 lines total for paired end entry
      push(@entry,$entry);
      $entryCount++;

      if($entryCount % $reportEvery == 0){
        $Q->enqueue(@entry);
        @entry=();
        if($verbose){
          my $numGood=sum(values(%threadStatus));
          my $freq_isClean=$numGood/($entryCount-$Q->pending);
          $freq_isClean=nearest(0.01,$freq_isClean);
          logmsg "Finished loading $entryCount pairs or singletons ($freq_isClean pass rate).";
          if($$settings{debug}){
            logmsg "DEBUG";
            last;
          }
        }
      }
    }
    close FQ;
    $Q->enqueue(@entry); # get the last pieces of it entered

  }
  logmsg "Done loading: $entryCount entries loaded.";
  
  # stop the threads
  $Q->enqueue(undef) for(1..$$settings{numcpus});
  queue_status_updater($Q,$entryCount,$settings);

  # grab all the results from the threads
  for(@t){
    my $tmp=$_->join;
  }
  return $entryCount;
}

sub trimCleanPolyWorker{
  my($Q,$printQueue,$singletonQueue,$settings)=@_;
  my $poly=$$settings{poly};
  my $notrim=$$settings{notrim};
  my $tid="TID".threads->tid;
  my(@entryOut);
  ENTRY:
  while(defined(my $entry=$Q->dequeue)){
    my $is_singleton=0; # this entry is not a set of singletons until proven otherwise
    my @read;
    my @entryLine=split(/\s*\n\s*/,$entry);
    for my $i (0..$poly-1){
      my $t={};
      ($$t{id},$$t{seq},undef,$$t{qual})=splice(@entryLine,0,4);
      dotsToNs($t,$settings);
      trimRead($t,$settings) if(!$notrim);
      $is_singleton=1 if(!read_is_good($t,$settings));
      $read[$i]=$t;
    }
    
    # see if these are singletons, and if so, add them to the singleton queue
    if($is_singleton){
      my @singleton;
      for my $read (@read){
        next if(!read_is_good($read,$settings));
        my $singletonOut="$$read{id}\n$$read{seq}\n+\n$$read{qual}\n";
        $singletonQueue->enqueue($singletonOut);
        $threadStatus{$tid}+=1/$poly;
      }
      # If these are singletons then they cannot be printed to the paired end output file. 
      # Move onto the next entry.
      next ENTRY; 
    }
    
    $threadStatus{$tid}++;
    my $entryOut="";
    for my $read (@read){
      $entryOut.="$$read{id}\n$$read{seq}\n+\n$$read{qual}\n";
    }
    $printQueue->enqueue($entryOut);
  }
  #return \@entryOut;
}

sub printWorker{
  my($outfile,$Q,$settings)=@_;
  if($$settings{outSuffix}=~/\.fastq$/){
    open(OUT,">",$outfile) or die "Error: could not open $outfile for writing because $!";
  }
  elsif($$settings{outSuffix}=~/\.fastq\.gz/){
    open (OUT, "| gzip -c > $outfile") or die "Error: could not open $outfile for writing because $!";
  }
  else{
    die "Could not determine the file type for writing based on your extension $$settings{outSuffix}";
  }
  while(defined(my $line=$Q->dequeue)){
    print OUT $line;
  }
  close OUT;
  logmsg "Reads are in $outfile";
}
sub singletonPrintWorker{
  my($outfile,$Q,$settings)=@_;

  # discard all singletons if the user does not want them
  if(!$$settings{singletons}){
    logmsg "User has specified no singletons, so I will be discarding singleton reads";
    while(defined(my $line=$Q->dequeue)){
      # do nothing except dequeue
    }
    return 0;
  }
  # END discarding singletons


  logmsg "User has specified singletons, so I will be placing singleton reads in $outfile";
  if($$settings{outSuffix}=~/\.fastq$/){
    open(OUT,">",$outfile) or die "Error: could not open $outfile for writing because $!";
  }
  elsif($$settings{outSuffix}=~/\.fastq\.gz/){
    open (OUT, "| gzip -c > $outfile") or die "Error: could not open $outfile for writing because $!";
  }
  else{
    die "Could not determine the file type for writing based on your extension $$settings{outSuffix}";
  }
  while(defined(my $line=$Q->dequeue)){
    print OUT $line;
  }
  close OUT;
  logmsg "Singleton reads are in $outfile";
}

# Trim a read using the trimming options
# The read is a hash of id, sequence, and qual
sub trimRead{
  my($read,$settings)=@_;
  # TODO just kill off this read right away if it is less than the min_length, to save on time
  my @qual=map(ord($_)-$$settings{qualOffset},split(//,$$read{qual}));
  my($numToTrim5,$numToTrim3)=numToTrim(\@qual,$read,$settings);
  if($numToTrim3 || $numToTrim5){
    $$read{seq}=substr($$read{seq},$numToTrim5,$$read{length}-$numToTrim5-$numToTrim3);
    $$read{qual}=substr($$read{qual},$numToTrim5,$$read{length}-$numToTrim5-$numToTrim3);
  }
}

# make any dot into an N in the read, and also convert any N's quality to 0
sub dotsToNs{
  my($read,$settings)=@_;
  return $read if($$read{seq} !~/\.|N|n/);
  my $zero=$$settings{'zero-quality'};

  # find which are naughty characters
  my @index=();
  while($$read{seq}=~/(\.|N|n)/g){
    push(@index,length($`));
  }

  # zero out these positions
  my @qual=split(//,$$read{qual});
  my @seq=split(//,$$read{seq});
  for my $i(@index){
    #next if($seq[$i]!~/\.nN/);
    $seq[$i]='N';
    $qual[$i]=$zero;
  }
  
  $$read{qual}=join("", @qual);
  $$read{seq}=join("",@seq);
}

# returns ($numToTrim5,$numToTrim3)
sub numToTrim{
  my($qual,$read,$settings)=@_;
  my @rQual=reverse(@$qual); # for 3' trimming
  my $length=$$read{length};
  # TODO set bases_to_trim equal to the read length if larger than
  #   ...or should I?  Might be a huge time cost.

  my($numToTrim3,$numToTrim5)=(0,0);
  if($$settings{'trim-on-average'}){
    for(my $i=0;$i<$$settings{bases_to_trim};$i++){
      if(sum(@$qual[0..$i])/($i+1)<$$settings{min_quality}){
        $numToTrim5++;
      } else {
        last if($i>5); # check minimum 5 bases before quitting
      }
    }
    for(my $i=0;$i<$$settings{bases_to_trim};$i++){
      if(sum(@rQual[0..$i])/($i+1)<$$settings{min_quality}){
        $numToTrim3++;
      } else {
        last if($i>5); # check minimum 5 bases before quitting
      }
    }
  } 
  # ELSE trimming until the first good quality
  else {
    for(my $i=0;$i<$$settings{bases_to_trim};$i++){
      if($$qual[$i]<$$settings{min_quality}){
        $numToTrim5++;
      } else {
        last;
      }
    }
    for(my $i=0;$i<$$settings{bases_to_trim};$i++){
      if($rQual[$i]<$$settings{min_quality}){
        $numToTrim3++;
      } else {
        last;
      }
    }
  }
  return($numToTrim5,$numToTrim3);
}

# Deterimine if a single read is good or bad (0 or 1)
# judging by the cleaning options
sub read_is_good{
  my($read,$settings)=@_;
  my $readLength=length($$read{seq});
  #die "short read: $$read{seq}\n$readLength < $$settings{min_length}\n" if($readLength<$$settings{min_length});
  return 0 if($readLength<$$settings{min_length});
  
  # avg quality
  my $qual_is_good=1;
  my @qual=map(ord($_)-$$settings{qualOffset},split(//,$$read{qual}));
  my $avgQual=sum(@qual)/$readLength;
  return 0 if($avgQual<$$settings{min_avg_quality});
  return 1;
}

sub queue_status_updater{
  my($Q,$entryCount,$settings)=@_;
  STATUS_UPDATE:
  while(my $p=$Q->pending){
    my $numGood=sum(values(%threadStatus));
    my $denominator=$entryCount-$Q->pending || 9e10;
    my $freq_isClean=$numGood/$denominator;
    $freq_isClean=nearest(0.01,$freq_isClean);
    logmsg "$p entries left to process. $freq_isClean pass rate.";
    # shave off some time if the queue empties while we are sleeping
    for(1..30){
      sleep 2;
      last STATUS_UPDATE if(!$Q->pending);
    }

  }
  return 1;
}

# TODO use this subroutine to call this script on the first X reads so that there is not any more code duplication
sub checkFirstXReads{
  my($infile,$numReads,$settings)=@_;

  (undef,undef,$$settings{inSuffix}) = fileparse($infile,@IOextensions);

  my $warning=0;
  die "Could not find the input file $infile" if(!-e $infile);
  if(-s $infile < 1000){
    warn "Warning: Your input file is pretty small";
    $warning++;
  }

  # figure out if the min_length is larger than any of the first X reads.
  # If so, spit out a warning. Do not change the value. The user should be allowed to make a mistake.
  logmsg "Checking minimum sequence length param on the first $numReads sequences.";
  if($$settings{inSuffix}=~/\.fastq$/){
    open(IN,'<',$infile) or die "Could not open $infile for reading: $!";
  }
  elsif($$settings{inSuffix}=~/\.fastq\.gz/){
    open(IN,"gunzip -c $infile | ") or die "Could not open $infile for reading: $!";
  }
  else{
    die "Could not determine the file type for reading based on your extension $$settings{inSuffix}";
  }
  my $i=1;
  while(my $line=<IN>){
    $line.=<IN> for (1..3);
    my ($id,$sequence)=(split("\n",$line))[0,1];
    my $length=length($sequence);
    if($length<$$settings{min_length}){
      warn "WARNING: Read $id has length less than the minimum sequence length $$settings{min_length}.\n" if(!$$settings{quieter});
      $warning++;
    }
    last if($i++>=$numReads);
  }
  close IN;
  return $warning;
}

################
## utility
################

sub autoChooseParameters{
  my($infile,$settings)=@_;
  my %newSettings=%$settings; # copy the hash
  
  # find the metrics and put them into a hash
  my $readMetrics=`run_assembly_readMetrics.pl '$infile' --fast`;
  chomp($readMetrics);
  my($header,$values)=split /\n/, $readMetrics;
  my @header=split /\t/,$header;
  my %metric;
  @metric{@header}=split /\t/,$values;
  
  # suggest some values in settings
  #$newSettings{bases_to_trim}=10;
  $newSettings{min_length}=$metric{avgReadLength} - $newSettings{bases_to_trim};
  $newSettings{min_avg_quality}=$metric{avgQuality} - 5;
  $newSettings{min_quality}=$metric{avgQuality} - 5;
  #print Dumper [\%metric,\%newSettings];die;

  my $newOptions; $newOptions.="$_: $newSettings{$_}, " for(qw(bases_to_trim min_length min_avg_quality min_quality));
  $newOptions=~s/,$//;
  logmsg "Auto-choose was specified. New options are now: $newOptions";
  
  return \%newSettings;
}

# returns 1 for SE, 2 for PE, and -1 is for internal error
sub checkPolyness{
  my ($infile,$settings)=@_;
  eval{
    require AKUtils;
  };
  if($@){
    warn "WARNING: I could not understand if the file is paired end.  I will assume not.";
    return 1;
  }

  my $poly=AKUtils::is_fastqPE($infile,{checkFirst=>20});
  $poly++;

  return $poly;
}

# until I find a better way to round
sub nearest{
  my ($nearestNumber,$num)=@_;
  $num=floor(($num+0.00)*100)/100;
  $num=sprintf("%.2f",$num);
  return $num;
}

sub getNumCPUs() {
  my $num_cpus;
  open(IN, '<', '/proc/cpuinfo'); while (<IN>) { /processor\s*\:\s*\d+/ or next; $num_cpus++; } close IN;
  return $num_cpus || 1;
}

sub usage{
  my ($settings)=@_;
  "Trim and clean a set of raw reads. Interleved reads if they are paired-end.
  Usage: $0 -i reads.fastq -o reads.filteredCleaned.fastq [-p 2]
    -i input file in fastq format
    -o output file in fastq format
  Additional options
  
  -p 1 or 2 (p for poly)
    1 for SE, 2 for paired end (PE). 0 for automatic detection (default)
  -quieter for somewhat quiet mode (use 1>/dev/null for totally quiet)
  --notrim to skip trimming of the reads. Useful for assemblers that require equal read lengths.
  --noremoveDots to skip the changing of dots to Ns and altering their qualities to zero
  -qual 64 to use an offset of 64 instead of 33(default).
  --numcpus 1 or --numcpus 2 for single or multithreaded. Currently: $$settings{numcpus}
  --nosingletons Do not output singleton reads, which are those whose pair has been filtered out.
  --trim-on-average Trim until an average quality is reached on the 5' and 3' ends. Ordinarily, trimming will stop when a base > min_quality is reached.
  Use phred scores (e.g. 20 or 30) or length in base pairs if it says P or L, respectively
  --auto                      # to choose the following values automatically based on run_assembly_readMetrics.pl (experimental)
  --min_quality P             # trimming
    currently: $$settings{min_quality}
  --bases_to_trim L           # trimming
    currently: $$settings{bases_to_trim}
  --min_avg_quality P         # cleaning
    currently: $$settings{min_avg_quality}
  --min_length L              # cleaning
    currently: $$settings{min_length}
  "
}
