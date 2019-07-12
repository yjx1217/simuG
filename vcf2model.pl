#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);

##############################################################
#  script: vcf2model.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  version: 1.0.0
#  last edited: 2019.07.12
#  description: vcf2model.pl will characterize mutational profile for SNPs and INDELs
#               based on user-provided vcf file.
##############################################################

my ($refseq, $vcf, $prefix, $qual_cutoff, $excluded_chr_list, $help, $man, $version);
$qual_cutoff = 0;

GetOptions(
  'help|h|?' => \$help,
  'man|m' => \$man,
  'version|v' => \$version,
  'vcf|i:s' => \$vcf,
  'excluded_chr_list|e:s' => \$excluded_chr_list,
  'prefix|p:s' => \$prefix,
  'qual|q:f' => \$qual_cutoff,
);

## Option check and print out Usage if needed.
# If usage (or version number) was explicitly requested, print usage (or version number).

if (defined $help) {
  pod2usage(-verbose => 1);
}

if (defined $man) {
  pod2usage(-verbose => 2);
}

if (defined $version) {
  pod2usage(-verbose => 99, -sections => qw(NAME|DESCRIPTION|VERSION));
}

if (not defined $vcf) {
  pod2usage(-message => "Mandatory argument '-vcf' is missing! Exit!",
    -exitval => 1,
	  -verbose => 1);
}

# define excluded chromosome(s) if any
my %excluded_chr_list = ();
if (defined $excluded_chr_list) {
  my $excluded_chr_list_fh = read_file($excluded_chr_list);
  %excluded_chr_list = parse_list_file($excluded_chr_list_fh);
}

my $vcf_fh = read_file($vcf);
my $snp_model_output = "$prefix.SNP_model.txt";
my $snp_model_output_fh = write_file($snp_model_output);
my $indel_model_output = "$prefix.INDEL_model.txt";
my $indel_model_output_fh = write_file($indel_model_output);

my %snp = ();
my $total_snp_count = 0;
my $total_transition_snp_count = 0;
my $total_transversion_snp_count = 0;

my %indels = ();
my %indel_freq_by_size = ();
my $total_deletion_count = 0;
my $total_insertion_count = 0;
my $total_homopolymer_indel_count = 0;

while (<$vcf_fh>) {
  chomp;
  /^#/ and next;
  /^\s*$/ and next;
  my @line = split /\t/, $_;
  my $chr = $line[0];
  my $pos = $line[1];
  my $ref_allele = $line[3];
  my $alt_allele = $line[4];
  my $qual = $line[5];
  my $filter = $line[6];
  if (not exists $excluded_chr_list{$chr}) {
	  if (($qual eq '.') or ($qual >= $qual_cutoff)) {
	    if (($ref_allele !~ /(N|n)/) and ($alt_allele !~ /(N|n)/)) {
		    if ($alt_allele !~ /,/) {
		      if ((length $ref_allele) ne (length $alt_allele)) {
			      # this is an deletion
			      $indels{$chr}{$pos}{'ref_allele'} = $ref_allele;
			      $indels{$chr}{$pos}{'alt_allele'} = $alt_allele;
			      $indels{$chr}{$pos}{'indel_size'} = (length $alt_allele) - (length $ref_allele);
			      my $indel_size_abs = abs($indels{$chr}{$pos}{'indel_size'});
			      if (exists $indel_freq_by_size{$indel_size_abs}) {
			        $indel_freq_by_size{$indel_size_abs}++;
			      } else {
			        $indel_freq_by_size{$indel_size_abs} = 1;
			      }
			      if ($indels{$chr}{$pos}{'indel_size'} > 0) {
			        $indels{$chr}{$pos}{'note'} = "insertion";
			        $total_insertion_count++;
			      } else {
			        $indels{$chr}{$pos}{'note'} = "deletion";
			        $total_deletion_count++;
			      }
		      } else {
			      # this is a SNP
			      $total_snp_count++;
			      if (exists $snp{$ref_allele}{$alt_allele}) {
			        $snp{$ref_allele}{$alt_allele}{'count'}++;
			      } else {
			        $snp{$ref_allele}{$alt_allele}{'count'} = 1;
			      }
			      # check transition vs. transversion
			      if (pupy($ref_allele) eq pupy($alt_allele)) {
			        $total_transition_snp_count++;
			      } else {
			        $total_transversion_snp_count++;
			      }
		      }
		    }
	    }
	  }
  }
}

my $titv_ratio;
if ($total_transversion_snp_count > 0) {
  $titv_ratio = $total_transition_snp_count/$total_transversion_snp_count;
} else {
  $titv_ratio = "Inf";
  print "Warning! total_transversion_snp_count = 0 >> set titv_ratio = \"Inf\"\n";
}
print $snp_model_output_fh "### SNP model ###\n";
print $snp_model_output_fh "titv_ratio=$titv_ratio\n";
print $snp_model_output_fh "#SNP_type\tfreq\n";

my @base = qw(A T G C);
foreach my $ref_allele (@base) {
  foreach my $alt_allele (@base) {
	  if ($ref_allele ne $alt_allele) {
	    if ($total_snp_count > 0) {
		    if (not exists $snp{$ref_allele}{$alt_allele}{'count'}) {
		      $snp{$ref_allele}{$alt_allele}{'freq'} = 0;
		    } else {
		      $snp{$ref_allele}{$alt_allele}{'freq'} = $snp{$ref_allele}{$alt_allele}{'count'}/$total_snp_count;
		    }
		    print $snp_model_output_fh "$ref_allele->$alt_allele\t$snp{$ref_allele}{$alt_allele}{'freq'}\n";
	    } else {
		    print "!!! Warning! total_snp_count = 0! Set the frequency of $ref_allele->$alt_allele substitution as NA!\n";
		    print $snp_model_output_fh "$ref_allele->$alt_allele\tNA\n";
	    }
	  }
  }
}

my $total_indel_count = $total_insertion_count + $total_deletion_count;
my $ins_del_ratio;
if ($total_deletion_count > 0) {
  $ins_del_ratio = $total_insertion_count/$total_deletion_count;
} else {
  $ins_del_ratio = "Inf";
  print "!!! Warning! total_deletion_count = 0! Set ins_del_ratio = \"Inf\"!\n";
}
my $accumulated_indel_freq = 0;
my %indel_freq = ();
for (my $indel_size = 50; $indel_size >= 0; $indel_size--) {
  my $size_specific_count = 0;
  if (exists $indel_freq_by_size{$indel_size}) {
	  $size_specific_count = $indel_freq_by_size{$indel_size};
  }
  if ($indel_size > 1) {
	  if ($total_indel_count > 0) {
	    $indel_freq{$indel_size} = $size_specific_count/$total_indel_count;
	    $accumulated_indel_freq += $indel_freq{$indel_size};
	  } else {
	    print "!!! Warning! total_indel_count = 0! Set the frequency of ${indel_size}-bp INDEL as NA!\n";
	    $indel_freq{$indel_size} = "NA";
	    $accumulated_indel_freq = "NA";
	  }
  } else {
	  if ($total_indel_count > 0) {
	    $indel_freq{$indel_size} = 1 - $accumulated_indel_freq;
	  } else {
	    print "!!! Warning! total_indel_count = 0! Set the frequency of ${indel_size}-bp INDEL as NA!\n";
	    $indel_freq{$indel_size} = "NA";
	  }
  }
}

print $indel_model_output_fh "### INDEL model ###\n";
print $indel_model_output_fh "# ins_del_ratio=total_number_of_insertions/total_number_of_deletions\n";
print $indel_model_output_fh "ins_del_ratio=$ins_del_ratio\n";
print $indel_model_output_fh "# indel_size(bp)\tindel_freq\n";

for (my $indel_size = 1; $indel_size <= 50; $indel_size++) {
  print $indel_model_output_fh "$indel_size\t$indel_freq{$indel_size}\n";
}


sub read_file {
  my $file = shift @_;
  my $fh;
  if ($file =~ /\.gz$/) {
    open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
  } else {
    open($fh, $file) or die "can't open $file";
  }
  return $fh;
}

sub write_file {
  my $file = shift @_;
  my $fh;
  if ($file =~ /\.gz$/) {
    open($fh, "| gzip -c >$file") or die "can't open $file\n";
  } else {
    open($fh, ">$file") or die "can't open $file\n";
  }
  return $fh;
}

sub parse_fasta_file {
  my ($fh, $input_hashref, $input_arrayref) = @_;
  my $seq_name = "";
  while (<$fh>) {
    chomp;
    if (/^\s*$/) {
      next;
    } elsif (/^\s*#/) {
	    next;
	  } elsif (/^>(\S+)/) {
	    $seq_name = $1;
	    push @$input_arrayref, $seq_name;
	    $$input_hashref{$seq_name} = "";
	  } else {
	    my $seq_line = uc $_;
	    $$input_hashref{$seq_name} .= $seq_line;
	  }
  }
}

sub parse_list_file {
  my $fh = shift @_;
  my %list = ();
  while (<$fh>) {
    chomp;
	  /^\s*$/ and next;
	  /^#/ and next;
	  if (exists $list{$_}) {
	    $list{$_}++;
	  } else {
	    $list{$_} = 1;
	  }
  }
  return %list;
}

sub pupy {
  my $base = shift @_;
  my $result;
  if ($base =~ /(A|a)/) {
	  $result = "purine";
  } elsif ($base =~ /(G|g)/) {
    $result = "purine";
  } elsif ($base =~ /(T|t)/) {
	  $result = "pyrimidine";
  } elsif ($base =~ /(C|c)/) {
	  $result = "pyrimidine";
  } else {
	  $result = "unknown";
  }
  return $result;
}

# sub find_motif {
#   my ($seq, $motif_regexp) = @_ ;
#   my $i = 0;
#   my %match = ();
#   $seq = uc $seq;
#   while ($seq =~ /$motif_regexp/g) {
# 	  $i++;
# 	  $match{$i}{'case'} = $&;
# 	  my $start = $-[0] + 1;
# 	  my $end = $+[0];
# 	  $match{$i}{'start'} = $start;
# 	  $match{$i}{'end'} = $end;
# 	  $match{$i}{'length'} = $end - $start + 1;
#   }
#   return %match;
# }


#-----------------------------------------------------------------
#----------------  Documentation / Usage / Help ------------------

=head1 NAME

vcf2model.pl - summarize SNP and INDEL mutational profile based on provided vcf file.

=head1 SYNOPSIS

perl vcf2model.pl [options] [file ...]

=head1 OPTIONS

=over 8

=item B<-help> or B<-h>

Print help message. Example: -h.

=item B<-man> or B<-m>

Print more detailed help message. Example: -m.

=item B<-version> or B<-v>

Print version information. Example: -v.

=item B<-vcf> or B<-i>

Specify the input variant calling file (in vcf or vcf.gz format). Example: -vcf input.vcf(.gz).

=item B<-excluded_chr_list> or B<-e>

Specify the chromosome(s) to be excluded for variant profiling. Example: -excluded_chr_list excluded_chr_list.txt. Default = "".

=item B<-qual> or B<-q>

Specify the cutoff of the minimal variant quality to be considered. Example: -qual 30. Default = 0.

=item B<-prefix> or B<-p>

Specify the file name prefix for the output files. Example: -prefix test_prefix. Default = "output_prefix".

=back

=head1 DESCRIPTION

B<vcf2model.pl> will summarize SNP and INDEL mutational profile based on the input vcf file. The resulting output can be used as the input for B<simuG.pl>.

=head1 AUTHOR

B<Jia-Xing Yue> (GitHub ID: yjx1217)

=head1 VERSION

B<version> v1.0.0

=cut
