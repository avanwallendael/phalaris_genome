#!/bin/bash
# Generate Circos for IWGC Phalaris minor genome report
# Nicholas A. Johnson - Michigan State University
# Note: must not use paths in FASTA or GFF input names as this will cause confusion when naming the files since the $genome variable is used to name them

# Call variables from input arguements
genome=$1
genes=$2
repeats=$3
samtools_sif=$4
bedtools_sif=$5
quartet_sif=$6
window=${7:-300000} # (optional) bedtools window size for all tracks, defaults to 300Kbp

# Track 1: Karyotype
# Generate fasta indexes
singularity exec ${samtools_sif} samtools faidx ${genome}
grep -Ei 'Chr0?[1-9][0-9]?' ${genome}.fai | sort -t$'\t' -k1,1 > ${genome}_chrs.fai

# Generate Circos karyotype file from fai
cut -f1,2 ${genome}_chrs.fai |\
sort -t$'\t' |\
awk '{
    chr=$1;
    gsub(/^[Cc][Hh][Rr]0?/, "", chr);  # remove "chr" or "chr0", case-insensitive
    print "chr", "-", $1, chr, 0, $2, "chr"NR
}' > ${genome}_karyotype.circos

# Identify telomeres on chr ends and add them to the karyotype file
singularity exec ${quartet_sif} quartet.py te -i ${genome} -c plant -m 50 -p ${genome}

# Extract telomere info from quartet output and format for circos bands
# Define genome length first (outside the pipeline)
genome_length=$(awk '{sum += $2} END {print sum}' "${genome}_chrs.fai")

# Then run the full processing pipeline with that value passed into awk
grep -v '#' "${genome}.telo.info" | grep -Ei 'Chr0?[1-9][0-9]?' |  awk '{ print $1, $2, $3, $4, "+", $6, "-" }' | \
awk -v genome_len="$genome_length" '
{
    telolen = int(genome_len * 0.005);  # 0.5% of total genome length

    if ($3 == "both") {
        $8 = 0;
        $9 = telolen;
        print;
        $8 = $2 - telolen;
        $9 = $2;
        print;
    }
    else if ($3 == "right") {
        $8 = $2 - telolen;
        $9 = $2;
        print;
    }
    else if ($3 == "left") {
        $8 = 0;
        $9 = telolen;
        print;
    }
    else if ($3 == "no") {
        next;
    }
    else {
        exit 1;
    }
}' | \
awk '{ print "band", $1, $1"_T"NR, $1"_T"NR, $8, $9, "vdgrey" }' > "${genome}_telomere_bands.bed"

# Append telomere bands to the Circos karyotype file
cat ${genome}_telomere_bands.bed >> ${genome}_karyotype.circos



#Track 2: Gene Density
# Make genomic windows with bedtools
singularity exec ${bedtools_sif} bedtools makewindows -g ${genome}_chrs.fai -w ${window} > ${genome}_windows.bed

# Pull gene coords from annotation gff3 file
grep 'gene' ${genes} | grep -Ei 'Chr0?[1-9][0-9]?' | awk '{ print $1, $4, $5 }' OFS='\t' |\
sort -t$'\t' -k1,1 -k2,2n > ${genes}_genes_coords.bed

# Calculate coverage of genes across genome windows
singularity exec ${bedtools_sif} bedtools coverage -sorted -a ${genome}_windows.bed -b ${genes}_genes_coords.bed | awk '{ print $1, $2, $3, $7 }' | \
sort -k1,1 -k2,2n -k3,3n | \
awk '
{
  raw_log = sqrt($4);  # sqrt transform density values
  data[NR]=$1"\t"$2"\t"$3;       # store coords
  values[NR]=raw_log;            # store transformed value
  if (NR==1 || raw_log < min) min = raw_log;
  if (NR==1 || raw_log > max) max = raw_log;
}
END {
  for (i=1; i<=NR; i++) {
    norm = (max == min) ? 0 : (values[i] - min) / (max - min);
    print data[i], norm;
  }
}' > ${genome}_gene_coverage.circos



# Track 3: Repeat Density
# Generate repeat coords bed file
grep -v '#' ${repeats} | grep -Ei 'Chr0?[1-9][0-9]?' | awk '{ print $1, $4, $5 }' OFS='\t' |\
sort -t$'\t' -k1,1 -k2,2n > ${repeats}_repeat_coords.bed

# Calculate coverage of repeats across genome windows
singularity exec ${bedtools_sif} bedtools coverage -sorted -a ${genome}_windows.bed -b ${repeats}_repeat_coords.bed | awk '{ print $1, $2, $3, $7 }' | \
sort -k1,1 -k2,2n -k3,3n | \
awk '
{
  transformed = ($4^3);         # cube transform TE density, good for high value clustering
  data[NR]=$1"\t"$2"\t"$3;        # store coordinates
  values[NR]=transformed;         # store transformed values
  if (NR==1 || transformed < min) min = transformed;
  if (NR==1 || transformed > max) max = transformed;
}
END {
  for (i=1; i<=NR; i++) {
    norm = (max == min) ? 0 : (values[i] - min) / (max - min);
    print data[i], norm;
  }
}' > ${repeats}_repeat_coverage.circos



# Track 4: Command for GC windows
# Generate GC windows
singularity exec ${bedtools_sif} \
  bedtools nuc -fi ${genome} -bed ${genome}_windows.bed > gc_content.bed

awk '{ print $1, $2, $3, $5}' gc_content.bed | grep -v '#' > ${genome}_gc.circos
