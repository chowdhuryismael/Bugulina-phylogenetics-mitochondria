# Bugulina-phylogenetics-mitochondria


#!/usr/bin/env fish

# Activate Conda for Fish
source /projects/health_sciences/bms/biochemistry/kenny_group/Ismael_chow/mamba/etc/fish/conf.d/conda.fish

# Activate MITOS environment
conda activate mitos_env

# Set the reference directory for MITOS
set refdir /projects/health_sciences/bms/biochemistry/kenny_group/Ismael_chow/X201SC26030018-Z01-F001/01.RawData/mitos_refseq/refseq89m

# Base output directory
set base_outdir ./mitos_outputs

# Create the base output folder if it doesn't exist
if not test -d $base_outdir
    mkdir $base_outdir
end

# Loop over all fasta files
for f in *.fasta
    # Get the base name without extension
    set base (basename $f .fasta)
    
    # Create output directory for this file
    set outdir $base_outdir/$base
    if not test -d $outdir
        mkdir $outdir
    end

    echo "Running MITOS on $f ..."

    # Run MITOS with circular genome option
    $CONDA_PREFIX/bin/runmitos.py \
        -i $f \
        -c 5 \
        -o $outdir \
        --circrot 0 \
        -r $refdir
end




aa gene extractions

!/bin/bash

# Make sure output folder exists
mkdir -p genes/aa/by_gene

# Define the 13 protein-coding genes
genes=("atp6" "atp8" "cox1" "cox2" "cox3" "cob" "nad1" "nad2" "nad3" "nad4" "nad4l" "nad5" "nad6")

# Loop through each sample folder
for sample_dir in */; do
    sample=$(basename "$sample_dir")
    faa_file="${sample_dir}result.faa"

    if [[ -f "$faa_file" ]]; then
        echo "Processing $sample"

        # Extract each gene
        for gene in "${genes[@]}"; do
            awk -v sample="$sample" -v gene="$gene" '
#!/bin/bash

mkdir -p genes/nuc/by_gene

genes=("atp6" "atp8" "cox1" "cox2" "cox3" "cob" "nad1" "nad2" "nad3" "nad4" "nad4l" "nad5" "nad6")

for sample_dir in */; do
    sample=$(basename "$sample_dir")
    
    # Use result.fas (nucleotides), not result.faa
    nuc_file="${sample_dir}result.fas"
    
    if [[ -f "$nuc_file" ]]; then
        echo "Extracting nucleotides from $sample"
        
        for gene in "${genes[@]}"; do
            awk -v sample="$sample" -v gene="$gene" '
BEGIN {outfile="genes/nuc/by_gene/" gene ".fasta"}
/^>/ {
    # Get the gene name from the end of the header
    match($0, /; ([^;]+)$/)
    g = tolower(substr($0, RSTART+2))
    keep = (g == gene)
}
keep {
    if (/^>/) {
        # Create new header with sample name
        print ">" sample "|" substr($0, 2) >> outfile
    } else {
        print >> outfile
    }
}
' "$nuc_file"
        done
    else
        echo "No result.fas in $sample_dir"
    fi
done


## gene check
#!/bin/bash

echo -e "Gene\tSample\tPresent" > gene_presence_absence.tsv

# Loop over each gene file
for gene_file in genes/aa/by_gene/*.faa; do
    gene=$(basename "$gene_file" .faa)
    
    # Extract the sample names from headers
    awk -v gene="$gene" '/^>/ {split($0,a,"|"); print gene"\t"a[1]"\t1"}' "$gene_file" >> gene_presence_absence.tsv
done

##


/bin/bash

# Activate your mamba environment
source /projects/health_sciences/bms/biochemistry/kenny_group/Ismael_chow/mamba/bin/activate bugulidae_analysis

# Loop through each gene in current directory and align
for gene in *.fas; do
    echo "Aligning $gene ..."
    mafft --auto "$gene" > "${gene%.fas}_aligned.fas"
done

echo "All genes aligned!"

#!/bin/bash

# Activate your mamba environment

source /projects/health_sciences/bms/biochemistry/kenny_group/Ismael_chow/mamba/bin/activate bugulidae_analysis
for aln in *_aligned.fasta; do
    trimal -in "$aln" -out "${aln%.fasta}_trimmed.fasta" -automated1
done

#!bin/fish

# For each clean file, we need to standardize headers
for file in *_clean.fasta
    set base (basename $file .fasta)
    # Create a temporary file with standardized headers
    awk 'BEGIN {OFS="\t"} 
         NR==FNR {headers[$0]=1; next} 
         /^>/ {
            # Extract base name
            gsub(/^>/, "", $0)
            gsub(/,_complete_genome_KM983335/, "", $0)
            gsub(/,_AW1141/, "", $0)
            gsub(/\|.*/, "", $0)
            if ($0 in headers) {
                print ">" $0
            } else {
                print ">" $0
            }
         }
         !/^>/ {print}' \
         master_headers.txt "$file" > "${base}_standard.fasta"
end








./AMAS.py concat -f fasta -d dna -t final_supermatrix_amas.fasta -p partition_amas.txt -u fasta -i atp6_aligned_trimmed_clean.fasta atp8_aligned_trimmed_clean.fasta cox1_aligned_trimmed_clean.fasta cox2_aligned_trimmed_clean.fasta cox3_aligned_trimmed_clean.fasta cob_aligned_trimmed_clean.fasta nad1_aligned_trimmed_clean.fasta nad2_aligned_trimmed_clean.fasta nad3_aligned_trimmed_clean.fasta nad4_aligned_trimmed_clean.fasta nad4l_aligned_trimmed_clean.fasta nad5_aligned_trimmed_clean.fasta nad6_aligned_trimmed_clean.fasta rrnL_aligned_1_clean.fasta rrnS_aligned_1_clean.fasta


iqtree2 -s final_supermatrix_amas.fasta -p partition_iqtree.txt -m MFP+MERGE -B 1000 -alrt 1000 -nt 32

Partition file is not in NEXUS format, assuming RAxML-style partition file...
Subset	Type	Seqs	Sites	Infor	Invar	Model	Name
1	DNA	35	687	519	103	GTR+MERGE	atp6
2	DNA	31	88	58	9	GTR+MERGE	atp8
3	DNA	35	1500	760	618	GTR+MERGE	cox1
4	DNA	33	670	446	147	GTR+MERGE	cox2
5	DNA	35	777	462	249	GTR+MERGE	cox3
6	DNA	34	1113	679	313	GTR+MERGE	cob
7	DNA	33	923	629	207	GTR+MERGE	nad1
8	DNA	32	918	723	96	GTR+MERGE	nad2
9	DNA	35	353	245	74	GTR+MERGE	nad3
10	DNA	34	1038	735	184	GTR+MERGE	nad4
11	DNA	33	255	191	33	GTR+MERGE	nad4l
12	DNA	34	1671	1227	269	GTR+MERGE	nad5
13	DNA	35	349	274	28	GTR+MERGE	nad6
14	DNA	34	1237	778	287	GTR+MERGE	rrnL
15	DNA	34	828	539	161	GTR+MERGE	rrnS







<img width="1188" height="923" alt="bugulina_nz_tree_version7phylogram2_final6" src="https://github.com/user-attachments/assets/0cf0511a-8af1-4d4a-9a7d-c8e9e824567e" />
