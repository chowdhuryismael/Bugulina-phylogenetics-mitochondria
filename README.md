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
BEGIN {outfile="genes/aa/by_gene/" gene ".faa"}
/^>/ {
    # get the part after last semicolon and remove spaces
    split($0, a, ";")
    g = a[length(a)]
    gsub(/^[ \t]+|[ \t]+$/, "", g)   # trim spaces
    keep = (g == gene)
}
keep {
    if (/^>/) print ">" sample "|" $0 >> outfile
    else print >> outfile
}
' "$faa_file"
	done
    else
        echo "No result.faa in $sample_dir"
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




