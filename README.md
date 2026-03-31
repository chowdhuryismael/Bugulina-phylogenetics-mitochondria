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
