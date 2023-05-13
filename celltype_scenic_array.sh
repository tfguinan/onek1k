#PBS -N onek1k_celltype_scenic
#PBS -M tfguinan@utas.edu.au
#PBS -m ab
#PBS -J 1-28:1
#PBS -l select=1:ncpus=28:mem=256gb
#PBS -j oe
#PBS -o /data/menzies_projects/onek1k/share/TG/git_analysis/logs/

QSCRIPT_VERSION='1'

echo Job ${PBS_ARRAY_INDEX}
echo Version $QSCRIPT_VERSION

# Follows R script for Seurat
# Follows R script for Seurat visualisation
# Some celltypes may require single qsub job (or higher array job resources)



# Update PATH to include pySCENIC user install
PATH=$PATH:/u/tfguinan/.local/bin
nwork=14

# Fill the path of the analysis directory
ANALYSIS_DIRECTORY="/data/menzies_projects/onek1k/share/TG/git_analysis/"
# Fill the path of the github pull
TOOL_PATH="/data/menzies_projects/onek1k/share/gitinit/"

cd $ANALYSIS_DIRECTORY

RAW_DIRECTORY="${ANALYSIS_DIRECTORY}raw_data/"
OUTPUT_DIRECTORY="${ANALYSIS_DIRECTORY}output_data/"

SCENIC_DIRECTORY=${OUTPUT_DIRECTORY}scenic/prod/

# Use per celltype h5ad as input
QUERY_DIRECTORY=${OUTPUT_DIRECTORY}seurat/celltype_split/

DATA=job${PBS_ARRAY_INDEX}-*-onek1k.h5ad
DATA_PATH=${QUERY_DIRECTORY}${DATA}


# !Rename data !IN SEPERATE TERMINAL! for use with job array index!
# Eg. job1_cell_type_onek1k.h5ad
# Eg. job1-cell_type-onek1k.h5ad

# pushd ${QUERY_DIRECTORY}
# number=0
# for filename in *.h5ad;do
# mv $filename job$(( ++number ))-${filename}
# done; popd


# Format to get celltype info without job number
POOL=`echo $DATA_PATH`
POOL=${POOL##${QUERY_DIRECTORY}job${PBS_ARRAY_INDEX}_}
POOL=${POOL%%.h5ad}


module load Python

# For reference in logs
echo Array job $PBS_ARRAY_INDEX for $DATA_PATH
echo Pool name is $POOL

# Pool is now celltype
mainloom=${SCENIC_DIRECTORY}${POOL}-filtered-scenic.loom
tflist=${RAW_DIRECTORY}hg19_tfs.txt
adjtsv=${SCENIC_DIRECTORY}${POOL}-adj.tsv

# Create loom to pySCENIC specification
# Use SCENIC_DIRECTORY for output
# TODO update post for same approach
echo python ${TOOL_PATH}pre_scenic.py $PWD $POOL $DATA_PATH $SCENIC_DIRECTORY $nwork
python ${TOOL_PATH}pre_scenic.py $PWD/ $POOL $DATA_PATH $SCENIC_DIRECTORY $nwork

# Derive co-expression modules
# Garbage collection and cpu time warning
# Use higher memory for worker excess

# Relies on pyscenic port detection and reattempt. 
# Could implement fixed port with --dask-scheduler-port

# echo pyscenic grn $mainloom $tflist -o $adjtsv --num_workers $nwork
pyscenic grn $mainloom $tflist -o $adjtsv --num_workers $nwork

ranfea=${RAW_DIRECTORY}/aux_scenic/*.mc9nr.genes_vs_motifs.rankings.feather
mottbl=${RAW_DIRECTORY}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
regcsv=${SCENIC_DIRECTORY}${POOL}-reg.csv

# Find enriched motifs
# echo pyscenic ctx $adjtsv $ranfea --annotations_fname $mottbl --expression_mtx_fname $mainloom --output $regcsv --num_workers $nwork --mask_dropouts --all_modules
pyscenic ctx $adjtsv $ranfea --annotations_fname $mottbl --expression_mtx_fname $mainloom --output $regcsv --num_workers $nwork --mask_dropouts --all_modules

# Quantify activity (independent of normalisation)
# echo pyscenic aucell $mainloom $regcsv --output ${SCENIC_DIRECTORY}${POOL}_scenic_output.loom --num_workers $nwork
pyscenic aucell $mainloom $regcsv --output ${SCENIC_DIRECTORY}${POOL}-scenic-output.loom --num_workers $nwork

python ${TOOL_PATH}post_scenic.py $PWD/ $POOL data_example $nwork