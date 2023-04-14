#PBS -N seurat_data
#PBS -M tfguinan@utas.edu.au
#PBS -m abe
#PBS -l select=1:ncpus=28:mem=256gb
#PBS -l walltime=168:00:00
#PBS -o /data/menzies_projects/onek1k/share/TG/git_analysis/logs
#PBS -j oe
##### /usr/local/bin/jobtail.sh -j xxxxxx
date
cd /data/menzies_projects/onek1k/share/TG/git_analysis/
module load R/4.1.2-foss-2021b
export R_LIBS="/data/menzies_projects/onek1k/share/installs/packages"
#
name=onek1k
h5seurat=/data/menzies_projects/onek1k/share/share_raw/onek1k_seurat_v210819.h5seurat
pairs=/data/menzies_projects/onek1k/share/share_raw/full_pairs.csv
tobid=/data/menzies_projects/onek1k/share/share_raw/OneK1K_sample_IDs_210304.tsv
#
Rscript /data/menzies_projects/onek1k/share/gitinit/seurat_data.R $PWD $name $h5seurat $pairs $tobid
date
# dmesg -l info,notice,warn,err,crit,alert,emerg > /scratch/tfguinan/misc_logs/dmesg.log
