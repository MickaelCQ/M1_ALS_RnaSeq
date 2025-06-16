#!/bin/bash

#SBATCH --job-name=featureCounts_SLA_exons
#SBATCH --output=logs/featureCounts_SLA_%j.out
#SBATCH --error=logs/featureCounts_SLA_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-user=mickael.coquerelle@etu.Umontpellier.fr

# ================================================================
# Auteur : Mickael Coquerelle
# Date   : 2025-06-16
# Objet  : Quantification exonique avec featureCounts (paired-end)
# Inspiré de la regle SnakeMake du laboratoire Bio2M
# ================================================================

# Répertoires et fichiers
GTF_PATH="/data/annotations/human/ensembl/Homo_sapiens.GRCh37.gtf"
OUTPUT_PATH="./all-counts_exons.tsv"
THREADS=10

# Création du répertoire pour les logs et la sortie
mkdir -p "$(dirname "$OUTPUT_PATH")"
mkdir -p logs

# Exécution de featureCounts en mode exonique (-f)
featureCounts -p --countReadPairs -T "$THREADS" -f \
-a "$GTF_PATH" -o "$OUTPUT_PATH" \
"output/STAR/bam/202402-2401121737-SLA.bam" \
"output/STAR/bam/202402-C01P026-Control.bam" \
"output/STAR/bam/202402-C01P034-Control.bam" \
"output/STAR/bam/202402-C01P040-Control.bam" \
"output/STAR/bam/202402-C01P030-SLA.bam" \
"output/STAR/bam/202402-2309211582-SLA.bam" \
"output/STAR/bam/202402-C01P027-Control.bam" \
"output/STAR/bam/202402-C01P039-SLA.bam" \
"output/STAR/bam/202505S13-2411071506-SLA.bam" \
"output/STAR/bam/202505S13-C01P035-SLA.bam" \
"output/STAR/bam/202505S13-C01P148-SLA.bam" \
"output/STAR/bam/202505S13-C01P038-SLA.bam" \
"output/STAR/bam/202505S13-C01P041-SLA.bam" \
"output/STAR/bam/202505S13-C01P022-SLA.bam" \
"output/STAR/bam/202505S13-C01P086-Control.bam" \
"output/STAR/bam/202505S13-C01P181-SLA.bam" \
"output/STAR/bam/202505S15-C01P196-SLA.bam" \
"output/STAR/bam/202505S15-C01P222-SLA.bam" \
"output/STAR/bam/202505S15-C01P238-SLA.bam" \
"output/STAR/bam/202505S15-C01P195-SLA.bam" \
"output/STAR/bam/202505S15-C01P163-SLA.bam" \
"output/STAR/bam/202505S15-C01P188-SLA.bam" \
"output/STAR/bam/202505S15-C01P260-SLA.bam" \
"output/STAR/bam/202505S15-C01P216-SLA.bam" \
"output/STAR/bam/202408-C01P088-Control.bam" \
"output/STAR/bam/202408-C01P073-SLA.bam" \
"output/STAR/bam/202408-C01P075-SLA.bam" \
"output/STAR/bam/202408-C01P094-Control.bam" \
"output/STAR/bam/202408-C01P083-SLA.bam" \
"output/STAR/bam/202408-C01P070-Control.bam" \
"output/STAR/bam/202408-C01P090-SLA.bam" \
"output/STAR/bam/202408-C01P095-Control.bam" \
"output/STAR/bam/202505S14-C01P092-SLA.bam" \
"output/STAR/bam/202505S14-C01P120-SLA.bam" \
"output/STAR/bam/202505S14-C01P179-SLA.bam" \
"output/STAR/bam/202505S14-C01P047-SLA.bam" \
"output/STAR/bam/202505S14-C01P167-SLA.bam" \
"output/STAR/bam/202505S14-C01P145-SLA.bam" \
"output/STAR/bam/202505S14-C01P051-SLA.bam" \
"output/STAR/bam/202505S14-C01P116-SLA.bam" \
"output/STAR/bam/202401-C01P013-Control.bam" \
"output/STAR/bam/202401-C01P017-Control.bam" \
"output/STAR/bam/202401-C01P024-Control.bam" \
"output/STAR/bam/202401-C01P003-Control.bam" \
"output/STAR/bam/202401-C01P014-Control.bam" \
"output/STAR/bam/202401-C01P012-Control.bam" \
"output/STAR/bam/202401-C01P016-Control.bam" \
"output/STAR/bam/202401-C01P023-Control.bam" \
"output/STAR/bam/202404-C01P069-Control.bam" \
"output/STAR/bam/202404-C01P059-Control.bam" \
"output/STAR/bam/202404-C01P053-SLA.bam" \
"output/STAR/bam/202404-C01P061-Control.bam" \
"output/STAR/bam/202404-C01P068-Control.bam" \
"output/STAR/bam/202404-C01P056-SLA.bam" \
"output/STAR/bam/202404-C01P060-Control.bam" \
"output/STAR/bam/202404-C01P057-SLA.bam" \
"output/STAR/bam/202409-C01P089-SLA.bam" \
"output/STAR/bam/202409-C01P128-SLA.bam" \
"output/STAR/bam/202409-C01P006-SLA.bam" \
"output/STAR/bam/202409-C01P035-SLA.bam" \
"output/STAR/bam/202409-C01P082-SLA.bam" \
"output/STAR/bam/202409-C01P154-SLA.bam" \
"output/STAR/bam/202409-C01P135-SLA.bam" \
"output/STAR/bam/202409-C01P138-SLA.bam" \
"output/STAR/bam/202501-C01P138-SLA.bam"

echo "featureCounts à l'échelle de l'exon. "
