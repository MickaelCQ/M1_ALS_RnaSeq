#!/bin/bash

############################################################################################################################################
# Auteur   : Mickael Coquerelle
# Date     : 29/04/2025
# Projet   : Stats BAM Log --  Version : 1.3.0
# Licence : CC BY 4.0 (Creative Commons Attribution 4.0 International)
#
# Ce script est libre d’utilisation, de modification et de redistribution
# à condition de créditer l’auteur original.
# Pour en savoir plus : https://creativecommons.org/licenses/by/4.0/

#############
# Objectifs :
#############
 
# La vocation de ce script est d'extraire les informations du Log final pour chaque Bam pour le répertoire courant/ss-dossiers,pour CRAC et STAR, afin d'effectuer de la statistique descriptive ensuite à l'aide de R. 
# En outre, une telle extraction nous permettra d'évaluer :
# - La qualité de l'échantillon (un ARN de mauvaise qualité).
# - D'éventuels effets de batch (qui peuvent s'identifier avec de mauvais taux de mapping).
# - Nous aider à d'éventuelles inclusions/exclusions lors d'une analyse d'expression différentielle ultérieure.
# - Répondre à l'hypothèse de dispersion des données entre nos runs RNASeq en corrélany les informations d'alignements et les données expérimentales.
# - Permet d'évaluer la reproductibilité des run RNASeq (pour notre cas d'étude). 

#############################################################################################################################################
# OUTPUT Log Final Star type , on utilisera chaque label avec grep et le "|" pour cut sur -f2 , sur chaque ligne la stratégie reste la même : 
#############################################################################################################################################

#                             Started job on |       Apr 28 10:35:
#                         Started mapping on |       Apr 28 10:38:
#                                Finished on |       Apr 28 10:39:
#   Mapping speed, Million of reads per hour |       84.82
#
#                      Number of input reads |       801045
#                  Average input read length |       302
# .....


# Demande du répertoire des logs STAR
#echo "Entrez le répertoire contenant les logs STAR à traiter (Log.final.out) au format : ~/Dir1/Dir2/..."
#read STAR_DIR

# Demande du répertoire des fichiers summary CRAC
#echo "Entrez le répertoire contenant les fichiers summary CRAC (.summary) : ~/Dir1/Dir2/...."
#read CRAC_DIR

STAR_DIR="/home/mickael/M1_Stage/M1_ALS_RnaSeq/4.LOG_STAR_01052025_HG37";
CRAC_DIR="/home/mickael/M1_Stage/M1_ALS_RnaSeq/output/crac/summary";


###########################################################	
# ETAPE 1 : Extraction des fichiers de Log finaux STAR    : 
###########################################################

SBL="Stats_Log_star.csv";
SBC="Stats_Log_crac.csv";
SBF="Stats_Log_merge.csv";


# Nettoyage des anciens LOGS
if [ -f "$SBL" ] || [ -f "$SBC" ] || [ -f "$SBF" ] ; then
    rm "$SBL" "$SBC" "$SBF"
fi
    
echo STAR_Run,STAR_Patient,STAR_Type,STAR_Date_Mapping,STAR_Total_reads,\
STAR_Unique_reads,STAR_Unique_pct,\
STAR_Multi_reads,STAR_Multi_pct,\
STAR_No_map_reads,STAR_No_map_pct_sum,STAR_No_map_pct_mismatch,STAR_No_map_pct_tooshort,STAR_No_map_pct_other,\
STAR_Avg_read_len,STAR_Avg_read_map_len,\
STAR_Splices_total,STAR_Splices_GTAG,STAR_Splices_GCAG,STAR_Splices_ATAC,STAR_Splices_Noncanonical,\
STAR_Mismatch_rate,STAR_Deletion_rate,STAR_Insertion_rate,\
STAR_Chimeric_reads,STAR_Chimeric_reads_pct > "$SBL"

echo CRAC_Sample,CRAC_Date_Mapping,CRAC_Total_reads,\
CRAC_Unique_reads,CRAC_Unique_pct,\
CRAC_Multi_reads,CRAC_Multi_pct,\
CRAC_No_map_reads,CRAC_No_map_pct,CRAC_Dup_reads,CRAC_Dup_pct,CRAC_Explainable_reads,CRAC_Explainable_pct,\
CRAC_Rep_reads,CRAC_Rep_pct,CRAC_Normal_reads,CRAC_Normal_pct,CRAC_AlmostNormal_reads,CRAC_AlmostNormal_pct,\
CRAC_SeqErr_reads,CRAC_SeqErr_pct,CRAC_SNV_reads,CRAC_SNV_pct,CRAC_ShortIndel_reads,CRAC_ShortIndel_pct,\
CRAC_Splice_reads,CRAC_Splice_pct,CRAC_WeakSplice_reads,CRAC_WeakSplice_pct,\
CRAC_Chimera_reads,CRAC_Chimera_pct,CRAC_PE_Chimera_reads,CRAC_PE_Chimera_pct,\
CRAC_BioUndetermined_reads,CRAC_BioUndetermined_pct,CRAC_Undetermined_reads,CRAC_Undetermined_pct > "$SBC"


# Fonctions d'extractions :
Ext_STAR() { cut -d '|' -f2 | tr -d ' '; }
Ext_absolue_CRAC()   { cut -d ':' -f2 | cut -d '(' -f1 | tr -d ' ' | head -n1; }
Ext_relative_CRAC()  { cut -d ':' -f2 | cut -d '(' -f2 | tr -d ' %)' | head -n1; }

# Pour chaque fichier du répertoire courant et des sous-répertoires :
for STAR_bam in "$STAR_DIR"/*.Log.final.out;
  do
    echo "Traitement de : $STAR_bam"
    
    #  1. ID :
    STAR_BaseName=$(basename "$STAR_bam")
    STAR_Sample=$(echo "$STAR_BaseName" | sed -E "s/.*Align_(.*)_Log.final.out/\1/")
    STAR_Run=$(echo "$STAR_Sample" | cut -d '-' -f1)
    STAR_Patient=$(echo "$STAR_Sample" | cut -d '-' -f2)
    STAR_Type=$(echo "$STAR_Sample" | cut -d '-' -f3 | cut -d '.' -f1)

    STAR_Date_Mapping=$(grep "Finished on" "$STAR_bam" | Ext_STAR);

    # 2. UNIQUE READS 
    STAR_Total_reads=$(grep "Number of input reads" "$STAR_bam" | Ext_STAR);	
    STAR_Unique_reads=$(grep "Uniquely mapped reads number" "$STAR_bam" | Ext_STAR);
    STAR_Unique_pct=$(grep "Uniquely mapped reads %" "$STAR_bam" | Ext_STAR | tr -d '%');

    # 3. MULTI-MAPPING READS 
    STAR_Multi_reads=$(grep "Number of reads mapped to multiple loci" "$STAR_bam" | Ext_STAR);
    STAR_Multi_pct=$(grep "% of reads mapped to multiple loci" "$STAR_bam" | Ext_STAR);

    # 4. UNMAPPED READS 
    STAR_No_map_pct_mismatch=$(grep "% of reads unmapped: too many mismatches" "$STAR_bam" | Ext_STAR | tr -d '%' || echo "0");
    STAR_No_map_pct_tooshort=$(grep "% of reads unmapped: too short" "$STAR_bam" | Ext_STAR | tr -d '%' || echo "0");
    STAR_No_map_pct_other=$(grep "% of reads unmapped: other" "$STAR_bam" | Ext_STAR | tr -d '%' || echo "0");
    STAR_No_map_pct_sum=$(echo "scale=2; $STAR_No_map_pct_mismatch + $STAR_No_map_pct_tooshort + $STAR_No_map_pct_other" | bc);

    STAR_No_map_reads=$(( 
  	$(grep "Number of reads unmapped: too many mismatches" "$STAR_bam" | Ext_STAR || echo 0) + 
  	$(grep "Number of reads unmapped: too short" "$STAR_bam" | Ext_STAR || echo 0) + 
  	$(grep "Number of reads unmapped: other" "$STAR_bam" | Ext_STAR || echo 0)
	));
   
    # 5. AVERAGES :
    STAR_Avg_read_len=$(grep "Average input read length" "$STAR_bam" | Ext_STAR); #Longueur moyenne des reads en entrée
    STAR_Avg_read_map_len=$(grep "Average mapped length" "$STAR_bam" | Ext_STAR); #Longueur moyenne des segments après alignement 
    	
    # 6. SPLICES :
    STAR_Splices_total=$(grep "Number of splices: Total" "$STAR_bam" | Ext_STAR);
    STAR_Splices_GTAG=$(grep -F "Number of splices: GT" "$STAR_bam" | Ext_STAR || echo "NA");
    STAR_Splices_GCAG=$(grep -F "Number of splices: GC" "$STAR_bam" | Ext_STAR || echo "NA");
    STAR_Splices_ATAC=$(grep -F "Number of splices: AT" "$STAR_bam" | Ext_STAR || echo "NA");
    STAR_Splices_Noncanonical=$(grep "Number of splices: Non-canonical" "$STAR_bam" | Ext_STAR); 
    
    #7. SUBSTITUTIONS :
    STAR_Mismatch_rate=$(grep "Mismatch rate per base" "$STAR_bam" | Ext_STAR | tr -d '%');
    STAR_Deletion_rate=$(grep "Deletion rate per base" "$STAR_bam" | Ext_STAR);
    STAR_Insertion_rate=$(grep "Insertion rate per base" "$STAR_bam" | Ext_STAR);

    #8. CHIMERICS :
    STAR_Chimeric_reads=$(grep "Number of chimeric reads" "$STAR_bam" | Ext_STAR);
    STAR_Chimeric_reads_pct=$(grep "% of chimeric reads" "$STAR_bam" | Ext_STAR | tr -d '%');
  
    #9. Écriture dans le tabulé :
    echo "${STAR_Run},${STAR_Patient},${STAR_Type},${STAR_Date_Mapping},${STAR_Total_reads},${STAR_Unique_reads},${STAR_Unique_pct},${STAR_Multi_reads},${STAR_Multi_pct},${STAR_No_map_reads},${STAR_No_map_pct_sum},${STAR_No_map_pct_mismatch},${STAR_No_map_pct_tooshort},${STAR_No_map_pct_other},${STAR_Avg_read_len},${STAR_Avg_read_map_len},${STAR_Splices_total},${STAR_Splices_GTAG},${STAR_Splices_GCAG},${STAR_Splices_ATAC},${STAR_Splices_Noncanonical},${STAR_Mismatch_rate},${STAR_Deletion_rate},${STAR_Insertion_rate},${STAR_Chimeric_reads},${STAR_Chimeric_reads_pct}" >> "$SBL";

done

########################################################	
# ETAPE 2 : Extraction des fichiers summary de CRAC    : 
########################################################

for CRAC_summary in "$CRAC_DIR"/*.summary; 
  do
    echo "Traitement de : $CRAC_summary"

    # 1. ID :
    CRAC_Echantillon=$(basename "$CRAC_summary" .summary)
    CRAC_Date_Mapping=$(date "+%Y-%m-%d") 
    
    # 2. UNIQUE READS :
    CRAC_Total_reads=$(grep "Total number of reads analyzed" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Unique_reads=$(grep "Single" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Unique_pct=$(grep "Single" "$CRAC_summary" | Ext_relative_CRAC)

    # 3. MULTI-MAPPING READS : 
    CRAC_Multi_reads=$(grep "Multiple" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Multi_pct=$(grep "Multiple" "$CRAC_summary" | Ext_relative_CRAC)

    # 4. UNMAPPED READS :
    CRAC_No_map_reads=$(grep "None" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_No_map_pct=$(grep "None" "$CRAC_summary" | Ext_relative_CRAC)

    # 5. OTHER METRICS :
    CRAC_Dup_reads=$(grep "Duplication" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Dup_pct=$(grep "Duplication" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Explainable_reads=$(grep "Explainable" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Explainable_pct=$(grep "Explainable" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Rep_reads=$(grep "Repetition" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Rep_pct=$(grep "Repetition" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Normal_reads=$(grep "Normal" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Normal_pct=$(grep "Normal" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_AlmostNormal_reads=$(grep -F "Almost-Normal" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_AlmostNormal_pct=$(grep -F "Almost-Normal" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_SeqErr_reads=$(grep -F "Sequence-Errors" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_SeqErr_pct=$(grep -F "Sequence-Errors" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_SNV_reads=$(grep "SNV" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_SNV_pct=$(grep "SNV" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_ShortIndel_reads=$(grep -F "Short-Indel" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_ShortIndel_pct=$(grep -F "Short-Indel" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Splice_reads=$(grep "Splice" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Splice_pct=$(grep  "Splice" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_WeakSplice_reads=$(grep -F "Weak-Splice" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_WeakSplice_pct=$(grep "Weak-Splice" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Chimera_reads=$(grep "Chimera" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Chimera_pct=$(grep "Chimera" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_PE_Chimera_reads=$(grep "Paired-end Chimera" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_PE_Chimera_pct=$(grep "Paired-end Chimera" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_BioUndetermined_reads=$(grep "Bio-Undetermined" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_BioUndetermined_pct=$(grep "Bio-Undetermined" "$CRAC_summary" | Ext_relative_CRAC)

    CRAC_Undetermined_reads=$(grep "Undetermined" "$CRAC_summary" | Ext_absolue_CRAC)
    CRAC_Undetermined_pct=$(grep "Undetermined" "$CRAC_summary" | Ext_relative_CRAC)

    # READ LINE CSV: 
    echo "${CRAC_Echantillon},${CRAC_Date_Mapping},${CRAC_Total_reads},${CRAC_Unique_reads},${CRAC_Unique_pct},${CRAC_Multi_reads},${CRAC_Multi_pct},${CRAC_No_map_reads},${CRAC_No_map_pct},${CRAC_Dup_reads},${CRAC_Dup_pct},${CRAC_Explainable_reads},${CRAC_Explainable_pct},${CRAC_Rep_reads},${CRAC_Rep_pct},${CRAC_Normal_reads},${CRAC_Normal_pct},${CRAC_AlmostNormal_reads},${CRAC_AlmostNormal_pct},${CRAC_SeqErr_reads},${CRAC_SeqErr_pct},${CRAC_SNV_reads},${CRAC_SNV_pct},${CRAC_ShortIndel_reads},${CRAC_ShortIndel_pct},${CRAC_Splice_reads},${CRAC_Splice_pct},${CRAC_WeakSplice_reads},${CRAC_WeakSplice_pct},${CRAC_Chimera_reads},${CRAC_Chimera_pct},${CRAC_PE_Chimera_reads},${CRAC_PE_Chimera_pct},${CRAC_BioUndetermined_reads},${CRAC_BioUndetermined_pct},${CRAC_Undetermined_reads},${CRAC_Undetermined_pct}" >> "$SBC"
  done


########################################################	
# ETAPE 3 : Fusion des données de mapping STAR/CRAC    : 
########################################################

    # Stockage dans un tableau global des métriques d'alignements : 
if [[ -f "$SBL" && -f "$SBC" ]]; then
    echo "Fusion des résultats STAR et CRAC en un seul fichier : $SBF"

if [[ $(( $(wc -l < "$SBL") - 1 )) -ne $(( $(wc -l < "$SBC") - 1 )) ]]; then
    echo "Attention $SBL et $SBC ont un nombre de lignes différents."
fi
    # Fusion des en-têtes
    paste -d ',' <(head -n 1 "$SBL") <(head -n 1 "$SBC") > "$SBF"
    # Fusion des contenus
    paste -d ',' <(tail -n +2 "$SBL") <(tail -n +2 "$SBC") >> "$SBF"

    echo "Fusion terminée. Résultat final : $SBF";
else
    echo "Erreur : l'un des fichiers $SBL ou $SBC est manquant.";
fi
