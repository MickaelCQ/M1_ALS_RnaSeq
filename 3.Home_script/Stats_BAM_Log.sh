#!/bin/bash

############################################################################################################################################
# Auteur   : Mickael Coquerelle
# Date     : 29/04/2025
# Projet   : SBL : Stats BAM Log --  Version : 1.2.0
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
# OUTPUT Log Final Star type , on utilisera chaque label avec grep et le "|" pour cut sur -f2 sur chaque ligne, la stratégie reste la même : 
#############################################################################################################################################

#                             Started job on |       Apr 28 10:35:
#                         Started mapping on |       Apr 28 10:38:
#                                Finished on |       Apr 28 10:39:
#   Mapping speed, Million of reads per hour |       84.82
#
#                      Number of input reads |       801045
#                  Average input read length |       302
# .....


##############################	
# Extraction partie STAR    : 
##############################

SBL="Stats_Log_Final.csv";

# Nettoyage de l'ancien log : 
if [ -f "$SBL" ]; then
    rm "$SBL"
fi
    
    echo Run,Patient,Type,Date_Mapping,Total_reads,\
    Unique_reads,Unique_pct,\
    Multi_reads,Multi_pct,\
    No_map_reads,No_map_pct,No_map_miss,No_map_short,No_map_other,\
    Avg_read_len,Avg_read_map_len,\
    Splices_total,SplicesGTAG,SplicesGCAG,SplicesATAC,Splices_Noncanonical,\
    MisMatch_rate,Deletion_rate,Insertion_rate,\
    Chimeric_reads,Chimeric_read_pct > "$SBL";
    

# Fonction pour extraire le 2ᵉ champ et enlever les espaces
Extract() { cut -d '|' -f2 | tr -d ' '; }

# Pour chaque fichier du répertoire courant et des sous-répertoires :
for bam in *.Log.final.out;
do
    #  1. ID :
    
    BaseName=$(basename "$bam")
    Sample=$(echo "$BaseName" | sed -E "s/.*Align_(.*)_Log.final.out/\1/")
    Run=$(echo "$Sample" | cut -d '-' -f1)
    Patient=$(echo "$Sample" | cut -d '-' -f2)
    Type=$(echo "$Sample" | cut -d '-' -f3 | cut -d '.' -f1)

    Date_Mapping=$(grep "Finished on" "$bam" | Extract);

    # 2. UNIQUE READS 
    
    Total_reads=$(grep "Number of input reads" "$bam" | Extract);	
    Unique_reads=$(grep "Uniquely mapped reads number" "$bam" | Extract);
    Unique_pct=$(grep "Uniquely mapped reads %" "$bam" | Extract | tr -d '%');

    # 3. MULTI-MAPPING READS 
    
    Multi_reads=$(grep "Number of reads mapped to multiple loci" "$bam" | Extract);
    Multi_pct=$(grep "% of reads mapped to multiple loci" "$bam" | Extract);

    # 4. UNMAPPED READS 
    
    # 4.1 Extraction des valeurs numériques 
    No_map_pct_mismatch=$(grep "% of reads unmapped: too many mismatches" "$bam" | Extract | tr -d '%' || echo "0");
    No_map_pct_tooshort=$(grep "% of reads unmapped: too short" "$bam" | Extract | tr -d '%' || echo "0");
    No_map_pct_other=$(grep "% of reads unmapped: other" "$bam" | Extract | tr -d '%' || echo "0");

	# 4.2 Somme des pourcentages avec bc (résultat avec 2 décimales)
    No_map_pct_sum=$(echo "scale=2; $No_map_pct_mismatch + $No_map_pct_tooshort + $No_map_pct_other" | bc);

	# 4.3 Nombre total de lectures non mappées
    No_map_reads=$(( 
  	$(grep "Number of reads unmapped: too many mismatches" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: too short" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: other" "$bam" | Extract || echo 0)
	));
   
    # 5. AVERAGES :
    Avg_read_len=$(grep "Average input read length" "$bam" | Extract); #Longueur moyenne des reads en entrée : on va vérifier la cohérence de taille entre les runs grace à cette metrique (entre runs)
    Avg_read_map_len=$(grep "Average mapped length" "$bam" | Extract); #Longueur moyenne des segments après alignement 
    	
    	#Memo:
    	#Des reads courts ou des longueurs moyennes faibles → suspicion de dégradation d’ARN (ex. tissus anciens, conservation défaillante).
	#Un mapped length bien plus faible que le read length → suspicion de soft clipping, ARN fragmentés ou contamination.

    # 6. SPLICES :
    Splices_total=$(grep "Number of splices: Total" "$bam" | Extract);
    Splices_GTAG=$(grep -F "Number of splices: GT" "$bam" | Extract || echo "NA");
    Splices_GCAG=$(grep -F "Number of splices: GC" "$bam" | Extract || echo "NA");
    Splices_ATAC=$(grep -F "Number of splices: AT" "$bam" | Extract || echo "NA");
    Splices_Noncanonical=$(grep "Number of splices: Non-canonical" "$bam" | Extract); 
    
    #7. SUBSTITUTIONS : Pourcentage de bases dans les lectures alignées qui ne correspondent pas à la séquence de référence :
    Mismatch_rate=$(grep "Mismatch rate per base" "$bam" | Extract | tr -d '%');
    Deletion_rate=$(grep "Deletion rate per base" "$bam" | Extract);
    Insertion_rate=$(grep "Insertion rate per base" "$bam" | Extract);

    #8. CHIMERICS : Plutot utile en somatique cela, mais je l'extrait comme cela on à tout à l'éxécution du .sh, (pour identifier des transcrits de fusion, transloc ?)
    Chimeric_reads=$(grep "Number of chimeric reads" "$bam" | Extract);
    Chimeric_reads_pct=$(grep "% of chimeric reads" "$bam" | Extract | tr -d '%');
  
    #9. Écriture dans le tabulé :
    echo "${Run},${Patient},${Type},${Date_Mapping},${Total_reads},${Unique_reads},${Unique_pct},${Multi_reads},${Multi_pct},${No_map_reads},${No_map_pct_sum},${No_map_pct_mismatch},${No_map_pct_tooshort},${No_map_pct_other},${Avg_read_len},${Avg_read_map_len},${Splices_total},${Splices_GTAG},${Splices_GCAG},${Splices_ATAC},${Splices_Noncanonical},${Mismatch_rate},${Deletion_rate},${Insertion_rate},${Chimeric_reads},${Chimeric_reads_pct}" >> "$SBL";

done

##############################	
# Extraction partie CRAC    : 
##############################
SBC="Stats_Log_crac.csv"

# Suppression de l'ancien fichier si existant
[ -f "$SBC" ] && rm "$SBC"

# Écriture de l’en-tête CSV
echo "Echantillon,Date_Mapping,Total_reads,Unique_reads,Unique_pct,Multi_reads,Multi_pct,No_map_reads,No_map_pct,Dup_reads,Dup_pct,Explainable_reads,Explainable_pct,Rep_reads,Rep_pct,Normal_reads,Normal_pct,AlmostNormal_reads,AlmostNormal_pct,SeqErr_reads,SeqErr_pct,SNV_reads,SNV_pct,ShortIndel_reads,ShortIndel_pct,Splice_reads,Splice_pct,WeakSplice_reads,WeakSplice_pct,Chimera_reads,Chimera_pct,PE_Chimera_reads,PE_Chimera_pct,BioUndetermined_reads,BioUndetermined_pct,Undetermined_reads,Undetermined_pct" > "$SBC"

# ========================
# Fonctions d'extraction
# ========================
Extract_absolue()   { cut -d ':' -f2 | cut -d '(' -f1 | tr -d ' ' | head -n1; }
Extract_relative()  { cut -d ':' -f2 | cut -d '(' -f2 | tr -d ' %%)' | head -n1; }

# ========================
# Traitement des fichiers
# ========================
for summary in $(find . -type f -name "*.summary"); do
    echo "Traitement de : $summary"

    # ID & date
    Echantillon=$(basename "$summary" .summary)
    Date_Mapping=$(date "+%Y-%m-%d") 

    # Lecture des métriques principales
    Total_reads=$(grep "Total number of reads analyzed" "$summary" | Extract_absolue)

    Unique_reads=$(grep "Single" "$summary" | Extract_absolue)
    Unique_pct=$(grep "Single" "$summary" | Extract_relative)

    Multi_reads=$(grep "Multiple" "$summary" | Extract_absolue)
    Multi_pct=$(grep "Multiple" "$summary" | Extract_relative)

    No_map_reads=$(grep "None" "$summary" | Extract_absolue)
    No_map_pct=$(grep "None" "$summary" | Extract_relative)

    Dup_reads=$(grep "Duplication" "$summary" | Extract_absolue)
    Dup_pct=$(grep "Duplication" "$summary" | Extract_relative)

    Explainable_reads=$(grep "Explainable" "$summary" | Extract_absolue)
    Explainable_pct=$(grep "Explainable" "$summary" | Extract_relative)

    Rep_reads=$(grep "Repetition" "$summary" | Extract_absolue)
    Rep_pct=$(grep "Repetition" "$summary" | Extract_relative)

    Normal_reads=$(grep "Normal" "$summary" | Extract_absolue)
    Normal_pct=$(grep "Normal" "$summary" | Extract_relative)

    AlmostNormal_reads=$(grep -F "Almost-Normal" "$summary" | Extract_absolue)
    AlmostNormal_pct=$(grep -F "Almost-Normal" "$summary" | Extract_relative)

    SeqErr_reads=$(grep -F "Sequence-Errors" "$summary" | Extract_absolue)
    SeqErr_pct=$(grep -F "Sequence-Errors" "$summary" | Extract_relative)

    SNV_reads=$(grep "SNV" "$summary" | Extract_absolue)
    SNV_pct=$(grep "SNV" "$summary" | Extract_relative)

    ShortIndel_reads=$(grep -F "Short-Indel" "$summary" | Extract_absolue)
    ShortIndel_pct=$(grep -F "Short-Indel" "$summary" | Extract_relative)

    Splice_reads=$(grep "Splice" "$summary" | Extract_absolue)
    Splice_pct=$(grep  "Splice" "$summary" | Extract_relative)

    WeakSplice_reads=$(grep -F "Weak-Splice" "$summary" | Extract_absolue)
    WeakSplice_pct=$(grep "Weak-Splice" "$summary" | Extract_relative)

    Chimera_reads=$(grep "Chimera" "$summary" | Extract_absolue)
    Chimera_pct=$(grep "Chimera" "$summary" | Extract_relative)

    PE_Chimera_reads=$(grep "Paired-end Chimera" "$summary" | Extract_absolue)
    PE_Chimera_pct=$(grep "Paired-end Chimera" "$summary" | Extract_relative)

    BioUndetermined_reads=$(grep "Bio-Undetermined" "$summary" | Extract_absolue)
    BioUndetermined_pct=$(grep "Bio-Undetermined" "$summary" | Extract_relative)

    Undetermined_reads=$(grep "Undetermined" "$summary" | Extract_absolue)
    Undetermined_pct=$(grep "Undetermined" "$summary" | Extract_relative)

    # Écriture de la ligne CSV
    echo "${Echantillon},${Date_Mapping},${Total_reads},${Unique_reads},${Unique_pct},${Multi_reads},${Multi_pct},${No_map_reads},${No_map_pct},${Dup_reads},${Dup_pct},${Explainable_reads},${Explainable_pct},${Rep_reads},${Rep_pct},${Normal_reads},${Normal_pct},${AlmostNormal_reads},${AlmostNormal_pct},${SeqErr_reads},${SeqErr_pct},${SNV_reads},${SNV_pct},${ShortIndel_reads},${ShortIndel_pct},${Splice_reads},${Splice_pct},${WeakSplice_reads},${WeakSplice_pct},${Chimera_reads},${Chimera_pct},${PE_Chimera_reads},${PE_Chimera_pct},${BioUndetermined_reads},${BioUndetermined_pct},${Undetermined_reads},${Undetermined_pct}" >> "$SBC"

done

