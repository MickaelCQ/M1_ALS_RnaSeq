#!/bin/bash

################################################################################## Mickael -- le 29/05/2025 - ####################################################################################
#                                                                                   SBL : Stats BAM Log                                                                                          #
##################################################################################################################################################################################################

#############
# Objectifs :
#############
 
# La vocation de ce script est d'extraire les informations du Log final pour chaque Bam pour le répertoire courant/ss-dossiers, afin d'effectuer de la statistique descriptive ensuite à l'aide de R. 
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

#                                UNIQUE READS:
#               Uniquely mapped reads number |       509079
#                    Uniquely mapped reads % |       63.55%
#                      Average mapped length |       278.06
#                   Number of splices: Total |       432182
#        Number of splices: Annotated (sjdb) |       0
#                   Number of splices: GT/AG |       416696
#                   Number of splices: GC/AG |       14478
#                   Number of splices: AT/AC |       21
#           Number of splices: Non-canonical |       987
#                  Mismatch rate per base, % |       0.24%
#                     Deletion rate per base |       0.00%
#                    Deletion average length |       1.84
#                    Insertion rate per base |       0.00%
#                   Insertion average length |       2.06

#                         MULTI-MAPPING READS:
#    Number of reads mapped to multiple loci |       201001
#         % of reads mapped to multiple loci |       25.09%
#    Number of reads mapped to too many loci |       368
#         % of reads mapped to too many loci |       0.05%

#                              UNMAPPED READS:
#number of reads unmapped: too many mismatches |       0
#   % of reads unmapped: too many mismatches |       0.00%
#        Number of reads unmapped: too short |       90538   --> No_map_short
#             % of reads unmapped: too short |       11.30%  
#            Number of reads unmapped: other |       59
#                 % of reads unmapped: other |       0.01%

#                              CHIMERIC READS:
#                   Number of chimeric reads |       0
#                        % of chimeric reads |       0.00% 

##################	
# Programme      : 
##################

SBL="Stats_Log_Final.csv";

# Nettoyage de l'ancien log : 
if [ -f "$SBL" ]; then
    rm "$SBL"
fi
    
    echo Sample,Date_Mapping,Total_reads,\
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
    
    Sample=$(echo "$bam" | sed -E "s/.*Align_(.*)_Log.final.out/\1/")
    Date_Mapping=$(grep "Finished on" "$bam" | Extract)

    # 2. UNIQUE READS 
    
    Total_reads=$(grep "Number of input reads" "$bam" | Extract)	
    Unique_reads=$(grep "Uniquely mapped reads number" "$bam" | Extract)
    Unique_pct=$(grep "Uniquely mapped reads %" "$bam" | Extract | tr -d '%')

    # 3. MULTI-MAPPING READS 
    
    Multi_reads=$(grep "Number of reads mapped to multiple loci" "$bam" | Extract)
    Multi_pct=$(grep "% of reads mapped to multiple loci" "$bam" | Extract)

    # 4. UNMAPPED READS 

	# 4.1 Extraction des valeurs numériques 
    No_map_pct_mismatch=$(grep "% of reads unmapped: too many mismatches" "$bam" | Extract | tr -d '%' || echo "0")
    No_map_pct_tooshort=$(grep "% of reads unmapped: too short" "$bam" | Extract | tr -d '%' || echo "0")
    No_map_pct_other=$(grep "% of reads unmapped: other" "$bam" | Extract | tr -d '%' || echo "0")

	# 4.2 Somme des pourcentages avec bc (résultat avec 2 décimales)
    No_map_pct_sum=$(echo "scale=2; $No_map_pct_mismatch + $No_map_pct_tooshort + $No_map_pct_other" | bc)

	# 4.3 Nombre total de lectures non mappées
    No_map_reads=$(( 
  	$(grep "Number of reads unmapped: too many mismatches" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: too short" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: other" "$bam" | Extract || echo 0)
	))
   
    # 5. AVERAGES :
    Avg_read_len=$(grep "Average input read length" "$bam" | Extract) #Longueur moyenne des reads en entrée : on va vérifier la cohérence de taille entre les runs grace à cette metrique (entre runs)
    Avg_read_map_len=$(grep "Average mapped length" "$bam" | Extract) #Longueur moyenne des segments après alignement 
    	
    	#Memo:
    	#Des reads courts ou des longueurs moyennes faibles → suspicion de dégradation d’ARN (ex. tissus anciens, conservation défaillante).
	#Un mapped length bien plus faible que le read length → suspicion de soft clipping, ARN fragmentés ou contamination.

    # 6. SPLICES :
    Splices_total=$(grep "Number of splices: Total" "$bam" | Extract)
    Splices_GTAG=$(grep -F "Number of splices: GT" "$bam" | Extract || echo "NA")
    Splices_GCAG=$(grep -F "Number of splices: GC" "$bam" | Extract || echo "NA")
    Splices_ATAC=$(grep -F "Number of splices: AT" "$bam" | Extract || echo "NA")
    Splices_Noncanonical=$(grep "Number of splices: Non-canonical" "$bam" | Extract) 
    
    #7. SUBSTITUTIONS : Pourcentage de bases dans les lectures alignées qui ne correspondent pas à la séquence de référence :
    Mismatch_rate=$(grep "Mismatch rate per base" "$bam" | Extract | tr -d '%')
    Deletion_rate=$(grep "Deletion rate per base" "$bam" | Extract)
    Insertion_rate=$(grep "Insertion rate per base" "$bam" | Extract)

    #8. CHIMERICS : Plutot utile en somatique cela, mais je l'extrait comme cela on à tout à l'éxécution du .sh, (pour identifier des transcrits de fusion, transloc ?)
    Chimeric_reads=$(grep "Number of chimeric reads" "$bam" | Extract)
    Chimeric_reads_pct=$(grep "% of chimeric reads" "$bam" | Extract | tr -d '%')
  
    #9. Écriture dans le tabulé :
    echo "${Sample},${Date_Mapping},${Total_reads},${Unique_reads},${Unique_pct},${Multi_reads},${Multi_pct},${No_map_reads},${No_map_pct_sum},${No_map_pct_mismatch},${No_map_pct_tooshort},${No_map_pct_other},${Avg_read_len},${Avg_read_map_len},${Splices_total},${Splices_GTAG},${Splices_GCAG},${Splices_ATAC},${Splices_Noncanonical},${Mismatch_rate},${Deletion_rate},${Insertion_rate},${Chimeric_reads},${Chimeric_reads_pct}" >> "$SBL"


done
