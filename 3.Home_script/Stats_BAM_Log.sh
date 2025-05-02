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
# - Permet d'évaluer la reproductibilité des run RNASeq (pour notre cas d'étude). 

##############
# Input Star: 
##############

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

# La première partie du script fait des vérifications d'usage,  et génère le fichier tabulé csv avec les en-têtes d'intérêts ;

SBL="Stats_Log_Final.csv";

# En-tête de l'output CSV si non générée: 
if [ ! -f "$SBL" ]; then
    echo Sample,Date_Mapping,Total_reads,Unique_reads,Unique_pct,Multi_reads,Multi_pct,No_map_reads,No_map_pct > "$SBL";
fi

# Fonction pour extraire le 2ᵉ champ et enlever les espaces
Extract() { cut -d '|' -f2 | tr -d ' '; }

# Pour chaque fichier du répertoire courant et des sous-répertoires :
for bam in *.Log.final.out;
do
    # ─── 1. ID ──────────────────────────────────────────────
    
    Sample=$(echo "$bam" | sed -E "s/.*Align_(.*)_Log.final.out/\1/")
    Date_Mapping=$(grep "Finished on" "$bam" | Extract)

    # ─── 2. UNIQUE READS ────────────────────────────────────
    
    Total_reads=$(grep "Number of input reads" "$bam" | Extract)	
    Unique_reads=$(grep "Uniquely mapped reads number" "$bam" | Extract)
    Unique_pct=$(grep "Uniquely mapped reads %" "$bam" | Extract | tr -d '%')

    # ─── 3. MULTI-MAPPING READS ─────────────────────────────
    
    Multi_reads=$(grep "Number of reads mapped to multiple loci" "$bam" | Extract)
    Multi_pct=$(grep "% of reads mapped to multiple loci" "$bam" | Extract)

    # ─── 4. UNMAPPED READS ──────────────────────────────────

	# Extraction des valeurs numériques 
    No_map_pct_mismatch=$(grep "% of reads unmapped: too many mismatches" "$bam" | Extract | tr -d '%' || echo "0")
    No_map_pct_tooshort=$(grep "% of reads unmapped: too short" "$bam" | Extract | tr -d '%' || echo "0")
    No_map_pct_other=$(grep "% of reads unmapped: other" "$bam" | Extract | tr -d '%' || echo "0")

	# Somme des pourcentages avec bc (résultat avec 2 décimales)
    No_map_pct_sum=$(echo "scale=2; $No_map_pct_mismatch + $No_map_pct_tooshort + $No_map_pct_other" | bc)

	# Nombre total de lectures non mappées
     No_map_reads=$(( 
  	$(grep "Number of reads unmapped: too many mismatches" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: too short" "$bam" | Extract || echo 0) + 
  	$(grep "Number of reads unmapped: other" "$bam" | Extract || echo 0)
	))

     # ─── 5. Écriture dans le CSV :
    echo "${Sample},${Date_Mapping},${Total_reads},${Unique_reads},${Unique_pct},${Multi_reads},${Multi_pct},${No_map_reads},${No_map_pct_sum},${No_map_pct_mismatch},${No_map_pct_tooshort},${No_map_pct_other}" >> "$SBL"


done
