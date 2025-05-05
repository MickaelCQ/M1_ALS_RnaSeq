# ##########################################################################################################################################################################
# Ce fichier est sous la licence MIT (Massachusetts Institute of Technology):
# ##########################################################################################################################################################################
# Fichier Makefile COUTEAU-SUISSE : Copyright (c) [2025] [COQUERELLE Mickael]
# ##########################################################################################################################################################################
# Permission est accordée, gratuitement, à toute personne obtenant une copie de ce logiciel et des fichiers associés (le "Logiciel"), de traiter le Logiciel sans restriction,
# y compris sans limitation les droits d'utilisation, de copie, de modification, de fusion, de publication, de distribution, de sous-licence et/ou de vente de copies du Logiciel,
# et de permettre à d'autres personnes à qui le Logiciel est fourni de le faire, sous les conditions suivantes :
# L'avis de droit d'auteur ci-dessus et cet avis de permission doivent être inclus dans toutes les copies ou portions substantielles du Logiciel.
# LE LOGICIEL EST FOURNI "TEL QUEL", SANS GARANTIE D'AUCUNE SORTE, EXPRESSE OU IMPLICITE, Y COMPRIS MAIS SANS S'Y LIMITER LES GARANTIES DE QUALITÉ MARCHANDE, D'ADAPTATION À UN USAGE PARTICULIER ET DE NON-VIOLATION.
# EN AUCUN CAS, LES AUTEURS OU DÉTENTEURS DES DROITS D'AUTEUR NE SERONT RESPONSABLES DE TOUTE RÉCLAMATION, DOMMAGES OU AUTRES RESPONSABILITÉS, QUE CE SOIT DANS UNE ACTION EN VERTU DU CONTRAT, DÉLICTUELLE OU AUTRE,
# DÉCOULANT DE, OU EN RELATION AVEC, LE LOGICIEL OU L'UTILISATION DE CE DERNIER.
# ##########################################################################################################################################################################

# ##############################
#    VARIABLES PRINCIPALES
# ##############################
MSQ = @
MAIN = main


# Paramètres pour latex et LA documentation doxygen :
TEXFILE = 2.Report_latex/$(MAIN).tex
PDF = $(MAIN).pdf
LATEX = $(MSQ)pdflatex -interaction=nonstopmode -shell-escape
BIB = biber
Doxy = 1.Documentation/Doxygen
BUILDLATEX = 2.Report_latex/PDF

# ##############################
#    VARIABLES ALIGNEMENT
# ##############################
THREADS = 16
REF_STAR = /data/indexes/STAR/2.7.8a/GRCh37

DATE = $(shell date +%d%m%Y)

OUTBAM_STAR = BAM_STAR_$(DATE)
OUTLOG_STAR = LOGS_STAR_$(DATE)
TMPDIR = TMP_STAR_$(DATE)


REF_CRAC = /data/indexes/crac/GRCh37
OUTBAM_CRAC = BAM_CRAC_$(DATE)
OUTSAM_CRAC = SAM_CRAC_$(DATE)
OUTLOG_CRAC = LOGS_CRAC_$(DATE)
OUTERR_CRAC = ERR_CRAC_$(DATE)
TMPDIR_CRAC = TMP_CRAC_$(DATE)
KMER_CRAC= 22

################################
#       REPERTOIRE FASTQ       
################################
FASTQ_DIR ?= fastq

SAMPLES = $(shell \
    for R1 in $(FASTQ_DIR)/*_1.fastq.gz; do \
        R2=`echo $$R1 | sed 's/_1.fastq.gz/_2.fastq.gz/'`; \
        if [ -f $$R2 ]; then \
            basename $$R1 _1.fastq.gz; \
        fi; \
    done)

#############################################################################################################################################################
# 					Lancement de CRAC en paire-end sur le format de fastq fourni
#############################################################################################################################################################

#Trame :  --summary output/crac/summary/202304-1409241167-SLA.summary --nb-threads 16 -r fastq/202304-1409241167-SLA_1.fastq.gz -o - fastq/202304-1409241167-SLA_2.fastq.gz
# Détection des échantillons par présence des fichiers _1 et _2
BAMS_CRAC = $(addprefix $(OUTBAM_CRAC)/, $(addsuffix .bam, $(SAMPLES)))


#Voila la cible initial que je te parler il faut garder au maximum le nom des variables d'environnement :
	
# Cible principale pour l'alignement CRAC
Crac_Paire: $(BAMS_CRAC)

$(OUTBAM_CRAC)/%.bam:
	@echo "Nettoyage des sorties eventuelles précédentes : "
	@rm -rf BAM_CRAC* ERR_CRAC* LOGS_CRAC* TMP_CRAC* OUT
	
	@mkdir -p $(OUTSAM_CRAC) $(OUTBAM_CRAC) $(OUTLOG_CRAC)/$* $(OUTERR_CRAC)/$* $(TMPDIR_CRAC)/$*
	@echo ">> Alignement de l'échantillon $* avec CRAC"
	crac --sam --stranded --rf \
	     -i $(REF_CRAC) \
	     -k $(KMER_CRAC) \
	     --summary $(OUTLOG_CRAC)/$*/Log.final.out \
	     --nb-threads $(THREADS) \
	     -r $(FASTQ_DIR)/$*_1.fastq.gz $(FASTQ_DIR)/$*_2.fastq.gz \
	     > $(OUTSAM_CRAC)/$*.sam \
	     2> $(OUTERR_CRAC)/$*/crac.ERR.log
	@mv $(TMPDIR_CRAC)/$*/Log.final.out $(OUTLOG_CRAC)/Log.final.out
	
	# Conversion du SAM en BAM avec SAMTOOLS
	@echo ">> Conversion du SAM en BAM avec Samtools : "
	samtools sort $(OUTSAM_CRAC)/$*.sam -@ $(THREADS) -o $(OUTBAM_CRAC)/$*.bam
	
	# Indexaction du fichier bam : 
	@echo ">> Indexation du fichier BAM : "
	samtools index $(OUTBAM_CRAC)/$*.bam $(OUTBAM_CRAC)/$*.bam.bai
	
	# Nettoyage après conversion en BAM
	rm -rf $(OUTSAM_CRAC)/$*.sam
	
#############################################################################################################################################################
# 					Lancement de STAR en paire-end sur le format de fastq fourni
#############################################################################################################################################################

# Détection des échantillons par présence des fichiers _1 et _2
BAMS_STAR = $(addprefix $(OUTBAM_STAR)/, $(addsuffix .bam, $(SAMPLES)))

# Cible principale
Star_Paire: $(BAMS_STAR)

$(OUTBAM_STAR)/%.bam:
	@mkdir -p $(OUTBAM_STAR) $(OUTLOG_STAR) $(TMPDIR)/$* 
	@echo ">> Alignement de l'échantillon $* avec STAR"
	STAR --runThreadN $(THREADS) \
	     --genomeDir $(REF_STAR) \
	     --readFilesIn $(FASTQ_DIR)/$*_1.fastq.gz $(FASTQ_DIR)/$*_2.fastq.gz \
	     --readFilesCommand zcat \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix $(TMPDIR)/$*/ 
	@mv $(TMPDIR)/$*/Aligned.sortedByCoord.out.bam $(OUTBAM_STAR)/$*.bam
	@mv $(TMPDIR)/$*/Log.final.out $(OUTLOG_STAR)/$*.Log.final.out
	
#############################################################################################################################################################
# 							REGLE DE COMPILATION documentation Doxygen
#############################################################################################################################################################

GenDoc:
	$(MSG) "Génération de la documentation avec Doxygen..."
	$(MSQ)doxygen 2.Documentation/Doxygen/Doxyfile -output $(Doxy)
	$(MSG) "Documentation générée avec succès dans $(Doxy)"

ReadDoc:
	$(MSG) "Ouverture de la documentation dans Firefox..."
	$(MSQ)firefox ./html/index.html &
	$(MSG) "Documentation ouverte dans Firefox."

CleanDoc:
	$(MSG) "Suppression des fichiers de documentation..."
	$(RM) "$(Doxy)/html" "$(Doxy)/latex"
	$(MSG) "Nettoyage de la documentation terminé."

GenLatex:
	$(MAKE) $(Doxy)/latex/refman.pdf

$(Doxy)/latex/refman.pdf: GenDoc
	$(MSG) "Compilation du LaTeX en PDF..."
	cd $(Doxy)/latex && pdflatex	 refman.tex && pdflatex refman.tex
	$(MSG) "PDF généré : $(Doxy)/latex/refman.pdf"

ReadLatex: $(Doxy)/latex/refman.pdf
	$(MSG) "Ouverture du PDF avec Evince..."
	$(MSQ)evince $(Doxy)/latex/refman.pdf &

# ############################ #
# REGLES COMPILATION LATEX
# ############################ #

pdfLatex:
	mkdir -p $(BUILDLATEX)
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	$(BIB) $(BUILDLATEX)/$(MAIN)
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	find $(BUILDLATEX) -type f -name "$(MAIN).*" ! -name "$(MAIN).pdf" -delete
	rm -rf ./_minted-$(MAIN)
	xdg-open $(BUILDLATEX)/$(MAIN).pdf > /dev/null 2>&1 &

cleanLatex:
	find $(BUILDLATEX) -type f -name "$(MAIN).*" ! -name "$(MAIN).pdf" -delete
	rm -rf ./_minted-$(MAIN)

# ############################ #
.PHONY: pdfLatex cleanLatex ReadLatex GenLatex CleanDoc ReadDoc GenDoc Star_Paire Crac_Paire

