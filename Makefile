# ##########################################################################################################################################################################
# Ce fichier est sous la licence MIT (Massachusetts Institute of Technology):
# ##########################################################################################################################################################################
# Fichier Makefile COUTEAU-SUISSE : Copyright (c) [2025] [COQUERELLE Mickael]
# ##########################################################################################################################################################################
# Permission est accordée, gratuitement, à toute personne obtenant une copie de ce programme  et des fichiers associés (le "Logiciel"), de traiter le Logiciel sans restriction,
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
REF_STAR = /data/indexes/STAR/2.7.8a/GRCh37/ 
REF_CRAC = /data/indexes/crac/GRCh37/GRCh37newIndexSansMito

DATE = $(shell date +%d%m%Y)

OUTLOG_STAR = LOGS_STAR_$(DATE)
OUTBAM_STAR = BAM_STAR_$(DATE)
TMPDIR = TMP_STAR_$(DATE)


#REF_CRAC = /data/indexes/crac/GRCh37/GRCh37newIndexSansMito	
REF_CRAC = /data/indexes/crac/GRCh37_with_MT/GRCh37newIndexAvecMito

OUTBAM_CRAC = BAM_CRAC_$(DATE)
OUTSAM_CRAC = SAM_CRAC_$(DATE)
TMPDIR_CRAC = TMP_CRAC_$(DATE)
KMER_CRAC= 22
UNZIP_FASTQ = fastq_unzipped

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

# Cibles finales souhaitées
CRAC_SAMS = $(addprefix output/crac/bam/, $(addsuffix .sam, $(SAMPLES)))

# Cible principale
Crac_Paire: $(CRAC_SAMS)

# Règle pour chaque .sam
output/crac/bam/%.sam: $(FASTQ_DIR)/%_1.fastq.gz $(FASTQ_DIR)/%_2.fastq.gz
	@echo "Vérification de l'échantillon : $*"
	@if [ ! -f output/crac/bam/$*.bam ]; then \
	    echo "Lancement de l'alignement de l'échantillon : $* avec Crac"; \
	    mkdir -p output/crac/summary output/crac/log output/crac/bam; \
	    gunzip -c $(FASTQ_DIR)/$*_1.fastq.gz > $(FASTQ_DIR)/$*_1.fastq; \
	    gunzip -c $(FASTQ_DIR)/$*_2.fastq.gz > $(FASTQ_DIR)/$*_2.fastq; \
	    crac --nb-tags-info-stored 10000 --bam --stranded \
	         -i $(REF_CRAC) -k $(KMER_CRAC) \
	         --summary output/crac/summary/$*.summary \
	         --nb-threads $(THREADS) \
	         -r $(FASTQ_DIR)/$*_1.fastq $(FASTQ_DIR)/$*_2.fastq \
	         -o output/crac/bam/$*.sam \
	         2> output/crac/log/$*_crac.log; \
	    rm -f $(FASTQ_DIR)/$*_1.fastq $(FASTQ_DIR)/$*_2.fastq; \
	else \
	    echo "Fichier déjà aligné : $@"; \
	fi
	
#############################################################################################################################################################
#							Compression des BAM en SAM et indexation des BAM 
#############################################################################################################################################################	
# Génère une liste des fichiers SAM sans extension
SAM_FILES := $(basename $(notdir $(wildcard output/crac/bam/*.sam)))

# Règle principale pour convertir tous les SAM en BAM+BAI
Convert_bam_sam: $(addprefix output/crac/bam/, $(addsuffix .bam, $(SAM_FILES)))

output/crac/bam/%.bam: output/crac/bam/%.sam
	@echo "Conversion SAM -> BAM pour l’échantillon : $*"
	#Conversion (view), en binaire (-b), avec une entrée SAM (-S) ET et trie du BAM en position génomique (sort) :
	samtools view -@ $(THREADS) -bS $< | samtools sort -@ $(THREADS) -o $@
	samtools index $@
	@rm -f $<

#############################################################################################################################################################
# 					Lancement de STAR en paire-end sur le format de fastq fourni
#############################################################################################################################################################

# Détection des échantillons par présence des fichiers _1 et _2
BAMS_STAR = $(addprefix $(OUTBAM_STAR)/, $(addsuffix .bam, $(SAMPLES)))

# Cible principale
Star_Paire_All: $(BAMS_STAR)

$(OUTBAM_STAR)/%.bam:
	@mkdir -p $(OUTBAM_STAR) $(OUTLOG_STAR) $(TMPDIR)/$* 
	@echo ">>Lancement de l'alignement de l'échantillon $* avec STAR"
	STAR --runThreadN $(THREADS) \
	     --genomeDir $(REF_STAR) \
	     --readFilesIn $(FASTQ_DIR)/$*_1.fastq.gz $(FASTQ_DIR)/$*_2.fastq.gz \
	     --readFilesCommand zcat \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix $(TMPDIR)/$*/ 
	@mv $(TMPDIR)/$*/Aligned.sortedByCoord.out.bam $(OUTBAM_STAR)/$*.bam
	@mv $(TMPDIR)/$*/Log.final.out $(OUTLOG_STAR)/$*.Log.final.out

Star_Paire_Alone:

	@if [ -z "$(FQ1)" ] || [ -z "$(FQ2)" ] || [ -z "$(OUT)" ]; then \
		echo "Utilisation : make Star_Paire_Alone FQ1=... FQ2=... OUT=SampleName" && exit 1; \
	fi
	@echo ">> Lancement de l'alignement de $(OUT) avec les fichiers :"
	@echo "   - $(FQ1)"
	@echo "   - $(FQ2)"
	STAR --runThreadN $(THREADS) \
	     --genomeDir $(REF_STAR) \
	     --readFilesIn $(FQ1) $(FQ2) \
	     --readFilesCommand zcat \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix $(TMPDIR)/$(OUT)/
	
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
#  SYSTEME  CONFIGURATION
# ############################ #

Update_swap:
	@if [ -z "$(size)" ]; then \
		echo "Soucis dans la valeur du swap :fournir une valeur comme :  make update_swap size=16)"; \
		exit 1; \
	fi; \
	sudo swapoff -a && sudo rm -f /swapfile; 
	# Nouveau fichier swap avec la taille donnée
	sudo fallocate -l $(size)G /swapfile && sudo chmod 600 /swapfile; \
	# Initialisation du fichier swap
	sudo mkswap /swapfile && sudo swapon /swapfile; \
	# Vérification du fichier /etc/fstab pour le rendre persistant
	if ! grep -q "/swapfile" /etc/fstab; then \
		echo "/swapfile none swap sw 0 0" | sudo tee -a /etc/fstab; \
	fi; \
	# Vérification de l'espace swap
	free -h && echo "Le swap a été mis à jour avec une taille de $(size) Go"
	
# ############################ #
.PHONY: Update_swap pdfLatex cleanLatex ReadLatex GenLatex CleanDoc ReadDoc GenDoc Star_Paire_All Star_Paire_Alone Crac_Paire
Alone
