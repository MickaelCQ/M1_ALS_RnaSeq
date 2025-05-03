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
#          VARIABLES
# ############################ #
MSQ = @
MAIN = main

# Paramètres pour latex et LA documentation doxygen :
TEXFILE = 2.Report_latex/$(MAIN).tex
PDF = $(MAIN).pdf
LATEX = $(MSQ)pdflatex -interaction=nonstopmode -shell-escape
BIB = biber
Doxy = 1.Documentation/Doxygen
BUILDLATEX = 2.Report_latex/PDF

#############################################################################################################################################################
# Lancement de STAR en paire-end sur le format de fastq fourni
#############################################################################################################################################################

THREADS = 16
REF = /data/indexes/STAR/2.7.8a/GRCh37
FASTQ_DIR ?= fastq
DATE = $(shell date +%d%m%Y)
OUTBAM_DIR = BAM_STAR_$(DATE)
OUTLOG_DIR = LOGS_STAR_$(DATE)
TMPDIR = TMP_STAR_$(DATE)

# Détection des échantillons par présence des fichiers _1 et _2
SAMPLES = $(shell \
    for R1 in $(FASTQ_DIR)/*_1.fastq.gz; do \
        R2=`echo $$R1 | sed 's/_1.fastq.gz/_2.fastq.gz/'`; \
        if [ -f $$R2 ]; then \
            basename $$R1 _1.fastq.gz; \
        fi; \
    done)

BAMS = $(addprefix $(OUTBAM_DIR)/, $(addsuffix .bam, $(SAMPLES)))

# Cible principale
Star_Paire: $(BAMS)

$(OUTBAM_DIR)/%.bam:
	@mkdir -p $(OUTBAM_DIR) $(OUTLOG_DIR) $(TMPDIR)/$*
	@echo ">> Alignement de l'échantillon $*"
	STAR --runThreadN $(THREADS) \
	     --genomeDir $(REF) \
	     --readFilesIn $(FASTQ_DIR)/$*_1.fastq.gz $(FASTQ_DIR)/$*_2.fastq.gz \
	     --readFilesCommand zcat \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix $(TMPDIR)/$*/ 
	@mv $(TMPDIR)/$*/Aligned.sortedByCoord.out.bam $(OUTBAM_DIR)/$*.bam
	@mv $(TMPDIR)/$*/Log.final.out $(OUTLOG_DIR)/$*.Log.final.out

#############################################################################################################################################################
# Cibles pour la documentation Doxygen
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
	cd $(Doxy)/latex && pdflatex refman.tex && pdflatex refman.tex
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
.PHONY: pdfLatex

