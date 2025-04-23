	
# ##########################################################################################################################################################################
#																					   #																						 
#						 Ce fichier est sous la licence MIT (Massachusetts Institute of Technology):                                               #
#																					   #
# ##########################################################################################################################################################################
#																					   #
# 						Fichier Makefile COUTEAU-SUISSE : Copyright (c) [2025] [COQUERELLE Mickael]						   #
#																					   #
# ##########################################################################################################################################################################
																					   #
# Permission est accordée, gratuitement, à toute personne obtenant une copie de ce logiciel et des fichiers associés (le "Logiciel"), de traiter le Logiciel               #
# sans restriction, y compris sans limitation les droits d'utilisation, de copie, de modification, de fusion, de publication, de distribution, de sous-licence et/ou       #
# de vente de copies du Logiciel, et de permettre à d'autres personnes à qui le Logiciel est fourni de le faire, sous les conditions suivantes :                           #
# 																					   #
# L'avis de droit d'auteur ci-dessus et cet avis de permission doivent être inclus dans toutes les copies ou portions substantielles du Logiciel.                          #
#  																					   #
# LE LOGICIEL EST FOURNI "TEL QUEL", SANS GARANTIE D'AUCUNE SORTE, EXPRESSE OU IMPLICITE, Y COMPRIS MAIS SANS S'Y LIMITER LES GARANTIES DE QUALITÉ MARCHANDE, D'ADAPTATION # 
# À UN USAGE PARTICULIER ET DE NON-VIOLATION. EN AUCUN CAS, LES AUTEURS OU DÉTENTEURS DES DROITS D'AUTEUR NE SERONT RESPONSABLES DE TOUTE RÉCLAMATION, DOMMAGES OU AUTRES  #
# RESPONSABILITÉS, QUE CE SOIT DANS UNE ACTION EN VERTU DU CONTRAT, DÉLICTUELLE OU AUTRE, DÉCOULANT DE, OU EN RELATION AVEC, LE LOGICIEL OU L'UTILISATION DE CE DERNIER.   #
# ##########################################################################################################################################################################

# ##############################
#          VARIABLES	       #	
# ############################ #

# Préfixe @ pour exécuter des commandes sans afficher cette dernière à chaque éxécution (lisibilité.)
MSQ      = @ 
MAIN     = main

TEXFILE  = 2.Report_latex/$(MAIN).tex
PDF      = $(MAIN).pdf
LATEX    = $(MSQ)pdflatex -interaction=nonstopmode -shell-escape
BIB      = biber

Doxy = 1.Documentation/Doxygen
BUILDLATEX   = 2.Report_latex/PDF


#############################################################################################################################################################
#                                                              Cibles pour la documentation  Doxygen                                                               
#############################################################################################################################################################

# Cible pour générer la documentation
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

GenLatex: $(Doxy)/latex/refman.pdf
	$(Doxy)/latex/refman.pdf: GenDoc
	$(MSG) "Compilation du LaTeX en PDF..."
	cd $(Doxy)/latex && pdflatex refman.tex && pdflatex refman.tex
	$(MSG) "PDF généré : $(Doxy)/latex/refman.pdf"

ReadLatex: $(Doxy)/latex/refman.pdf
	$(MSG) "Ouverture du PDF avec Evince..."
	$(MSQ)evince $(Doxy)/latex/refman.pdf &

# ################################ #
#      REGLES COMPILATION LATEX    #	
# ################################ #

# Cible principale de compilation LaTeX avec nettoyage intégré
pdfLatex:
	# Crée le dossier de build s'il n'existe pas
	mkdir -p $(BUILDLATEX)
	
	# 1ère passe LaTeX
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	
	# Bibliographie
	$(BIB) $(BUILDLATEX)/$(MAIN)
	
	# 2e et 3e passe LaTeX
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	$(LATEX) -output-directory=$(BUILDLATEX) $(TEXFILE)
	
	# Nettoyage des fichiers auxiliaires (dans le dossier PDF)
	find $(BUILDLATEX) -type f -name "$(MAIN).*" ! -name "$(MAIN).pdf" -delete
	
	# Suppression explicite du dossier _minted-main à la racine
	rm -rf ./_minted-$(MAIN)
	
	# Ouveture du pdf sans blocage du terminal : 
	xdg-open $(BUILDLATEX)/$(MAIN).pdf > /dev/null 2>&1 &

cleanLatex:
	find $(BUILDLATEX) -type f -name "$(MAIN).*" ! -name "$(MAIN).pdf" -delete
	rm -rf ./_minted-$(MAIN)

# ############################ #
.PHONY: pdfLatex 
