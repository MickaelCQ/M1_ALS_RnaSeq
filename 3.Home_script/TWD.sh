#!/bin/bash

################################################################################## Mickael -- le 27/12/2024 - ####################################################################################
#                                                                                   TWD : Time Work DevOps                                                                                       #
##################################################################################################################################################################################################

# L'objectif de ce script est de quantifier mon temps de travail de manière flexible et précise pour un projet donné, cette idée de quantification à émerger des nombreuses remarques du
# Professeur Alban Mancheron concernant le temps théorique à alloué pour légitimer l'octroi des ECTS (nos fameux crédits européens), ce petit script en plus de nourrir ma curiosité, me permet de :
    #    Quantifier de manière quantitative et autonome mon temps afin d'avoir un tableau du travail réellement fournit.
    #    Quantifier de manière qualititative mon temps , ce script peux me permettre d'améliorer ma manière de travailler, en identifiant une incohérence du volume horaire pour une tache.
    #    Fournir au Prof un indicateur qui se veut (qui essaye) d'être au plus juste sur mon investissement réel dans sa matière (et mon intérêt pour la bio informatique).
    #    Ma permis de m'exercer dans bash et donc d'apprendre de nouveaux trucs.

# L'utilisateur fournit les arguments suivants :
# $1 : Nom du projet (sera stocké dans une colonne "Projet" du CSV).
# $2 : Nom du fichier CSV à utiliser (doit exister, pas d'auto-création).
# $3 : Commentaire décrivant la session (par ex., "Début", "Fin du job journalier", etc.).

# Vérification des arguments
if [ "$#" -lt 3 ]; then
    echo "Trois arguments sont attendus à la suite de l'appel de TWD.sh <Nom_du_projet> <Nom_du_fichier_CSV> <Commentaire>"
    exit 1
fi

# Définition des variables
nom_projet="$1"
fichier_csv="$2"
commentaire="$3"

# Vérification d'usage .... 
if [ ! -f "$fichier_csv" ]; then
    echo "Erreur : Le fichier CSV  renseigné en argument '$fichier_csv' n'existe pas. Veuillez le créer avant d'exécuter ce script."
    exit 2
fi



# Récupération de ma dernière ligne correspondant au projet en cours
derniere_ligne=$(grep "$nom_projet" "$fichier_csv" | tail -n 1)

# Extraction des données de la dernière ligne
heure_debut=$(echo "$derniere_ligne" | cut -d ',' -f 1 | xargs)
heure_fin=$(echo "$derniere_ligne" | cut -d ',' -f 2 | xargs)

if [ -z "$heure_debut" ] || [ -n "$heure_fin" ]; then
    # Si aucune session en cours ou si la dernière session est déjà clôturée
    heure_debut=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$heure_debut, , , $nom_projet, \"$commentaire\"" >> "$fichier_csv"
    echo "Nouvelle session de travail commencée pour le projet '$nom_projet'."
else
    # Si une session en cours
    heure_fin=$(date "+%Y-%m-%d %H:%M:%S")
    # Calcul du delta en secondes
    secondes_debut=$(date -d "$heure_debut" +%s)
    secondes_fin=$(date -d "$heure_fin" +%s)
    duree_en_secondes=$((secondes_fin - secondes_debut))

    # Conversion en heures et minutes
    heures=$((duree_en_secondes / 3600))
    minutes=$(((duree_en_secondes % 3600) / 60))
    duree=$(printf "%02d:%02d" $heures $minutes)

    # Mise à jour de la dernière ligne avec l'heure de fin et la durée
    sed -i "s|^$heure_debut, , , $nom_projet, .*|$heure_debut, $heure_fin, $duree, $nom_projet, \"$commentaire\"|" "$fichier_csv"
    echo "Session de travail clôturée pour le projet '$nom_projet'. Temps total : $duree."
fi
