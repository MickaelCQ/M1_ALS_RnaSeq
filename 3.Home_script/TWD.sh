#!/bin/bash

##################################################################################
# Auteur   : Mickael Coquerelle
# Date     : 03/05/2025
# Projet   : TWD - Time Work DevOps --  Version : 1.2.0
# Licence : CC BY 4.0 (Creative Commons Attribution 4.0 International)
#
# Ce script est libre d’utilisation, de modification et de redistribution
# à condition de créditer l’auteur original.
# Pour en savoir plus : https://creativecommons.org/licenses/by/4.0/
#
# Objectif du script :
#   - Quantifier le temps de travail de manière flexible et précise pour un projet donné.
#   - Produire un tableau du travail réellement fourni, utile pour l’auto-évaluation.
#   - Qualifier le volume horaire pour chaque tâche et identifier des axes d’amélioration.
#
# Utilisation :
#   ./TWD.sh <Nom_du_projet> <Nom_du_fichier_CSV> <Tache(s)> [Commentaire] [Langage]
#
#   $1 : Nom du projet (sera stocké dans la colonne "Projet")
#   $2 : Nom du fichier CSV de log (doit exister préalablement)
#   $3 : Description de la tâche ou activité en cours
#   $4 : (optionnel) Commentaire libre (ex: "Sprint 1", "Fix bug X"…)
#   $5 : (optionnel) Langage utilisé (ex: BASH, R, C++, PYTHON, LATEX, RMARKDOWN…)
##################################################################################

# Nom du fichier de log principal (structure standard si utilisé sans $2)
TWD="Log_DevOps.csv"

# Création de l'en-tête si le fichier principal n'existe pas
if [ ! -f "$TWD" ]; then
    echo "Debut_de_la_session,Fin_de_la_session,Duree,Projet,Tache(s),Langage(s),Commentaire(s)," > "$TWD"
fi

# Vérification des arguments obligatoires
if [ "$#" -lt 3 ]; then
    echo "Erreur : Trois arguments minimum sont requis."
    echo "Syntaxe : ./TWD.sh <Nom_du_projet> <Nom_du_fichier_CSV> <Tache(s)> [Langage] [Commentaire] "
    exit 1
fi

# Attribution des variables à partir des arguments
nom_projet="$1"
fichier_csv="$2"
tache="$3"
commentaire="${5:-Aucun}"
langage="${4:-NA}"

# Vérification de l’existence du fichier CSV donné en argument
if [ ! -f "$fichier_csv" ]; then
    echo "Erreur : Le fichier CSV '$fichier_csv' n'existe pas. Veuillez le créer manuellement."
    exit 2
fi

# Récupération de la dernière ligne correspondant au projet
derniere_ligne=$(grep "$nom_projet" "$fichier_csv" | tail -n 1)

# Extraction des heures de début et de fin
heure_debut=$(echo "$derniere_ligne" | cut -d ',' -f 1 | xargs)
heure_fin=$(echo "$derniere_ligne" | cut -d ',' -f 2 | xargs)

# Cas 1 : Nouvelle session (pas de session ouverte ou déjà clôturée)
if [ -z "$heure_debut" ] || [ -n "$heure_fin" ]; then
    heure_debut=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$heure_debut,,,$nom_projet,\"$tache\",\"$langage\",\"$commentaire\"" >> "$fichier_csv"
    echo "Nouvelle session de travail commencée pour le projet '$nom_projet'."
else
    # Cas 2 : Clôture d'une session existante
    heure_fin=$(date "+%Y-%m-%d %H:%M:%S")

    # Calcul de la durée en secondes
    secondes_debut=$(date -d "$heure_debut" +%s)
    secondes_fin=$(date -d "$heure_fin" +%s)
    duree_en_secondes=$((secondes_fin - secondes_debut))

    # Conversion de la durée en format HH:MM
    heures=$((duree_en_secondes / 3600))
    minutes=$(((duree_en_secondes % 3600) / 60))
    duree=$(printf "%02d:%02d" $heures $minutes)

    # Sécurisation des champs à insérer dans la ligne mise à jour pour exploiter les caractères spéciaux. 
    escaped_commentaire=$(printf "%s" "$commentaire" | sed 's/[&/\]/\\&/g')
    escaped_tache=$(printf "%s" "$tache" | sed 's/[&/\]/\\&/g')
    escaped_langage=$(printf "%s" "$langage" | sed 's/[&/\]/\\&/g')

    # Mise à jour de la ligne dans le fichier CSV avec sed
    sed -i "s|^$heure_debut,,,$nom_projet.*|$heure_debut,$heure_fin,$duree,$nom_projet,\"$escaped_tache\",\"$escaped_langage\",\"$escaped_commentaire\",|" "$fichier_csv"

    echo "Session de travail clôturée pour le projet '$nom_projet'. Temps total : $duree."
fi

