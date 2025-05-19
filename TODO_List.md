<span style="color:#1b2a49"><h1> Cahier de labo – Analyse Kmerator / Transipédia</h1></span>

## Bibliographie

- **Kmerator** : [https://pmc.ncbi.nlm.nih.gov/articles/PMC8221386/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8221386/)

---

<span style="color:#1b2a49"><h2>Réunion avec Thérèse Commes – Lundi 05/05/2025</h2></span>

### Objectifs & pistes de traitement – *Kmerator / KmerExplorer*

1. **Utilisation de KmerExplorer pour extraire des informations techniques et biologiques**
   - Se référer à la documentation officielle.
   - Une catégorie est dédiée à la **quantification des adaptateurs**.
   - Objectif : obtenir une **corrélation entre les reads contenant des adaptateurs** et **les reads non mappés après alignement**.

2. **Visualisation de la représentativité des gènes du panel ciblé**
   - Réaliser une **représentation graphique** de la couverture des 40 gènes d’intérêt avec *KmerExplorer*.

3. **Extraction d’événements biologiques spécifiques**
   - Identifier :
     - les **mutations** ;
     - les **événements d’épissage (splicing)**.

4. **Création d’un jeu de kmers dédié aux gènes cibles**
   - Utiliser *Kmerator* pour déterminer les **kmers spécifiques à chaque gène**.
   - Générer un fichier `.yml` en sortie pour alimenter *KmerExplorer* avec ces kmers.

5. **Analyse de couverture et d’épissage avec Kmerator**
   - Chaque gène est représenté par un ensemble de kmers spécifiques.
   - Pour chaque échantillon :
     - **Compter les occurrences de ces kmers** pour estimer la **profondeur de couverture**.
     - Développer un outil graphique simple pour **visualiser les événements d’épissage**.

6. **Quantification au niveau des kmers**
   - En raison de l’hétérogénéité des résultats globaux, se focaliser sur des **comparaisons à l’échelle du kmer**.
   - Objectif : pour un kmer donné, **mesurer la différence d’abondance entre deux échantillons RNA-Seq**.
   - Avantage : une **quantification locale et fiable**.

---

## À approfondir

- Rechercher dans la documentation de CRAC/STAR une chaîne d’options permettant **de lister les sites d’épissage détectés**.

---

<span style="color:#1b2a49"><h2>Atelier Transipédia / Kmerator – 14/05/2025</h2></span>

### Thématique : utilisations pour les biologistes & interrogation de Transipédia

#### Cas d’usage

- **TCGA** : problématiques liées au cancer  
- **GTEx** : analyse des tissus sains  
- **CCLE** (ou "CCME") : modèles de lignées cellulaires tumorales  

#### Concept général

- À partir de fichiers **FASTQ/FASTA**, quelques analyses de qualité sont réalisées.
- Une **indexation par kmers** est effectuée, produisant une **table de comptage (nombre de kmers par read)**.
- Les données sont stockées dans une matrice construite avec **Reindeer index**.
- Une interface web permet ensuite de :
  - rechercher à partir d’une **liste de kmers**,
  - extraire la **sous-matrice correspondante**,
  - enrichir les résultats avec des **métadonnées** si disponibles.

#### Application à l’expression d’un transcrit

- Pour étudier un transcrit, on le découpe en **kmers spécifiques**.
- Problème : certains kmers peuvent être **non spécifiques**, ce qui biaise le comptage.
- Solution : utiliser *Kmerator* pour générer des **contigs uniquement composés de kmers spécifiques**.
- Ainsi, la quantification porte uniquement sur ce qui est propre au transcrit.

#### Valeur biologique

- Méthode **quantitative**, comparable aux approches classiques (ex. *Kallisto*).
- Elle s’aligne également avec la logique PCR : ciblage de séquences courtes (~20 nt).
- Le **choix de la taille des kmers** est donc crucial et peut s’inspirer de cette logique.

---

### Fonctionnement de Transipédia

- L’utilisateur interroge Transipédia avec des fichiers FASTQ.
- En arrière-plan :
  - des **index spécifiques** sont disponibles selon les bases,
  - *Kmerator* est utilisé pour extraire les **kmers spécifiques d’un gène ou transcrit**.
- L’accès aux données est variable selon la ressource interrogée.

---

