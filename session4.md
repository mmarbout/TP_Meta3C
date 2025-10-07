# Session 4


## Validation des bins obtenus

Différentes approches permettent de valider les bins obtenus. Nous allons essayer de passer en revue différentes méthodes permettant de valider notre binning.


###	Marqueurs taxonomiques

![checkM](docs/images/checkm.png)

Différents programmes existent afin de valider les bins obtenus après partitionnement d'un métagénome (CheckM, Micomplete). Dans notre cas nous utiliserons Micomplete qui est moins complet mais beacoup moins gourmand que CheckM.

La validation des bins avec ce programme consiste à rechercher un set de gènes bactériens (via des modèles HMM), essentiels et présents en une seule copie dans plus de 97% des génomes bactériens connus.

L’absence/présence et la multiplicité de ces marqueurs permettent ainsi d’évaluer : 

i - la complétude (mesure reliée au nombre de marqueurs au sein d'un bin par rapport au nombre attendu.

ii - la contamination (mesure reliée au nombre de marqueurs en plusieurs copies).

![micomplete](docs/images/micomplete.png)


pour lancer micomplete, il faut d'abord changer les extensions des bins (.fa --> .fna)

```sh
for f in binning/metator/overlapping_bin/*.fa ; do mv $f  binning/metator/overlapping_bin/`basename $f .fa`.fna ;done
```


il faut ensuite construire un fichier nécessaire au fonctionnement de micomplete

```sh
find binning/metator/overlapping_bin/ -maxdepth 1 -type f -name "*.fna" | miCompletelist.sh > binning/metator/overlapping_bin/listbins.tab
```

on peut ensuite lancer l'analyse (5 - 10 min):

```sh
miComplete binning/metator/overlapping_bin/listbins.tab --threads 8 --hmms Bact105 -o binning/metator/miComplete.txt 
```

et jeter un oeil aux résultats

```sh
cat binning/metator/miComplete.txt | head
```

Nous considèrerons un génome complet quand :

o	sa complétude se situe au-delà de 90% (>0.9)

o	sa contamination se situe en deçà de 5% (<1.05)

un génome HQ (high-quality):

o	sa complétude se situe au-delà de 90% (>0.9)

o	sa contamination se situe en deçà de 10% (<1.1)

un génome MQ (medium-quality):

o	sa complétude se situe au-delà de 70% (>0.7)

o	sa contamination se situe en deçà de 10% (<1.1)

un génome LQ (low-quality):

o	sa complétude se situe au-delà de 50% (>0.5)

o	sa contamination se situe en deçà de 10% (<1.1)

un génome PQ (poor-quality):

o	sa complétude se situe en deçà de 50% (<0.5)

o	sa contamination se situe en deçà de 10% (<1.1)

un génome Contaminé:

o	sa contamination se situe au delà de 10% (>1.1)


Q : Combien de génome(s) reconstruit(s) et complet(s) avez-vous ? Quelle proportion en terme de séquence cela représente t il ?

![outMAG](docs/images/MAG1.png)

![outMAG](docs/images/MAG2.png)

![outMAG](docs/images/MAG3.png)


<details><summary>Solution</summary>
<p>

```sh
bash scripts/micomplete_analysis.sh binning/metator/contig_data_partition.txt binning/metator/miComplete.txt figures/
```
</p>
</details>

vous savez maintenant evaluer la qualité de vos bins en utilisant MiComplete... Félicitations, vous pouvez mainteant faire la même chose avec le deuxième output de metator pour lequel nous avons réalisé 40 itérations avec un seuil de 80% et faire une comparaison des résultats obtenus.

<details><summary>Solution</summary>
<p>

```sh
for f in binning/metator_v2/overlapping_bin/*.fa ; do mv $f  binning/metator_v2/overlapping_bin/`basename $f .fa`.fna ;done
find binning/metator_v2/overlapping_bin/ -maxdepth 1 -type f -name "*.fna" | miCompletelist.sh > binning/metator_v2/overlapping_bin/listbins.tab
miComplete binning/metator_v2/overlapping_bin/listbins.tab --threads 8 --hmms Bact105 -o binning/metator_v2/miComplete.txt
```
</p>
</details>


### Louvain recursif

Il demeure dans votre binning des MAGs très contaminés. avez vous une idée pour décontaminer ces MAGs?

En effet, l'algorithme de louvain présente des limites de résolution quand il est appliqué a de grandss graphs. L'un des possibilité est d'isoler chaque sous-réseau d'interactions correspondants à des bins contaminés et de refaire tourner l'algorithme uniquement sur ces sous-réseaux.

la commande de metator pour cette étape est la commande validation 

```sh
metator validation -h
```

vous pouvez lancer la commande pour la prodédure récursive de cette manière:


```sh
for f in binning/metator_v2/overlapping_bin/*.fna ; do mv $f  binning/metator_v2/overlapping_bin/`basename $f .fna`.fa ;done
metator validation -a assemblage/assembly_all.fa -c binning/metator_v2/contig_data_partition.txt -f binning/metator_v2/overlapping_bin/ -n binning/metator/network.txt -o binning/metator_v2/ -t 4
```
vous pouvez explorer le repertoire de sortie:

refaites tourner miComplete sur cet output et les fichiers FastA des MAGs [final_bin/] afin de répondre à la question suivante ;)

Q : faites une analyse comparative de vos différents binning.

# ATTENTION PARTIE EN CHANTIER !!!

## Analyse des bins obtenus

pour cette session, nous allons travailler avec les données obtenues après le pipeline complet de MetaTOR (partitionnement itératif et récursif et évaluation des bins par checkM et non par miComplete). Nous allons analyser les bins que nous avons obtenus. 

vous pouvez explorer le repertoire de sortie

```sh
ls -l binning/metator_final/
```
il contient les mêmes fichiers que ceux que vous avez obtenus plus les fichiers de sortie de checkM:
- bin_summary.txt
- contig_data_final.txt
- checkM_results_complete.txt
- checkM_taxonomy_output.txt

Nous allons utiliser ces fichiers pour analyser notre communauté.

### Répartition taxonomique des MAGs (Metagenomic Assembled Genomes)

les bins de grande taille et de bonne qualité sont dénommés des MAGs pour Metagenomic Assembled Genomes. A l'aide du fichier bin_summary.txt, vous allez analyser la répartition taxonomique des MAGs que nous avons recontruit. Pour cela, vous pourrez vous inspirer des graphs présentés ci-dessous. Faites la même analyse en prenant en compte l'abondance des MAGs dans la communauté.

![outMAG](docs/images/MAG4.png)

![outMAG](docs/images/MAG5.png)


<details><summary>Solution</summary>
<p>

```sh
bash scripts/MAGs_repartition.sh metator_final/bin_summary.txt metator_final/checkM_taxonomy_output.txt figures/
```
</p>
</details>

### Couverture et contenu en GC

Une autre façon d'analyser la diversité de notre communauté microbienne mais également de voir les relations entre complétion/contamination et couverture/GC est de regarder la distribution de leur couverture et de leur contenu en GC.

à l'aide des données du fichier contig_data_final.txt, générez des graphs similaires à ceux ci-dessous (il vous faudra utiliser la fonction boxplot de R)

dans la figure ci-dessous, les boxplot sont colorés en fonction de la qualité des MAGs ... si vous arrivez à le faire ... bravo !!

![outMAG](docs/images/outMAG3.png)

![outMAG](docs/images/outMAG4.png)

si vous avez un peu de mal ... vous pouvez jeter un oeil au script GC_cov_analysis dans le dossier scripts/.

en explorant le script vous devriez être en mesure de le lancer avec les bons arguments ;)



<details><summary>Solution</summary>
<p>

```sh
bash scripts/GC_cov_analysis.sh metator_final/contig_data_final.txt metator_final/checkM_results_complete.txt figures/
```
</p>
</details>


### Analyse de bin unique

Il est également possible de générer des « density plot » pour chaque bin afin de vérifier leur homogénéité ou au contraire voir si il y a différentes populations de contigs.

![outMAG](docs/images/outMAG8.png)

pour celui là , je vais vous filer un coup de pouce ... il va fallor lancer le scripts bin_analysis.sh qui se trouve dans le dossier scripts/

lancement du script bin_analysis.sh qui prends 3 arguments en entrée [1-targeted_bin; 2-output_directory; 3-contig_data_file from MetaTOR]

```sh
bash scripts/bin_analysis.sh MetaTOR_21_1 figures/ metator_final/contig_data_final.txt
```

# FIN DE CHANTIER

### Matrices d’interactions

A partir de n'importe quel réseau ou fichier d’alignement, il est possible de générer une matrice qui est une méthode de visualisation de graphe.

Pour cela, nous allons utiliser une fonction de notre programme MetaTOR qui permet de générer des matrices d'interactions pour différents "objets" (contig, core_bin, overlapping, bin, final_bin)

```sh
metator contactmap -h
```
on peut ainsi générer des matrices d'interactions pour différents "objets". Nous allons commencer par générer la matrice du MAG MetaTOR_21_1

nous allons commencer par créer un répertoire de sortie

```sh
mkdir -p matrices/
```

nous allons ensuite lancer la commande suivante:

```sh
metator contactmap -t 8 -a assemblage/assembly_all.fa -c metator_final/contig_data_final.txt -e DpnII,HinfI -n "metator_00020_00000" -f -o matrices/ -O "final_bin" metator_final/alignment_0_sorted.pairs.gz
```
ce script génère uniquement la matrice au format txt qui est ensuite utilisable via le programme hicstuff qui est notre logiciel de traitement des matrices d'interactions.

```sh
hicstuff -h
```

nous allons utiliser la fonction "view"

```sh
hicstuff view -h
```
vous pouvez maintenant lancer la commande suivante:

```sh
hicstuff view -n -b 10kb -f matrices/metator_00020_00000.frags.tsv -o matrices/metator_00020_00000_10kb_norm.pdf matrices/metator_00020_00000.mat.tsv
```

lorsque vous utilisez cette commande, faites bien attention à la taille de vos bins (taille d'un pixel, i.e. l'option -b)

vous devriez obtenir une figure comme celle ci: 

![outMAG](docs/images/outMAG9.png)


il est également possible de faire cela pour des contigs. En vous servant de l'aide de la fonction "contactmap" et de tout ce que vous avez appris, générez la matrice d'interactions (bin = 1kb) du contig NODE_21821_len_26639.

que pouvez vous me dire sur ce contig ?


<details><summary>Solution</summary>
<p>

```sh
metator contactmap -t 8 -a assemblage/assembly_all.fa -c metator_final/contig_data_final.txt -e DpnII,HinfI -n "NODE_21821_len_26639" -f -o matrices/ -O "contig" metator_final/alignment_0_sorted.pairs.gz
hicstuff view -n -b 1kb -f matrices/NODE_21821_len_26639.frags.tsv -o matrices/NODE_21821_len_26639_1kb_norm.pdf matrices/NODE_21821_len_26639.mat.tsv
```
</p>
</details>


Pendant le temps qu'il nous reste, faites moi les matrices du plus grand contig, du contig le plus abondant, du MAG le plus abondant, du MAG le plus contaminé ... 





