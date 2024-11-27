# Session 4 : journée du 12/12


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
bash scripts/binning_analysis.sh /pasteur/homes/mmarbout/Desktop/TP_Meta3C/
```
</p>
</details>

vous savez maintenant evaluer la qualité de vos bins en utilisant MiComplete... Félicitations, vous pouvez mainteant passer à la suite !!

il demeure dans votre binning des MAGs très contaminés. avez vous une idée pour décontaminer ces MAGs?

vous trouverez sur l'espace GAIA, le répertoire de sortie de MetaTOR après un cycle complet du programme, copiez le:

```sh
scp -r votrelogin@@sftpcampus.pasteur.fr:/pasteur/gaia/projets/p01/Enseignements/GAIA_ENSEIGNEMENTS/AdG_2023_2024/TP_Meta/metator_final/  binning/
```

vous pouvez explorer le repertoire de sortie:

```sh
ls -l binning/metator_final/
```

refaites tourner miComplete sur cet output et les fichiers FastA des MAGs [final_bin/] afin de répondre à la question 36 ;)

Q : faites une analyse comparative de votre binning avec celui obtenu après l'étape recursive.

###Analyse des bins obtenus

pour cette session, nous allons travailler avec les données obtenues après le pipeline complet de MetaTOR (partitionnement itératif et récursif). Nous allons analyser les bins que nous avons obtenus. 

vous pouvez explorer le repertoire de sortie

```sh
ls -l binning/metator_final/
```
il contient notamment 2 fichiers conteant les résultats finaux sur les contigs et les bins obtenus:


**Contig_data_final.txt**
 : il s'agit du fichier contenant l'ensemble des informations sur les contigs issues de l'assemblage et sur les bins auxquels ils appartiennent.
 
|ID|Name|Size|GC_content|Hit|Shotgun_coverage|Restriction_site|Core_bin_ID|Core_bin_contigs|Core_bin_size|Overlapping_bin_ID|Overlapping_bin_contigs|Overlapping_bin_size|Recursive_bin_ID|Recursive_bin_contigs|Recursive_bin_size|Final_bin|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|1|NODE_1|642311|38.6876450815882|3837|41.1565|2006|1|65|2175226|1|396|6322353|1|52|2158803|MetaTOR_1_1|
|2|NODE_2|576356|30.235826468363303|1724|24.509|1256|2|40|1735419|2|401|735419|0|-|-|MetaTOR_2_0|
|3|NODE_3|540571|42.305266098255366|2188|14.5855|3405|3|127|6409484|3|431|13615480|1|112|6385126|MetaTOR_3_1|

**Bin_summary.txt**
 : il s'agit du fichier contenant les informations sur les bins générés avec des données sur toutes les étapes du processus de partitionnement.(The HiC coverage is the number of contacts (intra and inter contigs) per kilobase in the whole bin. The Shotgun coverage is the mean coverage normalized by the size of the shotgun reads from the depth file.)

||lineage|completness|contamination|size|contigs|N50|longest_contig|GC|coding_density|taxonomy|Coverage|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaTOR_8_1|o__Clostridiales|68.29|2.46|1431612|15|116129|291620|26.36|87.97|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales|146.46719755483332|
|MetaTOR_8_2|o__Clostridiales|58.42|2.01|1396934|58|41290|174682|28.89|83.70|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales|22.252416224710686|
|MetaTOR_8_3|o__Clostridiales|49.37|0.94|1420821|82|33095|89964|30.29|83.24|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae_3;g__Clostridium_3|44.27369196532141|


Nous allons utiliser ces deux fichiers pour analyser notre communauté.

### Répartition taxonomique des MAGs (Metagenomic Assembled Genomes)

les bins de grande taille et de bonne qualité sont dénommés des MAGs pour Metagenomic Assembled Genomes. A l'aide du fichier bin_summary.txt, vous allez analyser la répartition taxonomique des MAGs que nous avons recontruit. Pour cela, vous pourrez vous inspirer des graphs présentés ci-dessous. Faites la même analyse en prenant en compte l'abondance des MAGs dans la communauté.

![outMAG](docs/images/MAG4.png)

![outMAG](docs/images/MAG5.png)

### Couverture et contenu en GC

Une autre façon d'analyser la diversité de notre communauté microbienne est de regarder la distribution de leur couverture et de leur contenu en GC.

à l'aide des données du fichier contig_data_final.txt, générez des graphs similaires à ceux ci-dessous (il vous faudra utiliser la fonction boxplot de R)

dans la figure ci-dessous, les boxplot sont colorés en fonction de la qualité des MAGs ... si vous arrivez à le faire ... bravo !!

![outMAG](docs/images/outMAG3.png)

![outMAG](docs/images/outMAG4.png)

si vous avez un peu de mal ... vous pouvez jeter un oeil au script GC_cov_analysis dans le dossier scripts/.

en explorant le script vous devriez être en mesure de le lancer avec les bons arguments ;)


### Analyse de bin unique

Il est également possible de générer des « density plot » pour chaque bin afin de vérifier leur homogénéité ou au contraire voir si il y a différentes populations de contigs.

![outMAG](docs/images/outMAG8.png)

pour celui là , je vais vous filer un coup de pouce ... il va fallor lancer le scripts bin_analysis.sh qui se trouve dans le dossier scripts/

lancement du script bin_analysis.sh qui prends 3 arguments en entrée [1-targeted_bin; 2-output_directory; 3-contig_data_file from MetaTOR]

```sh
bash scripts/bin_analysis MetaTOR_22_2 figures/ binning/metator_final/contig_data_final.txt
```

### Matrices d’interactions

A partir de n'importe quel réseau ou fichier d’alignement, il est possible de générer une matrice qui est une méthode de visualisation de graphe.

Pour cela, nous allons utiliser une fonction de notre programme MetaTOR qui permet de générer des matrices d'interactions pour différents "objets" (contig, core_bin, overlapping, bin, final_bin)

```sh
metator contactmap -h
```
on peut ainsi générer des matrices d'interactions pour différents "objets". Nous allons commencer par générer la matrice du MAG MetaTOR_22_2

nous allons commencer par créer un répertoire de sortie

```sh
mkdir -p matrices/
```

nous allons ensuite lancer la commande suivante:

```sh
metator contactmap -t 8 -a assemblage/assembly_all.fa -c binning/metator_final/contig_data_final.txt -e DpnII,HinfI -n "MetaTOR_22_2" -f -o matrices/MetaTOR_22_2/ -O "final_bin" binning/metator_final/alignment_sorted.pairs.gz
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
hicstuff view -n -b 10kb -f matrices/MetaTOR_22_2/MetaTOR_22_2.frags.tsv -o matrices/MetaTOR_22_2/mat_10kb_norm.pdf matrices/MetaTOR_22_2/MetaTOR_22_2.mat.tsv
```

lorsque vous utilisez cette commande, faites bien attention à la taille de vos bins (taille d'un pixel, i.e. l'option -b)

vous devriez obtenir une figure comme celle ci: 

![outMAG](docs/images/outMAG9.png)


il est également possible de faire cela pour des contigs. En vous servant de l'aide de la fonction "contactmap" et de tout ce que vous avez appris, générez la matrice d'interactions (bin = 5kb) du contig de plus de 100kb le plus couvert de vos données.






