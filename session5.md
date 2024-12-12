# Session 5 - Journée du 13/12


vous allez vous mettre par groupe de travail et me préparer une figure comme pour une publication sur un élément particulier de notre communauté artificielle.

A vous de jouer !!!!!

Soyez imaginatifs !

Racontez-nous une histoire sur cet échantillon et son contenu !!

vous avez plusieurs sujets que vous pouvez traiter :

- 1 - Phage PhiKZ

- 2 - Phage Spp1

- 3 - gènes AMR

- 4 - les plasmides 



Vous avez internet, vous pouvez chercher des publications sur ce sujet

Je peux également vous fournir des données ou faire tourner des analyses si besoin ;)


# reconstruire des genomes à partir des cartes de contact 

le HiC permet de scaffolder des génomes , c'est à dire d'orientier de reorganiser les contigs les uns par rapport aux autres.

il existe une commande de metator pour essayer ce genre d'approches:

```sh
metator scaffold -b MetaTOR_21_1 -i matrices/MetaTOR_21_1.fa -o matrices/MetaTOR_21_1_reassembly.fa -O matrices/MetaTOR_21_1_reassembly.frag.tsv -t 4 metator_final/alignment_sorted.pairs.gz
```

vous devrez ensuite reconstruire une matrice à partir de ce nouveau génome:

```sh
hicstuff pipeline --enzyme "DpnII,HinfI" --binning 1000 --threads 8 --genome matrices/MetaTOR_21_1_reassembly.fa -o test_reassembly/ fastq/lib7_3C_for.fastq.gz fastq/lib7_3C_rev.fastq.gz
```

vous aurez ensuite besoin des coordonées en format UCS du scaffold que vous souhaitez mettre sous forme de matrice.
pour cela il faut utiliser la commande suivante:

```sh
cooler dump -t chroms test_reassembly/abs_fragments_contacts_weighted.mcool::/resolutions/5000
```

vous pouvez ensuite afficher la matrice en utilisant le programme cooler:

```sh
cooler show -b -s log2 -o test.png test_reassembly/abs_fragments_contacts_weighted.mcool::/resolutions/10000 metator_scaffold_0002:1-4873044
```

