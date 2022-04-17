# April 15, 2022
# Upload ACC data


# target syn29442253

ssh dataxfer

DIR=/sc/arion/projects/CommonMind/data/symlinks_for_mondale/RNAseq/*ACC*/Processed/RAPiD/

TMP=/sc/arion/scratch/hoffmg01/tmp
mkdir $TMP/featureCounts $TMP/kallisto $TMP/rsem $TMP/qc_metrics $TMP/salmon $TMP//star

# syn29442509
\ls $DIR/bams/*bam | parallel -P1 synapse add --parentid syn29442509 {}

# syn29443110
\ls $DIR/featureCounts/* | parallel -P1 cp {} $TMP/featureCounts/
\ls $TMP/featureCounts/* | parallel gzip
\ls $TMP/featureCounts/* | parallel -P1 synapse add --parentid syn29443110 {}

# syn29442519
for FILE in $(\ls $DIR/kallisto/abundance.tsv);
do
	cp $FILE $TMP/kallisto/$(echo $FILE | cut -d'/' -f9).abundance.tsv
done
\ls $TMP/kallisto/* | parallel gzip
\ls $TMP/kallisto/* | parallel -P1 synapse add --parentid syn29442519 {}

# syn29442526
\ls $DIR/rsem/*.results | parallel -P1 cp {} $TMP/rsem/
\ls $TMP/rsem/* | parallel gzip
\ls $TMP/rsem/* | parallel -P1 synapse add --parentid syn29442526 {}

# syn29442528
\ls $DIR/qc_metrics/*.RNASeqMetrics | parallel -P1 cp {} $TMP/qc_metrics/
\ls $TMP/qc_metrics/* | parallel gzip
\ls $TMP/qc_metrics/* | parallel -P1 synapse add --parentid syn29442528 {}

# syn29442530
\ls $DIR/salmon/*.none.quant.sf | grep -v txrevise | parallel -P1 cp {} $TMP/salmon/
\ls $TMP/salmon/* | parallel gzip
\ls $TMP/salmon/* | parallel -P1 synapse add --parentid syn29442530 {}

# syn29442534
\ls $DIR/star/*.Log.final.out | grep -v txrevise | parallel -P1 cp {} $TMP/star/
\ls $TMP/star/* | parallel gzip
\ls $TMP/star/* | parallel -P1 synapse add --parentid syn29442534 {}




