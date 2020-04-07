# generate matched input abundance (ab) and sequence (seq) files if needed
perl matchInputs.pl ab seq

# generate the workflow scripts
dir=`pwd`

fa=$dir/seq.fa
ab=$dir/ab.tsv
out=$dir/outdir

perl ../bin/oba.pl \
-fas $fa \
-exp $ab \
-out $out \
-group A1/A2/A3,B1/B2/B3,C1/C2/C3,D1/D2/D3,E1/E2/E3,F1/F2/F3 \
-name A,B,C,D,E,F \
-comp 2:1,3:1,4:1,5:1,6:1 \
-fc 0.5 \
-fdr 0.2 \
-spe hsa

# run this workflow 
nohup sh $out/OneStep.sh &
