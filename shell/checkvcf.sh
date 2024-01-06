VCF=/home/bloodmark/workarea/Allsample_merged.vcf.gz
genome=/media/bloodmark/SS_genome_assembly/pubref_genome_files/GCF_000002335.3_Tcas5.2_genomic.fasta

python /home/bloodmark/Tools/checkVCF/checkVCF.py -r $genome -o Allsample_merged $VCF
