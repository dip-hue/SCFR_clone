#########################################################################################################################
#Create the required folders
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
mkdir -p chrs/"$species"
mkdir -p SCFR/"$species"
mkdir -p SCFR/gaps/"$species"
mkdir -p genes/"$species"
mkdir -p GC/"$species"
mkdir -p SCFR_all
mkdir -p scripts
done


#########################################################################################################################
##Download all the seven primate genome
#Download the T2T-CHM13v2.0 genome
for chr in `cut -f 9 genome_reports/GCA_009914755.4_human.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/human/"$chr".fasta
done
#Download the NHGRI_mPanTro3-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_028858775.3_chimpanzee.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/chimpanzee/"$chr".fasta
done
#Download the NHGRI_mPanPan1-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_029289425.3_bonobo.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/bonobo/"$chr".fasta
done
#Download the NHGRI_mSymSyn1-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_028878055.3_gibbon.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/gibbon/"$chr".fasta
done
#Download the NHGRI_mGorGor1-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_029281585.3_gorilla.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/gorilla/"$chr".fasta
done
#Download the NHGRI_mPonPyg2-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_028885625.3_Bornean_orangutan.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/borangutan/"$chr".fasta
done
#Download the NHGRI_mPonAbe1-v2.1_pri genome
for chr in `cut -f 9 genome_reports/GCA_028885655.3_sumatran_orangutan.tsv|grep -v "RefSeq"`
do
echo $chr
esearch -db nucleotide -query "$chr" | efetch -format fasta > chrs/sorangutan/"$chr".fasta
done
#########################################################################################################################
#List all the SCFRs in the genomes of the 7 primate species
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
for chr in `ls -1 chrs/$species/*.fasta|cut -f 3 -d '/'|sed 's/\.fasta//g'`
do
echo $chr
python scripts/find_stop_codon_free_regions_with_reverse_gap_report.py chrs/$species/"$chr".fasta --gaps SCFR/gaps/"$species"/gap_regions_"$chr".bed > SCFR/"$species"/"$chr".fasta.SCFRs.out
done
done
#########################################################
#collate all gaps
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
cat SCFR/gaps/"$species"/gap*.bed|sort -k1,1 -k2n,2 > SCFR_all/gaps_"$species".bed
done
#########################################################################################################################
#get summary of all SCFRs in the genomes of the 7 primate species
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
cat SCFR/"$species"/*.fasta.SCFRs.out|sort -k1,1 -k2n,2|bedtools intersect -v -a stdin -b SCFR_all/gaps_"$species".bed > SCFR_all/"$species"_SCFR_all.out
Rscript scripts/summarize_SCFR_bed_frames_all.R SCFR_all/"$species"_SCFR_all.out
done
#########################################################################################################################
#get GC content of the SCFRs
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
for chr in `ls -1 chrs/$species/*.fasta|cut -f 3 -d '/'|sed 's/\.fasta//g'`
do
echo $chr
bedtools nuc -fi chrs/$species/"$chr".fasta -bed SCFR/"$species"/"$chr".fasta.SCFRs.out > GC/"$species"/"$chr".fasta.SCFRs_GC.out
done
cat GC/"$species"/*.out|awk '$11<1{print $0}' > SCFR_all/"$species"_SCFR_GC_all.out
cat SCFR_all/"$species"_SCFR_GC_all.out|awk '$13>10000{print $0}'|grep -v "^#"|sed 's/:/\t/g'|cut -f 1-13|sort -k1,1 -k2n,2 > SCFR_all/"$species"_long_SCFRs.bed
done
#########################################################################################################################
##Download all the seven primate genome annotation files
cd genes/human
#Download the gene annotation file of T2T-CHM13v2.0 genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
cd ../../
cd genes/chimpanzee
#Download the NHGRI_mPanTro3-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/bonobo
#Download the NHGRI_mPanPan1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/gibbon
#Download the NHGRI_mSymSyn1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf.gz
cd ../../
cd genes/gorilla
#Download the NHGRI_mGorGor1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gtf.gz
cd ../../
cd genes/borangutan
#Download the NHGRI_mPonPyg2-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gtf.gz
cd ../../
cd genes/sorangutan
#Download the NHGRI_mPonAbe1-v2.1_pri genome annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gtf.gz
cd ../../
#########################################################################################################################
#Extract coding exons annotated in the GTF
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
gunzip genes/"$species"/*.gtf.gz
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
#Extract coding exon coordinates from GTF
python scripts/gtf_to_bed.py genes/"$species"/"$gtf".gtf genes/"$species"/"$gtf".bed
done
#########################################################################################################################
#Plot a 2 dimensional histogram of length vs AT content and label the SCFR longer than 10 Kb that overlap coding exons
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo $species
gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
bedtools intersect -a SCFR_all/"$species"_long_SCFRs.bed -b genes/"$species"/"$gtf".bed -wa -wb|cut -f 5,13,17|sort -u > SCFR_all/"$species"_genes_of_interest.txt
Rscript scripts/plot_2dhist_seq_len_vs_pctAT.R SCFR_all/"$species"_SCFR_GC_all.out SCFR_all/"$species"_genes_of_interest.txt
done
#########################################################################################################################
