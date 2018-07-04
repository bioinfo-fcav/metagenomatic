#!/bin/bash
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo 
#

### Assembly based directory
input="$1"

### Assembly name
asmname="$2"

### Threads
threads=14

### Memory
# 10000 MB
memlimit=10000

# Minimum contig length for gene finding process
mincontiglen=500

# as linhas que iniciam com cerquilha são comentários
        
if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit 
else    
        if [ ! -d ${input} ]
        then
                echo "Wrong input directory ${input}"
                exit
        fi
fi

if [ ! ${asmname} ]
then   
	asmname="testasm"
fi

declare -A taxonomy=( ["s"]="species" ["g"]="genus" ["f"]="family" ["o"]="order" ["c"]="class" ["p"]="phylum" ["k"]="kingdom" ["d"]="domain" )

echo "Mapping data for ${input}..."

echo "   Indexing contig \"bowtie2-build\" ..."

curdir=`pwd`

cd ${input}

mkdir -p ./mapped

cd ./mapped

### Criando atalho

if [ ! -e final-assembly.fa ]; then
	ln -s ../assembled/final-assembly.fa
fi

### Indexação da montagem

bowtie2-build -f final-assembly.fa metagenome > bowtie2-build.out.txt 2> bowtie2-build.err.txt

echo "   Alignment reads X contigs \"bowtie2\" ..."

samples=()
contig_counts=()

rm -f contig_names.txt
rm -f header.txt
rm -f count_matrix.tab

for i in ../processed/prinseq/*_1.fastq; do
	
	name=`basename ${i} _1.fastq`
	sampname=`echo ${name} | sed 's/\..*//'`
	dirname=`dirname ${i}`
	
	echo "      ${sampname}"
	
	Uparam=""
	samples=($(printf "%s\n" ${samples[@]} ${sampname} | sort -u ))
	if [[ ( -e "${dirname}/${name}_1_singletons.fastq") || ( -e "${dirname}/${name}_2_singletons.fastq") ]]; then
		if [ -e "${dirname}/${name}_1_singletons.fastq" ]; then
			Uparam=" -U ${dirname}/${name}_1_singletons.fastq"
		fi
		if [ -e "${dirname}/${name}_2_singletons.fastq" ]; then
			Uparam="${Uparam} -U ${dirname}/${name}_2_singletons.fastq"
		fi
	fi

	if (( $(echo "$(samtools --version | grep '^samtools' | sed 's/^samtools //' | cut -d . -f 1,2) >= 1.3" |bc -l) )); then
		bowtie2 --very-sensitive --all -p ${threads} -q -x metagenome -1 ${dirname}/${name}_1.fastq -2 ${dirname}/${name}_2.fastq ${Uparam} 2> ${sampname}.bowtie2.err.txt | samtools view -q 10 -Su - 2> /dev/null |  samtools sort -n -o ${sampname}.metagenome.bam - 2> /dev/null
	else 
		# LEGACY samtools sort
		bowtie2 --very-sensitive --all -p ${threads} -q -x metagenome -1 ${dirname}/${name}_1.fastq -2 ${dirname}/${name}_2.fastq ${Uparam} 2> ${sampname}.bowtie2.err.txt | samtools view -q 10 -Su - 2> /dev/null |  samtools sort - ${sampname}.metagenome 2> /dev/null
	fi

	express --output-dir ./${sampname} --max-read-len 500 final-assembly.fa ${sampname}.metagenome.bam > ${sampname}.express.log.out.txt 2> ${sampname}.express.log.err.txt 
	
	if [ ! -e contig_names.txt ]; then
		cut -f 2 ./${sampname}/results.xprs | sed '1d' | sort -k1,1 -t$'\t' > contig_names.txt
	fi
	
	cut -f 2,8 ./${sampname}/results.xprs | sed '1d' | sort -k1,1 -t$'\t' | cut -f 2 > ${sampname}.counts.txt
	contig_counts=($(printf "%s\n" ${contig_counts[@]} "${sampname}.counts.txt" | sort -u ))
done

echo "   Accounting contigs abundance in samples ..."

header=("contig_name" ${samples[@]})

IFS=$'\t';echo "${header[*]}" > header.txt;IFS=$' \t\n'
paste contig_names.txt ${contig_counts[*]} | cat header.txt - > count_matrix.tab

echo "   Annotating OTUs in contig data \"kraken\" ..."

kraken --preload ./final-assembly.fa > ./sequences.kraken.out.txt 2> ./sequences.kraken.err.txt
kraken-translate --mpa-format ./sequences.kraken.out.txt > ./sequences.kraken.out.labels.txt 2> ./sequences.kraken.err.labels.txt

mergeSimpleAnnot.pl -t count_matrix.tab -k contig_name -i sequences.kraken.out.labels.txt -c 1 -a 2 -n taxonomy > count_matrix_tax_kraken.tab 2> count_matrix_tax_kraken.err.txt

# rank is the relative level of a group of organisms (a taxon) in a taxonomic hierarchy
for rank in 's' 'g' 'f' 'o' 'c' 'p' 'k' 'd'; do
       echo "Count OTUs using rank ${rank} (${taxonomy[${rank}]}) ..."
       countOTUs.pl -i count_matrix_tax_kraken.tab -c taxonomy -r ${rank}  > count_matrix_OTU_${taxonomy[${rank}]}_number.tab 2> count_matrix_OTU_${taxonomy[${rank}]}_number.err.txt
       countOTUs.pl -i count_matrix_tax_kraken.tab -c taxonomy -r ${rank} -n 100 > count_matrix_OTU_${taxonomy[${rank}]}_percent.tab 2> count_matrix_OTU_${taxonomy[${rank}]}_percent.err.txt
done

cd ../

echo "   Annotating genes in contig data \"prokka\" ..."

mkdir -p ./annotated

prokka ./assembled/final-assembly.fa --mincontiglen ${mincontiglen} --force --cpus ${threads} --outdir ./annotated --prefix ${asmname} --metagenome > ./annotated/prokka.out.txt 2> ./annotated/prokka.err.txt

cd ./annotated


mkdir -p ./align

cd ./align

if [ ! -e ${asmname}.ffn ]; then
	ln -s ../${asmname}.ffn
fi

echo "   Clustering gene predictions \"cd-hit-est\" ..."

### Agrupamento de genes

cd-hit-est -i ${asmname}.ffn -o ${asmname}.cdhit.ffn -c 0.95 -n 3 -l 10 -aS 0.9 -d 0 -B 0 -p 1 -g 1 -T ${threads} -M ${memlimit}  > cd-hit-est.out.txt 2> cd-hit-est.err.txt

cdhit-cluster-consensus -clustfile=${asmname}.cdhit.ffn.clstr -fastafile=${asmname}.ffn maxlen=1 -output=${asmname}.cdhit.consensus


### Indexação da montagem

echo "   Indexing gene predictions \"bowtie2-build\" ..."

bowtie2-build -f ${asmname}.cdhit.consensus.fasta  ${asmname} > bowtie2-build.out.txt 2> bowtie2-build.err.txt

echo "   Alignment reads X gene predictions \"bowtie2\" ..."

samples=()
gene_counts=()

rm -f gene_names.txt
rm -f header.txt
rm -f count_matrix.tab

for i in ../../processed/prinseq/*_1.fastq; do
	
	name=`basename ${i} _1.fastq`
	sampname=`echo ${name} | sed 's/\..*//'`
	dirname=`dirname ${i}`
	
	echo "      ${sampname}"
	
	Uparam=""
	samples=($(printf "%s\n" ${samples[@]} ${sampname} | sort -u ))
	
	if [[ ( -e "${dirname}/${name}_1_singletons.fastq") || ( -e "${dirname}/${name}_2_singletons.fastq") ]]; then
		if [ -e "${dirname}/${name}_1_singletons.fastq" ]; then
			Uparam=" -U ${dirname}/${name}_1_singletons.fastq"
		fi
		if [ -e "${dirname}/${name}_2_singletons.fastq" ]; then
			Uparam="${Uparam} -U ${dirname}/${name}_2_singletons.fastq"
		fi
	fi
	if (( $(echo "$(samtools --version | grep '^samtools' | sed 's/^samtools //' | cut -d . -f 1,2) >= 1.3" |bc -l) )); then
		bowtie2 --very-sensitive --all -p ${threads} -q -x ${asmname} -1 ${dirname}/${name}_1.fastq -2 ${dirname}/${name}_2.fastq ${Uparam} 2> ${sampname}.bowtie2.err.txt | samtools view -q 10 -Su - 2> /dev/null |  samtools sort -n -o ${sampname}.${asmname}.bam - 2> /dev/null
	else
		bowtie2 --very-sensitive --all -p ${threads} -q -x ${asmname} -1 ${dirname}/${name}_1.fastq -2 ${dirname}/${name}_2.fastq ${Uparam} 2> ${sampname}.bowtie2.err.txt | samtools view -q 10 -Su - 2> /dev/null |  samtools sort - ${sampname}.${asmname} 2> /dev/null
	fi

	express --output-dir ./${sampname} --max-read-len 500 ${asmname}.cdhit.consensus.fasta ${sampname}.${asmname}.bam > ${sampname}.express.log.out.txt 2> ${sampname}.express.log.err.txt 
	
	if [ ! -e contig_names.txt ]; then
		cut -f 2 ./${sampname}/results.xprs | sed '1d' | sort -k1,1 -t$'\t' > gene_names.txt
	fi
	
	cut -f 2,8 ./${sampname}/results.xprs | sed '1d' | sort -k1,1 -t$'\t' | cut -f 2 > ${sampname}.counts.txt
	gene_counts=($(printf "%s\n" ${gene_counts[@]} "${sampname}.counts.txt" | sort -u ))
done

echo "   Accounting gene abundance in samples ..."

header=("gene_name" ${samples[@]})

IFS=$'\t';echo "${header[*]}" > header.txt;IFS=$' \t\n'
paste gene_names.txt ${gene_counts[*]} | cat header.txt - > count_matrix.tab

mergefaAnnot.pl -f ${asmname}.cdhit.consensus.fasta -t count_matrix.tab -c gene_name -n clstr.description > ./count_matrix_clstr_desc.tab

recvfaAnnot.pl -f ${asmname}.cdhit.consensus.fasta -c ${asmname}.cdhit.ffn.clstr -r ${asmname}.cdhit.ffn -o ${asmname}.cdhit.consensus.prokka.fasta

mergefaAnnot.pl -f ${asmname}.cdhit.consensus.prokka.fasta -t ./count_matrix_clstr_desc.tab -c gene_name -n ref.description > ../count_matrix_clstr_ref_desc.tab

cd ../
mkdir -p ./iprscan/
sed 's/\*//g' ${asmname}.faa > ${asmname}.cleaned.faa
tempdir=`mktemp -d --tmpdir=/dev/shm/`
interproscan.sh --applications Hamap,ProDom,PANTHER,SMART,PRINTS,PIRSF,Pfam --seqtype p --tempdir /dev/shm --pathways --goterms --iprlookup --output-dir ./iprscan/ --tempdir ${tempdir} --input ${asmname}.cleaned.faa
rm -fr ${tempdir}

mergeiprAnnot.pl -t ./count_matrix_clstr_ref_desc.tab -i ./iprscan/${asmname}.cleaned.faa.tsv -c gene_name -a IPR -x ./align/${asmname}.cdhit.ffn.clstr > count_matrix_desc_ipr.tab

cd ${curdir}
