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

# Diretório onde está o diretório ./processed/prinseq com os arquivos .fastq processados
input="$1"

# Expressão regular para distinguir grupos de arquivos para a montagem
re="$2"

# Nome para os contigs/scaffolds montados
asmname="$3"

# Número de threads disponíveis para os processos
threads=14
# Número de threads para o processamento dos arquivos pmap utilizando parallel
pmap_threads=2

# Número de sequências dentro de cada partição
partitions_n=100000

# http://khmer.readthedocs.org/en/v1.0/choosing-table-sizes.html
# Memória disponível = -N (n_tables) * -x (tablesize)
memlimitGB=21
memlimitMB=$((memlimitMB*1000))
# khmer parâmetro -N
# n_tables = Memória disponível em GB / tablesize
# Retirado da documentação: "Just use use -N 4, always, and vary the -x parameter."
khmer_N="4"
# khmer parâmetro -x 
# tablesize = Memória disponível em GB / n_tables * 1 bilhão
khmer_byte_x="1e9"
khmer_byte_graph_x="$((memlimitGB/${khmer_N}))e9"

# khmer parâmetro -k
khmer_k=23

# Tamanho mínimo para a avaliação das montagens, realizada pela soma dos tamanhos de todos os contigs maiors que este valor
cutoff=200
# Check read length:
# perl -lane 'INIT { our $size=0; my $c = 0; } if ($.%4==2) { last if ($c > 10); $size+=length($_); $c++; } END { print $size/$c; }' ./processed/prinseq/sampleA1.scythe.cutadapt5p.filtered.prinseq_1.fastq
read_len=150
cutoff_perc=$(echo "scale=2; (${cutoff}/${read_len})" | bc -l)

TMP_DIR="/state/partition1"

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

biosamples=()

echo "Assembling data from ${input} for ${asmname} ..."

if [ ! ${re} ]
then   
	biosamples=("*")
else    
	for i in ${input}/processed/prinseq/*_1.fastq; do
		name=`basename ${i} .fastq | sed 's/\..*//'`
		if [[ $name =~ ${re} ]]; then
			biosamples=($(printf "%s\n" ${biosamples[@]} "${BASH_REMATCH[1]}" | sort -u ))
		fi
	done
fi

echo "   Selected samples: ${biosamples[*]} "

curdir=`pwd`

cd ${input}

finalasm=()

for samp in "${biosamples[@]}"; do
	
	sampname=""
	
	if [ "${samp}" == "*" ]; then
		sampname="all"
	else
		sampname="${samp}"
	fi
	
	mkdir -p ./${sampname}
	
	cd ./${sampname}
	
	echo "   Processing sample ${sampname}"
	
	echo "   Link and interleave reads ..."
	
	mkdir -p ./reads
	
	cd ./reads
	
	anyse=0
	# s1_pe s1_se s2_pe s2_se
	for i in ../../processed/prinseq/${samp}*_1.fastq; do
		
		bn=`basename ${i} _1.fastq`
		
		simplename=`echo ${bn} | sed 's/.scythe.cutadapt5p.filtered.prinseq//'` 
		
		ln -f -s ../../processed/prinseq/${bn}_1.fastq s1_pe.fq
		ln -f -s ../../processed/prinseq/${bn}_2.fastq s2_pe.fq
		
		interleave-reads.py s?_pe.fq > ${simplename}_combined.pe.fq 2> ${simplename}_combined.pe.interleave-reads.err.txt
		
		gzip -9c ${simplename}_combined.pe.fq > ${simplename}_combined.pe.fq.gz
		
		if [ -e ../../processed/prinseq/${bn}_1_singletons.fastq ]; then
			ln -f -s ../../processed/prinseq/${bn}_1_singletons.fastq s1_se.fq
			cat s1_se.fq > ${simplename}_combined.se.fq
		fi
		if [ -e ../../processed/prinseq/${bn}_2_singletons.fastq ]; then
			ln -f -s ../../processed/prinseq/${bn}_2_singletons.fastq s2_se.fq
			cat s2_se.fq >> ${simplename}_combined.se.fq
		fi
		
		if [ -e ${simplename}_combined.se.fq ]; then
			anyse=1
			gzip -9c ${simplename}_combined.se.fq > ${simplename}_combined.se.fq.gz
		fi
		
		rm -f *.fq
	done
	
	cd ../
	
	echo "   Digital normalization - \"diginorm\" ..."
	
	mkdir -p diginorm
	
	cd ./diginorm
	
	### Normaliza tudo para uma cobertura de 20, considerando um k-mer de 20
	# PE (-p)
	normalize-by-median.py -k ${khmer_k} -C 20 -N ${khmer_N} -x ${khmer_byte_x} -p --savetable normC20k20.kh ../reads/*.pe.fq.gz > normC20k20.pe.out.txt 2> normC20k20.pe.err.txt
	
	# SE
	if [ ${anyse} ]; then
		normalize-by-median.py -k ${khmer_k} -C 20 -N ${khmer_N} -x ${khmer_byte_x} --savetable normC20k20.kh --loadtable normC20k20.kh ../reads/*.se.fq.gz > normC20k20.se.out.txt 2> normC20k20.se.err.txt
	fi
	
	### Poda leituras em k-mers que são pouco abundantes em leituras de alta cobertura. A opção -V é usada para datasets com cobertura variável.
	filter-abund.py -V normC20k20.kh *.keep > filter-abund.out.txt 2> filter-abund.err.txta
	
	### Extração de arquivos PE (final .pe) e SE (final .se)
	for i in *.pe.fq.gz.keep.abundfilt; do
		extract-paired-reads.py ${i} > ${i}.extract-paired-reads.out.txt 2> ${i}.extract-paired-reads.err.txt
	done
	
	### Após eliminar k-mers errôneos, vamos abandonar mais alguns dados de alta cobertura.
	# PE (-p)
	normalize-by-median.py -C 5 -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_x} -p --savetable normC5k20.kh *.pe.fq.gz.keep.abundfilt.pe > normC5k20.pe.out.txt 2> normC5k20.pe.err.txt
	
	# SE
	if [ ${anyse} ]; then
		normalize-by-median.py -C 5 -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_x} --savetable normC5k20.kh --loadtable normC5k20.kh *.pe.fq.gz.keep.abundfilt.se *.se.fq.gz.keep.abundfilt > normC5k20.se.out.txt 2> normC5k20.se.err.txt
	fi
	
	# Compactar cada amostra (gzip) em arquivos com nome base mais curto
	for i in `ls *.pe.fq.gz.keep.abundfilt.pe.keep`; do bn=`basename ${i} .pe.fq.gz.keep.abundfilt.pe.keep`; gzip -9c ${i} > ${bn}.pe.kak.fq.gz; done
	if [ ${anyse} ]; then
		for i in `ls *.pe.fq.gz.keep.abundfilt.se.keep`; do bn=`basename ${i} .pe.fq.gz.keep.abundfilt.se.keep`; gzip -9c ${i} > ${bn}.se.kak.fq.gz; done
		for i in `ls *.se.fq.gz.keep.abundfilt.keep`; do bn=`basename ${i} .se.fq.gz.keep.abundfilt.keep`; gzip -9c ${i} >> ${bn}.se.kak.fq.gz; done
	fi
	
	# Remover arquivos desnecessários
	rm -f normC20k20.kh *.keep *.abundfilt *.pe *.se
	
	readstats.py *.kak.fq.gz ../../processed/prinseq/${samp}*.fastq > diginorm.out.txt 2> diginorm.err.txt
	
	cd ../
	
	echo "   Partitioning ..."
	
	mkdir -p ./partitioned
	
	cd ./partitioned
	
	### Eliminação dos k-mers altamente repetitivos que podem juntar múltiplas espécies (quimeras) e renomear de forma apropriada:
	
	filter-below-abund.py ../diginorm/normC5k20.kh ../diginorm/*.fq.gz > filter-below-abund.out.txt 2> filter-below-abund.err.txt
	
	# Renomear para .below.fq
	for i in *.below; do mv ${i} ${i}.fq; done
	
	## Carrega grafo de k-mers
	load-graph.py -k ${khmer_k} -T ${threads} -N ${khmer_N} -x ${khmer_byte_graph_x} lump *.below.fq > load-graph.out.txt 2> load-graph.err.txt
	
	## Encontrando k-mers altamente conectados iniciais (possíveis artefatos)
	make-initial-stoptags.py -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_graph_x} lump > make-initial-stoptags.out.txt 2> make-initial-stoptags.err.txt
	
	## Particiona o grafo de acordo com a sua conectividade
	partition-graph.py --threads ${threads}  --stoptags lump.stoptags lump > partition-graph.out.txt 2> partition-graph.err.txt
	
	pmap_count=`ls -l lump.*.pmap | wc -l`
	
	# Evitar divisão por zero
	if [ ${pmap_count} -lt ${pmap_threads} ]; then
		pmap_count=${pmap_count}
	fi
	
	## Encontrando k-mers altamente conectados (possíveis artefatos)
	
	pmap_count=`ls -l lump.*.pmap | wc -l`
	
	
	if [ ${pmap_count} -lt ${pmap_threads} ]; then
	        # 1 em cada thread
	        pmap_limit=1
	else    
	        # (pmap_count / pmap_thread) em cada thread)
	        pmap_limit=$((pmap_count / pmap_threads))
	fi
	
	echo "      PMAP count: ${pmap_count} ${pmap_limit}"
	
	pmap_c=0
	pmap_dir=0
	rm -f ./run-find-knots.sh
	for i in lump.*.pmap; do
		if [ $((pmap_c % pmap_limit)) == 0 ]; then
			pmap_dir=$((pmap_dir + 1))
			mkdir -p ./fk.${pmap_dir}
			ln -s $(readlink -f lump.pt) ./fk.${pmap_dir}/lump.pt
			ln -s $(readlink -f lump.tagset) ./fk.${pmap_dir}/lump.tagset
			echo "(cd ./fk.${pmap_dir} && find-knots.py lump > find-knots.out.txt 2> find-knots.err.txt)" >> ./run-find-knots.sh
		fi
		echo "         Move ${i} to fk.${pmap_dir}"
		mv ${i} ./fk.${pmap_dir}/
	done
	
	echo "      Running find-knots.sh using ${pmap_threads} parallel processes ..."
	
	parallel --gnu -j ${pmap_threads} < ./run-find-knots.sh
	
	echo "      Merging all .stoptags files to merge.stoptags ..."
	
	merge-stoptags.py -k ${khmer_k} ./fk > merge-stoptags.out.txt 2> merge-stoptags.err.txt
	
	echo "      Filtering stoptags ..."
	
	### Poda sequências nos k-mers considerados artefatos
	filter-stoptags.py -k ${khmer_k} merge.stoptags *.kak.fq.gz.below.fq > filter-stoptags.out.txt 2> filter-stoptags.err.txt
	
	
	echo "      do-partitions ..."
	
	### Particionamento, gera arquivos que contêm as anotações das partições
	do-partition.py -N ${khmer_N} -k ${khmer_k} -x ${khmer_byte_graph_x} -T ${threads} kak *.kak.fq.gz.below.fq.stopfilt > do-partition.out.txt 2> do-partition.err.txt
	
	echo "      extract-partitions ..."
	
	### Extraindo partições em grupos
	extract-partitions.py -m 0 -X ${partitions_n} kak *.part > extract-partitions.out.txt 2> extract-partitions.err.txt
	
	# Extraindo arquivos PE e SE (dn - digital normalization)
	for i in kak*.fq; do 
		extract-paired-reads.py ${i} > ${i}.extract-paired-reads.out.txt 2> ${i}.extract-paired-reads.err.txt
		name=$(basename ${i} .fq)
		mv ${name}.fq.pe ${name}.dn.pe.fq
		mv ${name}.fq.se ${name}.dn.se.fq
	done
	
	mkdir -p ../input
	
	# Cópia	
	cp *.dn.?e.fq ../input
	
	# Compressão
	gzip *.dn.?e.fq
	
	echo "      sweep-reads ..."
	
	### Alocar os dados processados (sem normalização) em partições - Recupera reads baseadas no compartilhamento de k-mers
	
	sweep-files.py -k ${khmer_k} -N ${khmer_N} -x ${khmer_byte_graph_x} --db kak.group*.fq --query ../reads/*.?e.fq.gz > sweep-files.out.txt 2> sweep-files.err.txt
	
	# Extraindo arquivos PE e SE  (nodn - no digital normalization)
	for i in *kak*.sweep; do
		sweep=$(basename $i .fq)
		mv $i $sweep.fq
		extract-paired-reads.py ${sweep}.fq
		mv $sweep.fq.pe ${sweep}.nodn.pe.fq
		mv $sweep.fq.se ${sweep}.nodn.se.fq
	done
	
	# Cópia
	cp *.nodn.?e.fq ../input

	# Compressão
	gzip *.nodn.se.fq *.nodn.pe.fq
	
	# Removendo arquivos desnecessários
	#rm -f *.sh *.part *.sweep.fq *.below.fq *.stopfilt *.fq
	#rm -f ../diginorm/normC5k20.kh
	
	
	cd ../
	
	echo "   Assembling ..."
	
	mkdir -p ./assembled
	
	cd ./assembled
	
	groups=()
	
	infiles=`ls ../input/*.*.pe.fq`;
	
	for i in ${infiles[@]}; do
		
		name=$(basename $i .pe.fq);
		
		indir=$(dirname $i)
		
		echo "   [${name}] ..."
		
		if [ ! -s ${indir}/${name}.se.fq ]; then
			echo "      Put a fake read to single end data ..."
			echo -e "@SAMPLE0000000000/1\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" > ${indir}/${name}.se.fq
			echo -e "@SAMPLE0000000000/2\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" >> ${indir}/${name}.se.fq
		fi	

		echo "      Change input files' format"
		
		# minia
		cat ${indir}/${name}.se.fq | fastq2fasta.pl -n 100 > ${indir}/${name}.fa
		cat ${indir}/${name}.pe.fq | fastq2fasta.pl -n 100 >> ${indir}/${name}.fa
		
		# NEWBLER
		cat ${indir}/${name}.pe.fq             | deinterleave_pairs -o ${indir}/${name}.pe_1.fq ${indir}/${name}.pe_2.fq
		
		# idba_ud
		fastq-se2ipe.pl ${indir}/${name}.se.fq | fastq2fasta.pl -n 100 > ${indir}/${name}.se.pe.fa
		cat ${indir}/${name}.pe.fq | fastq2fasta.pl -n 100             > ${indir}/${name}.pe.fa
	
		# ABYSS - nomes idênticos em R1 e R2 com a única diferença de possuir um /1 ou /2, respectivamente (não tenho certeza se necessário).
		cat ${indir}/${name}.pe_1.fq | renameSeqs.pl -p ABYSSPE | sed 's/\(^@ABYSSPE[^ ]*\).*/\1\/1/' > ${indir}/${name}.abyss.pe_1.fq
		cat ${indir}/${name}.pe_2.fq | renameSeqs.pl -p ABYSSPE | sed 's/\(^@ABYSSPE[^ ]*\).*/\1\/2/' > ${indir}/${name}.abyss.pe_2.fq
		cat ${indir}/${name}.se.fq   | renameSeqs.pl -p ABYSSSE > ${indir}/${name}.abyss.se.fq
		
		echo "      IDBA-UD ..."
		# scaffold.fa
		idba_ud --num_threads ${threads} --read ${indir}/${name}.pe.fa --read_level_2 ${indir}/${name}.se.pe.fa --min_contig ${cutoff} --maxk 51 --step 16 --mink 19 --min_pairs 0 -o ${name}.idba_ud.0.d > ${name}.idba_ud.0.log.out.txt 2> ${name}.idba_ud.0.log.err.txt
		   
		echo "      Newbler ..."
		# 454AllContigs.fna
		newAssembly ${name}.newbler.0.d > ${name}.newbler.0.log.out.txt 2> ${name}.newbler.0.log.err.txt
		addRun -lib PE ${name}.newbler.0.d ${indir}/${name}.pe_1.fq >> ${name}.newbler.0.log.out.txt 2>> ${name}.newbler.0.log.err.txt
		addRun -lib PE ${name}.newbler.0.d ${indir}/${name}.pe_2.fq >> ${name}.newbler.0.log.out.txt 2>> ${name}.newbler.0.log.err.txt
		addRun -lib SE ${name}.newbler.0.d ${indir}/${name}.se.fq >> ${name}.newbler.0.log.out.txt 2>> ${name}.newbler.0.log.err.txt
		runProject -mi 95 -ml 20 -cpu ${threads} ${name}.newbler.0.d >> ${name}.newbler.0.log.out.txt 2>> ${name}.newbler.0.log.err.txt
		
		echo "      MEGAHIT ..."
		# final.contigs.fa
		megahit --mem-flag 2  -t ${threads} --min-contig-len ${cutoff} -o ${name}.megahit.0.d -m 0.9 --presets meta-large --12 ${indir}/${name}.pe.fq -r ${indir}/${name}.se.fq --min-count 1 --k-min 19 --k-max 51 --k-step 16 > ${name}.megahit.0.log.out.txt 2> ${name}.megahit.0.log.err.txt
		
		echo "      SPAdes ..."
		# scaffolds.fasta
		spades.py --phred-offset 33 --memory ${memlimitGB} -t ${threads} --cov-cutoff auto --12 ${indir}/${name}.pe.fq -s ${indir}/${name}.se.fq -o ${name}.spades.0.d > ${name}.spades.0.log.out.txt 2> ${name}.spades.0.log.err.txt
		
		echo "      metaSPAdes ..."
		# scaffolds.fasta
		metaspades.py --phred-offset 33 --memory ${memlimitGB} -t ${threads} --12 ${indir}/${name}.pe.fq -s ${indir}/${name}.se.fq -o ${name}.metaspades.0.d > ${name}.metaspades.0.log.out.txt 2> ${name}.metaspades.0.log.err.txt
		
		echo "      Minia ..."
		# é necessário criar o diretório antes de rodar o minia (https://www.biostars.org/p/168676/)
		# kak.contigs.fa
		
		for k in {19..51..16}; do
			echo "         ${k} ..."
			mkdir -p ./${name}.minia.$k.d
			minia -nb-cores ${threads} -kmer-size $k -max-memory ${memlimitGB} -out-dir ${name}.minia.$k.d -out ${name}.minia.$k.d/kak -in ${indir}/${name}.fa > ${name}.minia.$k.log.out.txt 2> ${name}.minia.$k.log.err.txt
		done
		
		echo "      ABySS ..."
		# kak-scaffolds.fa
		for k in {19..51..16}; do
			echo "         ${k} ..."
			mkdir -p ./${name}.abyss.$k.d
			rsep=`readlink -f ${indir}/${name}.abyss.se.fq`
			rpep1=`readlink -f ${indir}/${name}.abyss.pe_1.fq`
			rpep2=`readlink -f ${indir}/${name}.abyss.pe_2.fq`
			cd ./${name}.abyss.$k.d
			# não usar lib='pe' (pe é um argumento de abyss-pe)
			abyss-pe k=${k} name=kak j=${threads} lib='pe1' pe1="${rpep1} ${rpep2}" se="${rsep}" e=2 c=2 s=100 N=1 n=1 S=${cutoff} scaffolds > ../${name}.abyss.$k.log.out.txt 2> ../${name}.abyss.$k.log.err.txt
			cd ../
		done
		
		echo "      Velvet ..."
		# contigs.fa
		for k in {19..51..16}; do
			echo "         ${k} ..."
			velveth ${name}.velvet.$k.d $k -fastq -long ${indir}/${name}.se.fq -fastq -longPaired -interleaved ${indir}/${name}.pe.fq > ${name}.velveth.$k.log.out.txt 2> ${name}.velveth.$k.log.err.txt
			velvetg ${name}.velvet.$k.d -exp_cov auto -cov_cutoff auto -scaffolding yes -conserveLong yes -min_contig_lgth ${cutoff} > ${name}.velvetg.$k.log.out.txt 2> ${name}.velvetg.$k.log.err.txt
		done
		
		echo "      metaVelvet ..."
		# meta-velvetg.contigs.fa
		for k in {19..51..16}; do
			echo "         ${k} ..."
			mkdir -p ${name}.metavelvet.$k.d
			cd ${name}.metavelvet.$k.d/
			find ../${name}.velvet.$k.d/ -maxdepth 1 -type f -exec ln -f -s {} \;
			cd ../
			meta-velvetg ${name}.metavelvet.$k.d -exp_cov auto -cov_cutoff auto -scaffolding yes -min_contig_lgth ${cutoff} > ${name}.metavelvetg.$k.log.out.txt 2> ${name}.metavelvetg.$k.log.err.txt
		done
		
		echo "      Ray ..."
		# Scaffolds.fasta
		for k in {19..51..16}; do
			echo "         ${k} ..."
			mpiexec -n ${threads} -num-cores ${threads} Ray -minimum-seed-length 50 -k ${k} -i ${indir}/${name}.pe.fq -s ${indir}/${name}.se.pq -o ${name}.ray.$k.d > ${name}.ray.$k.log.out.txt 2> ${name}.ray.$k.log.err.txt
		done
		
		echo "      MetaPlatanus ..."
		
		mkdir -p ${name}.meta_platanus.0.d/
		# kak_finalClusters_all.fa
		# The parameters -c and -C must be increased with large datasets add 1 for each 50000000
		meta_platanus_n=$(wc -l ${indir}/${name}.pe_1.fq)
		meta_platanus_cov=$(( $(echo "scale=0; (${meta_platanus_n}/50000000)" | bc -l)+1 ))
		meta_platanus assemble -tmp ${TMP_DIR} -k 0.15 -K 0.5 -c ${meta_platanus_cov} -C ${meta_platanus_cov} -l ${cutoff_perc} -t ${threads} -m ${memlimitGB} -f ${indir}/${name}.pe_1.fq ${indir}/${name}.pe_2.fq ${indir}/${name}.se.fq -o ${name}.meta_platanus.0.d/kak > ${name}.meta_platanus-assemble.0.log.out.txt 2> ${name}.meta_platanus-assemble.0.log.err.txt
		meta_platanus scaffold -tmp ${TMP_DIR} -k ${name}.meta_platanus.0.d/kak_kmer_occ.bin -t ${threads} -c ${name}.meta_platanus.0.d/kak_contig.fa -IP1 ${indir}/${name}.pe_1.fq ${indir}/${name}.pe_2.fq -o ${name}.meta_platanus.0.d/kak > ${name}.meta_platanus-scaffold.0.log.out.txt 2> ${name}.meta_platanus-scaffold.0.log.err.txt
		cd ${name}.meta_platanus.0.d/
		meta_platanus iterate -tmp ${TMP_DIR} -m ${memlimitGB} -t ${threads} -k kak_kmer_occ.bin -c kak_scaffold.fa -IP1 ${indir}/../${name}.pe_1.fq ${indir}/../${name}.pe_2.fq -o kak > ../${name}.meta_platanus-iterate.0.log.out.txt 2> ../${name}.meta_platanus-iterate.0.log.err.txt
		meta_platanus cluster_scaffold -tmp ${TMP_DIR} -t ${threads} -c out_iterativeAssembly.fa -IP1 ${indir}/../${name}.pe_1.fq ${indir}/../${name}.pe_2.fq -o kak > ../${name}.meta_platanus-cluster_scaffold.0.log.out.txt 2> ../${name}.meta_platanus-cluster_scaffold.0.log.err.txt
		
		if [ ! -e "kak_finalClusters_all.fa" ]; then
			if [ ! -e "kak_contig.fa" ]; then
				touch kak_contig.fa
			fi
			ln -f -s kak_contig.fa kak_finalClusters_all.fa
		fi
		
		cd ../
		
		groups=($(printf "%s\n" ${groups[@]} ${name} | sort -u ))
	
	done
	
	assemstats3.py ${cutoff} *.idba_ud.*.d/scaffold.fa *.newbler.*.d/assembly/454AllContigs.fna *.megahit.*.d/final.contigs.fa *.spades.*.d/scaffolds.fasta *.metaspades.*.d/scaffolds.fasta *.minia.*.d/kak.contigs.fa *.velvet.*.d/contigs.fa *.metavelvet.*.d/meta-velvetg.contigs.fa *.ray.*.d/Scaffolds.fasta *.abyss.*.d/kak-scaffolds.fa *.meta_platanus.*.d/kak_finalClusters_all.fa > assemstats3.out.txt 2> assemstats3.err.txt
	
	for g in ${groups[@]}; do
		echo "   Evaluating assembly of group ${g} ..."
		calc-best-assembly.py -C ${cutoff} -q ${g}.{idba_ud.*.d/scaffold.fa,newbler.*.d/assembly/454AllContigs.fna,megahit.*.d/final.contigs.fa,spades.*.d/scaffolds.fasta,metaspades.*.d/scaffolds.fasta,minia.*.d/kak.contigs.fa,velvet.*.d/contigs.fa,metavelvet.*.d/meta-velvetg.contigs.fa,ray.*.d/Scaffolds.fasta,abyss.*.d/kak-scaffolds.fa,*.meta_platanus.*.d/kak_finalClusters_all.fa} -o ${g}.best.fa > calc-best-assembly.out.txt 2> calc-best-assembly.err.txt
	done
	
	multi-rename.py ${asmname} *.best.fa > final-assembly.fa
	
	finalasm=($(printf "%s\n" ${finalasm[@]} "../${sampname}/assembled/final-assembly.fa" | sort -u ))
	
	assemblathon_stats.pl final-assembly.fa > assemblathon_stats.out.txt 2> assemblathon_stats.err.txt
	
	echo "Finish ${sampname} assembly"
	
	cd ../../
done

# Fundindo montagens de amostras
echo "Merging assemblies ..."

mkdir -p ./assembled

cd ./assembled

### MeGAMerge
#MeGAMerge-1.1.pl -overlap=${cutoff} -minID=99 -cpu=${threads} -force -o=final-assembly.fa . ${finalasm[*]}

### usearch
cat ${finalasm[*]} > finalasm.tmp
usearch81 -cluster_fast finalasm.tmp -id 0.99 -sort size -threads ${threads} -centroids nr.fasta -uc clusters.uc -consout final-assembly.fa > usearch.out.txt 2> usearch.err.txt
rm -f finalasm.tmp

assemblathon_stats.pl MergedContigs.fasta > assemblathon_stats.out.txt 2> assemblathon_stats.err.txt

cd ../

echo "Finish Final Assembly (${asmname})"

cd ${curdir}
