#!/bin/csh
umask 002

set seq454_orig_sff = $1
set sanger_orig_fasta = $2
set solexa_orig_fastq = $3
set solexa_orig_trimpts = $4
set db_name = $5

set ROOT_DIR = "/usr/local/VHTNGS"
set PROJECT_DIR = "${ROOT_DIR}/project"
set REF_DIR = "${ROOT_DIR}/references"
set TOOLS_DIR = "${ROOT_DIR}/tools"
set TOOLS_BINARIES_DIR = "${TOOLS_DIR}/BINARIES"
set TOOLS_PERL_DIR = "${TOOLS_DIR}/PERL"
set TOOLS_RUBYDIR_DIR = "${TOOLS_DIR}/RUBYLIB"
set TOOLS_SFF_DIR = "${TOOLS_DIR}/SFF"
set TOOLS_FASTX_DIR = "${TOOLS_DIR}/FASTX"

cd ${PROJECT_DIR}

set seq454_orig_sff_dir = "${PROJECT_DIR}/input_sff"
set seq454_orig_sff_filename = `basename ${seq454_orig_sff}`
set seq454_orig_sff_file = "${seq454_orig_sff_dir}/${seq454_orig_sff_filename}"
mkdir -p ${seq454_orig_sff_dir}
cp ${seq454_orig_sff} ${seq454_orig_sff_file}

set sanger_orig_fasta_dir = "${PROJECT_DIR}/input_sanger"
set sanger_orig_fasta_filename = `basename ${sanger_orig_fasta}`
set sanger_orig_fasta_file = "${sanger_orig_fasta_dir}/${sanger_orig_fasta_filename}"
mkdir -p ${sanger_orig_fasta_dir}
cp ${sanger_orig_fasta} ${sanger_orig_fasta_file}

set solexa_orig_fastq_dir = "${PROJECT_DIR}/input_solexa"
set solexa_orig_fastq_filename = `basename ${solexa_orig_fastq}`
set solexa_orig_fastq_file = "${solexa_orig_fastq_dir}/${solexa_orig_fastq_filename}"
mkdir -p ${solexa_orig_fastq_dir}
cp ${solexa_orig_fastq} ${solexa_orig_fastq_file}

set solexa_orig_trimpts_filename = `basename ${solexa_orig_trimpts}`
set solexa_orig_trimpts_file = "${solexa_orig_fastq_dir}/${solexa_orig_trimpts_filename}"
cp ${solexa_orig_trimpts} ${solexa_orig_trimpts_file}


set flu_a = 0
switch ($db_name)
  case barda:
    echo "Using Influenza A reference data for database [${db_name}]"
    set ref_dir = ${REF_DIR}/influenza_a_virus
    set blast_db_dir = ${REF_DIR}/influenza_a_virus_full_length_NT
    set segments = "HA MP NA NP NS PA PB1 PB2"
    set seg_cov = "HA:175000 MP:100000 NA:145000 NP:155000 NS:89000 PA:220000 PB1:235000 PB2:235000"
    set flu_a = 1
  breaksw
  case giv:
    echo "Using Influenza A reference data for database [${db_name}]"
    set ref_dir = ${REF_DIR}/influenza_a_virus
    set blast_db_dir = ${REF_DIR}/influenza_a_virus_full_length_NT
    set segments = "HA MP NA NP NS PA PB1 PB2"
    set seg_cov = "HA:175000 MP:100000 NA:145000 NP:155000 NS:89000 PA:220000 PB1:235000 PB2:235000"
    set flu_a = 1
  breaksw
  case giv3:
    echo "Using Influenza A reference data for database [${db_name}]"
    set ref_dir = ${REF_DIR}/influenza_a_virus
    set blast_db_dir = ${REF_DIR}/influenza_a_virus_full_length_NT
    set segments = "HA MP NA NP NS PA PB1 PB2"
    set seg_cov = "HA:175000 MP:100000 NA:145000 NP:155000 NS:89000 PA:220000 PB1:235000 PB2:235000"
    set flu_a = 1
  breaksw
  case piv:
    echo "Using Influenza A reference data for database [${db_name}]"
    set ref_dir = ${REF_DIR}/influenza_a_virus
    set blast_db_dir = ${REF_DIR}/influenza_a_virus_full_length_NT
    set segments = "HA MP NA NP NS PA PB1 PB2"
    set seg_cov = "HA:175000 MP:100000 NA:145000 NP:155000 NS:89000 PA:220000 PB1:235000 PB2:235000"
    set flu_a = 1
  breaksw
  case swiv:
    echo "Using Influenza A reference data for database [${db_name}]"
    set ref_dir = ${REF_DIR}/influenza_a_virus
    set blast_db_dir = ${REF_DIR}/influenza_a_virus_full_length_NT
    set segments = "HA MP NA NP NS PA PB1 PB2"
    set seg_cov = "HA:175000 MP:100000 NA:145000 NP:155000 NS:89000 PA:220000 PB1:235000 PB2:235000"
    set flu_a = 1
  breaksw
  case rtv:
    echo "Using Rotavirus reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/rota_virus
    set blast_db_dir = ${REF_DIR}/rota_virus_full_length_NT
    set segments = "VP1 VP2 VP3 VP4 NSP1 VP6 NSP3 NSP2 VP7 NSP4 NSP5"
    set seg_cov = "VP1:326700 VP2:268600 VP3:255000 VP4:232400 NSP1:151800 VP6:132300 NSP3:104100 NSP2:102200 VP7:103000 NSP4:70800 NSP5:62900"
    set flu_a = 0
  breaksw
  case gcv:
    echo "Using Coronavirus reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/corona_virus
    set blast_db_dir = ${REF_DIR}/corona_virus_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:3000000"
    set flu_a = 0
  breaksw
  case veev:
    echo "Using VEEV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/veev
    set blast_db_dir = ${REF_DIR}/veev_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:1200000"
    set flu_a = 0
  breaksw
  case hadv:
    echo "Using HADV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/hadv
    set blast_db_dir = ${REF_DIR}/hadv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:4500000"
    set flu_a = 0
  breaksw
  case mpv:
    echo "Using MPV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/mpv
    set blast_db_dir = ${REF_DIR}/mpv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:1335000"
    set flu_a = 0
  breaksw
  case norv:
    echo "Using NORV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/norv
    set blast_db_dir = ${REF_DIR}/norv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:774600"
    set flu_a = 0
  breaksw
  case vzv:
    echo "Using VZV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/vzv
    set blast_db_dir = ${REF_DIR}/vzv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:12500000"
    set flu_a = 0
  breaksw
  case rsv:
    echo "Using RSV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/rsv
    set blast_db_dir = ${REF_DIR}/rsv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:1530000"
    set flu_a = 0
  breaksw
  case jev:
    echo "Using JEV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/jev
    set blast_db_dir = ${REF_DIR}/jev_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:1100000"
    set flu_a = 0
  breaksw
  case yfv:
    echo "Using YFV reference data for database [${db_name}]"
    set ref_dir      = ${REF_DIR}/yfv
    set blast_db_dir = ${REF_DIR}/yfv_full_length_NT
    set segments = "MAIN"
    set seg_cov = "MAIN:1090000"
    set flu_a = 0
  breaksw
  default:
    echo "Using no reference data for database [${db_name}]"
    set ref_dir      = ""
    set blast_db_dir = ""
    set segments = ""
    set seg_cov = ""
    set flu_a = 0
  breaksw
endsw

echo ${seq454_orig_sff_file}
echo ${sanger_orig_fasta_file}
echo ${solexa_orig_fastq_file}
echo ${db_name}
echo ${ref_dir}
echo ${blast_db_dir}
echo ${segments}
echo ${seg_cov}
echo ${flu_a}

echo "INFO: converting sff to fasta"
touch ${seq454_orig_sff_file}.fna
${TOOLS_SFF_DIR}/sffinfo -s ${seq454_orig_sff_file} | grep -v " length=0 " >> ${seq454_orig_sff_file}.fna 

if ( `cat ${seq454_orig_sff_file}.fna | wc -l` > 0 ) then
  echo "INFO: formatdb of SFF fasta"
  ${TOOLS_BINARIES_DIR}/formatdb -p F -i ${seq454_orig_sff_file}.fna
endif

if ( `cat ${sanger_orig_fasta_file} | wc -l` > 0 ) then
  echo "INFO: formatdb of Sanger fasta"
  ${TOOLS_BINARIES_DIR}/formatdb -p F -i ${sanger_orig_fasta_file}
endif

set tblastx_outdir = ${PROJECT_DIR}/tblastx_output
mkdir -p ${tblastx_outdir}

foreach seg ( `echo ${segments} | tr ' ' '\n' ` )
  echo "INFO: tblastx segment data [${seg}] against SFF and Sanger reads databases for [${db_name}]"
  set ref_fna = ${ref_dir}/${seg}.fasta

  if ( `cat ${seq454_orig_sff_file}.fna | wc -l` > 0 ) then
    set blastdb = ${seq454_orig_sff_file}.fna
    set blastout = ${tblastx_outdir}/${seg}.out
    ${TOOLS_BINARIES_DIR}/blastall \
      -p tblastx \
      -d ${blastdb} \
      -i ${ref_fna} \
      -m 9 \
      -b 100000 \
      -v 100000 \
      -o ${blastout}
  else
    touch ${blastout}
  endif

  if ( `cat ${sanger_orig_fasta_file} | wc -l` > 0 ) then
    set blastdb = ${sanger_orig_fasta_file}
    set blastout = ${tblastx_outdir}/${seg}_sanger.out
    ${TOOLS_BINARIES_DIR}/blastall \
      -p tblastx \
      -d ${blastdb} \
      -i ${ref_fna} \
      -m 9 \
      -b 100000 \
      -v 100000 \
      -o ${blastout}
  else
    touch ${blastout}
  endif
end

set noninter_chimera_list = ${tblastx_outdir}/noninter_chimera_reads.uaccno_list
set inter_chimera_list = ${tblastx_outdir}/inter_chimera_reads.uaccno_list
foreach seg ( `echo ${segments} | tr ' ' '\n' ` )
  echo "INFO: parsing tblastx output for segment [${seg}] against [${db_name}]"
  set blastout = ${tblastx_outdir}/${seg}.out
  set blastout_sanger = ${tblastx_outdir}/${seg}_sanger.out
  set nonintra_chimera_list = ${tblastx_outdir}/${seg}_nonintra_chimera_reads.uaccno_list
  set intra_chimera_list = ${tblastx_outdir}/${seg}_intra_chimera_reads.uaccno_list
  if ( ${flu_a} > 0 ) then
    cat ${blastout} ${blastout_sanger} | \
      gawk '{if($0 !~ "#" && $3>60 && $4 > 25 ) {if(($7-$8)*($9-$10)>0){o=1;d=$8-$10;}else{o=0;d=$8+$10;}{printf("%s\t%s\t%s\t%06d\n",$1,$2, o,d);}}}' | \
      sort -g | \
      uniq | \
      gawk '{if( $1==q && $2==s && $3==o && sqrt(($4-d)^2) < 6){d=$4;}else {print $0;q=$1;s=$2;o=$3;d=$4;}}' | \
      cut -f 1,2 | \
      uniq -c | \
      grep " 1 " | \
      gawk '{print $3}' | \
      sort | \
      uniq > ${nonintra_chimera_list}
    cat ${blastout} ${blastout_sanger} | \
      gawk '{if($0 !~ "#" && $3>60 && $4 > 25 ) {if(($7-$8)*($9-$10)>0){o=1;d=$8-$10;}else{o=0;d=$8+$10;}{printf("%s\t%s\t%s\t%06d\n",$1,$2, o,d);}}}' | \
      sort -g | \
      uniq | \
      gawk '{if( $1==q && $2==s && $3==o && sqrt(($4-d)^2) < 6){d=$4;}else {print $0;q=$1;s=$2;o=$3;d=$4;}}' | \
      cut -f 1,2 | \
      uniq -c | \
      grep -v " 1 " | \
      gawk '{print $3}' | \
      sort | \
      uniq > ${intra_chimera_list}
  else
    cat ${blastout} ${blastout_sanger} | \
      gawk '{if($0 !~ "#" && $3>60 && $4 > 25 ) {print $2;}}' | \
      sort | \
      uniq > ${nonintra_chimera_list}
    echo "" > ${intra_chimera_list}     
  endif
end
cat ${tblastx_outdir}/*_nonintra_chimera_reads.uaccno_list | \
  sort | \
  uniq -c | \
  tr '\t' ' ' | \
  grep " 1 " | \
  gawk '{print $2}' > ${noninter_chimera_list}
cat ${tblastx_outdir}/*_nonintra_chimera_reads.uaccno_list | \
  sort | \
  uniq -c | \
  tr '\t' ' ' | \
  grep -v " 1 " | \
  gawk '{print $2}' > ${inter_chimera_list}

foreach seg ( `echo ${segments} | tr ' ' '\n' ` )
  set nonintra_chimera_list = ${tblastx_outdir}/${seg}_nonintra_chimera_reads.uaccno_list
  set non_chimera_list = ${tblastx_outdir}/${seg}_nonchimera_reads.uaccno_list
  join -1 1 -2 1 \
    ${noninter_chimera_list} \
    ${nonintra_chimera_list} > \
    ${non_chimera_list}
  echo "INFO: creating sff of non_chimeric reads from reads matching segment [${seg}] for [${db_name}]"

  set sample_seg_sff_file = ${seq454_orig_sff_dir}/sample_nonchimera_${seg}.sff
  ${TOOLS_SFF_DIR}/sfffile \
    -i ${non_chimera_list} \
    -o ${sample_seg_sff_file} \
    ${seq454_orig_sff_file}

  set sample_seg_sanger_file = ${sanger_orig_fasta_dir}/sample_nonchimera_${seg}.fasta
  if ( `cat ${sanger_orig_fasta_file} | wc -l` > 0 ) then
    echo "INFO: creating Sanger fasta of non_chimeric reads from reads matching segment [${seg}]"
    ${TOOLS_SFF_DIR}/fnafile \
      -i ${non_chimera_list} \
      -o ${sample_seg_sanger_file} \
      ${sanger_orig_fasta_file}
  endif

  echo "INFO: creating 100x max coverage sff of non_chimeric reads from reads matching segment [${seg}]"
  set bps = `echo ${seg_cov} | tr ' ' '\n' | grep ${seg} | cut -d ':' -f 2`
  set sample_seg_100x_sff_file = ${seq454_orig_sff_dir}/sample_nonchimera_${seg}_100x.sff
  ${TOOLS_SFF_DIR}/sfffile \
    -pick ${bps} \
    -o ${sample_seg_100x_sff_file} \
    ${sample_seg_sff_file}
  
  set seg_assembly_dir = ${PROJECT_DIR}/assembly_by_segment/${seg}
  mkdir -p ${seg_assembly_dir}

  pushd ${seg_assembly_dir} >& /dev/null
    echo "INFO: performing de novo assembly of 100x coverage for nonchimera reads from segment [${seg}] for [${db_name}]"
    if ( -e sample_nonchimera_${seg}_100x.fasta ) then
      rm sample_nonchimera_${seg}_100x.fasta
    endif
    touch sample_nonchimera_${seg}_100x.fasta
    ${TOOLS_SFF_DIR}/sffinfo \
      -s ${sample_seg_100x_sff_file} | \
      grep -v " length=0 " \
      >> sample_nonchimera_${seg}_100x.fasta
    if ( -e ${sample_seg_sanger_file} ) then
      cat ${sample_seg_sanger_file} >> sample_nonchimera_${seg}_100x.fasta
    endif
    
    /usr/local/bin/cap3 sample_nonchimera_${seg}_100x.fasta
    mv sample_nonchimera_${seg}_100x.fasta.cap.contigs ${seg}_100x_contigs.fasta
    rm sample_nonchimera_${seg}_100x.fasta*
    
    set blast_db = ${blast_db_dir}/${seg}_full_length_NT_complete.fa
    set best_reference = ${seg_assembly_dir}/${seg}_best_reference.fna
    if ( -e ${seg}_100x_contigs.fasta ) then
      echo "INFO: finding best FL reference for segment [${seg}] for [${db_name}]"
      set best_hit = `${TOOLS_BINARIES_DIR}/blastall -p blastn -d ${blast_db} -b 1 -v 1 -m 8 -i ${seg}_100x_contigs.fasta | gawk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12);}' | sort -nrk12 | head -n 1 | gawk '{printf("%s\n", $2);}'`
      touch ${best_reference}
      chmod 755 ${best_reference}
      ${TOOLS_BINARIES_DIR}/fastacmd -d ${blast_db} -p F -s "${best_hit}"  -o ${best_reference}
      grep "^>" ${best_reference} | cut -c 2- | gawk -v s=${seg} '{printf(">%s %s\n", s, $0);}' > ${best_reference}_mod
      grep -v "^>" ${best_reference} >> ${best_reference}_mod
      mv ${best_reference}_mod ${best_reference}
    else
      echo "ERROR: missing de novo assembly of 100x coverage for nonchimera reads from segment [${seg}] for [${db_name}]"
    endif
  popd >& /dev/null
end

echo "INFO: consolidating best FL reference sequences for [${db_name}]"
set seg_best_ref_dir = ${PROJECT_DIR}/reference_fasta
mkdir -p ${seg_best_ref_dir}
set best_refs_file = ${seg_best_ref_dir}/reference.fasta
cat ${PROJECT_DIR}/assembly_by_segment/*/*_best_reference.fna > ${best_refs_file}

echo "INFO: consolidating nonchimera 454 reads for [${db_name}]"
set non_chimera_list = ${tblastx_outdir}/nonchimera_reads.uaccno_list
cat ${tblastx_outdir}/*_nonchimera_reads.uaccno_list > ${non_chimera_list}

set deconvolved_sff = ${seq454_orig_sff_dir}/sample_nonchimera.sff
${TOOLS_SFF_DIR}/sfffile -i ${non_chimera_list} -o ${deconvolved_sff} ${seq454_orig_sff_file}

echo "INFO: mapping viral sequences for [${db_name}]"
set sample_mapping_dir = ${PROJECT_DIR}/mapping
mkdir -p ${sample_mapping_dir}

pushd ${sample_mapping_dir} >& /dev/null
  set deconvolved_sff_fna = ${deconvolved_sff}.fna
  set solexa_orig_fastq_fna_file = ${solexa_orig_fastq_file}.fna
  
  echo "INFO: converting non-chimeric sff to fasta"
  touch ${deconvolved_sff_fna}
  ${TOOLS_SFF_DIR}/sffinfo -s ${deconvolved_sff} | grep -v " length=0 " >> ${deconvolved_sff_fna}
  
  echo "INFO: bowtie-build on best references fasta"
  /usr/bin/bowtie-build ${best_refs_file} BEST_REFS

  echo "INFO: using BOWTIE and SAMTOOLS to find sff SNPs for [${db_name}]"
  /usr/bin/bowtie -S -f BEST_REFS ${deconvolved_sff_fna} sample_454_only_gb_refs.sam
  /usr/bin/samtools view -bS -o sample_454_only_gb_refs.bam sample_454_only_gb_refs.sam
  /usr/bin/samtools sort sample_454_only_gb_refs.bam sample_454_only_gb_refs.sorted
  /usr/bin/samtools mpileup -ugf ${best_refs_file} sample_454_only_gb_refs.sorted.bam | /usr/bin/bcftools view -bvcg - > sample_454_only_gb_refs.raw.bcf
  /usr/bin/bcftools view sample_454_only_gb_refs.raw.bcf | ${TOOLS_PERL_DIR}/vcfutils.pl varFilter -D 500 > sample_454_only_gb_refs.SNPs.txt
  grep -v "^#" sample_454_only_gb_refs.SNPs.txt | gawk -F'\t' '{split($10,a,":"); printf("%s:%s:%s:%s\n",$1,$2,$5,a[2]);}' | gawk -F':' '{if(index($3,",")>0) {split($4,s,",") split($3,b,","); if(s[3]>s[6]) printf("%s:%s:%s\n",$1,$2,b[2]); else printf("%s:%s:%s\n",$1,$2,b[1]);} else printf("%s:%s:%s\n",$1,$2,$3);}' > sample_454_only_gb_refs.SNPs.reduced.txt

  echo "INFO: converting solexa to fasta"
  touch ${solexa_orig_fastq_fna_file}
  ${TOOLS_FASTX_DIR}/fastq_to_fasta -Q 33 -i ${solexa_orig_fastq_file} -o ${solexa_orig_fastq_fna_file}
  
  echo "INFO: using BOWTIE and SAMTOOLS to find fastq SNPs for [${db_name}]"
  /usr/bin/bowtie -S -f BEST_REFS ${solexa_orig_fastq_fna_file} sample_solexa_only_gb_refs.sam
  /usr/bin/samtools view -bS -o sample_solexa_only_gb_refs.bam sample_solexa_only_gb_refs.sam
  /usr/bin/samtools sort sample_solexa_only_gb_refs.bam sample_solexa_only_gb_refs.sorted
  /usr/bin/samtools mpileup -ugf ${best_refs_file} sample_solexa_only_gb_refs.sorted.bam | /usr/bin/bcftools view -bvcg - > sample_solexa_only_gb_refs.raw.bcf
  /usr/bin/bcftools view sample_solexa_only_gb_refs.raw.bcf | ${TOOLS_PERL_DIR}/vcfutils.pl varFilter -D 500 > sample_solexa_only_gb_refs.SNPs.txt
  grep -v "^#" sample_solexa_only_gb_refs.SNPs.txt | gawk -F'\t' '{split($10,a,":"); printf("%s:%s:%s:%s\n",$1,$2,$5,a[2]);}' | gawk -F':' '{if(index($3,",")>0) {split($4,s,",") split($3,b,","); if(s[3]>s[6]) printf("%s:%s:%s\n",$1,$2,b[2]); else printf("%s:%s:%s\n",$1,$2,b[1]);} else printf("%s:%s:%s\n",$1,$2,$3);}' > sample_solexa_only_gb_refs.SNPs.reduced.txt
  
  if ( `cat sample_solexa_only_gb_refs.SNPs.reduced.txt | wc -l` > 0 ) then
      sdiff \
      sample_454_only_gb_refs.SNPs.reduced.txt \
      sample_solexa_only_gb_refs.SNPs.reduced.txt | \
      grep -v "[<|>]" | \
      cut -f 1 > \
      sample_454_solexa_common_gb_refs.SNPs.reduced.txt
  else
      cp \
      sample_454_only_gb_refs.SNPs.reduced.txt \
      sample_454_solexa_common_gb_refs.SNPs.reduced.txt
  endif

  echo "INFO: building edited references based on common sff and fastq SNPs for [${db_name}]"
  foreach seg ( `grep "^>" ${best_refs_file} | cut -d ' ' -f 1 | cut -c 2-` )
    nthseq -sequence ${best_refs_file} \
      -number `grep "^>" ${best_refs_file} | cut -d ' ' -f 1 | cut -c 2- | grep -n ${seg} | cut -d ':' -f 1` \
      -outseq sample_${seg}.extracted >& /dev/null
    cat sample_454_solexa_common_gb_refs.SNPs.reduced.txt | \
      grep ${seg} | \
      cut -d ':' -f 2-3 | \
      tr '\n ' ' ' > sample_${seg}.edits
    ${TOOLS_PERL_DIR}/delta2seq.pl \
      -r sample_${seg}.extracted \
      -f sample_${seg}.edits \
      -q sample_${seg}.extracted.edited
    grep "^>" sample_${seg}.extracted > \
      sample_${seg}.extracted.edited.fasta
    grep -v "^>" sample_${seg}.extracted.edited >> \
      sample_${seg}.extracted.edited.fasta
  end
  set best_edited_refs_file = sample_reference_edited.fasta
  cat sample_*.extracted.edited.fasta > \
    ${best_edited_refs_file}
  
  echo "INFO: using 454 mapper for final chimera check for [${db_name}]"

  if ( -d 454_mapping_best_refs_chimera_check ) then
    rm -Rf 454_mapping_best_refs_chimera_check
  endif
  ${TOOLS_SFF_DIR}/newMapping 454_mapping_best_refs_chimera_check
  ${TOOLS_SFF_DIR}/setRef 454_mapping_best_refs_chimera_check ${best_edited_refs_file}
  ${TOOLS_SFF_DIR}/addRun 454_mapping_best_refs_chimera_check ${deconvolved_sff}

  if ( `cat ${solexa_orig_fastq_file} | wc -l` > 0 ) then
    ${TOOLS_SFF_DIR}/addRun 454_mapping_best_refs_chimera_check ${solexa_orig_fastq_file}
  endif

#  ${TOOLS_SFF_DIR}/runProject -no 454_mapping_best_refs_chimera_check >& runProject_454_mapping_best_refs_chimera_check.log
#  grep "Chimeric" 454_mapping_best_refs_chimera_check/mapping/454ReadStatus.txt | \
#    gawk '{print $1}' > exclude_list.txt

  cat ${solexa_orig_trimpts_file} | gawk -F'\t' '{if($2!=29){print $1;}}' >> exclude_list.txt
  
  set final_sff_reads = sample_final.sff
  set final_fastq_reads = sample_final.fastq
  set final_fastq_fna_reads = sample_final.fastq.fna
  set final_fasta_reads = sample_final.fasta

  if ( `cat ${sanger_orig_fasta_file} | wc -l` > 0 ) then
    cp ${sanger_orig_fasta_file} ${final_fasta_reads}
  endif
  
  ${TOOLS_SFF_DIR}/sfffile -o ${final_sff_reads} -e exclude_list.txt ${deconvolved_sff}
  
  touch ${final_fastq_reads}
  
  if ( `cat ${solexa_orig_fastq_file} | wc -l` > 0 ) then
    ${TOOLS_PERL_DIR}/fastqfile.pl \
      -o ${final_fastq_reads} \
      -e exclude_list.txt \
      -f ${solexa_orig_fastq_file}
    
    cat ${final_fastq_reads} | \
      gawk '{t=NR % 4;\
             if(t==1){\
               if(length(sid) > 0 ) {printf("%s\t%s\t%s\t%s\n", sid,s,qid,q)};\
               sid=substr($0,2);\
             }\
             else if (t==2){s=$0;}\
             else if (t==3){qid=substr($0,2);}\
             else if (t==0){q=$0;}\
            }\
            END {\
              if(length(sid) > 0 ) {printf("%s\t%s\t%s\t%s\n", sid,s,qid,q)};\
              sid=substr($0,2);\
            }' | \
      sort | \
      gawk -F'\t' '{printf("@%s\n%s\n+%s\n%s\n", $1, $2, $3, $4);}' > ${final_fastq_reads}.sorted
    mv ${final_fastq_reads} ${final_fastq_reads}.unsorted
    mv ${final_fastq_reads}.sorted ${final_fastq_reads}
  endif
  
  echo "INFO: running BOWTIE and SAMTOOLS for [${db_name}]"

  set input_read_files = ""
  set input_read_files = `echo "${input_read_files} -q ${final_sff_reads}"`
  
  if ( `cat ${final_fasta_reads} | wc -l` > 0 ) then
    set input_read_files = `echo "${input_read_files} -q ${final_fasta_reads}"`
  endif
  
  if ( `cat ${final_fastq_reads} | wc -l` > 0 ) then
    echo "INFO: converting solexa to fasta"
    touch ${final_fastq_fna_reads}
    ${TOOLS_FASTX_DIR}/fastq_to_fasta -Q 33 -i ${final_fastq_reads} -o ${final_fastq_fna_reads}
    set input_read_files = `echo "${input_read_files} -q ${final_fastq_fna_reads}"`
  endif
  
  echo "INFO: bowtie-build on best edited references fasta"
  /usr/bin/bowtie-build ${best_edited_refs_file} BEST_EDITED_REFS
  
  echo "INFO: using BOWTIE and SAMTOOLS to acquire final consensus for [${db_name}]"
  /usr/bin/bowtie -S -f BEST_EDITED_REFS ${final_sff_reads},${final_fasta_reads},${final_fastq_fna_reads} sample_hybrid_edited_refs.sam
  /usr/bin/samtools view -bS -o sample_hybrid_edited_refs.bam sample_hybrid_edited_refs.sam
  /usr/bin/samtools mpileup -uf ${best_edited_refs_file} sample_hybrid_edited_refs.bam | /usr/bin/bcftools view -cg - | ${TOOLS_PERL_DIR}/vcfutils.pl vcf2fq > sample_hybrid_edited_refs_consensus.fastq

  ${TOOLS_FASTX_DIR}/fastq_to_fasta -Q 33 -i sample_hybrid_edited_refs_consensus.fastq -o sample_hybrid_edited_refs_consensus.fasta
    
popd >& /dev/null
