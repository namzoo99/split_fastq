workflow split_fastq {

		R_ch = Channel.value(file("${params.nf_home}/R/nextflow"))


		raw_fq_samples = make_fastq_samples( params.dna_legacy_file_sheet )

		read_fq(raw_fq_samples)

		count = read_fq.out.map{aliquot, RGPU_fq1, count, n_read, cmds -> [ count ].collect()}
						   .unique()

		make_manifest_for_single_rgid(count, R_ch)

		gdc_fastq_splitter(raw_fq_samples.combine(read_fq.out.map{aliquot, RGPU_fq1, count, n_read, cmds -> [ aliquot, count ]}, by:0))

		validate_splitted_fq(gdc_fastq_splitter.out.map{aliquot, splitted_FQs, splitted_FQpath, cmds -> [ aliquot, splitted_FQs ]},
						 	 R_ch)
	 
		flag = validate_splitted_fq.out.map{aliquot, n_read, RGPU, flag, cmds -> [ flag ].collect()}
									   .unique()
									    
		make_manifest(gdc_fastq_splitter.out.map{aliquot, splitted_FQs, splitted_FQpath, cmds -> [ splitted_FQpath ]}.collect(),
					  flag,
					  R_ch)

}

/*
================================================================================

                               Processes

================================================================================
*/
process read_fq {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/fq_header/${aliquot_barcode}", pattern: "*", mode: 'copy'

	input: 
		tuple val(aliquot_barcode), path(fq1), path(fq2)

	output:
		tuple val(aliquot_barcode), path("${aliquot_barcode}_RGPU_fq1.txt"), env(count), path("*_fq1_n_of_read.txt"), path("*{command,exitcode}*", hidden:true)

	script:

		"""
		#!/bin/bash

		zcat ${fq1} | grep '^@' | awk '{split(\$0, a, ":"); print a[1] ":" a[2] ":" a[3] ":" a[4]}' | sort -u > ./${aliquot_barcode}_RGPU_fq1.txt
		#zcat ${fq2} | grep '@' | cut -c 1-23 | uniq > ./${aliquot_barcode}_RGPU_fq2.txt

		zcat ${fq1} | wc -l | awk '{print \$1/4}' > ./${aliquot_barcode}_fq1_n_of_read.txt
		#zcat ${fq2} | wc -l | awk '{print \$1/4}' > ./${aliquot_barcode}_fq2_n_of_read.txt

		count=\$(awk 'END {print NR}' ${aliquot_barcode}_RGPU_fq1.txt)

		"""

}

process make_manifest_for_single_rgid { 

	tag "all files"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/final_manifest/no_multiplexed", pattern: "*", mode: 'copy'

	label "r_for_split_fastq"

	input: val(count)
		   file(R_dir)

	output:
		path("*.csv")

	when:
		count = "1"

	script:

		"""
		#!/bin/bash

		Rscript ${R_dir}/make_final_manifest_for_single_rgid.R \
		    --raw_manifest ${params.dna_legacy_file_sheet} \
		    --previous_output_path ${params.scratch_dir}/results/${params.workflow_name} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> env.txt

		"""

}

process gdc_fastq_splitter { 

	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/gdc_fastq_splitter/${aliquot_barcode}", pattern: "*", mode: 'copy'

	label "gdc_fastq_splitter"

	input: 
		tuple val(aliquot_barcode), path(fq1), path(fq2), val(count)

	output:
		tuple val(aliquot_barcode), path("*fq.gz"), path("*splitted_fq_path.txt"), path("*{command,exitcode}*", hidden:true) 

	when:
		count != "1"

	script:

		"""
		#!/bin/bash

		gdc-fastq-splitter \
		-o ${aliquot_barcode}_ \
		${fq1} \
		${fq2} && \
		find . -name "*.fq.gz" -exec realpath {} \\; > ${aliquot_barcode}_splitted_fq_path.txt
		"""

}

process validate_splitted_fq { 

	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/validation/${aliquot_barcode}", pattern: "*", mode: 'copy'

	label "r_for_split_fastq"
	
	input:
		tuple val(aliquot_barcode), path(splitted_FQs)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("*_n_read.txt"), path("*_RGPU.txt"), path("*.flag"), path("*{command,exitcode}*", hidden:true) 

	script:

		"""
		#!/bin/bash

		for f in ${splitted_FQs}; do
  			result1=\$(zcat "\$f" | wc -l | awk '{print \$1/4}')
 			output_file1="\${f%.gz}_n_read.txt"
  			echo "\$result1" > "\$output_file1"
  			result2=\$(zcat "\$f" | grep '@' | cut -c 1-23 | uniq)
 			output_file2="\${f%.gz}_RGPU.txt"
  			echo "\$result2" > "\$output_file2"		
  		done && \
  		Rscript ${R_dir}/validate_splitted_fqs.R \
		    --previous_output_path ${params.scratch_dir}/results/${params.workflow_name}/fq_header/${aliquot_barcode}

		"""

}

process make_manifest { 

	tag "all files"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/final_manifest/multiplexed", pattern: "*", mode: 'copy'

	label "r_for_split_fastq"

	input:
		path(splitted_FQpath)
		val(flag)
		file(R_dir)

	output:
		path("*.csv")

	when:
		flag = "success.flag" 

	script:

		"""
		#!/bin/bash

		Rscript ${R_dir}/make_final_manifest.R \
		    --raw_manifest ${params.dna_legacy_file_sheet} \
		    --previous_output_path ${params.scratch_dir}/results/${params.workflow_name} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> env.txt

		"""

}


/*
================================================================================

                               Functions

================================================================================
*/

def make_fastq_samples( samplesFile )  {

    Channel
        .fromPath(samplesFile)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aliquot_barcode,
                            row.sample_barcode,
                            row.FQ1,
                            row.FQ2,
                            row.action) }  
        .filter{
                aliquot, sample, FQ1, FQ2, action  ->
                action == 'run'
            }
        .map {  aliquot, sample, FQ1, FQ2, action ->
                                    def fq1 = file("${FQ1}", checkIfExists: true)
                                    def fq2 = file("${FQ2}", checkIfExists: true)
        return([    aliquot, fq1, fq2    ] ) }
        .unique()

    }
