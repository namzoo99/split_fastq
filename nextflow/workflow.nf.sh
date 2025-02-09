// ######################################################################
//
// workflow_name: split_fastq
//
// ######################################################################
include { split_fastq } from './split_fastq/split_fastq.nf.sh' addParams(workflow_name: "split_fastq")

  // you should not add 'split_fastq' in step2List !!
    
    if (['split_fastq'].contains(params.step)) split_fastq()
