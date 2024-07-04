process index_ref {

    tag { sample_id + ' / ' + ref_filename }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "ref.fa"

    input:
    tuple val(sample_id), path(ref)

    output:
    tuple val(sample_id), path('ref.fa*'), emit: ref

    script:
    ref_filename = ref.getName()
    """
    cp ${ref} ref.fa
    bwa index ref.fa
    """
}


process filter_residual_adapters {

    tag { sample_id }

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("*1_posttrim_filter.fq.gz"), path("*2_posttrim_filter.fq.gz")

    script:
    """
    filter_residual_adapters.py --input_R1 $forward --input_R2 $reverse
    """
}


process bwa_mem {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}.{bam,bam.bai}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}{.bam,.bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_bwa_mem_provenance.yml"), emit: provenance
    
    script:
    bwa_threads = task.cpus - 8
    short_long = "short"
    samtools_view_filter_flags = params.skip_alignment_cleaning ? "0" : "1548"
    samtools_fixmate_remove_secondary_and_unmapped = params.skip_alignment_cleaning ? "" : "-r"
    samtools_markdup_remove_duplicates = params.skip_alignment_cleaning ? "" : "-r"
    """
    printf -- "- process_name: bwa_mem\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: bwa\\n"        >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(bwa 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "    - tool_name: samtools\\n"   >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n" >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      subcommand: sort\\n"      >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "      parameters:\\n"           >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "        - parameter: -l\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: 0\\n"          >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "        - parameter: -m\\n"     >> ${sample_id}_bwa_mem_provenance.yml
    printf -- "          value: 1000M\\n"      >> ${sample_id}_bwa_mem_provenance.yml

    bwa mem \
	-t ${bwa_threads} \
	-R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
	${ref[0]} \
	${reads_1} \
	${reads_2} \
	| samtools view -@ 2 -h \
	| samtools sort -@ 2 -l 0 -m 1000M \
	> ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}


process trim_primer_sequences {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.mapped.primertrimmed.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(bedfile)

    output:
    tuple val(sample_id), path("${sample_id}.mapped.primertrimmed.sorted{.bam,.bam.bai}"), emit: primer_trimmed_alignment
    tuple val(sample_id), path("${sample_id}_trim_primer_sequences_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: trim_primer_sequences\\n"     >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "    - tool_name: ivar\\n"        >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "      tool_version: \$(ivar version 2>&1 | head -n 1 | cut -d ' ' -f 3)\\n"      >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "      subcommand: trim\\n"      >> ${sample_id}_trim_primer_sequences_provenance.yml

    # Filter out unmapped reads
    samtools view -F4 -o ${sample_id}.mapped.bam ${alignment[0]}
    samtools index ${sample_id}.mapped.bam
    
    ivar trim \
	-e \
	-i ${sample_id}.mapped.bam \
	-b ${bedfile} \
	-p ivar.out
    
    samtools sort -o ${sample_id}.mapped.primertrimmed.sorted.bam ivar.out.bam
    samtools index ${sample_id}.mapped.primertrimmed.sorted.bam
    """
}


process make_consensus {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.primertrimmed.consensus.fa", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment)

    output:
    tuple val(sample_id), path("${sample_id}.primertrimmed.consensus.fa"), emit: consensus
    tuple val(sample_id), path("${sample_id}_make_consensus_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: make_consensus\\n"     >> ${sample_id}_make_consensus_provenance.yml
    printf -- "  tools:\\n"                           >> ${sample_id}_make_consensus_provenance.yml
    printf -- "    - tool_name: ivar\\n"              >> ${sample_id}_make_consensus_provenance.yml
    printf -- "      tool_version: \$(ivar version 2>&1 | head -n 1 | cut -d ' ' -f 3)\\n" >> ${sample_id}_make_consensus_provenance.yml
    printf -- "      subcommand: consensus\\n"        >> ${sample_id}_make_consensus_provenance.yml

    samtools mpileup -aa -A -B -d ${params.max_depth} -Q0 ${alignment[0]} \
	| ivar consensus -t ${params.unambiguous_allele_freq_threshold} -m ${params.min_depth} \
        -n N -p ${sample_id}.primertrimmed.consensus
    """
}

process align_consensus_to_ref {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_vs_ref.aln.fa", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_vs_ref.aln.fa"), emit: alignment
    tuple val(sample_id), path("${sample_id}_align_consensus_to_ref_provenance.yml"), emit: provenance

    script:
    awk_string = '/^>/ {printf("\\n%s\\n", $0); next; } { printf("%s", $0); }  END { printf("\\n"); }'
    """
    printf -- "- process_name: align_consensus_to_ref\\n" >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "  tools:\\n"                               >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "    - tool_name: mafft\\n"                  >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "      tool_version: \$(mafft --version 2>&1 | cut -d ' ' -f 1)\\n" >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "      parameters:\\n"                       >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "        - parameter: --preservecase\\n"     >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "          value: null\\n"                   >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "        - parameter: --keeplength\\n"       >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "          value: null\\n"                   >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "        - parameter: --add\\n"              >> ${sample_id}_align_consensus_to_ref_provenance.yml
    printf -- "          value: null\\n"                   >> ${sample_id}_align_consensus_to_ref_provenance.yml

    mafft \
        --preservecase \
        --keeplength \
        --add \
        ${consensus} \
        ${ref[0]} \
        > ${sample_id}_vs_ref_multi_line.aln.fa
    awk '${awk_string}' ${sample_id}_vs_ref_multi_line.aln.fa > ${sample_id}_vs_ref.aln.fa
    """
}


process minimap2 {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${short_long}.{bam,bam.bai}"

    input:
    tuple val(sample_id), path(reads), path(ref), path(ref_index)

    output:
    tuple val(sample_id), val(short_long), path("${sample_id}_${short_long}.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_minimap2_provenance.yml"), emit: provenance
    
    script:
    short_long = "long"
    minimap2_threads = task.cpus - 4
    samtools_view_filter_flags = params.skip_alignment_cleaning ? "0" : "1540"
    """
    printf -- "- process_name: \"minimap2\"\\n" >> ${sample_id}_minimap2_provenance.yml
    printf -- "  tools:\\n"                     >> ${sample_id}_minimap2_provenance.yml
    printf -- "    - tool_name: minimap2\\n"    >> ${sample_id}_minimap2_provenance.yml
    printf -- "      tool_version: \$(minimap2 --version)\\n"  >> ${sample_id}_minimap2_provenance.yml
    printf -- "      parameters:\\n"            >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -a\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: null\\n"        >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -x\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: map-ont\\n"     >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -MD\\n"     >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: null\\n"        >> ${sample_id}_minimap2_provenance.yml
    printf -- "    - tool_name: samtools\\n"    >> ${sample_id}_minimap2_provenance.yml
    printf -- "      tool_version: \$(samtools 2>&1 | grep 'Version' | cut -d ' ' -f 2)\\n"  >> ${sample_id}_minimap2_provenance.yml
    printf -- "      subcommand: view\\n"       >> ${sample_id}_minimap2_provenance.yml
    printf -- "      parameters:\\n"            >> ${sample_id}_minimap2_provenance.yml
    printf -- "        - parameter: -F\\n"      >> ${sample_id}_minimap2_provenance.yml
    printf -- "          value: ${samtools_view_filter_flags}\\n"        >> ${sample_id}_minimap2_provenance.yml

    minimap2 \
	-t ${minimap2_threads} \
	-ax map-ont \
	-R "@RG\\tID:${sample_id}-ONT\\tSM:${sample_id}\\tPL:ONT" \
	-MD \
	${ref} \
	${reads} \
	| samtools view -@ 2 -h -F ${samtools_view_filter_flags} \
	| samtools sort -@ 2 -l 0 -m 1000M -O bam \
	> ${sample_id}_${short_long}.bam

    samtools index ${sample_id}_${short_long}.bam
    """
}



