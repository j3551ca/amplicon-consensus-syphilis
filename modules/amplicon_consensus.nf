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
    bwa_threads = task.cpus - 4
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
    printf -- "  tools:\\n"                                  >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "    - tool_name: ivar\\n"                     >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "      tool_version: \$(ivar version 2>&1 | head -n 1 | cut -d ' ' -f 3)\\n"  >> ${sample_id}_trim_primer_sequences_provenance.yml
    printf -- "      subcommand: trim\\n"                    >> ${sample_id}_trim_primer_sequences_provenance.yml

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


process qualimap_bamqc {

    tag { sample_id }

    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_alignment_qc.csv"
    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_genome_results.txt"
    publishDir  "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_report.pdf"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_alignment_qc.csv"), emit: alignment_qc
    tuple val(sample_id), path("${sample_id}_qualimap_report.pdf"), emit: report, optional: true
    tuple val(sample_id), path("${sample_id}_qualimap_genome_results.txt"), emit: genome_results, optional: true
    tuple val(sample_id), path("${sample_id}_qualimap_bamqc_provenance.yml"), emit: provenance
    
    script:
    """
    printf -- "- process_name: \"qualimap_bamqc\"\\n"  >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "  tools:\\n"                            >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "    - tool_name: qualimap\\n"           >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "      tool_version: \$(qualimap bamqc | head | grep QualiMap | cut -d ' ' -f 2)\\n" >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "      parameters:\\n"                   >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "        - parameter: --collect-overlap-pairs\\n" >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "          value: null\\n"               >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "        - parameter: --cov-hist-lim\\n" >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "          value: ${params.qualimap_coverage_histogram_limit}\\n" >> ${sample_id}_qualimap_bamqc_provenance.yml

    # Assume qualimap exits successfully
    # If it fails we will re-assign the exit code
    # and generate an empty qualimap alignment qc
    qualimap_exit_code=0

    qualimap \
	--java-mem-size=${params.qualimap_memory} \
	bamqc \
	--paint-chromosome-limits \
	--collect-overlap-pairs \
	--cov-hist-lim ${params.qualimap_coverage_histogram_limit} \
	--output-genome-coverage ${sample_id}_genome_coverage.txt \
	-nt ${task.cpus} \
	-bam ${alignment[0]} \
	-outformat PDF \
	--outdir ${sample_id}_bamqc \
	|| qualimap_exit_code=\$?

    if [ \${qualimap_exit_code} -ne 0 ]; then
    echo "Qualimap failed with exit code \${qualimap_exit_code}"
        echo "Generating empty qualimap alignment qc"
        qualimap_bamqc_genome_results_to_csv.py \
	-s ${sample_id} \
	--failed \
	> ${sample_id}_qualimap_alignment_qc.csv
    else
	qualimap_bamqc_genome_results_to_csv.py \
	-s ${sample_id} \
	--qualimap-bamqc-genome-results ${sample_id}_bamqc/genome_results.txt \
	> ${sample_id}_qualimap_alignment_qc.csv
        cp ${sample_id}_bamqc/report.pdf ${sample_id}_qualimap_report.pdf
        cp ${sample_id}_bamqc/genome_results.txt ${sample_id}_qualimap_genome_results.txt
    fi
    """
}


process samtools_stats {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_samtools_stats*.{txt,tsv,csv}"

    input:
    tuple val(sample_id), path(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_samtools_stats.txt"), emit: stats
    tuple val(sample_id), path("${sample_id}_samtools_stats_summary.txt"), emit: stats_summary
    tuple val(sample_id), path("${sample_id}_samtools_stats_summary.csv"), emit: stats_summary_csv
    tuple val(sample_id), path("${sample_id}_samtools_stats_insert_sizes.tsv"), emit: insert_sizes
    tuple val(sample_id), path("${sample_id}_samtools_stats_coverage_distribution.tsv"), emit: coverage_distribution
    tuple val(sample_id), path("${sample_id}_samtools_stats_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: samtools_stats\\n" >> ${sample_id}_samtools_stats_provenance.yml
    printf -- "  tools:\\n"                       >> ${sample_id}_samtools_stats_provenance.yml
    printf -- "    - tool_name: samtools\\n"      >> ${sample_id}_samtools_stats_provenance.yml
    printf -- "      tool_version: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_samtools_stats_provenance.yml
    printf -- "      subcommand: stats\\n"        >> ${sample_id}_samtools_stats_provenance.yml

    samtools stats \
	--threads ${task.cpus} \
	${alignment[0]} > ${sample_id}_samtools_stats.txt

    grep '^SN' ${sample_id}_samtools_stats.txt | cut -f 2-  > ${sample_id}_samtools_stats_summary.txt

    parse_samtools_stats_summary.py -i ${sample_id}_samtools_stats_summary.txt -s ${sample_id} > ${sample_id}_samtools_stats_summary.csv

    echo "insert_size,pairs_total,inward_oriented_pairs,outward_oriented_pairs,other_pairs" | tr ',' '\t' > ${sample_id}_samtools_stats_insert_sizes.tsv
    grep '^IS' ${sample_id}_samtools_stats.txt | cut -f 2-  >> ${sample_id}_samtools_stats_insert_sizes.tsv

    echo "coverage,depth" | tr ',' '\t' > ${sample_id}_samtools_stats_coverage_distribution.tsv
    grep '^COV' ${sample_id}_samtools_stats.txt | cut -f 2- >> ${sample_id}_samtools_stats_coverage_distribution.tsv	
    """
}


process samtools_mpileup {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_depths.tsv"

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_depths.tsv"), emit: depths
    tuple val(sample_id), path("${sample_id}_samtools_mpileup_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: samtools_mpileup\\n" >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "  tools:\\n"                         >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "    - tool_name: samtools\\n"        >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "      tool_version: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "      subcommand: mpileup\\n"        >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "      parameters:\\n"                >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: -a\\n"          >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "          value: null\\n"            >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: --min-BQ\\n"    >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "          value: 0\\n"               >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "        - parameter: --count-orphans\\n" >> ${sample_id}_samtools_mpileup_provenance.yml
    printf -- "          value: null\\n"                >> ${sample_id}_samtools_mpileup_provenance.yml

    samtools faidx ${ref}

    printf "chrom\tpos\tref\tdepth\n" > ${sample_id}_depths.tsv

    samtools mpileup -a \
	--fasta-ref ${ref} \
	--min-BQ 0 \
	--count-orphans \
	${alignment[0]} | cut -f 1-4 >> ${sample_id}_depths.tsv
    """
}


process amplicon_coverage {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_amplicon_coverage.tsv"

    input:
    tuple val(sample_id), path(alignment), path(bedfile)

    output:
    tuple val(sample_id), path("${sample_id}_amplicon_coverage.tsv"), emit: depths

    script:
    """
    make_amplicon_bed.py --primer-bed ${bedfile} > amplicons.bed

    echo -e "reference_name\tstart\tend\tamplicon_id\tmean_depth" > ${sample_id}_amplicon_coverage.tsv
    
    bedtools coverage \
	-mean \
	-a amplicons.bed \
	-b ${alignment[0]} \
	>> ${sample_id}_amplicon_coverage.tsv
    """
}


process call_variants {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.variants.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.variants.tsv")

    script:
    """
    samtools faidx ${ref}

    samtools mpileup -aa -A -d ${params.max_depth} -B -Q 0 --reference ${ref} ${alignment[0]} \
	| ivar variants \
	-r ${ref} \
	-m ${params.min_depth}  \
	-q ${params.min_qual_for_variant_calling} \
	-t ${params.ambiguous_allele_freq_threshold} \
	-p ${sample_id}.variants
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

    samtools mpileup -aa -A -B -d ${params.max_depth} -Q 0 ${alignment[0]} \
	| ivar consensus \
	-t ${params.unambiguous_allele_freq_threshold} \
	-m ${params.min_depth} \
        -n N \
	-p ${sample_id}.primertrimmed.consensus
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

process plot_coverage {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_coverage.png"

    input:
    tuple val(sample_id), path(depths), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_coverage.png"), optional: true

    script:
    log_scale = params.coverage_plot_log_scale ? "--log-scale" : ""
    """
    plot-coverage.py \
	--sample-id ${sample_id} \
	--ref ${ref} \
	--depths ${depths} \
	--threshold ${params.min_depth} \
	--y-limit ${params.coverage_plot_y_limit} \
	--width-inches-per-mb ${params.coverage_plot_width_inches_per_mb} \
	--height-inches-per-chrom ${params.coverage_plot_height_inches_per_chrom} \
	--window ${params.coverage_plot_window_size} \
	${log_scale} \
	--output ${sample_id}_coverage.png
    """
}


process plot_amplicon_coverage {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_amplicon_coverage.png"

    input:
    tuple val(sample_id), path(depths)

    output:
    tuple val(sample_id), path("${sample_id}_amplicon_coverage.png"), optional: true

    script:
    log_scale = params.coverage_plot_log_scale ? "--log-scale" : ""
    """
    plot-amplicon-coverage.py \
	--sample-id ${sample_id} \
	--depths ${depths} \
	--threshold ${params.min_depth} \
	--y-limit ${params.coverage_plot_y_limit} \
	${log_scale} \
	--output ${sample_id}_amplicon_coverage.png
    """
}
