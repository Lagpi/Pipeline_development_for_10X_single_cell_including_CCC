#!/usr/bin/env nextflow

/*
 * Nextflow pipeline for STARsolo single-cell RNA-seq processing
 */

nextflow.enable.dsl = 2

// Help message
if (params.help) {
    log.info """
    STARsolo Nextflow Pipeline (Complete)
    ====================================
    
    Usage:
        nextflow run main.nf [options]
    
    Options:
        --fastq_dir         Input directory with FASTQ files [${params.fastq_dir}]
        --output_dir        Output directory [${params.output_dir}]
        --genome_dir        STAR genome index directory [${params.genome_dir}]
        --whitelist         Cell barcode whitelist file [${params.whitelist}]
        --threads           Number of threads [${params.threads}]
        --conditions        Conditions to process as comma-separated list [${params.conditions}]
        --help              Show this help message
    
    Example:
        nextflow run main.nf --conditions "condition1,condition2" --threads 16
    """.stripIndent()
    
    exit 0
}

// =============================================================================
// PROCESS DEFINITIONS  
// =============================================================================

process MAKE_SAMPLESHEET {
    publishDir "${params.output_dir}", mode: 'copy', pattern: "sample_manifest.txt"
    cache 'lenient'
    
    cpus 1
    memory '4 GB'
    time '10m'
    
    output:
    path "sample_manifest.txt", emit: samplesheet
    path "samplesheet_creation.log", emit: error_report, optional: true
    
    script:
    """
    #!/bin/bash
    echo "Creating samplesheet from input directory structure..."
    printf "sample_name\\tcondition\\toriginal_sample\\tfastq_R1\\tfastq_R2\\n" > sample_manifest.txt
    
    # Scan input directory structure: condition/sample
    for condition_dir in ${params.fastq_dir}/*/; do
        if [ -d "\$condition_dir" ]; then
            condition=\$(basename "\$condition_dir")
            for sample_dir in "\$condition_dir"/*/; do
                if [ -d "\$sample_dir" ]; then
                    original_sample=\$(basename "\$sample_dir")
                    sample_name="\${condition}_\${original_sample}"
                    
                    # Find FASTQ files
                    fastq_r1=\$(ls "\$sample_dir"/*_R1_*.fastq.gz 2>/dev/null | head -n1)
                    fastq_r2=\$(ls "\$sample_dir"/*_R2_*.fastq.gz 2>/dev/null | head -n1)
                    
                    if [ -n "\$fastq_r1" ] && [ -n "\$fastq_r2" ]; then
                        printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$sample_name" "\$condition" "\$original_sample" "\$fastq_r1" "\$fastq_r2" >> sample_manifest.txt
                    fi
                fi
            done
        fi
    done
    
    n_samples=\$(tail -n +2 sample_manifest.txt | wc -l)
    echo "Samplesheet created with \$n_samples samples" > samplesheet_creation.log
    echo "Samplesheet created successfully with \$n_samples samples"
    """
}

process STARSOLO_SAMPLE {
    tag "STAR_${sample_name}"
    conda params.conda_starsolo_env
    cache 'deep'

    publishDir "${params.output_dir}/00_starsolo/${condition}/${sample}", mode: 'copy', pattern: "*_Solo.out/**"
    publishDir "${params.output_dir}/00_starsolo/logs", mode: 'copy', pattern: "*.log"

    cpus params.starsolo_threads
    memory '120 GB'
    time '12h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 2

    input:
    tuple val(sample_name), val(condition), val(sample), path(fastq_R1), path(fastq_R2)

    output:
    tuple val(condition), val(sample_name), path("${sample_name}_Solo.out/Gene/filtered/*"), emit: results
    path "timing_${sample_name}.log", emit: timing
    path "${sample_name}_Solo.out/**", emit: all_outputs

    script:
    """
    # Use the STARsolo module with nextflow config parameters
    ${workflow.projectDir}/Scripts/utils/00_starsolo.sh \\
        "${params.genome_dir}" \\
        "${fastq_R1}" \\
        "${fastq_R2}" \\
        "${sample_name}" \\
        "." \\
        "${params.whitelist}" \\
        "${task.cpus}" \\
        "${params.soloUMIlen}" \\
        "${params.soloCBlen}" \\
        "${params.soloUMIstart}" \\
        "${params.soloMemoryLimit}"
    
    # Create timing log
    if [ -f ${sample_name}_starsolo.log ]; then
        cp ${sample_name}_starsolo.log timing_${sample_name}.log
    else
        echo "STARsolo completed for ${sample_name}" > timing_${sample_name}.log
    fi
    """
}

process QC_INTEGRATION {
    conda params.conda_seurat_env
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 4
    memory '32 GB'
    time '2h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path samplesheet
    val starsolo_ready
    path error_report

    output:
    path "01_qc_integration/data/seurat_qc.rds", emit: seurat_qc
    path "01_qc_integration/data/seurat_normalized.rds", emit: seurat_normalized
    path "01_qc_integration/data/seurat_normalized_hvg.rds", emit: seurat_normalized_hvg
    path "01_qc_integration/data/integrated_seurat.rds", emit: integrated_seurat
    path "01_qc_integration/data/integrated_stacas.rds", emit: integrated_stacas
    path "01_qc_integration/data/integrated_primary.rds", emit: integrated_primary
    path "01_qc_integration/**", emit: all_outputs

    script:
    """
    echo "Starting QC integration after STARsolo completion..."
    echo "All STARsolo tasks completed: ${starsolo_ready}"
    
    Rscript -e "library(Seurat); cat('Seurat version:', as.character(packageVersion('Seurat')), '\\n')"
    
    Rscript ${projectDir}/Scripts/01_qc_integration.R \
        --samplesheet ${samplesheet} \
        --starsolo_dir ${params.output_dir}/00_starsolo \
        --output_dir 01_qc_integration
    """
}

process CLUSTERING_ANNOTATION {
    conda params.conda_annotation_env  // Use new environment with SingleR 2.x
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 8
    memory '240 GB'  // Increased for marker gene finding with ~78K cells
    time '8h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 || task.exitStatus == 137 || task.exitStatus == 140 ? 'ignore' : 'retry' }
    maxRetries 2

    input:
    path integrated_primary
    path singler_ref

    output:
    path "02_clustering/clustered_annotated_seurat.rds", emit: clustered_rds
    path "02_clustering/**", emit: all_outputs

    script:
    """
    echo "Starting clustering and annotation..."
    
    Rscript ${projectDir}/Scripts/02_clustering_annotation.R \
        --input_rds ${integrated_primary} \
        --output_dir 02_clustering \
        --singler_ref ${singler_ref} \\
        --pca_search_space ${params.pca_search_space} \\
        --resolutions "${params.resolutions}" \\
        --umap_neighbors ${params.umap_neighbors} \\
        --species ${params.species}
    """
}

process DIFFERENTIAL_ANALYSIS {
    tag "Differential"
    conda params.conda_seurat_env
    cache 'lenient'
    
    publishDir "${params.output_dir}/03_differential/plots", mode: 'copy', pattern: "03_differential/plots/*"
    publishDir "${params.output_dir}/03_differential/interactive", mode: 'copy', pattern: "03_differential/interactive/*"
    publishDir "${params.output_dir}/03_differential/data", mode: 'copy', pattern: "03_differential/data/*"

    cpus 8
    memory '250 GB'
    time '12h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path clustered_rds

    output:
    path "03_differential/data/analysis_summary.rds", emit: diff_results
    path "03_differential/interactive/*.html", optional: true, emit: diff_report
    path "03_differential/plots/*", optional: true, emit: diff_plots
    path "03_differential/**", emit: all_outputs

    script:
    """
    echo "Starting differential analysis..."
    
    Rscript ${projectDir}/Scripts/03_differential_analysis.R \
        --seurat ${clustered_rds} \
        --out 03_differential \
        --milo_fdr ${params.milo.fdr_threshold} \\
        --milo_fdr_weighting ${params.milo.fdr_weighting} \\
        --seurat_test_use ${params.seurat.findmarkers.test_use} \\
        --seurat_logfc_threshold ${params.seurat.findmarkers.logfc_threshold} \\
        --seurat_min_pct ${params.seurat.findmarkers.min_pct} \\
        --seurat_plot_fdr ${params.seurat.findmarkers.plot_fdr} \\
        --seurat_plot_logfc ${params.seurat.findmarkers.plot_logfc} \\
        --muscat_fdr ${params.muscat.fdr} \\
        --muscat_logfc ${params.muscat.logfc}
    """
}

process CELLCHAT_ANALYSIS {
    conda params.conda_cellchat_v2
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 8  
    memory '120 GB'
    time '6h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path clustered_rds

    output:
    path "04_cellcomm/data/cellchat_summary.rds", emit: cellchat_results
    path "04_cellcomm/plots/*", emit: cellchat_plots, optional: true
    path "04_cellcomm/rds/*", emit: cellchat_objects, optional: true
    path "04_cellcomm/**", emit: all_outputs

    script:
    """
    echo "Starting CellChat analysis..."
    
    Rscript ${projectDir}/Scripts/04a_cellchat.R \
        --seurat ${clustered_rds} \
        --out 04_cellcomm \
        --species ${params.species ?: 'mouse'}
    """
}

process NICHENET_ANALYSIS {
    conda params.conda_cellcomm_multi
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 8  
    memory '80 GB'
    time '4h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path clustered_rds
    path normalized_full_rds
    path cellchat_summary
    val diff_done  // Wait for differential analysis to complete

    output:
    path "04_cellcomm/data/nichenet_summary.rds", emit: nichenet_results
    path "04_cellcomm/data/analysis_summary.rds", emit: analysis_summary
    path "04_cellcomm/**", emit: all_outputs

    script:
    """
    echo "Starting NicheNet analysis (official vignette approach)..."
    echo "CellChat completed, now running NicheNet..."
    echo "Using normalized Seurat with all genes for expression data..."
    
    # Copy CellChat summary to expected location
    mkdir -p 04_cellcomm/data
    cp ${cellchat_summary} 04_cellcomm/data/cellchat_summary.rds
    
    Rscript ${projectDir}/Scripts/04b_nichenet.R \
        --seurat ${clustered_rds} \
        --normalized ${normalized_full_rds} \
        --out 04_cellcomm \
        --species ${params.species ?: 'mouse'} \
        --nichenet_dir ${params.nichenet_dir}
    
    # Create final analysis summary combining both
    echo "Creating combined analysis summary..."
    Rscript -e "
    cellchat_sum <- readRDS('04_cellcomm/data/cellchat_summary.rds')
    nichenet_sum <- readRDS('04_cellcomm/data/nichenet_summary.rds')
    final_sum <- list(
      status = 'completed',
      time = Sys.time(),
      cellchat = cellchat_sum,
      nichenet = nichenet_sum
    )
    saveRDS(final_sum, '04_cellcomm/data/analysis_summary.rds')
    "
    """
}

process GENE_PROGRAM_DISCOVERY {
    conda params.conda_nmf_env
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    cpus 8
    memory = '200 GB'  // Increased memory for NMF with 80k cells

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path normalized_hvg
    path clustered_rds

    output:
    path "05_gene_programs/data/seurat_with_gene_programs.rds", emit: seurat_with_gp
    path "05_gene_programs/data/nmf_result.rds", emit: nmf_results
    path "05_gene_programs/**", emit: all_outputs

    script:
    """
    echo "Starting gene program discovery..."
    
    Rscript ${projectDir}/Scripts/05_gene_program_discovery.R \
        --normalized ${normalized_hvg} \
        --annotations ${clustered_rds} \
        --out 05_gene_programs \
        --min_rank ${params.nmf_min_rank} \\
        --max_rank ${params.nmf_max_rank}
    """
}

process TRAJECTORY_ANALYSIS {
    conda params.conda_trajectory_env
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    errorStrategy 'ignore'
    maxRetries 0

    input:
    path clustered_rds

    output:
    path "06_trajectory/data/*.rds", emit: trajectory_results, optional: true
    path "06_trajectory/**", emit: all_outputs

    script:
    """
    echo "Starting trajectory analysis..."
    
    Rscript ${projectDir}/Scripts/06_trajectory.R \
        --seurat ${clustered_rds} \
        --out 06_trajectory
    """
}

process BENCHMARKING {
    conda params.conda_benchmarking_env
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 4
    memory '150 GB'
    time '4h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 || task.exitStatus == 137 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    val ready  // Just a signal that previous steps are done

    output:
    path "07_benchmarking/**", emit: all_outputs

    script:
    """
    echo "Starting benchmarking analysis..."
    echo "Previous steps completed: ${ready}"
    
    Rscript ${projectDir}/Scripts/07_benchmarking.R \\
        --input_dir ${params.output_dir} \\
        --out 07_benchmarking
    """
}

process FINAL_DASHBOARD {
    conda params.conda_seurat_env
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode: 'copy'

    cpus 2
    memory '16 GB'
    time '1h'

    errorStrategy { task.exitStatus == 104 || task.exitStatus == 105 ? 'ignore' : 'retry' }
    maxRetries 1

    input:
    path "qc_outputs/**"
    path "clustering_outputs/**"
    path "diff_outputs/**"
    path "all_other_outputs/**"

    output:
    path "08_reports/comprehensive_analysis_dashboard.html", emit: dashboard
    path "08_reports/dashboard_data.rds", emit: dashboard_data

    script:
    """
    echo "Creating comprehensive dashboard..."
    
    Rscript ${projectDir}/Scripts/08_dashboard.R \\
        --out 08_reports
    """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {
    
    // Parse --conditions parameter into list (or empty = all)
    def cond_list = params.conditions && params.conditions != 'auto' ? params.conditions.tokenize(',').collect { it.trim() } : []
    def allowed = cond_list as Set

    log.info """
    STARsolo Nextflow Pipeline (Complete)
    ====================================
    fastq_dir    : ${params.fastq_dir}
    output_dir   : ${params.output_dir}
    genome_dir   : ${params.genome_dir}
    whitelist    : ${params.whitelist}
    conditions   : ${params.conditions ?: 'ALL'} (${cond_list.size() == 0 ? 'all' : cond_list.size()} selected)
    threads      : ${params.threads}
    """.stripIndent()

    // Create samplesheet from input directory structure
    SAMPLESHEET = MAKE_SAMPLESHEET()

    // Discover samples and create a channel
    samples_ch = Channel.fromPath("${params.fastq_dir}/*/*", type: 'dir')
                    .map { p -> 
                        def condition = p.getParent().getFileName().toString()
                        def sample = p.getFileName().toString()
                        def sample_name = "${condition}_${sample}"
                        log.info "Found sample: ${sample_name} (condition: ${condition}, sample: ${sample})"
                        tuple(sample_name, condition, sample, file("${p}/*_R1_*.fastq.gz"), file("${p}/*_R2_*.fastq.gz"))
                    }
                    .branch {
                        selected: allowed.size() == 0 || allowed.contains(it[1])
                        other: true
                    }
                    .selected

    // Run STARsolo per-sample in parallel
    STAR_TASKS = STARSOLO_SAMPLE(samples_ch)

    // Run QC and Integration after all STAR tasks finished
    QC_INTEGRATION(
        SAMPLESHEET.samplesheet,
        STAR_TASKS.results.count().map { "ready_${it}_samples" },
        SAMPLESHEET.error_report
    )


    // Run Clustering & Annotation (Step 02) using QC output
    singler_ref = Channel.fromPath(params.singler_ref)
    CLUST_OUT = CLUSTERING_ANNOTATION(
        QC_INTEGRATION.out.integrated_primary,
        singler_ref
    )

    // Run Differential Analysis (Step 03) using clustered Seurat object
    DIFF_OUT = DIFFERENTIAL_ANALYSIS(
        CLUST_OUT.clustered_rds
    )

    // Run Intercellular Communication analysis (Step 04) - split into CellChat and NicheNet
    CELLCHAT_OUT = CELLCHAT_ANALYSIS(
        CLUST_OUT.clustered_rds
    )
    
    NICHENET_OUT = NICHENET_ANALYSIS(
        CLUST_OUT.clustered_rds,
        QC_INTEGRATION.out.seurat_normalized,
        CELLCHAT_OUT.cellchat_results,
        DIFF_OUT.all_outputs.collect().map { "ready" }
    )
    
    // Combine outputs for downstream compatibility
    CELLCOMM_OUT = NICHENET_OUT

    // Run Gene Program Discovery (Step 05) - waits for DIFF
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        QC_INTEGRATION.out.seurat_normalized_hvg,
        CLUST_OUT.clustered_rds
    )

    // Run Trajectory Analysis (Step 06) - uses Gene Programs AND waits for DIFF
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )

    // Run Benchmarking (Step 07) - waits for all other analyses to complete
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard ONLY after ALL analysis steps complete
    FINAL_DASHBOARD(
        QC_INTEGRATION.out.all_outputs,
        CLUST_OUT.all_outputs,
        DIFF_OUT.all_outputs,
        CELLCOMM_OUT.all_outputs.concat(
            GENE_PROGRAM_OUT.all_outputs,
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

// Alternative workflow: Start from QC (skip STARsolo)
workflow from_qc {
    // Load existing samplesheet  
    samplesheet = Channel.fromPath("${params.output_dir}/sample_manifest.txt").first()
    error_report = Channel.fromPath("${params.output_dir}/sample_manifest.txt").map { file("${params.output_dir}/samplesheet_creation.log") }.ifEmpty(file('NO_FILE'))
    
    // Run QC Integration
    QC_INTEGRATION(
        samplesheet,
        "starsolo_completed",
        error_report
    )
    
    // Run Clustering & Annotation
    singler_ref = Channel.fromPath(params.singler_ref)
    CLUST_OUT = CLUSTERING_ANNOTATION(
        QC_INTEGRATION.out.integrated_primary,
        singler_ref
    )
    
    // Run Differential Analysis
    DIFF_OUT = DIFFERENTIAL_ANALYSIS(
        CLUST_OUT.clustered_rds
    )
    
    // Run Intercellular Communication (Step 04) - split into CellChat and NicheNet
    CELLCHAT_OUT = CELLCHAT_ANALYSIS(
        CLUST_OUT.clustered_rds
    )
    
    NICHENET_OUT = NICHENET_ANALYSIS(
        CLUST_OUT.clustered_rds,
        CELLCHAT_OUT.all_outputs.concat(DIFF_OUT.all_outputs).collect().map { "ready" }
    )
    
    // Combine outputs for downstream compatibility
    CELLCOMM_OUT = NICHENET_OUT
    
    // Run Gene Program Discovery (Step 05) - waits for CELLCOMM
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        CELLCOMM_OUT.all_outputs.collect().map { QC_INTEGRATION.out.seurat_normalized_hvg },
        CELLCOMM_OUT.all_outputs.collect().map { CLUST_OUT.clustered_rds }
    )
    
    // Run Trajectory Analysis (Step 06) - waits for GENE_PROGRAM
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )
    
    // Run Benchmarking
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )
    
    // Generate dashboard
    FINAL_DASHBOARD(
        QC_INTEGRATION.out.all_outputs,
        CLUST_OUT.all_outputs,
        DIFF_OUT.all_outputs,
        CELLCOMM_OUT.all_outputs.concat(
            GENE_PROGRAM_OUT.all_outputs,
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

// Alternative workflow: Start from clustering (skip STARsolo and QC)
workflow from_clustering {
    // Load existing QC output
    qc_seurat = Channel.fromPath("${params.output_dir}/01_qc_integration/data/integrated_primary.rds")
    singler_ref = Channel.fromPath(params.singler_ref)
    
    // Run Clustering & Annotation (Step 02)
    CLUST_OUT = CLUSTERING_ANNOTATION(
        qc_seurat,
        singler_ref
    )

    // Run Differential Analysis (Step 03)
    DIFF_OUT = DIFFERENTIAL_ANALYSIS(
        CLUST_OUT.clustered_rds
    )

    // Run Intercellular Communication analysis (Step 04) - split into CellChat and NicheNet
    CELLCHAT_OUT = CELLCHAT_ANALYSIS(
        CLUST_OUT.clustered_rds
    )
    
    NICHENET_OUT = NICHENET_ANALYSIS(
        CLUST_OUT.clustered_rds,
        CELLCHAT_OUT.cellchat_results,
        DIFF_OUT.all_outputs.collect().map { "ready" }
    )
    
    // Combine outputs for downstream compatibility
    CELLCOMM_OUT = NICHENET_OUT

    // Run Gene Program Discovery (Step 05) - waits for CELLCOMM
    // Load normalized HVG from QC output since QC_INTEGRATION is not in this workflow
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        CELLCOMM_OUT.all_outputs.collect().map { file("${params.output_dir}/01_qc_integration/data/seurat_normalized_hvg.rds") },
        CELLCOMM_OUT.all_outputs.collect().map { CLUST_OUT.clustered_rds }
    )

    // Run Trajectory Analysis (Step 06) - uses Seurat with Gene Programs
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )

    // Run Benchmarking (Step 07)
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.empty(),  // Skip QC outputs
        CLUST_OUT.all_outputs,
        DIFF_OUT.all_outputs,
        CELLCOMM_OUT.all_outputs.concat(
            GENE_PROGRAM_OUT.all_outputs,
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

// Alternative workflow: Start from differential analysis
workflow from_differential {
    // Load existing clustering output and create a value channel that can be reused
    clustered_file = file("${params.output_dir}/02_clustering/clustered_annotated_seurat.rds")
    
    // Run Differential Analysis (Step 03)
    DIFF_OUT = DIFFERENTIAL_ANALYSIS(
        Channel.fromPath(clustered_file)
    )

    // Run Intercellular Communication analysis (Step 04) - split into CellChat and NicheNet
    CELLCHAT_OUT = CELLCHAT_ANALYSIS(
        Channel.fromPath(clustered_file)
    )
    
    NICHENET_OUT = NICHENET_ANALYSIS(
        Channel.fromPath(clustered_file),
        CELLCHAT_OUT.cellchat_results,
        DIFF_OUT.all_outputs.collect().map { "ready" }
    )
    
    // Combine outputs for downstream compatibility
    CELLCOMM_OUT = NICHENET_OUT

    // Run Gene Program Discovery (Step 05) - waits for CELLCOMM
    // Load normalized HVG and clustered files
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        CELLCOMM_OUT.all_outputs.collect().map { file("${params.output_dir}/01_qc_integration/data/seurat_normalized_hvg.rds") },
        CELLCOMM_OUT.all_outputs.collect().map { clustered_file }
    )

    // Run Trajectory Analysis (Step 06) - uses Seurat with Gene Programs
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )

    // Run Benchmarking (Step 07)
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.empty(),  // Skip QC outputs
        Channel.empty(),  // Skip clustering outputs
        DIFF_OUT.all_outputs,
        CELLCOMM_OUT.all_outputs.concat(
            GENE_PROGRAM_OUT.all_outputs,
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

// Alternative workflow: Start from cell communication analysis
workflow from_cellcomm {
    // Load existing clustering and differential analysis outputs
    clustered_file = file("${params.output_dir}/02_clustering/clustered_annotated_seurat.rds")
    normalized_full_file = file("${params.output_dir}/01_qc_integration/data/seurat_normalized.rds")
    diff_summary = file("${params.output_dir}/03_differential")
    normalized_hvg_file = file("${params.output_dir}/01_qc_integration/data/seurat_normalized_hvg.rds")
    
    // Run Intercellular Communication analysis (Step 04) - split into CellChat and NicheNet
    CELLCHAT_OUT = CELLCHAT_ANALYSIS(
        Channel.fromPath(clustered_file)
    )
    
    NICHENET_OUT = NICHENET_ANALYSIS(
        Channel.fromPath(clustered_file),
        Channel.fromPath(normalized_full_file),
        CELLCHAT_OUT.cellchat_results,
        Channel.fromPath(diff_summary).collect().map { "ready" }
    )
    
    // Combine outputs for downstream compatibility
    CELLCOMM_OUT = NICHENET_OUT

    // Run Gene Program Discovery (Step 05) - waits for CELLCOMM completion
    // Load both required files from disk
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        CELLCOMM_OUT.all_outputs.collect().map { normalized_hvg_file },
        CELLCOMM_OUT.all_outputs.collect().map { clustered_file }
    )

    // Run Trajectory Analysis (Step 06) - uses Seurat with Gene Programs
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )

    // Run Benchmarking (Step 07)
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.fromPath("${params.output_dir}/01_qc_integration").collect().ifEmpty([]),
        Channel.fromPath("${params.output_dir}/02_clustering").collect().ifEmpty([]),
        Channel.fromPath("${params.output_dir}/03_differential").collect().ifEmpty([]),
        CELLCOMM_OUT.all_outputs.concat(
            GENE_PROGRAM_OUT.all_outputs,
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

// Alternative workflow: Start from NicheNet analysis only
workflow from_nichenet {
    // Load existing files
    clustered_file = file("${params.output_dir}/02_clustering/clustered_annotated_seurat.rds")
    normalized_full_file = file("${params.output_dir}/01_qc_integration/data/seurat_normalized.rds")
    diff_summary = file("${params.output_dir}/03_differential")
    cellchat_summary = file("${params.output_dir}/04_cellcomm/data/cellchat_summary.rds")
    
    // Run only NicheNet analysis
    NICHENET_OUT = NICHENET_ANALYSIS(
        Channel.fromPath(clustered_file),
        Channel.fromPath(normalized_full_file),
        Channel.fromPath(cellchat_summary),
        Channel.fromPath(diff_summary).collect().map { "ready" }
    )
}

workflow from_geneprog {
    // Load existing files from disk
    normalized_hvg = file("${params.output_dir}/01_qc_integration/data/seurat_normalized_hvg.rds")
    clustered_file = file("${params.output_dir}/02_clustering/clustered_annotated_seurat.rds")
    
    // Run Gene Program Discovery (Step 05)
    GENE_PROGRAM_OUT = GENE_PROGRAM_DISCOVERY(
        Channel.fromPath(normalized_hvg),
        Channel.fromPath(clustered_file)
    )

    // Run Trajectory Analysis (Step 06) - uses Seurat with Gene Programs
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        GENE_PROGRAM_OUT.seurat_with_gp
    )

    // Run Benchmarking (Step 07)
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.empty(),  // Skip QC outputs
        Channel.empty(),  // Skip clustering outputs
        Channel.fromPath("${params.output_dir}/03_differential").collect().ifEmpty([]),
        GENE_PROGRAM_OUT.all_outputs.concat(
            TRAJECTORY_OUT.all_outputs,
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

workflow from_trajectory {
    // Load existing clustering output
    clustered_file = file("${params.output_dir}/02_clustering/clustered_annotated_seurat.rds")
    
    // Run Trajectory Analysis (Step 06)
    TRAJECTORY_OUT = TRAJECTORY_ANALYSIS(
        Channel.fromPath(clustered_file)
    )

    // Run Benchmarking (Step 07)
    BENCHMARK_OUT = BENCHMARKING(
        TRAJECTORY_OUT.all_outputs.collect().ifEmpty('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.empty(),  // Skip QC outputs
        Channel.empty(),  // Skip clustering outputs
        Channel.fromPath("${params.output_dir}/03_differential").collect().ifEmpty([]),
        TRAJECTORY_OUT.all_outputs.concat(
            BENCHMARK_OUT.all_outputs
        ).collect()
    )
}

workflow from_benchmark {
    // Run Benchmarking only (Step 07) - loads all previous outputs from disk
    BENCHMARK_OUT = BENCHMARKING(
        Channel.of('ready')
    )

    // Generate comprehensive dashboard
    FINAL_DASHBOARD(
        Channel.fromPath("${params.output_dir}/01_qc_integration").collect().ifEmpty([]),
        Channel.fromPath("${params.output_dir}/02_clustering").collect().ifEmpty([]),
        Channel.fromPath("${params.output_dir}/03_differential").collect().ifEmpty([]),
        Channel.fromPath("${params.output_dir}/04_cellcomm").collect()
            .concat(
                Channel.fromPath("${params.output_dir}/05_gene_programs").collect().ifEmpty([]),
                Channel.fromPath("${params.output_dir}/06_trajectory").collect().ifEmpty([]),
                BENCHMARK_OUT.all_outputs
            ).collect()
    )
}

// Workflow completion
workflow.onComplete {
    log.info """
    🎉 Pipeline completed successfully!
    Output directory: ${params.output_dir}
    
    📊 Key outputs:
    - STARsolo alignment: ${params.output_dir}/00_starsolo/
    - QC & Integration: ${params.output_dir}/01_qc_integration/
    - Clustering & Annotation: ${params.output_dir}/02_clustering/
    - Differential Analysis: ${params.output_dir}/03_differential/
    - Cell Communication: ${params.output_dir}/04_cellcomm/
    - Gene Program Discovery: ${params.output_dir}/05_gene_programs/
    - Trajectory Analysis: ${params.output_dir}/06_trajectory/
    - Benchmarking: ${params.output_dir}/07_benchmarking/
    - Final Dashboard: ${params.output_dir}/08_reports/
    
    📈 Comprehensive Analysis Report:
    - HTML Dashboard: ${params.output_dir}/08_reports/comprehensive_analysis_dashboard.html
    
    Open the HTML report in your browser to explore all results!
    """.stripIndent()
}