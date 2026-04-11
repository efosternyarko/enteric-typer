#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * enteric-typer: species-gated genotyping workflow for Enterobacteriaceae
 *
 * Supports:  Escherichia coli  |  Salmonella enterica  |  Shigella spp.
 *
 * Phase 1 – species gate (Mash) against reference sketch
 * Phase 2 – species-specific typing in parallel
 * Phase 2c– K-locus typing (E. coli only): G2/G3 → G1/G4 on untypeables
 * Phase 3 – SKA2 core-SNP alignment + IQ-TREE per species (skip with --skip_local_phylo)
 * Phase 4 – Aggregation (one TSV per species)
 * Phase 5 – Summary plots + tree annotation + SNP distance heatmap
 */

include { SPECIES_CHECK                             } from './modules/species_check'
include { MLST          as MLST_ECOLI               } from './modules/mlst'
include { MLST          as MLST_SALMONELLA          } from './modules/mlst'
include { AMRFINDER     as AMRFINDER_ECOLI          } from './modules/amrfinder'
include { AMRFINDER     as AMRFINDER_SALMONELLA     } from './modules/amrfinder'
include { ECTYPER                                   } from './modules/ectyper'
include { KLEBORATE                                 } from './modules/kleborate'
include { EZCLERMONT                               } from './modules/ezclermont'
include { MYKROBE                                   } from './modules/mykrobe'
include { PARSE_MYKROBE                             } from './modules/mykrobe'
include { SHIGEIFINDER                              } from './modules/shigeifinder'
include { PINV_SCREEN                               } from './modules/pinv_screen'
include { IS_SCREEN                                 } from './modules/is_screen'
include { KAPTIVE_G2G3                              } from './modules/kaptive'
include { KAPTIVE_G1G4                              } from './modules/kaptive'
include { PARSE_KAPTIVE                             } from './modules/kaptive'
include { SISTR                                     } from './modules/sistr'
include { PLASMIDFINDER as PLASMIDFINDER_ECOLI      } from './modules/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SALMONELLA } from './modules/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHIGELLA   } from './modules/plasmidfinder'
include { SKA2_BUILD    as SKA2_ECOLI               } from './modules/ska2'
include { SKA2_BUILD    as SKA2_SALMONELLA          } from './modules/ska2'
include { SKA2_BUILD    as SKA2_SHIGELLA            } from './modules/ska2'
include { IQTREE        as IQTREE_ECOLI             } from './modules/iqtree'
include { IQTREE        as IQTREE_SALMONELLA        } from './modules/iqtree'
include { IQTREE        as IQTREE_SHIGELLA          } from './modules/iqtree'
include { AGGREGATE     as AGGREGATE_ECOLI          } from './modules/aggregate'
include { AGGREGATE     as AGGREGATE_SALMONELLA     } from './modules/aggregate'
include { AGGREGATE     as AGGREGATE_SHIGELLA       } from './modules/aggregate'
include { MLST          as MLST_SHIGELLA            } from './modules/mlst'
include { AMRFINDER     as AMRFINDER_SHIGELLA       } from './modules/amrfinder'
include { PLOT_SUMMARY      as PLOT_SUMMARY_ECOLI      } from './modules/plot_summary'
include { PLOT_SUMMARY      as PLOT_SUMMARY_SALMONELLA } from './modules/plot_summary'
include { PLOT_SUMMARY      as PLOT_SUMMARY_SHIGELLA   } from './modules/plot_summary'
include { TREE_ANNOTATION   as TREE_ANNOTATION_ECOLI      } from './modules/tree_annotation'
include { TREE_ANNOTATION   as TREE_ANNOTATION_SALMONELLA } from './modules/tree_annotation'
include { TREE_ANNOTATION   as TREE_ANNOTATION_SHIGELLA   } from './modules/tree_annotation'
include { SNP_HEATMAP       as SNP_HEATMAP_ECOLI          } from './modules/snp_heatmap'
include { SNP_HEATMAP       as SNP_HEATMAP_SALMONELLA      } from './modules/snp_heatmap'
include { SNP_HEATMAP       as SNP_HEATMAP_SHIGELLA        } from './modules/snp_heatmap'
include { ASSEMBLY_QC       as ASSEMBLY_QC_ECOLI            } from './modules/assembly_qc'
include { ASSEMBLY_QC       as ASSEMBLY_QC_SALMONELLA       } from './modules/assembly_qc'
include { ASSEMBLY_QC       as ASSEMBLY_QC_SHIGELLA         } from './modules/assembly_qc'
include { PLOT_ASSEMBLY_METRICS as PLOT_ASSEMBLY_METRICS_ECOLI      } from './modules/plot_assembly_metrics'
include { PLOT_ASSEMBLY_METRICS as PLOT_ASSEMBLY_METRICS_SALMONELLA } from './modules/plot_assembly_metrics'
include { PLOT_ASSEMBLY_METRICS as PLOT_ASSEMBLY_METRICS_SHIGELLA   } from './modules/plot_assembly_metrics'

// ── Help ──────────────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║                          enteric-typer                                  ║
    ║   Species-gated genotyping for E. coli, Salmonella & Shigella           ║
    ╚══════════════════════════════════════════════════════════════════════════╝

    USAGE
      nextflow run main.nf -profile <profile> [options]

    INPUT (one required)
      --input_dir   <dir>   Directory of FASTA assemblies (.fasta/.fa/.fna/.fas)
      --samplesheet <csv>   CSV with columns: id,fasta

    COMMON OPTIONS
      --outdir            <dir>    Output directory (default: results)
      --skip_local_phylo           Skip SKA2 + IQ-TREE phylogenetics
      --iqtree_model      <model>  IQ-TREE substitution model (default: GTR+G)
      --iqtree_bootstraps <n>      Ultrafast bootstrap replicates (default: 1000)
      --max_cpus          <n>      Max CPUs per process (default: 16)

    EXAMPLES
      # Linux / Intel Mac
      nextflow run main.nf -profile conda \\
          --input_dir /path/to/assemblies/

      # macOS Apple Silicon (M1/M2/M3/M4)
      CONDA_SUBDIR=osx-64 nextflow run main.nf -profile conda,arm64 \\
          --input_dir /path/to/assemblies/

      # HPC (SLURM + Singularity)
      nextflow run main.nf -profile singularity,slurm \\
          --input_dir /path/to/assemblies/

      # Skip phylogenetics (faster)
      nextflow run main.nf -profile conda \\
          --input_dir /path/to/assemblies/ \\
          --skip_local_phylo

      # Resume after interruption
      nextflow run main.nf -profile conda \\
          --input_dir /path/to/assemblies/ -resume

    PROFILES
      conda           Local — conda environments
      mamba           Local — mamba (faster env solving)
      arm64           Add on Apple Silicon (combine with conda/mamba)
      docker          Local — Docker
      singularity     HPC — Singularity/Apptainer
      slurm           HPC — SLURM executor
      pbs             HPC — PBS/Torque executor

    Full documentation: README.md
    """.stripIndent()
    exit 0
}

// ── Parameter validation ──────────────────────────────────────────────────────
if (!params.samplesheet && !params.input_dir) {
    error "ERROR: Provide --samplesheet <csv> or --input_dir <folder>.\nSee README.md for details."
}

// ── Pre-flight checks ─────────────────────────────────────────────────────────

// 1. Reference sketch
def sketchFile = file("${projectDir}/assets/enteric_species_refs.msh")
if (!sketchFile.exists()) {
    error """
    ════════════════════════════════════════════════════════════════
    ERROR: Mash reference sketch not found:
      ${sketchFile}

    Build it once before running the pipeline:
      bash assets/build_references.sh

    This downloads 7 reference genomes from NCBI (~5 min) and
    produces assets/enteric_species_refs.msh
    ════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// 2. Apple Silicon + conda without arm64 profile
def isArm64      = System.properties['os.arch'] == 'aarch64'
def activeProfiles = workflow.profile?.tokenize(',')*.trim() ?: []
def usingConda   = activeProfiles.any { it in ['conda', 'mamba'] }
def usingArm64   = activeProfiles.contains('arm64')
if (isArm64 && usingConda && !usingArm64) {
    log.warn """
    ════════════════════════════════════════════════════════════════
    WARNING: Apple Silicon detected without the arm64 profile.

    Some Bioconda packages have no native arm64 build and will fail
    to install, killing the pipeline.

    Rerun with the arm64 profile to install osx-64 envs via Rosetta:
      -profile ${workflow.profile},arm64
    ════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// 3. Input FASTA files exist (input_dir path check)
if (params.input_dir && !file(params.input_dir).exists()) {
    error """
    ════════════════════════════════════════════════════════════════
    ERROR: --input_dir path does not exist:
      ${params.input_dir}

    Check the path and try again.
    ════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// 4. Samplesheet exists
if (params.samplesheet && !file(params.samplesheet).exists()) {
    error """
    ════════════════════════════════════════════════════════════════
    ERROR: --samplesheet file not found:
      ${params.samplesheet}

    Generate one with:
      python bin/make_samplesheet.py --input /path/to/assemblies/ --output samples.csv
    ════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // ── Build input channel ───────────────────────────────────────────────────
    if (params.input_dir) {
        ch_samples = Channel
            .fromPath("${params.input_dir}/*.{fasta,fa,fna,fas,fsa}", checkIfExists: true)
            .map { fasta -> tuple(fasta.baseName.replaceAll(/\.(fasta|fa|fna|fas|fsa)$/, ''), fasta) }
    } else {
        ch_samples = Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def fasta = file(row.fasta)
                if (!fasta.exists()) error "FASTA not found for '${row.id}': ${row.fasta}"
                tuple(row.id, fasta)
            }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 1: Species identification (Mash)
    // ─────────────────────────────────────────────────────────────────────────
    ch_ref_sketch = file("${projectDir}/assets/enteric_species_refs.msh",
                          checkIfExists: true)
    SPECIES_CHECK(ch_samples, ch_ref_sketch)

    ch_classified = ch_samples
        .join(SPECIES_CHECK.out.species)
        .map { id, fasta, sp_file ->
            def cols    = sp_file.text.trim().split('\t')
            def species = cols[0]?.trim() ?: 'Unknown'
            def dist    = cols.size() > 1 ? cols[1].trim().toFloat() : 1.0f
            tuple(id, fasta, species, dist)
        }

    ch_classified.branch {
        ecoli:      it[2] == 'E_coli'              && it[3] < 0.05f
        salmonella: it[2] == 'Salmonella_enterica'  && it[3] < 0.05f
        shigella:   it[2] == 'Shigella'            && it[3] < 0.05f
        other:      true
    }.set { ch_branched }

    ch_branched.other
        .map { id, fasta, species, dist ->
            "  • ${id}: closest match = ${species} (Mash dist = ${String.format('%.4f', dist)})"
        }
        .collect()
        .subscribe { lines ->
            if (lines) {
                log.warn """
════════════════════════════════════════════════════════════════
WARNING: ${lines.size()} sample(s) did not pass the species gate
         and will NOT be typed (no results produced):

${lines.join('\n')}

Reason: Mash distance ≥ 0.05 to all reference genomes.
If any of these are non-sonnei Shigella (flexneri / boydii /
dysenteriae), the reference sketch needs rebuilding:
  bash assets/build_references.sh
════════════════════════════════════════════════════════════════""".stripIndent()
            }
        }

    ch_ecoli      = ch_branched.ecoli.map      { id, fasta, sp, d -> tuple(id, fasta) }
    ch_salmonella = ch_branched.salmonella.map { id, fasta, sp, d -> tuple(id, fasta) }
    ch_shigella   = ch_branched.shigella.map   { id, fasta, sp, d -> tuple(id, fasta) }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 1b: Assembly QC — per-species metrics figure
    // 8-panel figure: genome length, N50, contig count, GC% (hist + boxplot each)
    // ─────────────────────────────────────────────────────────────────────────
    ASSEMBLY_QC_ECOLI(ch_ecoli)
    ASSEMBLY_QC_SALMONELLA(ch_salmonella)
    ASSEMBLY_QC_SHIGELLA(ch_shigella)

    PLOT_ASSEMBLY_METRICS_ECOLI(
        ASSEMBLY_QC_ECOLI.out.stats.map { id, f -> f }.collect().ifEmpty([]),
        'ecoli'
    )
    PLOT_ASSEMBLY_METRICS_SALMONELLA(
        ASSEMBLY_QC_SALMONELLA.out.stats.map { id, f -> f }.collect().ifEmpty([]),
        'salmonella'
    )
    PLOT_ASSEMBLY_METRICS_SHIGELLA(
        ASSEMBLY_QC_SHIGELLA.out.stats.map { id, f -> f }.collect().ifEmpty([]),
        'shigella'
    )

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2a: E. coli typing (parallel per sample)
    // ─────────────────────────────────────────────────────────────────────────
    MLST_ECOLI(ch_ecoli,          'ecoli_achtman_4')
    AMRFINDER_ECOLI(ch_ecoli,     'Escherichia')
    ECTYPER(ch_ecoli)
    KLEBORATE(ch_ecoli)
    EZCLERMONT(ch_ecoli)
    PLASMIDFINDER_ECOLI(ch_ecoli, 'enterobacteriaceae')

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2c: E. coli K-locus typing
    // ─────────────────────────────────────────────────────────────────────────
    ch_g2g3_db        = file("${projectDir}/assets/kaptive_dbs/EC-K-typing_group2and3_v3.0.0.gbk",
                              checkIfExists: true)
    ch_g1g4_db        = file("${projectDir}/assets/kaptive_dbs/EC-K-typing_group1and4_v0.9.gbk",
                              checkIfExists: true)
    ch_normalise_script = file("${projectDir}/bin/normalise_kaptive_scores.py",
                              checkIfExists: true)

    KAPTIVE_G2G3(ch_ecoli, ch_g2g3_db)

    KAPTIVE_G2G3.out.results
        .branch {
            typed:     it[2].text.contains('Perfect') || it[2].text.contains('Very High') ||
                       it[2].text.contains('High')    || it[2].text.contains('Good') ||
                       it[2].text.contains('Low')     || it[2].text.contains('\tTypeable')
            untypeable: true
        }
        .set { ch_g2g3_branched }

    ch_for_g1g4 = ch_g2g3_branched.untypeable.map { id, fasta, g2g3_tsv -> tuple(id, fasta) }
    KAPTIVE_G1G4(ch_for_g1g4, ch_g1g4_db, ch_normalise_script)

    ch_g2g3_for_merge = KAPTIVE_G2G3.out.results.map { id, fasta, tsv -> tuple(id, tsv) }
    ch_g1g4_for_merge = KAPTIVE_G1G4.out.results
        .mix(
            ch_g2g3_branched.typed.map { id, fasta, g2g3_tsv ->
                tuple(id, file("${workDir}/NO_G1G4_${id}.tsv"))
            }
        )

    PARSE_KAPTIVE(
        ch_g2g3_for_merge.join(ch_g1g4_for_merge)
    )

    ch_ecoli_ktype = PARSE_KAPTIVE.out.ktype.map { id, f -> f }.collect().ifEmpty([])

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2b: Salmonella typing (parallel per sample)
    // ─────────────────────────────────────────────────────────────────────────
    MLST_SALMONELLA(ch_salmonella,          'senterica_achtman_2')
    AMRFINDER_SALMONELLA(ch_salmonella,     'Salmonella')
    SISTR(ch_salmonella)
    PLASMIDFINDER_SALMONELLA(ch_salmonella, 'enterobacteriaceae')

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2d: Shigella typing (parallel per sample)
    // ─────────────────────────────────────────────────────────────────────────
    ch_pinv_db          = file("${projectDir}/assets/pinv_markers.fasta",
                               checkIfExists: true)
    ch_is_db            = file("${projectDir}/assets/shigella_is_elements.fasta",
                               checkIfExists: true)
    ch_parse_mykrobe_script = file("${projectDir}/bin/parse_mykrobe.py",
                               checkIfExists: true)

    MLST_SHIGELLA(ch_shigella,             'ecoli_achtman_4')
    AMRFINDER_SHIGELLA(ch_shigella,        'Escherichia')
    SHIGEIFINDER(ch_shigella)
    PLASMIDFINDER_SHIGELLA(ch_shigella,    'enterobacteriaceae')
    MYKROBE(ch_shigella)
    PARSE_MYKROBE(MYKROBE.out, ch_parse_mykrobe_script)
    PINV_SCREEN(ch_shigella, ch_pinv_db)
    IS_SCREEN(ch_shigella, ch_is_db)

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 3: Core-SNP phylogenetics (SKA2 + IQ-TREE) per species
    //          Skip entirely when --skip_local_phylo is set.
    // ─────────────────────────────────────────────────────────────────────────
    ch_ecoli_tree            = Channel.value(file('NO_FILE'))
    ch_salmonella_tree       = Channel.value(file('NO_FILE'))
    ch_shigella_tree         = Channel.value(file('NO_FILE'))
    ch_ecoli_snp_matrix      = Channel.value(file('NO_FILE'))
    ch_salmonella_snp_matrix = Channel.value(file('NO_FILE'))
    ch_shigella_snp_matrix   = Channel.value(file('NO_FILE'))

    if (!params.skip_local_phylo) {
        ch_ecoli_fastas      = ch_ecoli.map      { id, fasta -> fasta }.collect().filter { it.size() > 0 }
        ch_salmonella_fastas = ch_salmonella.map { id, fasta -> fasta }.collect().filter { it.size() > 0 }
        ch_shigella_fastas   = ch_shigella.map   { id, fasta -> fasta }.collect().filter { it.size() > 0 }

        SKA2_ECOLI(ch_ecoli_fastas)
        SKA2_SALMONELLA(ch_salmonella_fastas)
        SKA2_SHIGELLA(ch_shigella_fastas)

        IQTREE_ECOLI(SKA2_ECOLI.out.alignment)
        IQTREE_SALMONELLA(SKA2_SALMONELLA.out.alignment)
        IQTREE_SHIGELLA(SKA2_SHIGELLA.out.alignment)

        ch_ecoli_tree      = IQTREE_ECOLI.out.treefile.ifEmpty(file('NO_FILE'))
        ch_salmonella_tree = IQTREE_SALMONELLA.out.treefile.ifEmpty(file('NO_FILE'))
        ch_shigella_tree   = IQTREE_SHIGELLA.out.treefile.ifEmpty(file('NO_FILE'))

        ch_ecoli_snp_matrix      = SKA2_ECOLI.out.snp_matrix.ifEmpty(file('NO_FILE'))
        ch_salmonella_snp_matrix = SKA2_SALMONELLA.out.snp_matrix.ifEmpty(file('NO_FILE'))
        ch_shigella_snp_matrix   = SKA2_SHIGELLA.out.snp_matrix.ifEmpty(file('NO_FILE'))
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 4: Aggregate per species
    // ─────────────────────────────────────────────────────────────────────────
    AGGREGATE_ECOLI(
        MLST_ECOLI.out.map         { id, f -> f }.collect().ifEmpty([]),
        AMRFINDER_ECOLI.out.map    { id, f -> f }.collect().ifEmpty([]),
        ECTYPER.out.map            { id, f -> f }.collect().ifEmpty([]),
        PLASMIDFINDER_ECOLI.out.map{ id, f -> f }.collect().ifEmpty([]),
        ch_ecoli_ktype,
        'NO_FILE',
        file("${projectDir}/assets/ecoli_st_complexes.tsv"),
        file("${projectDir}/assets/amrrules/Escherichia_coli.tsv"),
        'ecoli',
        KLEBORATE.out.map          { id, f -> f }.collect().ifEmpty([]),
        EZCLERMONT.out.map         { id, f -> f }.collect().ifEmpty([]),
        [], [], []                 // mykrobe, pinv, is (shigella only)
    )

    AGGREGATE_SALMONELLA(
        MLST_SALMONELLA.out.map         { id, f -> f }.collect().ifEmpty([]),
        AMRFINDER_SALMONELLA.out.map    { id, f -> f }.collect().ifEmpty([]),
        SISTR.out.map                   { id, f -> f }.collect().ifEmpty([]),
        PLASMIDFINDER_SALMONELLA.out.map{ id, f -> f }.collect().ifEmpty([]),
        [],
        'NO_FILE',
        file("${projectDir}/assets/salmonella_st_complexes.tsv"),
        file("${projectDir}/assets/amrrules/Salmonella_enterica.tsv"),
        'salmonella',
        [], [],                    // kleborate, clermont (ecoli only)
        [], [], []                 // mykrobe, pinv, is (shigella only)
    )

    AGGREGATE_SHIGELLA(
        MLST_SHIGELLA.out.map            { id, f -> f }.collect().ifEmpty([]),
        AMRFINDER_SHIGELLA.out.map       { id, f -> f }.collect().ifEmpty([]),
        SHIGEIFINDER.out.map             { id, f -> f }.collect().ifEmpty([]),
        PLASMIDFINDER_SHIGELLA.out.map   { id, f -> f }.collect().ifEmpty([]),
        [],
        'NO_FILE',
        file("${projectDir}/assets/sonnei_st_complexes.tsv"),
        file("${projectDir}/assets/amrrules/Escherichia_coli.tsv"),
        'shigella',
        [], [],                    // kleborate, clermont (ecoli only)
        PARSE_MYKROBE.out.map        { id, f -> f }.collect().ifEmpty([]),
        PINV_SCREEN.out.map          { id, f -> f }.collect().ifEmpty([]),
        IS_SCREEN.out.map            { id, f -> f }.collect().ifEmpty([])
    )

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 5: Summary plots (always runs — no API keys required)
    // ─────────────────────────────────────────────────────────────────────────
    PLOT_SUMMARY_ECOLI(AGGREGATE_ECOLI.out.results,         'ecoli')
    PLOT_SUMMARY_SALMONELLA(AGGREGATE_SALMONELLA.out.results, 'salmonella')
    PLOT_SUMMARY_SHIGELLA(AGGREGATE_SHIGELLA.out.results,   'shigella')

    TREE_ANNOTATION_ECOLI(
        ch_ecoli_tree,
        AGGREGATE_ECOLI.out.results,
        'ecoli'
    )
    TREE_ANNOTATION_SALMONELLA(
        ch_salmonella_tree,
        AGGREGATE_SALMONELLA.out.results,
        'salmonella'
    )
    TREE_ANNOTATION_SHIGELLA(
        ch_shigella_tree,
        AGGREGATE_SHIGELLA.out.results,
        'shigella'
    )

    // SNP distance heatmaps — .filter skips species with no samples (NO_FILE sentinel)
    if (!params.skip_local_phylo) {
        SNP_HEATMAP_ECOLI(
            ch_ecoli_snp_matrix.filter      { it.name != 'NO_FILE' },
            'ecoli'
        )
        SNP_HEATMAP_SALMONELLA(
            ch_salmonella_snp_matrix.filter { it.name != 'NO_FILE' },
            'salmonella'
        )
        SNP_HEATMAP_SHIGELLA(
            ch_shigella_snp_matrix.filter   { it.name != 'NO_FILE' },
            'shigella'
        )
    }

}

// ── Completion summary ────────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ============================================================
    enteric-typer  ${workflow.success ? 'completed successfully' : 'finished with errors'}
    Results  : ${params.outdir}
    Duration : ${workflow.duration}
    ============================================================
    """.stripIndent()
}

// ── Error handler: match known failure patterns and suggest fixes ─────────────
workflow.onError {
    def err = (workflow.errorMessage ?: '') + '\n' + (workflow.errorReport ?: '')

    // Known error patterns → fix messages
    def KNOWN_ERRORS = [
        [
            pattern: /Failed to create Conda environment/,
            title:   "Conda environment creation failed",
            fix:     """\
                A conda environment could not be built. Common causes:

                1. Network issue — check internet connection and retry
                2. Apple Silicon without arm64 profile — rerun with:
                     -profile conda,arm64
                3. Corrupt conda cache — clear it and retry:
                     rm -rf work/conda/
                4. Conda not initialised — run:
                     conda init && source ~/.zshrc  (or ~/.bashrc)
                   then retry"""
        ],
        [
            pattern: /No such file.*enteric_species_refs\.msh|checkIfExists.*enteric_species_refs/,
            title:   "Mash reference sketch missing",
            fix:     """\
                Build the reference sketch before running the pipeline:

                  bash assets/build_references.sh

                This downloads 7 reference genomes from NCBI and produces
                assets/enteric_species_refs.msh  (~1-5 min depending on bandwidth)."""
        ],
        [
            pattern: /amrfinder.*update|AMRFinder.*database|No database/,
            title:   "AMRFinder database missing or outdated",
            fix:     """\
                Update the AMRFinder database:

                  amrfinder --update

                Or inside the conda environment:

                  conda run -n <amrfinder_env> amrfinder --update"""
        ],
        [
            pattern: /No files match.*fasta|No such file.*\.fasta|checkIfExists.*fasta/,
            title:   "No FASTA files found in --input_dir",
            fix:     """\
                No assemblies were found. Check:

                1. The path is correct:        --input_dir /path/to/assemblies/
                2. Files end in a recognised extension:  .fasta .fa .fna .fas
                3. Alternatively use a samplesheet:
                     python bin/make_samplesheet.py --input /path/to/assemblies/ --output samples.csv
                     nextflow run main.nf --samplesheet samples.csv ..."""
        ],
        [
            pattern: /FASTA not found for/,
            title:   "Samplesheet contains missing FASTA paths",
            fix:     """\
                One or more FASTA files listed in your samplesheet do not exist.
                Regenerate the samplesheet from your assembly folder:

                  python bin/make_samplesheet.py --input /path/to/assemblies/ --output samples.csv"""
        ],
        [
            pattern: /OutOfMemoryError|java\.lang\.OutOfMemory|Insufficient memory/,
            title:   "Java / process ran out of memory",
            fix:     """\
                A process exceeded its memory allocation. Options:

                1. Increase process limits in a custom config:
                     process { withLabel: 'high' { memory = '64 GB' } }
                   then run with:  -c custom.config
                2. On HPC, request more RAM via the SLURM/PBS profile.
                3. For IQ-TREE specifically, reduce bootstrap replicates:
                     --iqtree_bootstraps 100"""
        ],
        [
            pattern: /iqtree.*error|IQ-TREE.*failed|IQTREE.*error/,
            title:   "IQ-TREE failed",
            fix:     """\
                IQ-TREE failed to build a tree. Common causes:

                1. Too few variable sites in the SKA2 alignment — this can happen
                   with very closely related samples. Try without phylogenetics:
                     --skip_local_phylo
                2. Insufficient memory — see memory fix above.
                3. Check the IQ-TREE log:
                     cat work/<hash>/iqtree.log"""
        ],
    ]

    def matched = false
    for (e in KNOWN_ERRORS) {
        if (err =~ e.pattern) {
            log.error """
════════════════════════════════════════════════════════════════
KNOWN ERROR: ${e.title}
────────────────────────────────────────────────────────────────
${e.fix.stripIndent()}
════════════════════════════════════════════════════════════════
            """.stripIndent()
            matched = true
            break
        }
    }

    if (!matched) {
        log.error """
════════════════════════════════════════════════════════════════
Pipeline failed — no specific fix matched this error.
Check the full log for details:
  cat .nextflow.log
Or the failed task's working directory (shown in the error above):
  cat work/<ab>/<hash>/.command.err
════════════════════════════════════════════════════════════════
        """.stripIndent()
    }
}
