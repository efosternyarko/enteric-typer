#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * enteric-typer: species-gated genotyping workflow for Enterobacteriaceae
 *
 * Supports:  Escherichia coli  |  Salmonella enterica
 *
 * Phase 1 – species gate (Mash) against reference sketch
 * Phase 2 – species-specific typing in parallel
 * Phase 2c– K-locus typing (E. coli only): G2/G3 → G1/G4 on untypeables
 * Phase 3 – SKA2 core-SNP alignment + IQ-TREE per species
 * Phase 4 – Pathogenwatch upload + multi-threshold cgMLST clustering (optional)
 * Phase 5 – Aggregation + Microreact project creation (optional)
 */

// DSL2 allows each module to be aliased so we can call the same process
// logic with different parameters for E. coli vs Salmonella
include { SPECIES_CHECK                             } from './modules/species_check'
include { MLST          as MLST_ECOLI               } from './modules/mlst'
include { MLST          as MLST_SALMONELLA          } from './modules/mlst'
include { AMRFINDER     as AMRFINDER_ECOLI          } from './modules/amrfinder'
include { AMRFINDER     as AMRFINDER_SALMONELLA     } from './modules/amrfinder'
include { ECTYPER                                   } from './modules/ectyper'
include { KAPTIVE_G2G3                              } from './modules/kaptive'
include { KAPTIVE_G1G4                              } from './modules/kaptive'
include { PARSE_KAPTIVE                             } from './modules/kaptive'
include { SISTR                                     } from './modules/sistr'
include { PLASMIDFINDER as PLASMIDFINDER_ECOLI      } from './modules/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SALMONELLA } from './modules/plasmidfinder'
include { SKA2_BUILD    as SKA2_ECOLI               } from './modules/ska2'
include { SKA2_BUILD    as SKA2_SALMONELLA          } from './modules/ska2'
include { IQTREE        as IQTREE_ECOLI             } from './modules/iqtree'
include { IQTREE        as IQTREE_SALMONELLA        } from './modules/iqtree'
include { PATHOGENWATCH as PATHOGENWATCH_ECOLI      } from './modules/pathogenwatch'
include { PATHOGENWATCH as PATHOGENWATCH_SALMONELLA } from './modules/pathogenwatch'
include { AGGREGATE     as AGGREGATE_ECOLI          } from './modules/aggregate'
include { AGGREGATE     as AGGREGATE_SALMONELLA     } from './modules/aggregate'
include { MICROREACT_UPLOAD as MICROREACT_ECOLI     } from './modules/microreact_upload'
include { MICROREACT_UPLOAD as MICROREACT_SALMONELLA} from './modules/microreact_upload'
include { PLOT_SUMMARY      as PLOT_SUMMARY_ECOLI      } from './modules/plot_summary'
include { PLOT_SUMMARY      as PLOT_SUMMARY_SALMONELLA } from './modules/plot_summary'
include { TREE_ANNOTATION   as TREE_ANNOTATION_ECOLI      } from './modules/tree_annotation'
include { TREE_ANNOTATION   as TREE_ANNOTATION_SALMONELLA } from './modules/tree_annotation'

// ── Parameter validation ──────────────────────────────────────────────────────
if (!params.samplesheet && !params.input_dir) {
    error "ERROR: Provide --samplesheet <csv> or --input_dir <folder>.\nSee README.md for details."
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

    // Join species label back and branch into species-specific channels
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
        other:      true
    }.set { ch_branched }

    ch_branched.other.subscribe { id, fasta, species, dist ->
        log.warn "SKIPPED '${id}': best match = ${species} (Mash dist = ${String.format('%.4f', dist)}) — not E. coli or Salmonella within threshold"
    }

    ch_ecoli      = ch_branched.ecoli.map      { id, fasta, sp, d -> tuple(id, fasta) }
    ch_salmonella = ch_branched.salmonella.map { id, fasta, sp, d -> tuple(id, fasta) }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2a: E. coli typing (parallel per sample)
    // ─────────────────────────────────────────────────────────────────────────
    MLST_ECOLI(ch_ecoli,          'ecoli_achtman_4')
    AMRFINDER_ECOLI(ch_ecoli,     'Escherichia')
    ECTYPER(ch_ecoli)
    PLASMIDFINDER_ECOLI(ch_ecoli, 'enterobacteriaceae')

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2c: E. coli K-locus typing
    //   Step 1: Run G2/G3 database on ALL E. coli samples
    //   Step 2: Run G1/G4 (scores + normalise) on G2/G3-untypeable samples only
    //   Step 3: Merge results — G1/G4 and G2/G3 are mutually exclusive
    // ─────────────────────────────────────────────────────────────────────────
    ch_g2g3_db        = file("${projectDir}/assets/kaptive_dbs/EC-K-typing_group2and3_v3.0.0.gbk",
                              checkIfExists: true)
    ch_g1g4_db        = file("${projectDir}/assets/kaptive_dbs/EC-K-typing_group1and4_v0.9.gbk",
                              checkIfExists: true)
    ch_normalise_script = file("${projectDir}/bin/normalise_kaptive_scores.py",
                              checkIfExists: true)

    KAPTIVE_G2G3(ch_ecoli, ch_g2g3_db)

    // Route samples: G2/G3 typeable → skip G1/G4; untypeable → run G1/G4
    KAPTIVE_G2G3.out.results
        .branch {
            typed:     it[2].text.contains('Perfect') || it[2].text.contains('Very High') ||
                       it[2].text.contains('High')    || it[2].text.contains('Good') ||
                       it[2].text.contains('Low')
            untypeable: true
        }
        .set { ch_g2g3_branched }

    ch_for_g1g4 = ch_g2g3_branched.untypeable.map { id, fasta, g2g3_tsv -> tuple(id, fasta) }
    KAPTIVE_G1G4(ch_for_g1g4, ch_g1g4_db, ch_normalise_script)

    // Pair G2/G3 + G1/G4 results (G1/G4 result is absent for G2/G3-typed samples)
    ch_g2g3_for_merge = KAPTIVE_G2G3.out.results.map { id, fasta, tsv -> tuple(id, tsv) }
    ch_g1g4_for_merge = KAPTIVE_G1G4.out.results
        .mix(
            // For G2/G3-typed samples, emit a no-op placeholder for G1/G4
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
    MLST_SALMONELLA(ch_salmonella,          'salmonella')
    AMRFINDER_SALMONELLA(ch_salmonella,     'Salmonella')
    SISTR(ch_salmonella)
    PLASMIDFINDER_SALMONELLA(ch_salmonella, 'enterobacteriaceae')

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 3: Core-SNP phylogenetics (SKA2 + IQ-TREE) per species
    // ─────────────────────────────────────────────────────────────────────────
    ch_ecoli_fastas      = ch_ecoli.map      { id, fasta -> fasta }.collect().ifEmpty(['NO_FASTAS'])
    ch_salmonella_fastas = ch_salmonella.map { id, fasta -> fasta }.collect().ifEmpty(['NO_FASTAS'])

    SKA2_ECOLI(ch_ecoli_fastas)
    SKA2_SALMONELLA(ch_salmonella_fastas)

    IQTREE_ECOLI(SKA2_ECOLI.out.alignment)
    IQTREE_SALMONELLA(SKA2_SALMONELLA.out.alignment)

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 4: Pathogenwatch (optional — requires PW_API_KEY)
    // ─────────────────────────────────────────────────────────────────────────
    ch_pw_ecoli_out      = Channel.value('NO_FILE')
    ch_pw_salmonella_out = Channel.value('NO_FILE')

    if (params.run_pathogenwatch) {
        // Collect (id, fasta) tuples for batch upload
        ch_ecoli_batch      = ch_ecoli.collect()
        ch_salmonella_batch = ch_salmonella.collect()

        def ecoli_collection_name = params.pathogenwatch_ecoli_collection_name
            ?: "enteric-typer-ecoli-${workflow.runName}"
        def salm_collection_name  = params.pathogenwatch_salmonella_collection_name
            ?: "enteric-typer-salmonella-${workflow.runName}"

        PATHOGENWATCH_ECOLI(
            ch_ecoli_batch,
            'ecoli',
            ecoli_collection_name
        )
        PATHOGENWATCH_SALMONELLA(
            ch_salmonella_batch,
            'salmonella',
            salm_collection_name
        )

        ch_pw_ecoli_out      = PATHOGENWATCH_ECOLI.out.results.map      { it.toString() }
        ch_pw_salmonella_out = PATHOGENWATCH_SALMONELLA.out.results.map { it.toString() }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 5: Aggregate per species
    // ─────────────────────────────────────────────────────────────────────────
    AGGREGATE_ECOLI(
        MLST_ECOLI.out.map         { id, f -> f }.collect().ifEmpty([]),
        AMRFINDER_ECOLI.out.map    { id, f -> f }.collect().ifEmpty([]),
        ECTYPER.out.map            { id, f -> f }.collect().ifEmpty([]),
        PLASMIDFINDER_ECOLI.out.map{ id, f -> f }.collect().ifEmpty([]),
        ch_ecoli_ktype,
        ch_pw_ecoli_out,
        file("${projectDir}/assets/ecoli_st_complexes.tsv"),
        'ecoli'
    )

    AGGREGATE_SALMONELLA(
        MLST_SALMONELLA.out.map         { id, f -> f }.collect().ifEmpty([]),
        AMRFINDER_SALMONELLA.out.map    { id, f -> f }.collect().ifEmpty([]),
        SISTR.out.map                   { id, f -> f }.collect().ifEmpty([]),
        PLASMIDFINDER_SALMONELLA.out.map{ id, f -> f }.collect().ifEmpty([]),
        [],   // no K-typing for Salmonella
        ch_pw_salmonella_out,
        file("${projectDir}/assets/salmonella_st_complexes.tsv"),
        'salmonella'
    )

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 6: Microreact upload (optional — requires MICROREACT_TOKEN)
    // ─────────────────────────────────────────────────────────────────────────
    if (params.upload_microreact) {
        ch_ecoli_tree_for_mr      = IQTREE_ECOLI.out.treefile.ifEmpty(file('NO_FILE'))
        ch_salmonella_tree_for_mr = IQTREE_SALMONELLA.out.treefile.ifEmpty(file('NO_FILE'))

        MICROREACT_ECOLI(
            AGGREGATE_ECOLI.out.results,
            ch_ecoli_tree_for_mr,
            "${params.microreact_project} – E. coli"
        )
        MICROREACT_SALMONELLA(
            AGGREGATE_SALMONELLA.out.results,
            ch_salmonella_tree_for_mr,
            "${params.microreact_project} – Salmonella"
        )
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 7: Summary plots (default — no API keys required)
    // ─────────────────────────────────────────────────────────────────────────
    PLOT_SUMMARY_ECOLI(AGGREGATE_ECOLI.out.results,      'ecoli')
    PLOT_SUMMARY_SALMONELLA(AGGREGATE_SALMONELLA.out.results, 'salmonella')

    ch_ecoli_tree_for_annot      = IQTREE_ECOLI.out.treefile.ifEmpty(file('NO_FILE'))
    ch_salmonella_tree_for_annot = IQTREE_SALMONELLA.out.treefile.ifEmpty(file('NO_FILE'))

    TREE_ANNOTATION_ECOLI(
        ch_ecoli_tree_for_annot,
        AGGREGATE_ECOLI.out.results,
        'ecoli'
    )
    TREE_ANNOTATION_SALMONELLA(
        ch_salmonella_tree_for_annot,
        AGGREGATE_SALMONELLA.out.results,
        'salmonella'
    )
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
