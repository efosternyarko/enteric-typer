// ── KAPTIVE: E. coli capsular K-locus typing ─────────────────────────────────
//
// Implements the two-step workflow described in efosternyarko/EC-K-typing-G1G4:
//
//   Step 1 — KAPTIVE_G2G3:
//     Run Kaptive v3 on all E. coli assemblies using the Group 2/3 database
//     (Gladstone et al. 2024). Samples typed here are Group 2 or 3.
//     Untypeable samples (Match confidence == "Untypeable") are forwarded to
//     Step 2.
//
//   Step 2 — KAPTIVE_G1G4:
//     Run Kaptive v3 in --scores mode on G2/G3-untypeable assemblies using
//     the Group 1/4 database (651 loci). Raw scores are then normalised by
//     expected locus CDS length (normalise_kaptive_scores.py) to achieve
//     100% typeability on G1/G4 strains.
//
// NOTE: G1/G4 and G2/G3 loci are MUTUALLY EXCLUSIVE in E. coli — running both
//       databases simultaneously causes wzy-interference misclassification.
//       The sequential approach (G2/G3 first, G1/G4 on untypeables) is correct.

// ── Step 1: G2/G3 typing ─────────────────────────────────────────────────────

process KAPTIVE_G2G3 {
    tag "${sample_id}"
    label 'medium'

    conda     "${projectDir}/envs/kaptive.yml"
    container 'quay.io/biocontainers/kaptive:3.0.0b5--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(fasta)
    path(g2g3_db)

    output:
    tuple val(sample_id), path(fasta), path("${sample_id}_kaptive_g2g3.tsv"), emit: results

    script:
    """
    kaptive assembly \\
        ${g2g3_db} \\
        ${fasta} \\
        -o ${sample_id}_kaptive_g2g3.tsv \\
        -t ${task.cpus} \\
        2>${sample_id}_kaptive_g2g3.log \\
    || true

    # Fallback: write empty result
    if [ ! -s "${sample_id}_kaptive_g2g3.tsv" ]; then
        printf 'Assembly\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\n' \\
            > "${sample_id}_kaptive_g2g3.tsv"
        printf '%s\\tnone\\tnone\\tUntypeable\\t\\n' "${sample_id}" \\
            >> "${sample_id}_kaptive_g2g3.tsv"
    fi
    """
}

// ── Step 2: G1/G4 typing (scores mode + normalisation) ───────────────────────

process KAPTIVE_G1G4 {
    tag "${sample_id}"
    label 'medium'

    conda     "${projectDir}/envs/kaptive.yml"
    container 'quay.io/biocontainers/kaptive:3.0.0b5--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(fasta)   // only G2/G3-untypeable samples
    path(g1g4_db)
    path(normalise_script)

    output:
    tuple val(sample_id), path("${sample_id}_kaptive_g1g4.tsv"), emit: results

    script:
    """
    # Run Kaptive in --scores mode (gene-by-gene alignment — achieves 100% typeability)
    kaptive assembly \\
        ${g1g4_db} \\
        ${fasta} \\
        -s ${sample_id}_kaptive_g1g4_raw_scores.tsv \\
        -t ${task.cpus} \\
        2>${sample_id}_kaptive_g1g4.log \\
    || true

    # Normalise scores by expected locus CDS length
    if [ -s "${sample_id}_kaptive_g1g4_raw_scores.tsv" ]; then
        python3 ${normalise_script} \\
            --db  ${g1g4_db} \\
            --in  ${sample_id}_kaptive_g1g4_raw_scores.tsv \\
            --out ${sample_id}_kaptive_g1g4.tsv \\
            --min-coverage 0.50 \\
            2>>${sample_id}_kaptive_g1g4.log
    else
        # Kaptive produced no scores — write minimal untypeable result
        printf 'Assembly\\tBest match locus\\tBest match confidence\\tGenes found\\tGenes expected\\tGene coverage\\tRaw AS\\tNorm AS\\n' \\
            > "${sample_id}_kaptive_g1g4.tsv"
        printf '%s\\tnone\\tUntypeable\\t0\\t0\\t0.0%%\\t0\\t0.0\\n' "${sample_id}" \\
            >> "${sample_id}_kaptive_g1g4.tsv"
    fi
    """
}

// ── PARSE_KAPTIVE: merge G2/G3 + G1/G4 into one K-type per sample ───────────

process PARSE_KAPTIVE {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11-slim'

    input:
    tuple val(sample_id), path(g2g3_tsv), path(g1g4_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_ktype.tsv"), emit: ktype

    script:
    """
    parse_kaptive.py \\
        --sample    ${sample_id} \\
        --g2g3      ${g2g3_tsv} \\
        --g1g4      ${g1g4_tsv} \\
        --output    ${sample_id}_ktype.tsv
    """
}
