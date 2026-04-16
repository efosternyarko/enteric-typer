// ── AMRFINDER_UPDATE: download AMRFinder Plus database once ───────────────────
// AMRFinder Plus conda package does not bundle its database. This process
// runs once before any per-sample AMRFINDER calls to ensure the database
// is present. All AMRFINDER processes receive the output path as a value
// channel so they depend on (and reuse) the same database directory.

process AMRFINDER_UPDATE {
    label 'low'

    conda     "${projectDir}/envs/amrfinder.yml"
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0'

    // Cache indefinitely — only reruns if the conda env changes
    storeDir "${projectDir}/assets/amrfinderplus_db"

    output:
    path('db_ready.txt'), emit: ready

    script:
    """
    amrfinder -u
    amrfinder --database_version > db_ready.txt
    """
}
