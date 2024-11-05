// Function for platform and paired-end or single-end
def get_meta(LinkedHashMap row) {
    platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]
    if (row.platform in platforms) {
        meta = tuple(row.id, tuple(file(row.read1), file(row.read2)), row.platform)
        if (row.read2) {
            reads = tuple(file(row.read1), file(row.read2))
        } else {
            reads = tuple(file(row.read1))
        }
    } else {
        exit 1, "ERROR: Please check input samplesheet -> Platform is not one of the following:!\n-${platforms.join('\n-')}"
    }
    return meta, reads
}
