// Function for platform and paired-end or single-end
def get_seqplat_meta(LinkedHashMap row) {
    platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]
    if (row.platform in platforms) {
        meta = tuple(row.id, row.platform)
    } else {
        exit 1, "ERROR: Please check input samplesheet -> Platform is not one of the following:!\n-${platforms.join('\n-')}"
    }
    return meta
}

def get_reads(LinkedHashMap row) {
    platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]
    if (row.platform in platforms) {
        if (row.read2) {
            reads = tuple(row.id, tuple(file(row.read1), file(row.read2)))
        } else {
            reads = tuple(row.id, tuple(file(row.read1)))
        }
    } else {
        exit 1, "ERROR: Please check input samplesheet -> Platform is not one of the following:!\n-${platforms.join('\n-')}"
    }
    return reads
}
