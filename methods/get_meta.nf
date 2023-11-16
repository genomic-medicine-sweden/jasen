// Function for platform and paired-end or single-end
def get_meta(LinkedHashMap row) {
    platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]
    if (row.platform in platforms) {
        if (row.read2) {
            meta = tuple(row.id, tuple(file(row.read1), file(row.read2)), row.platform)
        } else {
            meta = tuple(row.id, tuple(file(row.read1)), row.platform)
        }
    } else {
        exit 1, "ERROR: Please check input samplesheet -> Platform is not one of the following:!\n-${platforms.join('\n-')}"
    }
    return meta
}
