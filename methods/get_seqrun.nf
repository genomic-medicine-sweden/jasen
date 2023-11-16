// Function for platform and paired-end or single-end
def get_seqrun(LinkedHashMap row) {
    if (row.sequencing_run) {
        meta = tuple(row.id, row.sequencing_run)
    } else {
        exit 1, "ERROR: Please check input csv -> sequencing_run is missing!"
    }
    return meta
}