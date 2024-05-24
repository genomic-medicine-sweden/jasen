// Function for sequencing_run & clarity_sample_id
def get_seqrun_meta(LinkedHashMap row) {
    if (row.sequencing_run && row.clarity_sample_id) {
        meta = tuple(row.id, row.sequencing_run, row.clarity_sample_id)
    } else {
        exit 1, "ERROR: Please check input csv -> sequencing_run and/or clarity_sample_id is missing!"
    }
    return meta
}