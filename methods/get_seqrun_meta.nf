// Function for sequencing_run & clarity_sample_id
def get_seqrun_meta(LinkedHashMap row) {
    if (row.sequencing_run) {
        meta = [row.id, row.sequencing_run]
        if (row.clarity_sample_id) {
            meta << row.clarity_sample_id
        } else {
            row.clarity_sample_id = row.id
            meta << row.clarity_sample_id
        }
        if (row.sample_name) {
            meta << row.sample_name
        } else {
            row.sample_name = row.id
            meta << row.sample_name
        }
    } else {
        exit 1, "ERROR: Please check input csv -> sequencing_run is missing!"
    }
    return tuple(meta)
}