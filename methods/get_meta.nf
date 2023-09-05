idList = []     // Global variable to save IDs from samplesheet

// Function for platform and paired-end or single-end
def get_meta(LinkedHashMap row) {
  platforms = ["illumina", "nanopore", "pacbio", "iontorrent"]

  // Error messages
  idErrorMessage = "\nERROR: Please check input samplesheet -> ID '${row.id}' is not a unique name in the samplesheet."
  platformErrorMessage = "\nERROR: Please check input samplesheet -> Platform is not one of the following:\n-${platforms.join('\n-')}"
  charErrorMessage = "\nERROR in ID '${row.id}': Please check input samplesheet -> ID must only consist of ASCII characters."
  lengthErrorMessage = "\nERROR in ID '${row.id}': Please check input samplesheet -> ID must only be between 5-50 characters long."
  errorMessages = []

  // Check if rows fullfill requirements
  idList.add(row.id)
  identicalIds=idList.clone().unique().size() != idList.size()
  correctPlatform = (row.platform in platforms)
  correctCharacters = ((String) row.id) =~ /^[\x00-\x7F]+$/
  correctLength = (row.id.length() >= 5 && row.id.length() <= 50)
  
  // Append error messages if false
  if (identicalIds){
    errorMessages+=idErrorMessage
  }
  if (!correctPlatform){
    errorMessages+=platformErrorMessage
  }
  if (!correctCharacters){
    errorMessages+=charErrorMessage
  }
  if (!correctLength){
    errorMessages+=lengthErrorMessage
  }

  if (correctPlatform && !identicalIds && correctCharacters && correctLength) {
    if (row.read2) {
      meta = tuple(row.id, tuple(file(row.read1), file(row.read2)), row.platform)
    } else {
      meta = tuple(row.id, tuple(file(row.read1)), row.platform)
    }
  } else {
    exit 1, errorMessages.each { println it }
  }
  return meta
}
