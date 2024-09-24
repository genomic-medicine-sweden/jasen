def getSpeciesTaxonName(fullName) {
  "Convert the full name to the abbreviated version"
  names = fullName.split(' ')
  if (fullName == "escherichia coli") {
    return names[0].capitalize()
  } else if (fullName == "streptococcus") {
    return null
  }
  return names[0].capitalize() + "_" + names[1]
}

def checkTaxon(species) {
  "Check if organism name exists"
  if (species == "streptococcus") {
    return "Other"
  }
  return species
}