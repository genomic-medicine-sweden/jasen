def get_species_taxon_name(fullName) {
  "Convert the full name to the abbreviated version"
  names = fullName.split(' ')
  if (fullName == "escherichia coli") {
    return names[0].capitalize()
  } else if (fullName in ["streptococcus", "staphylococcus"]) {
    return null
  }
  return names[0].capitalize() + "_" + names[1]
}
