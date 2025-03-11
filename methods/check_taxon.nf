def check_taxon(species) {
  "Check if organism name exists"
  if (species == "streptococcus" || species == "streptococcus pyogenes") {
    return "Other"
  }
  return species
}
