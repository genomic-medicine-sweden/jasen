def check_taxon(species) {
  "Check if organism name exists"
  if (species in ["streptococcus", "streptococcus pyogenes"]) {
    return "Other"
  }
  return species
}
