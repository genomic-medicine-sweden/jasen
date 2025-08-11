def check_taxon(species) {
  "Check if organism name exists"
  if (species in ["staphylococcus", "streptococcus", "streptococcus pyogenes"]) {
    return "Other"
  }
  return species
}
