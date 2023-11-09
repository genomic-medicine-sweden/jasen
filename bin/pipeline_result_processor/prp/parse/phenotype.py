"""Parse output of various phenotype predicting tools."""

import re
import json
import logging
import pandas as pd
from typing import Any, Dict, Tuple

from ..models.metadata import SoupVersion, SoupVersions
from ..models.phenotype import (
    ElementTypeResult,
    ElementType,
    ResistanceGene,
    ResistanceVariant,
    VirulenceGene,
)
from ..models.phenotype import PredictionSoftware as Software
from ..models.sample import MethodIndex

LOG = logging.getLogger(__name__)


def _get_resfinder_amr_sr_profie(resfinder_result, limit_to_phenotypes=None):
    """Get resfinder susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()
    for phenotype in resfinder_result["phenotypes"].values():
        # skip phenotype if its not part of the desired category
        if (
            limit_to_phenotypes is not None
            and phenotype["key"] not in limit_to_phenotypes
        ):
            continue

        if "amr_resistant" in phenotype.keys():
            if phenotype["amr_resistant"]:
                resistant.add(phenotype["amr_resistance"])
            else:
                susceptible.add(phenotype["amr_resistance"])
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_resfinder_amr_genes(
    resfinder_result, limit_to_phenotypes=None
) -> Tuple[ResistanceGene, ...]:
    """Get resistance genes from resfinder result."""
    results = []

    if not "seq_regions" in resfinder_result:
        results = _default_resistance().genes
        return results

    for info in resfinder_result["seq_regions"].values():
        # Get only acquired resistance genes
        if not info["ref_database"][0].startswith("Res"):
            continue

        # Get only genes of desired phenotype
        if limit_to_phenotypes is not None:
            intersect = set(info["phenotypes"]) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue

        # store results
        gene = ResistanceGene(
            gene_symbol=info["name"],
            accession=info["ref_acc"],
            depth=info["depth"],
            identity=info["identity"],
            coverage=info["coverage"],
            ref_start_pos=info["ref_start_pos"],
            ref_end_pos=info["ref_end_pos"],
            ref_gene_length=info["ref_seq_lenght"],
            alignment_length=info["alignment_length"],
            phenotypes=info["phenotypes"],
            ref_database=info["ref_database"][0],
            ref_id=info["ref_id"],
            contig_id=None,
            sequence_name=None,
            ass_start_pos=None,
            ass_end_pos=None,
            strand=None,
            element_type=None,
            element_subtype=None,
            target_length=None,
            res_class=None,
            res_subclass=None,
            method=None,
            close_seq_name=None,
        )
        results.append(gene)
    return results


def _parse_resfinder_amr_variants(
    resfinder_result, limit_to_phenotypes=None
) -> Tuple[ResistanceVariant, ...]:
    """Get resistance genes from resfinder result."""
    results = []
    igenes = []
    for info in resfinder_result["seq_variations"].values():
        # Get only variants from desired phenotypes
        if limit_to_phenotypes is not None:
            intersect = set(info["phenotypes"]) & set(limit_to_phenotypes)
            if len(intersect) == 0:
                continue
        # get gene depth
        if "seq_regions" in resfinder_result:
            info["depth"] = resfinder_result["seq_regions"][info["seq_regions"][0]][
                "depth"
            ]
        else:
            info["depth"] = 0
        # translate variation type bools into classifier
        if info["substitution"]:
            var_type = "substitution"
        elif info["insertion"]:
            var_type = "insertion"
        elif info["deletion"]:
            var_type = "deletion"
        else:
            raise ValueError("Output has no known mutation type")
        if not "seq_regions" in info:
            # igenes = _default_resistance().genes
            igenes = [""]
        variant = ResistanceVariant(
            variant_type=var_type,
            genes=igenes,
            phenotypes=info["phenotypes"],
            position=info["ref_start_pos"],
            ref_nt=info["ref_codon"],
            alt_nt=info["var_codon"],
            depth=info["depth"],
            ref_database=info["ref_database"],
            ref_id=info["ref_id"],
        )
        results.append(variant)
    return results


def _parse_amrfinder_amr_results(predictions: dict) -> Tuple[ResistanceGene, ...]:
    """Parse amrfinder prediction results from amrfinderplus."""
    genes = []
    for prediction in predictions:
        gene = ResistanceGene(
            virulence_category=None,
            accession=prediction["close_seq_accn"],
            depth=None,
            identity=prediction["ref_seq_identity"],
            coverage=prediction["ref_seq_cov"],
            ref_start_pos=None,
            ref_end_pos=None,
            ref_gene_length=prediction["ref_seq_len"],
            alignment_length=prediction["align_len"],
            ref_database=None,
            phenotypes=[],
            ref_id=None,
            contig_id=prediction["contig_id"],
            gene_symbol=prediction["gene_symbol"],
            sequence_name=prediction["sequence_name"],
            ass_start_pos=prediction["Start"],
            ass_end_pos=prediction["Stop"],
            strand=prediction["Strand"],
            element_type=prediction["element_type"],
            element_subtype=prediction["element_subtype"],
            target_length=prediction["target_length"],
            res_class=prediction["Class"],
            res_subclass=prediction["Subclass"],
            method=prediction["Method"],
            close_seq_name=prediction["close_seq_name"],
        )
        genes.append(gene)
    return ElementTypeResult(phenotypes=[], genes=genes, mutations=[])


def _get_mykrobe_amr_sr_profie(mykrobe_result):
    """Get mykrobe susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()

    if not mykrobe_result:
        return {}

    for element_type in mykrobe_result:
        if mykrobe_result[element_type]["predict"].upper() == "R":
            resistant.add(element_type)
        else:
            susceptible.add(element_type)
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_mykrobe_amr_genes(mykrobe_result) -> Tuple[ResistanceGene, ...]:
    """Get resistance genes from mykrobe result."""
    results = []

    if not mykrobe_result:
        results = _default_resistance().genes
        return results
    
    for element_type in mykrobe_result:
        if mykrobe_result[element_type]["predict"].upper() == "R":
            hits = mykrobe_result[element_type]["called_by"]
            for hit in hits:
                gene = ResistanceGene(
                    gene_symbol=hit.split("_")[0],
                    accession=None,
                    depth=hits[hit]["info"]["coverage"]["alternate"]["median_depth"],
                    identity=None,
                    coverage=hits[hit]["info"]["coverage"]["alternate"]["percent_coverage"],
                    ref_start_pos=None,
                    ref_end_pos=None,
                    ref_gene_length=None,
                    alignment_length=None,
                    phenotypes=element_type,
                    ref_database=None,
                    ref_id=None,
                    contig_id=None,
                    sequence_name=None,
                    ass_start_pos=None,
                    ass_end_pos=None,
                    strand=None,
                    element_type=None,
                    element_subtype=None,
                    target_length=None,
                    res_class=None,
                    res_subclass=None,
                    method=None,
                    close_seq_name=None,
                )
                results.append(gene)
    return results


def _parse_mykrobe_amr_variants(mykrobe_result) -> Tuple[ResistanceVariant, ...]:
    """Get resistance genes from mykrobe result."""
    results = []
    def get_mutation_type(var_nom):
        try:
            ref_idx = re.search(r"\d", var_nom, 1).start()
            alt_idx = re.search(r"\d(?=[^\d]*$)", var_nom).start()+1
        except AttributeError:
            return [None]*4
        ref = var_nom[:ref_idx]
        alt=var_nom[alt_idx:]
        position = int(var_nom[ref_idx:alt_idx])
        if len(ref) > len(alt):
            mut_type = "deletion"
        elif len(ref) < len(alt):
            mut_type = "insertion"
        else:
            mut_type = "substitution"
        return mut_type, ref, alt, position

    for element_type in mykrobe_result:
        if mykrobe_result[element_type]["predict"].upper() == "R":
            hits = mykrobe_result[element_type]["called_by"]
            for hit in hits:
                if hits[hit]["variant"] == None:
                    var_info = hit.split("-")[1]
                    _, ref_nt, alt_nt, position = get_mutation_type(var_info)
                    var_nom = hit.split("-")[0].split("_")[1]
                    var_type, _, _, _ = get_mutation_type(var_nom)
                    variant = ResistanceVariant(
                        variant_type=var_type,
                        genes=[hit.split("_")[0]],
                        phenotypes=[element_type],
                        position=position,
                        ref_nt=ref_nt,
                        alt_nt=alt_nt,
                        depth=hits[hit]["info"]["coverage"]["alternate"]["median_depth"],
                        ref_database=None,
                        ref_id=None,
                        type=None,
                        change=var_nom,
                        nucleotide_change=None,
                        protein_change=None,
                        annotation=None,
                        drugs=None,
                    )
                    results.append(variant)
    if not results:
        results = _default_variant().mutations
        return results

    return results


def _get_tbprofiler_amr_sr_profie(tbprofiler_result):
    """Get tbprofiler susceptibility/resistance profile."""
    susceptible = set()
    resistant = set()
    drugs = ["ofloxacin", "moxifloxacin", "isoniazid", "delamanid", 
             "kanamycin", "amikacin", "ethambutol", "ethionamide", 
             "streptomycin", "ciprofloxacin", "levofloxacin", "pyrazinamide", 
             "linezolid", "rifampicin", "capreomycin"]

    if not tbprofiler_result:
        return {}

    for hit in tbprofiler_result["dr_variants"]:
        for drug in hit["gene_associated_drugs"]:
            resistant.add(drug)
    susceptible = [drug for drug in drugs if drug not in resistant]
    return {"susceptible": list(susceptible), "resistant": list(resistant)}


def _parse_tbprofiler_amr_variants(tbprofiler_result) -> Tuple[ResistanceVariant, ...]:
    """Get resistance genes from tbprofiler result."""
    results = []

    for hit in tbprofiler_result["dr_variants"]:
        var_type = "substitution"
        variant = ResistanceVariant(
            variant_type=var_type,
            genes=[hit["gene"]],
            phenotypes=hit["gene_associated_drugs"],
            position=int(hit["genome_pos"]),
            ref_nt=hit["ref"],
            alt_nt=hit["alt"],
            depth=hit["depth"],
            ref_database=tbprofiler_result["db_version"]["name"],
            ref_id=None,
            type=hit["type"],
            change=hit["change"],
            nucleotide_change=hit["nucleotide_change"],
            protein_change=hit["protein_change"],
            annotation=hit["annotation"],
            drugs=hit["drugs"],
        )
        results.append(variant)

    if not results:
        results = _default_variant().mutations
        return results
    
    return results


def parse_resfinder_amr_pred(
    prediction: Dict[str, Any], resistance_category
) -> Tuple[SoupVersions, ElementTypeResult]:
    """Parse resfinder resistance prediction results."""
    # resfinder missclassifies resistance the param amr_category by setting all to amr
    LOG.info("Parsing resistance prediction")
    meta = [
        SoupVersion(
            **{
                "name": metrics["database_name"],
                "version": metrics["database_version"],
                "type": metrics["type"],
            }
        )
        for metrics in prediction["databases"].values()
    ]
    # parse resistance based on the category
    categories = {
        ElementType.BIOCIDE: [
            "formaldehyde",
            "benzylkonium chloride",
            "ethidium bromide",
            "chlorhexidine",
            "cetylpyridinium chloride",
            "hydrogen peroxide",
        ],
        ElementType.HEAT: ["temperature"],
    }
    categories[ElementType.AMR] = list(
        {k for k in prediction["phenotypes"].keys()}
        - set(categories[ElementType.BIOCIDE] + categories[ElementType.HEAT])
    )

    # parse resistance
    resistance = ElementTypeResult(
        phenotypes=_get_resfinder_amr_sr_profie(
            prediction, categories[resistance_category]
        ),
        genes=_parse_resfinder_amr_genes(prediction, categories[resistance_category]),
        mutations=_parse_resfinder_amr_variants(
            prediction, categories[resistance_category]
        ),
    )
    return MethodIndex(
        type=resistance_category, software=Software.RESFINDER, result=resistance
    )


def parse_amrfinder_amr_pred(file, element_type: str) -> ElementTypeResult:
    """Parse amrfinder resistance prediction results."""
    LOG.info("Parsing amrfinder amr prediction")
    with open(file, "rb") as tsvfile:
        hits = pd.read_csv(tsvfile, delimiter="\t")
        hits = hits.rename(
            columns={
                "Contig id": "contig_id",
                "Gene symbol": "gene_symbol",
                "Sequence name": "sequence_name",
                "Element type": "element_type",
                "Element subtype": "element_subtype",
                "Target length": "target_length",
                "Reference sequence length": "ref_seq_len",
                "% Coverage of reference sequence": "ref_seq_cov",
                "% Identity to reference sequence": "ref_seq_identity",
                "Alignment length": "align_len",
                "Accession of closest sequence": "close_seq_accn",
                "Name of closest sequence": "close_seq_name",
            }
        )
        hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
        hits = hits.where(pd.notnull(hits), None)
        if element_type == ElementType.AMR:
            predictions = hits[hits["element_type"] == "AMR"].to_dict(orient="records")
            results: ElementTypeResult = _parse_amrfinder_amr_results(predictions)
        elif element_type == ElementType.HEAT:
            predictions = hits[(hits["element_subtype"] == "HEAT")].to_dict(orient="records")
            results: ElementTypeResult = _parse_amrfinder_amr_results(predictions)
        elif element_type == ElementType.BIOCIDE:
            predictions = hits[
                (hits["element_subtype"] == "ACID")
                & (hits["element_subtype"] == "BIOCIDE")
            ].to_dict(orient="records")
            results: ElementTypeResult = _parse_amrfinder_amr_results(predictions)
        elif element_type == ElementType.METAL:
            predictions = hits[hits["element_subtype"] == "METAL"].to_dict(
                orient="records"
            )
            results: ElementTypeResult = _parse_amrfinder_amr_results(predictions)
        else:
            results = _default_resistance()
    return MethodIndex(type=element_type, result=results, software=Software.AMRFINDER)


def _parse_virulencefinder_vir_results(pred: str) -> ElementTypeResult:
    """Parse virulence prediction results from ARIBA."""
    results = {}
    # parse virulence finder results
    species = [k for k in pred["virulencefinder"]["results"].keys()]
    for key, genes in pred["virulencefinder"]["results"][species[0]].items():
        virulence_category = key.split("_")[1]
        vir_genes = []
        if not genes == "No hit found":
            for gn in genes.values():
                start_pos, end_pos = map(int, gn["position_in_ref"].split(".."))
                gene = VirulenceGene(
                    name=gn["virulence_gene"],
                    virulence_category=virulence_category,
                    accession=gn["accession"],
                    depth=None,
                    identity=gn["identity"],
                    coverage=gn["coverage"],
                    ref_start_pos=start_pos,
                    ref_end_pos=end_pos,
                    ref_gene_length=gn["template_length"],
                    alignment_length=gn["HSP_length"],
                    ref_database="virulenceFinder",
                    ref_id=gn["hit_id"],
                    contig_id=None,
                    gene_symbol=None,
                    sequence_name=None,
                    ass_start_pos=int(None),
                    ass_end_pos=int(None),
                    strand=None,
                    element_type=None,
                    element_subtype=None,
                    target_length=int(None),
                    res_class=None,
                    res_subclass=None,
                    method=None,
                    close_seq_name=None,
                )
            vir_genes.append(gene)
        results[virulence_category] = vir_genes

    return ElementTypeResult(results)


def _parse_amrfinder_vir_results(predictions: dict) -> ElementTypeResult:
    """Parse amrfinder prediction results from amrfinderplus."""
    genes = []
    for prediction in predictions:
        gene = VirulenceGene(
            name=None,
            virulence_category=None,
            accession=prediction["close_seq_accn"],
            depth=None,
            identity=prediction["ref_seq_identity"],
            coverage=prediction["ref_seq_cov"],
            ref_start_pos=None,
            ref_end_pos=None,
            ref_gene_length=prediction["ref_seq_len"],
            alignment_length=prediction["align_len"],
            ref_database=None,
            phenotypes=[],
            ref_id=None,
            contig_id=prediction["contig_id"],
            gene_symbol=prediction["gene_symbol"],
            sequence_name=prediction["sequence_name"],
            ass_start_pos=int(prediction["Start"]),
            ass_end_pos=int(prediction["Stop"]),
            strand=prediction["Strand"],
            element_type=prediction["element_type"],
            element_subtype=prediction["element_subtype"],
            target_length=int(prediction["target_length"]),
            res_class=prediction["Class"],
            res_subclass=prediction["Subclass"],
            method=prediction["Method"],
            close_seq_name=prediction["close_seq_name"],
        )
        genes.append(gene)
    return ElementTypeResult(phenotypes=[], genes=genes, mutations=[])


def _default_virulence() -> ElementTypeResult:
    gene = VirulenceGene(
        name=None,
        virulence_category=None,
        accession="",
        depth=None,
        identity=0,
        coverage=0,
        ref_start_pos=0,
        ref_end_pos=0,
        ref_gene_length=0,
        alignment_length=0,
        ref_database="",
        ref_id=0,
        contig_id=None,
        gene_symbol=None,
        sequence_name=None,
        ass_start_pos=None,
        ass_end_pos=None,
        strand=None,
        element_type=None,
        element_subtype=None,
        target_length=None,
        res_class=None,
        res_subclass=None,
        method=None,
        close_seq_name=None,
    )
    genes = list()
    genes.append(gene)
    return ElementTypeResult(phenotypes=[], genes=genes, mutations=[])


def _default_resistance() -> ElementTypeResult:
    gene = ResistanceGene(
        name=None,
        virulence_category=None,
        accession=None,
        depth=None,
        identity=None,
        coverage=None,
        ref_start_pos=None,
        ref_end_pos=None,
        ref_gene_length=None,
        alignment_length=None,
        ref_database=None,
        phenotypes=[],
        ref_id=None,
        contig_id=None,
        sequence_name=None,
        ass_start_pos=None,
        ass_end_pos=None,
        strand=None,
        element_type=None,
        element_subtype=None,
        target_length=None,
        res_class=None,
        res_subclass=None,
        method=None,
        close_seq_name=None,
    )
    genes = list()
    genes.append(gene)
    return ElementTypeResult(phenotypes=[], genes=genes, mutations=[])

def _default_variant() -> ElementTypeResult:
    mutation = ResistanceGene(
        variant_type=None,
        genes=None,
        phenotypes=[],
        position=None,
        ref_nt=None,
        alt_nt=None,
        depth=None,
    )
    mutations = list()
    mutations.append(mutation)
    return ElementTypeResult(phenotypes=[], genes=[], mutations=mutations)


def parse_virulencefinder_vir_pred(file: str) -> ElementTypeResult:
    """Parse virulencefinder virulence prediction results."""
    LOG.info("Parsing virulencefinder virulence prediction")
    pred = json.load(file)
    if "not virulencefinder" in pred:
        results: ElementTypeResult = _parse_virulencefinder_vir_results(pred)
    else:
        results: ElementTypeResult = _default_virulence()
    return MethodIndex(
        type=ElementType.VIR, software=Software.VIRFINDER, result=results
    )


def parse_amrfinder_vir_pred(file: str):
    """Parse amrfinder virulence prediction results."""
    LOG.info("Parsing amrfinder virulence prediction")
    with open(file, "rb") as tsvfile:
        hits = pd.read_csv(tsvfile, delimiter="\t")
        hits = hits.rename(
            columns={
                "Contig id": "contig_id",
                "Gene symbol": "gene_symbol",
                "Sequence name": "sequence_name",
                "Element type": "element_type",
                "Element subtype": "element_subtype",
                "Target length": "target_length",
                "Reference sequence length": "ref_seq_len",
                "% Coverage of reference sequence": "ref_seq_cov",
                "% Identity to reference sequence": "ref_seq_identity",
                "Alignment length": "align_len",
                "Accession of closest sequence": "close_seq_accn",
                "Name of closest sequence": "close_seq_name",
            }
        )
        hits = hits.drop(columns=["Protein identifier", "HMM id", "HMM description"])
        hits = hits.where(pd.notnull(hits), None)
        predictions = hits[hits["element_type"] == "VIRULENCE"].to_dict(orient="records")
        results: ElementTypeResult = _parse_amrfinder_vir_results(predictions)
    return MethodIndex(
        type=ElementType.VIR, software=Software.AMRFINDER, result=results
    )


def parse_mykrobe_amr_pred(prediction: Dict[str, Any], resistance_category) -> Tuple[SoupVersions, ElementTypeResult]:
    """Parse mykrobe resistance prediction results."""
    LOG.info("Parsing mykrobe prediction")
    meta = [
        SoupVersion(
            **{
                "name": "mykrobe-predictor",
                "version": prediction["version"]["mykrobe-predictor"],
                "type": "database",
            }
        )
    ]
    pred = prediction["susceptibility"]
    resistance = ElementTypeResult(
        phenotypes=_get_mykrobe_amr_sr_profie(pred),
        genes=[], #_parse_mykrobe_amr_genes(pred),
        mutations=_parse_mykrobe_amr_variants(pred),
    )
    return MethodIndex(type=resistance_category, software=Software.MYKROBE, result=resistance)


def parse_tbprofiler_amr_pred(prediction: Dict[str, Any], resistance_category) -> Tuple[SoupVersions, ElementTypeResult]:
    """Parse tbprofiler resistance prediction results."""
    LOG.info("Parsing tbprofiler prediction")
    db_info = [
            SoupVersion(
                **{
                    "name": prediction["db_version"]["name"],
                    "version": prediction["db_version"]["commit"],
                    "type": "database",
                }
            )
        ]
    resistance = ElementTypeResult(
        phenotypes=_get_tbprofiler_amr_sr_profie(prediction),
        genes=[],
        mutations=_parse_tbprofiler_amr_variants(prediction),
    )
    return MethodIndex(type=resistance_category, software=Software.TBPROFILER, result=resistance)
