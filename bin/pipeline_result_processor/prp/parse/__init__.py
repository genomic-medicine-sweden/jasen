"""Parse output of softwares in pipeline."""

from .phenotype import parse_resistance_pred, parse_virulence_pred
from .qc import parse_quast_results
from .specie import parse_kraken_result
from .typing import parse_cgmlst_results, parse_mlst_results
