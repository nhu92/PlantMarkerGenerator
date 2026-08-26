#!/usr/bin/env python3
"""Shared functions for the chloroplast annotation cleanup pipeline."""

from __future__ import annotations

import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable, TextIO

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


DEFAULT_CONFIG: dict[str, Any] = {
    "annotator_priority": {
        "Chloe, blatX; merged": 3,
        "Chloe": 2,
        "blatX": 1,
    },
    "valid_start_codons": ["ATG", "GTG", "TTG"],
    "valid_stop_codons": ["TAA", "TAG", "TGA"],
    "rna_edit_start_genes": ["ndhD", "psbL", "rpl2"],
    "rna_edit_stop_genes": ["petD"],
    "editable_stop_codons": ["CAA", "CAG", "CGA"],
    "trans_spliced_genes": ["rps12"],
    "maximum_stop_extension_codons": 20,
    "maximum_length_difference_aa": 2,
    "minimum_cds_aa": 10,
    "minimum_cds_fraction_of_median": 0.5,
    "trna_minimum_nt": 55,
    "trna_maximum_nt": 120,
    "trna_expected_nt": 75,
    "expected_rrna_lengths": {
        "rrn4.5": 103,
        "rrn5": 121,
        "rrn16": 1491,
        "rrn23": 2810,
    },
}


def load_config(path: Path | None) -> dict[str, Any]:
    """Load JSON configuration and merge it over conservative defaults."""
    config = deepcopy(DEFAULT_CONFIG)
    if path is None:
        return config
    supplied = json.loads(path.read_text(encoding="utf-8"))
    for key, value in supplied.items():
        if isinstance(value, dict) and isinstance(config.get(key), dict):
            config[key].update(value)
        else:
            config[key] = value
    return config


def first_qualifier(feature: SeqFeature | None, *keys: str, default: str = "") -> str:
    """Return the first non-empty qualifier value from a feature."""
    if feature is None:
        return default
    for key in keys:
        values = feature.qualifiers.get(key)
        if values:
            return str(values[0])
    return default


def feature_gene(feature: SeqFeature) -> str:
    """Return the biological gene symbol when one is available."""
    return first_qualifier(feature, "gene")


def sequence_digest(sequence) -> str:
    """Return a case-insensitive SHA-256 digest for sequence matching."""
    return hashlib.sha256(bytes(sequence.upper())).hexdigest()


def read_records(path: Path, fmt: str) -> list[SeqRecord]:
    """Read all records and fail early on an empty input file."""
    records = list(SeqIO.parse(path, fmt))
    if not records:
        raise ValueError(f"No {fmt} records were found in {path}")
    return records


def match_fasta_to_genbank(
    fasta_records: list[SeqRecord], genbank_records: list[SeqRecord]
) -> list[tuple[SeqRecord, SeqRecord]]:
    """Match FASTA records to GenBank records by exact sequence content."""
    by_digest: dict[str, list[SeqRecord]] = defaultdict(list)
    for record in genbank_records:
        by_digest[sequence_digest(record.seq)].append(record)

    matched: list[tuple[SeqRecord, SeqRecord]] = []
    used_objects: set[int] = set()
    for fasta_record in fasta_records:
        digest = sequence_digest(fasta_record.seq)
        candidates = [r for r in by_digest.get(digest, []) if id(r) not in used_objects]
        if not candidates:
            raise ValueError(
                f"FASTA record {fasta_record.id!r} has no exact sequence match in the GenBank file"
            )
        source = candidates[0]
        used_objects.add(id(source))
        matched.append((fasta_record, source))

    if len(matched) != len(genbank_records):
        raise ValueError(
            "FASTA and GenBank record counts differ or one GenBank sequence was not matched"
        )
    return matched


def infer_binomial(text: str) -> str:
    """Infer a simple scientific name from an identifier or description."""
    tokens = [token for token in re.split(r"[_\s]+", text.strip()) if token]
    for index in range(len(tokens) - 1):
        genus, species = tokens[index], tokens[index + 1]
        if (
            genus[:1].isupper()
            and genus[1:].islower()
            and genus.isalpha()
            and species[:1].islower()
            and re.fullmatch(r"[a-z-]+", species)
        ):
            name = f"{genus} {species}"
            if index + 3 < len(tokens) and tokens[index + 2] in {"var", "subsp", "f"}:
                name += f" {tokens[index + 2]}. {tokens[index + 3]}"
            return name
    return ""


def source_metadata(fasta_record: SeqRecord, genbank_record: SeqRecord) -> dict[str, str]:
    """Collect source metadata without performing taxonomy normalization."""
    source = next((f for f in genbank_record.features if f.type == "source"), None)
    organism = first_qualifier(source, "organism")
    if (
        not organism
        or " " not in organism
        or organism.lower() in {"unknown", "unidentified plant", "synthetic construct"}
    ):
        organism = infer_binomial(fasta_record.description)
    if not organism:
        organism = infer_binomial(genbank_record.description) or "unidentified plant"

    isolate = first_qualifier(source, "isolate", "strain", "cultivar")
    return {
        "organism": organism,
        "isolate": isolate,
        "collection_date": first_qualifier(source, "collection_date"),
        "geo_loc_name": first_qualifier(source, "geo_loc_name", "country"),
        "specimen_voucher": first_qualifier(source, "specimen_voucher"),
    }


def coding_state(record: SeqRecord, feature: SeqFeature) -> dict[str, Any]:
    """Calculate the unedited table-11 translation state of a CDS."""
    nucleotide = feature.extract(record.seq).upper()
    try:
        codon_start = int(first_qualifier(feature, "codon_start", default="1"))
    except ValueError:
        codon_start = 1
    nucleotide = nucleotide[codon_start - 1 :]
    usable = nucleotide[: len(nucleotide) - len(nucleotide) % 3]
    protein = str(usable.translate(table=11, cds=False))
    internal_stops = protein[:-1].count("*") if protein.endswith("*") else protein.count("*")
    return {
        "nucleotide": nucleotide,
        "protein": protein,
        "start": str(nucleotide[:3]),
        "stop": str(nucleotide[-3:]),
        "modulo_three": len(nucleotide) % 3,
        "internal_stops": internal_stops,
        "amino_acids": len(protein.rstrip("*")),
    }


def locations_overlap(left: SeqFeature, right: SeqFeature) -> bool:
    """Return True when any pair of feature intervals overlaps."""
    for left_part in left.location.parts:
        for right_part in right.location.parts:
            if int(left_part.start) < int(right_part.end) and int(right_part.start) < int(left_part.end):
                return True
    return False


def overlap_clusters(features: list[SeqFeature]) -> list[list[SeqFeature]]:
    """Group features connected by direct or transitive interval overlap."""
    clusters: list[list[SeqFeature]] = []
    unused = set(range(len(features)))
    while unused:
        seed = unused.pop()
        group = [seed]
        queue = [seed]
        while queue:
            left_index = queue.pop()
            for right_index in list(unused):
                if locations_overlap(features[left_index], features[right_index]):
                    unused.remove(right_index)
                    group.append(right_index)
                    queue.append(right_index)
        clusters.append([features[index] for index in group])
    return clusters


def build_gene_length_medians(records: Iterable[SeqRecord], config: dict[str, Any]) -> dict[str, float]:
    """Estimate robust protein-length expectations from intact batch models."""
    starts = set(config["valid_start_codons"])
    stops = set(config["valid_stop_codons"])
    lengths: dict[str, list[int]] = defaultdict(list)
    for record in records:
        for feature in record.features:
            gene = feature_gene(feature)
            if feature.type != "CDS" or not gene or not first_qualifier(feature, "product"):
                continue
            state = coding_state(record, feature)
            if (
                not state["modulo_three"]
                and not state["internal_stops"]
                and state["start"] in starts
                and state["stop"] in stops
            ):
                lengths[gene].append(int(state["amino_acids"]))
    return {gene: float(statistics.median(values)) for gene, values in lengths.items() if values}


def build_gene_product_consensus(records: Iterable[SeqRecord]) -> dict[str, str]:
    """Choose the most frequent non-empty product name for each gene in a batch."""
    products: dict[str, Counter[str]] = defaultdict(Counter)
    for record in records:
        for feature in record.features:
            gene = feature_gene(feature)
            product = first_qualifier(feature, "product")
            if feature.type == "CDS" and gene and product and "fragment" not in gene.lower():
                products[gene][product] += 1
    return {
        gene: sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
        for gene, counts in products.items()
    }


def feature_score(
    record: SeqRecord, feature: SeqFeature, median_length: float | None, config: dict[str, Any]
) -> tuple[Any, ...]:
    """Rank alternative annotation models using frame, evidence, and length."""
    state = coding_state(record, feature)
    gene = feature_gene(feature)
    start_ok = state["start"] in set(config["valid_start_codons"])
    start_ok = start_ok or (
        state["start"] == "ACG" and gene in set(config["rna_edit_start_genes"])
    )
    stop_ok = state["stop"] in set(config["valid_stop_codons"])
    stop_ok = stop_ok or (
        state["stop"] in set(config["editable_stop_codons"])
        and gene in set(config["rna_edit_stop_genes"])
    )
    distance = abs(float(state["amino_acids"]) - median_length) if median_length is not None else 0.0
    annotator = first_qualifier(feature, "annotator")
    rank = int(config["annotator_priority"].get(annotator, 0))
    return (
        int(not state["modulo_three"] and not state["internal_stops"]),
        int(bool(first_qualifier(feature, "product"))),
        int(bool(first_qualifier(feature, "translation"))),
        int(start_ok),
        int(stop_ok),
        rank,
        -distance,
    )


def extend_terminal(record: SeqRecord, feature: SeqFeature, codons: int) -> SeqFeature | None:
    """Extend the terminal CDS interval without crossing sequence boundaries."""
    result = deepcopy(feature)
    parts = list(result.location.parts)
    terminal = parts[-1]
    if terminal.strand == 1:
        new_end = int(terminal.end) + 3 * codons
        if new_end > len(record.seq):
            return None
        parts[-1] = FeatureLocation(terminal.start, new_end, strand=1)
    else:
        new_start = int(terminal.start) - 3 * codons
        if new_start < 0:
            return None
        parts[-1] = FeatureLocation(new_start, terminal.end, strand=-1)
    result.location = parts[0] if len(parts) == 1 else CompoundLocation(parts, operator="join")
    return result


def first_downstream_stop(
    record: SeqRecord, feature: SeqFeature, stop_codons: set[str], maximum: int
) -> int | None:
    """Find the first in-frame downstream stop within a bounded window."""
    terminal = list(feature.location.parts)[-1]
    for number in range(1, maximum + 1):
        if terminal.strand == 1:
            start = int(terminal.end) + 3 * (number - 1)
            codon = record.seq[start : start + 3]
        else:
            end = int(terminal.start) - 3 * (number - 1)
            codon = record.seq[end - 3 : end].reverse_complement() if end >= 3 else ""
        if len(codon) == 3 and str(codon.upper()) in stop_codons:
            return number
    return None


def codon_location(location, first: bool) -> str:
    """Format the first or terminal codon for an NCBI transl_except qualifier."""
    part = list(location.parts)[0 if first else -1]
    if part.strand == 1:
        start = int(part.start) + 1 if first else int(part.end) - 2
        return f"{start}..{start + 2}"
    start = int(part.end) - 2 if first else int(part.start) + 1
    return f"complement({start}..{start + 2})"


def clean_cds(
    record: SeqRecord,
    feature: SeqFeature,
    median_length: float | None,
    config: dict[str, Any],
) -> tuple[SeqFeature | None, list[str]]:
    """Return a conservative submission-safe CDS or a documented rejection."""
    result = deepcopy(feature)
    actions: list[str] = []
    gene = feature_gene(result)
    starts = set(config["valid_start_codons"])
    stops = set(config["valid_stop_codons"])
    editable_stops = set(config["editable_stop_codons"])

    if not gene:
        return None, ["dropped: missing biological gene symbol"]
    if not first_qualifier(result, "product"):
        return None, ["dropped: missing product name"]

    codon_start = first_qualifier(result, "codon_start", default="1")
    if codon_start not in {"1", "2", "3"}:
        return None, [f"dropped: invalid codon_start {codon_start!r}"]

    state = coding_state(record, result)
    if state["modulo_three"] or state["internal_stops"]:
        return None, ["dropped: coding frame is not intact"]

    transl_except: list[str] = []
    if state["start"] not in starts:
        if state["start"] == "ACG" and gene in set(config["rna_edit_start_genes"]):
            transl_except.append(f"(pos:{codon_location(result.location, True)},aa:Met)")
            actions.append("RNA-edited start ACG->AUG")
        else:
            return None, [f"dropped: unsupported start codon {state['start']}"]

    terminal_edit = False
    if state["stop"] not in stops:
        current_length = float(state["amino_acids"])
        next_stop = first_downstream_stop(
            record,
            result,
            stops,
            int(config["maximum_stop_extension_codons"]),
        )
        current_distance = (
            abs((current_length - 1) - median_length) if median_length is not None else float("inf")
        )
        extension_distance = (
            abs((current_length + next_stop - 1) - median_length)
            if median_length is not None and next_stop is not None
            else float("inf")
        )
        tolerance = float(config["maximum_length_difference_aa"])
        if (
            gene in set(config["rna_edit_stop_genes"])
            and state["stop"] in editable_stops
            and current_distance <= tolerance
            and current_distance < extension_distance
        ):
            transl_except.append(f"(pos:{codon_location(result.location, False)},aa:TERM)")
            terminal_edit = True
            actions.append(f"RNA-edited stop {state['stop']}->stop")
        elif next_stop is not None and (next_stop <= 2 or extension_distance <= tolerance):
            extended = extend_terminal(record, result, next_stop)
            if extended is None:
                return None, ["dropped: stop extension would cross a sequence boundary"]
            result = extended
            actions.append(f"extended CDS by {3 * next_stop} nt to stop")
        else:
            return None, [f"dropped: no credible stop near terminal codon {state['stop']}"]

    state = coding_state(record, result)
    minimum_length = float(config["minimum_cds_aa"])
    if median_length is not None:
        minimum_length = max(
            minimum_length,
            float(config["minimum_cds_fraction_of_median"]) * median_length,
        )
    if float(state["amino_acids"]) < minimum_length:
        return None, [
            f"dropped: implausibly short CDS ({state['amino_acids']} aa; expected at least {minimum_length:g} aa)"
        ]

    protein = state["protein"]
    if protein.endswith("*"):
        protein = protein[:-1]
    if terminal_edit and protein:
        protein = protein[:-1]
    if transl_except and any("aa:Met" in value for value in transl_except) and protein:
        protein = "M" + protein[1:]

    qualifiers: dict[str, list[str]] = {
        "gene": [gene],
        "product": [first_qualifier(result, "product")],
        "codon_start": [codon_start],
        "transl_table": ["11"],
        "translation": [protein],
    }
    if transl_except:
        qualifiers["transl_except"] = transl_except
    if gene in set(config["trans_spliced_genes"]) and len(result.location.parts) > 1:
        qualifiers["exception"] = ["trans-splicing"]
        actions.append("marked trans-splicing")
    result.qualifiers = qualifiers
    return result, actions


def choose_rnas(features: list[SeqFeature], config: dict[str, Any]) -> list[SeqFeature]:
    """Deduplicate RNA models and reject implausible feature lengths."""
    by_key: dict[tuple[str, str], list[SeqFeature]] = defaultdict(list)
    for feature in features:
        gene = feature_gene(feature)
        product = first_qualifier(feature, "product")
        if not gene or not product or "fragment" in gene.lower():
            continue
        by_key[(feature.type, gene)].append(feature)

    selected: list[SeqFeature] = []
    annotator_rank = config["annotator_priority"]
    for (feature_type, gene), items in by_key.items():
        for cluster in overlap_clusters(items):
            if feature_type == "tRNA":
                plausible = [
                    feature
                    for feature in cluster
                    if int(config["trna_minimum_nt"])
                    <= len(feature.location)
                    <= int(config["trna_maximum_nt"])
                ]
                if not plausible:
                    continue
                target = float(config["trna_expected_nt"])
            else:
                expected = config["expected_rrna_lengths"].get(gene)
                target = float(expected) if expected is not None else float(
                    statistics.median(len(feature.location) for feature in cluster)
                )
                plausible = [
                    feature
                    for feature in cluster
                    if 0.65 * target <= len(feature.location) <= 1.35 * target
                ]
                if not plausible:
                    continue
            winner = min(
                plausible,
                key=lambda feature: (
                    abs(len(feature.location) - target),
                    -int(annotator_rank.get(first_qualifier(feature, "annotator"), 0)),
                ),
            )
            winner = deepcopy(winner)
            winner.qualifiers = {
                "gene": [gene],
                "product": [first_qualifier(winner, "product")],
            }
            selected.append(winner)
    return selected


def is_origin_spanning(feature: SeqFeature, record_length: int) -> bool:
    """Identify compound features that join the two ends of a circular record."""
    return (
        len(feature.location.parts) > 1
        and int(feature.location.end) - int(feature.location.start) > record_length / 2
    )


def gene_location(child: SeqFeature, record_length: int, config: dict[str, Any]):
    """Build a gene feature that contains its CDS or RNA child."""
    parts = list(child.location.parts)
    gene = feature_gene(child)
    if gene in set(config["trans_spliced_genes"]) and len(parts) > 1:
        return deepcopy(child.location)
    if is_origin_spanning(child, record_length):
        return deepcopy(child.location)
    start = min(int(part.start) for part in parts)
    end = max(int(part.end) for part in parts)
    return FeatureLocation(start, end, strand=child.location.strand)


def feature_intervals(location) -> list[tuple[int, int]]:
    """Convert Biopython zero-based intervals to one-based feature-table rows."""
    intervals: list[tuple[int, int]] = []
    for part in location.parts:
        if part.strand == -1:
            intervals.append((int(part.end), int(part.start) + 1))
        else:
            intervals.append((int(part.start) + 1, int(part.end)))
    return intervals


def write_feature_table_feature(handle: TextIO, feature: SeqFeature) -> None:
    """Write one feature in the NCBI five-column table format."""
    for index, (start, stop) in enumerate(feature_intervals(feature.location)):
        handle.write(f"{start}\t{stop}\t{feature.type if index == 0 else ''}\n")
    qualifier_order = ["gene", "product", "codon_start", "transl_table", "exception", "transl_except"]
    for key in qualifier_order:
        # Child features inherit their overlapping gene. Omitting duplicate
        # child /gene values avoids cross-linking identically named IR copies.
        if key == "gene" and feature.type != "gene":
            continue
        for value in feature.qualifiers.get(key, []):
            handle.write(f"\t\t\t{key}\t{value}\n")


def safe_sequence_id(prefix: str, number: int, width: int) -> str:
    """Create a short identifier accepted by NCBI FASTA parsers."""
    clean_prefix = re.sub(r"[^A-Za-z0-9_.:#*-]", "_", prefix)
    return f"{clean_prefix}{number:0{width}d}"
