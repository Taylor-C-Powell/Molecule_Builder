#!/usr/bin/env python3
"""Build molbuilder/data/ord_conditions.json from Open Reaction Database files.

Requires: pip install molbuilder[ord]   (installs ord-schema)

Usage:
    python scripts/build_ord_data.py /path/to/ord/data/

The script walks the given directory for .pb.gz protobuf files, extracts
reaction conditions and outcomes, classifies each reaction by named
reaction type, and writes aggregated statistics to
molbuilder/data/ord_conditions.json.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
#  Named reaction keyword mapping
# ---------------------------------------------------------------------------

# Maps keyword substrings (lowercased) found in ORD reaction identifiers or
# procedure text to the MolBuilder ORD key used in ord_conditions.json.
NAMED_REACTION_KEYWORDS: dict[str, str] = {
    "suzuki": "suzuki_coupling",
    "heck": "heck_coupling",
    "sonogashira": "sonogashira_coupling",
    "buchwald": "buchwald_hartwig_amination",
    "hartwig": "buchwald_hartwig_amination",
    "grignard": "grignard_reaction",
    "wittig": "wittig_olefination",
    "swern": "swern_oxidation",
    "dess-martin": "dess_martin_oxidation",
    "dess martin": "dess_martin_oxidation",
    "tempo": "tempo_oxidation",
    "jones": "jones_oxidation",
    "aldol": "aldol_condensation",
    "friedel-crafts": "friedel_crafts_alkylation",
    "friedel crafts": "friedel_crafts_alkylation",
    "diels-alder": "diels_alder",
    "diels alder": "diels_alder",
    "mitsunobu": "mitsunobu_reaction",
    "appel": "appel_reaction",
    "stille": "stille_coupling",
    "negishi": "negishi_coupling",
    "boc protection": "boc_protection",
    "boc deprotection": "boc_deprotection",
    "hydrogenation": "catalytic_hydrogenation",
    "fischer esterif": "fischer_esterification",
    "williamson": "williamson_ether_synthesis",
    "reductive amination": "n_alkylation",
    "chan-lam": "chan_lam_coupling",
    "chan lam": "chan_lam_coupling",
    "ullmann": "ullmann_coupling",
    "kumada": "kumada_coupling",
    "hiyama": "hiyama_coupling",
    "miyaura borylation": "miyaura_borylation",
    "tsuji": "tsuji_trost_reaction",
    "trost": "tsuji_trost_reaction",
    "click": "cuaac_click",
    "huisgen": "cuaac_click",
    "metathesis": "ring_closing_metathesis",
    "olefin metathesis": "cross_metathesis",
    "wacker": "wacker_oxidation",
    "sharpless epox": "sharpless_epoxidation",
    "sharpless dihydrox": "sharpless_dihydroxylation",
    "sharpless ad": "sharpless_dihydroxylation",
    "jacobsen": "jacobsen_epoxidation",
    "shi epox": "shi_epoxidation",
    "baeyer-villiger": "baeyer_villiger_oxidation",
    "baeyer villiger": "baeyer_villiger_oxidation",
    "pinnick": "pinnick_oxidation",
    "ibx": "ibx_oxidation",
    "tpap": "tpap_oxidation",
    "pcc": "pcc_oxidation",
    "permanganate": "permanganate_oxidation",
    "ozonolysis": "ozonolysis",
    "nabh4": "nabh4_reduction",
    "sodium borohydride": "nabh4_reduction",
    "lialh4": "lialh4_reduction",
    "lithium aluminum hydride": "lialh4_reduction",
    "lithium aluminium hydride": "lialh4_reduction",
    "dibal": "dibal_reduction",
    "birch reduction": "birch_reduction",
    "wolff-kishner": "wolff_kishner_reduction",
    "wolff kishner": "wolff_kishner_reduction",
    "clemmensen": "clemmensen_reduction",
    "luche": "luche_reduction",
    "lindlar": "lindlar_hydrogenation",
    "cbs reduction": "cbs_reduction",
    "noyori": "noyori_asymmetric_hydrogenation",
    "evans aldol": "evans_aldol",
    "michael addition": "michael_addition",
    "robinson annulation": "robinson_annulation",
    "mannich": "mannich_reaction",
    "strecker": "strecker_synthesis",
    "ugi": "ugi_reaction",
    "passerini": "passerini_reaction",
    "beckmann": "beckmann_rearrangement",
    "curtius": "curtius_rearrangement",
    "claisen rearrangement": "claisen_rearrangement",
    "cope rearrangement": "cope_rearrangement",
    "horner-wadsworth": "horner_wadsworth_emmons",
    "horner wadsworth": "horner_wadsworth_emmons",
    "knoevenagel": "knoevenagel_condensation",
    "weinreb": "weinreb_amide",
    "reformatsky": "reformatsky_reaction",
    "henry reaction": "henry_reaction",
    "baylis-hillman": "baylis_hillman_reaction",
    "baylis hillman": "baylis_hillman_reaction",
    "sandmeyer": "sandmeyer_reaction",
    "finkelstein": "finkelstein_reaction",
    "gabriel synthesis": "gabriel_synthesis",
    "wohl-ziegler": "wohl_ziegler_bromination",
    "wohl ziegler": "wohl_ziegler_bromination",
    "minisci": "minisci_reaction",
    "barton decarboxylation": "barton_decarboxylation",
    "amide coupling": "amide_coupling",
    "edc": "amide_coupling",
    "hatu": "amide_coupling",
    "fischer indole": "fischer_indole_synthesis",
    "paal-knorr": "paal_knorr_synthesis",
    "paal knorr": "paal_knorr_synthesis",
    "gewald": "gewald_reaction",
    "hantzsch": "hantzsch_pyridine_synthesis",
    "vilsmeier": "vilsmeier_haack_reaction",
    "romp": "romp",
    "staudinger": "staudinger_reduction",
    "hydroboration": "hydroboration_oxidation",
    "simmons-smith": "simmons_smith_reaction",
    "simmons smith": "simmons_smith_reaction",
    "oppenauer": "oppenauer_oxidation",
    "fmoc protection": "fmoc_protection",
    "fmoc deprotection": "fmoc_deprotection",
    "tbs protection": "tbs_protection",
    "tbs deprotection": "tbs_deprotection",
    "tbaf": "tbs_deprotection",
    "benzyl protection": "benzyl_protection",
    "hydrogenolysis": "hydrogenolysis",
    "acetylation": "acetylation",
}


def _classify_reaction(reaction) -> str | None:
    """Classify an ORD Reaction protobuf by named reaction.

    Returns an ORD key string or None if unclassifiable.
    """
    # 1. Check reaction identifiers for named reaction labels
    for ident in reaction.identifiers:
        if ident.type == ident.REACTION_TYPE or ident.type == ident.NAME:
            text = ident.value.lower()
            for keyword, key in NAMED_REACTION_KEYWORDS.items():
                if keyword in text:
                    return key

    # 2. Check procedure details text
    proc = reaction.notes.procedure_details.lower() if reaction.notes else ""
    if proc:
        for keyword, key in NAMED_REACTION_KEYWORDS.items():
            if keyword in proc:
                return key

    return None


def _extract_temperature_c(conditions) -> float | None:
    """Extract temperature in Celsius from ORD TemperatureConditions."""
    tc = conditions.temperature
    if not tc.HasField("setpoint"):
        return None
    sp = tc.setpoint
    if sp.units == sp.CELSIUS:
        return sp.value
    if sp.units == sp.KELVIN:
        return sp.value - 273.15
    if sp.units == sp.FAHRENHEIT:
        return (sp.value - 32.0) * 5.0 / 9.0
    return None


def _extract_pressure_atm(conditions) -> float | None:
    """Extract pressure in atm from ORD PressureConditions."""
    pc = conditions.pressure
    if not pc.HasField("setpoint"):
        return None
    sp = pc.setpoint
    if sp.units == sp.ATMOSPHERE:
        return sp.value
    if sp.units == sp.BAR:
        return sp.value * 0.986923
    if sp.units == sp.PASCAL:
        return sp.value / 101325.0
    if sp.units == sp.TORR or sp.units == sp.MMHG:
        return sp.value / 760.0
    if sp.units == sp.PSI:
        return sp.value / 14.696
    return None


def _extract_atmosphere(conditions) -> str | None:
    """Extract atmosphere string from ORD conditions."""
    atm = conditions.pressure.atmosphere
    if not atm:
        return None
    type_map = {1: "air", 2: "N2", 3: "Ar", 4: "O2", 5: "H2"}
    return type_map.get(atm.type, None)


def _extract_solvents(reaction) -> list[str]:
    """Extract solvent names from ORD reaction inputs."""
    solvents = []
    for inp in reaction.inputs.values():
        for comp in inp.components:
            if comp.reaction_role == comp.SOLVENT:
                for ident in comp.identifiers:
                    if ident.type == ident.NAME:
                        solvents.append(ident.value)
                        break
    return solvents


def _extract_catalysts(reaction) -> list[str]:
    """Extract catalyst names from ORD reaction inputs."""
    catalysts = []
    for inp in reaction.inputs.values():
        for comp in inp.components:
            if comp.reaction_role == comp.CATALYST:
                for ident in comp.identifiers:
                    if ident.type == ident.NAME:
                        catalysts.append(ident.value)
                        break
    return catalysts


def _extract_yield(reaction) -> float | None:
    """Extract yield percentage from ORD reaction outcomes."""
    for outcome in reaction.outcomes:
        for product in outcome.products:
            for meas in product.measurements:
                if meas.type == meas.YIELD:
                    if meas.HasField("percentage"):
                        return meas.percentage.value
    return None


def _extract_time_hours(conditions) -> float | None:
    """Extract reaction time in hours from ORD conditions."""
    tc = conditions.conditions_are_dynamic
    # Check for reaction time in the main conditions
    rt = conditions.time
    if rt and rt.HasField("value"):
        val = rt.value
        if rt.units == rt.HOUR:
            return val
        if rt.units == rt.MINUTE:
            return val / 60.0
        if rt.units == rt.SECOND:
            return val / 3600.0
    return None


def _percentile(data: list[float], p: float) -> float:
    """Compute percentile (0-100) of a sorted list."""
    if not data:
        return 0.0
    data_sorted = sorted(data)
    k = (len(data_sorted) - 1) * p / 100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return data_sorted[int(k)]
    return data_sorted[f] * (c - k) + data_sorted[c] * (k - f)


def _aggregate(records: list[dict]) -> dict[str, Any]:
    """Aggregate a list of per-reaction records into summary statistics."""
    temps = [r["temperature_C"] for r in records if r["temperature_C"] is not None]
    pressures = [r["pressure_atm"] for r in records if r["pressure_atm"] is not None]
    yields = [r["yield_pct"] for r in records if r["yield_pct"] is not None]
    times = [r["time_hours"] for r in records if r["time_hours"] is not None]

    # Solvent counts
    solvent_counts: dict[str, int] = defaultdict(int)
    for r in records:
        for s in r["solvents"]:
            solvent_counts[s] += 1
    total_solvent = sum(solvent_counts.values()) or 1
    top_solvents = sorted(solvent_counts.items(), key=lambda x: x[1], reverse=True)[:5]

    # Catalyst counts
    catalyst_counts: dict[str, int] = defaultdict(int)
    for r in records:
        for c in r["catalysts"]:
            catalyst_counts[c] += 1
    total_catalyst = sum(catalyst_counts.values()) or 1
    top_catalysts = sorted(catalyst_counts.items(), key=lambda x: x[1], reverse=True)[:3]

    # Atmosphere counts
    atmo_counts: dict[str, int] = defaultdict(int)
    for r in records:
        if r["atmosphere"]:
            atmo_counts[r["atmosphere"]] += 1
    total_atmo = sum(atmo_counts.values()) or 1

    result: dict[str, Any] = {"n": len(records)}

    if temps:
        result["temperature_C"] = {
            "median": round(statistics.median(temps), 1),
            "p25": round(_percentile(temps, 25), 1),
            "p75": round(_percentile(temps, 75), 1),
        }
    if pressures:
        result["pressure_atm"] = {
            "median": round(statistics.median(pressures), 2),
            "p25": round(_percentile(pressures, 25), 2),
            "p75": round(_percentile(pressures, 75), 2),
        }
    if top_solvents:
        result["solvents"] = [
            {
                "name": name,
                "count": count,
                "fraction": round(count / total_solvent, 2),
            }
            for name, count in top_solvents
        ]
    if top_catalysts:
        result["catalysts"] = [
            {
                "name": name,
                "count": count,
                "fraction": round(count / total_catalyst, 2),
            }
            for name, count in top_catalysts
        ]
    if atmo_counts:
        result["atmosphere"] = {
            k: round(v / total_atmo, 2) for k, v in atmo_counts.items()
        }
    if yields:
        result["yield_pct"] = {
            "median": round(statistics.median(yields), 1),
            "p25": round(_percentile(yields, 25), 1),
            "p75": round(_percentile(yields, 75), 1),
        }
    if times:
        result["reaction_time_hours"] = {
            "median": round(statistics.median(times), 1),
            "p25": round(_percentile(times, 25), 1),
            "p75": round(_percentile(times, 75), 1),
        }

    return result


def process_ord_directory(data_dir: str) -> dict[str, Any]:
    """Walk *data_dir* for .pb.gz files and aggregate reaction statistics."""
    try:
        from ord_schema.proto import reaction_pb2
    except ImportError:
        print(
            "ERROR: ord-schema is not installed.\n"
            "Install it with: pip install molbuilder[ord]",
            file=sys.stderr,
        )
        sys.exit(1)

    classified: dict[str, list[dict]] = defaultdict(list)
    total_processed = 0
    total_classified = 0

    data_path = Path(data_dir)
    pb_files = list(data_path.rglob("*.pb.gz")) + list(data_path.rglob("*.pb"))

    if not pb_files:
        print(f"WARNING: No .pb.gz or .pb files found in {data_dir}", file=sys.stderr)

    for pb_file in pb_files:
        try:
            if str(pb_file).endswith(".gz"):
                with gzip.open(pb_file, "rb") as f:
                    raw = f.read()
            else:
                with open(pb_file, "rb") as f:
                    raw = f.read()

            dataset = reaction_pb2.Dataset()
            dataset.ParseFromString(raw)

            for reaction in dataset.reactions:
                total_processed += 1
                rxn_key = _classify_reaction(reaction)
                if not rxn_key:
                    continue

                total_classified += 1
                cond = reaction.conditions if reaction.HasField("conditions") else None

                record = {
                    "temperature_C": _extract_temperature_c(cond) if cond else None,
                    "pressure_atm": _extract_pressure_atm(cond) if cond else None,
                    "atmosphere": _extract_atmosphere(cond) if cond else None,
                    "solvents": _extract_solvents(reaction),
                    "catalysts": _extract_catalysts(reaction),
                    "yield_pct": _extract_yield(reaction),
                    "time_hours": _extract_time_hours(cond) if cond else None,
                }
                classified[rxn_key].append(record)

        except Exception as e:
            print(f"WARNING: Failed to process {pb_file}: {e}", file=sys.stderr)
            continue

    # Aggregate
    reactions: dict[str, Any] = {}
    for key in sorted(classified.keys()):
        records = classified[key]
        if len(records) >= 3:  # require at least 3 data points
            reactions[key] = _aggregate(records)

    output = {
        "_meta": {
            "ord_version": "0.3.99",
            "generated": __import__("datetime").date.today().isoformat(),
            "total_reactions_processed": total_processed,
            "matched_reactions": total_classified,
            "license": "CC-BY-SA-4.0",
        },
        "reactions": reactions,
    }

    print(f"Processed {total_processed} reactions, classified {total_classified}")
    print(f"Generated stats for {len(reactions)} reaction types")

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Build ORD condition statistics for MolBuilder"
    )
    parser.add_argument(
        "data_dir",
        help="Path to directory containing ORD .pb.gz files",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output JSON path (default: molbuilder/data/ord_conditions.json)",
    )
    args = parser.parse_args()

    if args.output is None:
        # Default to the package data directory
        script_dir = Path(__file__).resolve().parent
        output_path = script_dir.parent / "molbuilder" / "data" / "ord_conditions.json"
    else:
        output_path = Path(args.output)

    result = process_ord_directory(args.data_dir)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=True)

    print(f"Wrote {output_path} ({output_path.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
