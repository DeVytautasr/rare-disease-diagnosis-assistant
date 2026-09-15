"""Derive the evidence scoring tiers from bam_tools' own source.

Requirement (Phase 7): the ceiling shown in the UI must not be a literal. If
the scoring function is revised the displayed ceiling has to move with it, so
the tiers are read out of the live source rather than copied. If the source
stops matching the shape this parser understands, derivation FAILS LOUDLY and
the UI shows "not derivable" instead of a stale number -- never a figure that
cannot be recomputed.
"""
import ast
import inspect

from stage1_igv_assistant.tools import bam_tools

# variable name in bam_tools -> (layer key, the observed field the tier reads)
_LAYERS = {
    "discordant_pair_score": ("discordant_pairs", "discordant_fraction"),
    "soft_clip_score":       ("soft_clipped_reads", "max_clips_at_position"),
    "split_read_score":      ("split_reads", "split_read_fraction"),
    # the depth ladder assigns to `candidate_score`, which is then gated on
    # dip_is_at_focus before it becomes depth_score
    "candidate_score":       ("read_depth", "depth_ratio_min_to_mean"),
}

# layers whose tier result is additionally gated before it can score
_EXTRA_GATE = {"read_depth": "dip_is_at_focus must be true; a dip not localised "
                             "to the queried position is suppressed to 0"}


class TierDerivationError(RuntimeError):
    pass


def _const(node):
    """Literal, or a module-level constant in bam_tools referenced by name."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.Name):
        v = getattr(bam_tools, node.id, None)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return float(v)
    return None


def _scan(fn_node):
    """Collect [(score_var, op, threshold, score_value)] from if/elif chains."""
    found = {}
    for node in ast.walk(fn_node):
        if not isinstance(node, ast.If):
            continue
        cur = node
        while isinstance(cur, ast.If):
            test = cur.test
            comps = [test] if isinstance(test, ast.Compare) else (
                [v for v in test.values if isinstance(v, ast.Compare)]
                if isinstance(test, ast.BoolOp) else [])
            assigned = None
            for stmt in cur.body:
                if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1 \
                        and isinstance(stmt.targets[0], ast.Name) \
                        and stmt.targets[0].id in _LAYERS:
                    val = _const(stmt.value)
                    if val is not None:
                        assigned = (stmt.targets[0].id, val)
            if assigned and comps:
                c = comps[0]
                thr = _const(c.comparators[0])
                if thr is not None:
                    op = type(c.ops[0]).__name__
                    found.setdefault(assigned[0], []).append(
                        {"op": {"GtE": ">=", "Gt": ">", "LtE": "<=", "Lt": "<"}.get(op, op),
                         "threshold": thr, "score": assigned[1]})
            cur = cur.orelse[0] if len(cur.orelse) == 1 and isinstance(cur.orelse[0], ast.If) else None
    return found


def _scan_bands(fn_node):
    """Collect [(strength, op, threshold)] from the evidence_strength if/elif
    chain. Only branches whose test compares `evidence_score` to a constant are
    taken: the QUALITY-LIMITED branch tests the MAPQ gate, not the score, and is
    deliberately not a band."""
    found = {}
    for node in ast.walk(fn_node):
        if not isinstance(node, ast.If) or not isinstance(node.test, ast.Compare):
            continue
        c = node.test
        if not (isinstance(c.left, ast.Name) and c.left.id == "evidence_score"):
            continue
        thr = _const(c.comparators[0])
        if thr is None:
            continue
        for stmt in node.body:
            if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1 \
                    and isinstance(stmt.targets[0], ast.Name) \
                    and stmt.targets[0].id == "evidence_strength" \
                    and isinstance(stmt.value, ast.Constant) \
                    and isinstance(stmt.value.value, str):
                op = type(c.ops[0]).__name__
                found[stmt.value.value] = {
                    "strength": stmt.value.value,
                    "op": {"GtE": ">=", "Gt": ">", "LtE": "<=", "Lt": "<"}.get(op, op),
                    "threshold": thr}
    return found


def _parsed_summary_fn():
    src = inspect.getsource(bam_tools.summarize_breakpoint_evidence)
    tree = ast.parse(inspect.cleandoc(src) if src.startswith("def") else src)
    fn = next((n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name == "summarize_breakpoint_evidence"), None)
    if fn is None:
        raise TierDerivationError("summarize_breakpoint_evidence not found in parsed source")
    return fn


def derive_bands():
    """The evidence_score -> evidence_strength ladder, read from the live source.

    Phase 9 measured that no model at any tier could state the ceiling argument,
    because the score at which "strong" begins appears in no tool return and no
    tool description -- so "this cannot reach strong" was unstatable. This makes
    the boundary derivable instead of a literal, on the same terms as the tiers:
    if the bands are ever revised the reported ceiling moves with them, and if
    the chain stops matching this parser it FAILS LOUDLY rather than reporting a
    stale number.

    Returned highest-threshold-first: [{"strength","op","threshold"}, ...]
    """
    found = _scan_bands(_parsed_summary_fn())
    if "strong" not in found:
        raise TierDerivationError(
            "could not derive the 'strong' band boundary — the evidence_strength "
            "chain structure changed")
    return sorted(found.values(), key=lambda b: -b["threshold"])


def strong_band(bands):
    """The score at which the top band begins."""
    for b in bands:
        if b["strength"] == "strong":
            return b["threshold"]
    return None


# Which of the four layers is structurally suppressed by a balanced rearrangement,
# and which is structurally capped by heterozygosity. Both are properties of the
# event, not of the run: no deeper sequencing or wider window moves either.
_BALANCED_ZERO_LAYER = "read_depth"
_HET_CAPPED_LAYER = "discordant_pairs"


def ceiling_from_observed(tiers, bands, observed):
    """Attainable-score analysis, derived from the live tiers and bands.

    Shared by the UI panel and the MCP summary tool so there is one derivation,
    not two that can drift. `observed` is {layer: observed_value_or_None}.
    """
    if tiers is None:
        return {"derivable": False, "reason": "tiers not derivable"}
    per, flat = {}, 0.0
    for layer, t in tiers.items():
        v = observed.get(layer)
        per[layer] = {
            "observed_field": t["observed_field"], "observed": v,
            "score_now": score_for(t, v) if v is not None else None,
            "max_score": t["max_score"],
            "next_tier": next_tier_up(t, v) if v is not None else None,
            "tiers": t["tiers"], "extra_gate": t.get("extra_gate"),
        }
        if layer != _BALANCED_ZERO_LAYER:
            flat += t["max_score"]
    # Locus-specific ceiling: hold the paired-read layer at the tier its OBSERVED
    # fraction actually reaches (heterozygosity, not run quality, is what pins it
    # there), let the two read-level layers reach their top tier, and let depth
    # contribute nothing because a balanced event leaves copy number unchanged.
    disc_now = per.get(_HET_CAPPED_LAYER, {}).get("score_now")
    attainable = None
    if disc_now is not None:
        attainable = disc_now
        for k, t in tiers.items():
            if k not in (_HET_CAPPED_LAYER, _BALANCED_ZERO_LAYER):
                attainable += t["max_score"]
    sb = strong_band(bands) if bands else None
    disc_top = None
    if _HET_CAPPED_LAYER in tiers:
        disc_top = max(t["threshold"] for t in tiers[_HET_CAPPED_LAYER]["tiers"])
    basis = (
        "the paired-read measurement held at the band its measured value reaches — in a "
        "HETEROZYGOUS rearrangement roughly half the reads crossing the breakpoint come from "
        "the intact homolog, so the discordant fraction is capped well below the "
        + (f"{disc_top} its top band requires" if disc_top is not None else "top band's threshold")
        + " and no amount of extra coverage or a wider window moves it; the soft-clip and "
        "split-read measurements at their best possible band; and the depth measurement "
        "contributing nothing, because a balanced rearrangement gains and loses no DNA.")
    return {
        "derivable": True, "per_layer": per,
        "max_with_flat_depth": flat,
        "max_all_layers": sum(t["max_score"] for t in tiers.values()),
        "attainable_here": attainable,
        "attainable_basis": basis,
        "bands": bands,
        "strong_band": sb,
        "strong_band_reachable_here": (None if (attainable is None or sb is None)
                                       else attainable >= sb),
        "note": ("Calculated from the scoring rules in force right now, so this stays correct "
                 "if those rules are changed. The flat-depth figure is the highest score "
                 "reachable when the depth measurement contributes nothing, which is what a "
                 "balanced rearrangement should look like."),
    }


def ceiling_sentence(ceiling):
    """One plain sentence stating the consequence, for the observation channel.

    Phase 8c and Phase 9 both measured that models quote a tool's sentences far
    more readily than they recompute from its fields, so the conclusion is
    offered as prose as well as as numbers.
    """
    if not ceiling.get("derivable"):
        return None
    a, sb = ceiling.get("attainable_here"), ceiling.get("strong_band")
    if a is None or sb is None:
        return None
    if a >= sb:
        return (f"Attainable-score note: for a heterozygous balanced rearrangement the highest "
                f"score reachable at this locus is {a}/100, which does reach the {sb}/100 at "
                f"which the top band begins.")
    return (f"Attainable-score note: for a heterozygous balanced rearrangement the highest score "
            f"reachable at this locus is {a}/100, below the {sb}/100 at which the top band "
            f"begins — so the top band is UNREACHABLE here by arithmetic, and a middling band is "
            f"the ceiling rather than a weak result. Two of the four layers are structurally "
            f"limited for such an event: the depth layer contributes nothing because no DNA is "
            f"gained or lost, and the paired-read layer is capped because roughly half the reads "
            f"crossing the breakpoint come from the intact homolog. Judge a balanced event on the "
            f"four measurements, not on the band it lands in.")


def derive_tiers():
    """{layer: {"observed_field":..., "tiers":[{op,threshold,score}...], "max_score":float}}"""
    found = _scan(_parsed_summary_fn())
    missing = set(_LAYERS) - set(found)
    if missing:
        raise TierDerivationError(
            f"could not derive tiers for {sorted(missing)} — scoring function structure changed")
    out = {}
    for var, (layer, field) in _LAYERS.items():
        seen, tiers = set(), []
        for t in sorted(found[var], key=lambda t: -t["score"]):
            key = (t["op"], t["threshold"], t["score"])
            if key not in seen:
                seen.add(key); tiers.append(t)
        out[layer] = {"observed_field": field, "tiers": tiers,
                      "max_score": max(t["score"] for t in tiers),
                      "extra_gate": _EXTRA_GATE.get(layer)}
    return out


def score_for(layer_tiers, observed):
    """Replays the derived tier ladder for one observed value."""
    if observed is None:
        return None
    for t in layer_tiers["tiers"]:
        ok = (observed >= t["threshold"] if t["op"] == ">=" else
              observed > t["threshold"] if t["op"] == ">" else
              observed <= t["threshold"] if t["op"] == "<=" else
              observed < t["threshold"])
        if ok:
            return t["score"]
    return 0.0


def next_tier_up(layer_tiers, observed):
    """The lowest tier strictly above what `observed` currently attains."""
    cur = score_for(layer_tiers, observed) or 0.0
    above = [t for t in layer_tiers["tiers"] if t["score"] > cur]
    return min(above, key=lambda t: t["score"]) if above else None
