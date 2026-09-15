"""
chat.py — local-model panel for the breakpoint instrument.

The model is a SECOND CONSUMER of ToolRecorder, never a second path to the
tools. It receives schemas generated from list_tools(), requests calls, and
those calls execute through the same recorder the manual UI uses. Nothing here
can put a number on screen that did not come back from a recorded call.

Ollama is called directly at /api/chat with its native `tools` parameter. No
agent framework: reports of qwen3.5 emitting malformed tool calls trace to
wrapper layers rather than to the model, and a wrapper would also put a
translation step between the schema and what the model actually sees.

The model never receives a filesystem path. Path parameters are rewritten to
dataset LABELS with an enum of the registered ones, and the label is resolved
to a path only inside the executor. The model therefore cannot name a file it
was not given, cannot read one, and cannot leak one into its prose.
"""
import json
import os
import re
import time
import urllib.error
import urllib.request

# Phase 11: configurable, because the default assumes ollama runs on the same
# machine as the UI. SV_OLLAMA_URL lets it live elsewhere, and lets the
# unreachable-backend path be tested by pointing at a port nothing listens on.
OLLAMA = os.environ.get("SV_OLLAMA_URL", "http://127.0.0.1:11434").rstrip("/")

# tool parameter name -> ("label kind", replacement parameter name)
_PATH_PARAMS = {
    "bam_path":   ("dataset", "dataset"),
    "bam_paths":  ("dataset", "datasets"),
    "path":       ("candidates", "candidates"),
}


def _post(path, payload, timeout=600):
    req = urllib.request.Request(
        OLLAMA + path, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


_CAPS_CACHE = {}


def model_capabilities(model):
    """Capabilities Ollama reports for a model, e.g. {'completion','tools','thinking'}.
    Cached: /api/show is not free and the answer cannot change while loaded."""
    if model in _CAPS_CACHE:
        return _CAPS_CACHE[model]
    try:
        caps = set(_post("/api/show", {"model": model}, timeout=60).get("capabilities") or [])
    except Exception:
        caps = set()
    _CAPS_CACHE[model] = caps
    return caps


def supports_thinking(model):
    return "thinking" in model_capabilities(model)


def list_models():
    try:
        with urllib.request.urlopen(OLLAMA + "/api/tags", timeout=10) as r:
            return sorted(m["name"] for m in json.loads(r.read()).get("models", []))
    except Exception:
        return []


# ── description trimming (model path only) ──────────────────────────────────
# The 15 schemas cost 5,114 prompt tokens on EVERY call, which is what starved
# qwen3.5 at an 8k context. Trimming happens HERE and not in server.py, because
# the MCP descriptions have a second consumer: benchmark/mcp_client.py forwards
# tool.description to models, and six recorded run sets in benchmark/runs/ were
# produced against the current text. Editing the docstrings would silently make
# those non-comparable, and test_server.py additionally asserts no description
# is empty.
#
# Two things are never trimmed away:
#   * the leading summary, which is what a model needs to call the tool at all
#   * safety-relevant sentences deliberately placed in the description, above
#     all the split_reads min_mapq 0-versus-20 divergence, because that text is
#     what the model reads before quoting a split-read count
_SAFETY_KEEP = [
    re.compile(r"[^.]*min_mapq DEFAULTS TO 0[^.]*\.", re.I),
    re.compile(r"[^.]*(?:must not describe|cannot see this image|never describe)[^.]*\.", re.I),
    re.compile(r"[^.]*do NOT receive[^.]*\.", re.I),
    re.compile(r"[^.]*(?:QUALITY-LIMITED|withheld)[^.]*\.", re.I),
    re.compile(r"[^.]*not a calibrated[^.]*\.", re.I),
]


# MEASURED AND REJECTED AS A DEFAULT. Trimming cut the schema cost 5,114 -> 3,430
# tokens, but at matched settings (ctx=32768, think=ON, max_iters=40) it nearly
# doubled the schema-invalid argument rate, 16.4% -> 31.8%, and completion fell
# 16/18 -> 15/18. The tokens it buys are only worth having at a small context,
# and the model that needs a small context (qwen3.5:9b, 2,048 tokens GPU-resident)
# cannot hold the schemas even with every description deleted -- the floor is
# 2,690 tokens. So the trim has no beneficiary and is off by default. Kept
# because the measurement is the finding, and because a future tool surface with
# fewer tools would change the arithmetic.
def trim_description(desc, summary_chars=240):
    """Leading summary + any safety sentence, deduplicated, order preserved."""
    d = " ".join((desc or "").split())
    if not d:
        return "."
    head = d[:summary_chars]
    if len(d) > summary_chars:                    # do not cut mid-sentence
        cut = head.rfind(". ")
        head = head[:cut + 1] if cut > 60 else head.rstrip() + "..."
    keep = [head]
    for pat in _SAFETY_KEEP:
        for m in pat.finditer(d):
            snip = m.group(0).strip()
            if snip and snip not in " ".join(keep):
                keep.append(snip)
    return " ".join(keep)


def build_tools(mcp_tools, datasets, candidate_files, trim=False):
    """Ollama `tools` array built FROM list_tools(), with path parameters
    rewritten to label enums. Returns (tools, name->server map)."""
    out, where = [], {}
    for server_name, tools in mcp_tools.items():
        for t in tools:
            schema = json.loads(json.dumps(t.parameters or {"type": "object", "properties": {}}))
            props = schema.get("properties", {})
            required = list(schema.get("required", []))
            for pname, (kind, newname) in _PATH_PARAMS.items():
                if pname not in props:
                    continue
                enum = sorted(datasets) if kind == "dataset" else sorted(candidate_files)
                is_list = props[pname].get("type") == "array"
                props.pop(pname)
                required = [newname if r == pname else r for r in required]
                props[newname] = ({"type": "array", "items": {"type": "string", "enum": enum},
                                   "description": f"registered {kind} label(s)"}
                                  if is_list else
                                  {"type": "string", "enum": enum,
                                   "description": f"registered {kind} label"})
            if "mask_path" in props:      # a reference file, not something to name
                props.pop("mask_path")
                props["exclude_masked"] = {"type": "boolean",
                                           "description": "drop junctions in the caller's exclude regions"}
            schema["properties"] = props
            schema["required"] = [r for r in required if r in props]
            desc = (t.description or "").strip()[:900]
            if trim:
                desc = trim_description(desc)
            out.append({"type": "function", "function": {
                "name": t.name, "description": desc,
                "parameters": schema}})
            where[t.name] = server_name
    return out, where


def resolve_args(name, args, datasets, candidate_files, mask_path):
    """Label -> path, inside the executor only. Returns (args, error)."""
    a = dict(args or {})
    if "dataset" in a:
        lbl = a.pop("dataset")
        if lbl not in datasets:
            return None, f"unknown dataset label {lbl!r}; registered: {sorted(datasets)}"
        a["bam_path"] = datasets[lbl]
    if "datasets" in a:
        lbls = a.pop("datasets") or []
        if isinstance(lbls, str):
            lbls = [lbls]
        bad = [l for l in lbls if l not in datasets]
        if bad:
            return None, f"unknown dataset label(s) {bad}; registered: {sorted(datasets)}"
        a["bam_paths"] = [datasets[l] for l in lbls]
    if "candidates" in a:
        lbl = a.pop("candidates")
        if lbl not in candidate_files:
            return None, f"unknown candidates label {lbl!r}; registered: {sorted(candidate_files)}"
        a["path"] = candidate_files[lbl]
    if a.pop("exclude_masked", False):
        a["mask_path"] = mask_path
    return a, None


# ── tool-call-as-text detection ─────────────────────────────────────────────
# Ollama has an open issue where qwen3.5 sometimes PRINTS a tool call instead of
# emitting a structured one. That is a failure, not an answer, so it must be
# detected rather than rendered as prose.
_TEXT_CALL_PATTERNS = [
    re.compile(r"<tool_call>", re.I),
    re.compile(r"<\|tool_call\|>", re.I),
    re.compile(r'"(?:name|function)"\s*:\s*"[a-z_]+"\s*,\s*"(?:arguments|parameters)"\s*:', re.I),
    re.compile(r'\{\s*"(?:arguments|parameters)"\s*:\s*\{', re.I),
    re.compile(r"^\s*```(?:json|tool_code)?\s*\{\s*\"name\"", re.I | re.M),
]


def detect_text_tool_call(content, tool_names):
    """A tool call printed as prose. Returns a reason string, or None."""
    if not content:
        return None
    for pat in _TEXT_CALL_PATTERNS:
        m = pat.search(content)
        if m:
            return f"matched {pat.pattern[:44]!r} at offset {m.start()}"
    # a bare function-call-looking line naming a real tool
    for n in tool_names:
        if re.search(rf"(?<![\w.]){re.escape(n)}\s*\(\s*[a-z_]+\s*=", content):
            return f"prose contains a call-shaped invocation of {n}"
    return None


# ── number verification ─────────────────────────────────────────────────────
_NUM = re.compile(r"(?<![\w.\-])(\d+(?:,\d{3})*(?:\.\d+)?)(?![\w.]*[A-Za-z_])")

# Text that contains digits but asserts no measurement. Stripped before
# extraction, because a validator that manufactures violations is worse than
# none: it trains the reader to ignore it. Each pattern below was added after
# seeing it produce a false positive on real model output, not in advance.
_NOT_A_CLAIM = [
    re.compile(r"^\s{0,3}\d+[.)]\s", re.M),          # markdown list ordinals "1. " "2) "
    re.compile(r"\bchr[0-9XYM]+\b", re.I),            # chr20, chr21
    re.compile(r"\bchromosomes?\s+\d+(?:\s*/\s*\d+)*", re.I),   # "chromosome 20/21"
    re.compile(r"\bIMP\d+\b"),                       # dataset labels
    re.compile(r"\bq?\d+_K_M\b", re.I),              # quantisation tags
]


def _strip_non_claims(text):
    out = text or ""
    for pat in _NOT_A_CLAIM:
        out = pat.sub(" ", out)
    return out


def _collect(obj, nums, strings):
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        nums.add(float(obj))
        return
    if isinstance(obj, str):
        strings.append(obj)
        return
    if isinstance(obj, dict):
        for v in obj.values():
            _collect(v, nums, strings)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _collect(v, nums, strings)


def supported_values(results):
    """Every number a tool actually returned, plus the forms prose renders them
    in (percent, rounded), plus the literal text of every returned string."""
    nums, strings = set(), []
    for r in results:
        _collect(r, nums, strings)
    expanded = set()
    for n in nums:
        expanded.add(n)
        expanded.add(round(n))
        for d in (1, 2, 3):
            expanded.add(round(n, d))
        if 0.0 <= n <= 1.0:                    # fractions are quoted as percentages
            for d in (0, 1, 2):
                expanded.add(round(n * 100, d))
    return expanded, "\n".join(strings)


def verify_numbers(prose, results, tolerance=0.011):
    """Every number in the model's prose, checked against the trace."""
    supported, text = supported_values(results)
    found, unsupported = [], []
    scanned = _strip_non_claims(prose)
    for m in _NUM.finditer(scanned):
        raw = m.group(1)
        val = float(raw.replace(",", ""))
        ok = any(abs(val - s) <= tolerance for s in supported) or raw in text \
            or raw.replace(",", "") in text
        found.append({"text": raw, "value": val, "supported": ok})
        if not ok:
            unsupported.append(raw)
    return {"numbers_in_prose": len(found), "unsupported": unsupported,
            "no_tool_calls": not results,
            "details": found}


# ── what the model is shown of a tool return ────────────────────────────────
# Prose first. The measured failure pattern in this project is that the model
# quotes a tool's sentence far more readily than it recomputes from the fields,
# so the observation strings must survive truncation while bulk does not.
# breakpoint_evidence_summary nests full copies of the four layer dicts, which
# the model has already been shown as separate returns; those get dropped.
_BULKY = ("discordant_pairs", "soft_clips", "split_reads", "depth_profile", "locus_stats")


def shrink_for_model(name, result, budget=2600):
    if not isinstance(result, dict):
        return result
    out = dict(result)
    if name == "breakpoint_evidence_summary":
        for k in _BULKY:
            if isinstance(out.get(k), dict):
                out[k] = "<already returned as its own tool call>"
    if name == "read_depth_profile" and isinstance(out.get("windows"), list):
        out["windows"] = f"<{len(out['windows'])} bins omitted; see summary>"
    if name == "list_candidates" and isinstance(out.get("candidates"), list):
        out["candidates"] = out["candidates"][:12]
    s = json.dumps(out)
    if len(s) <= budget:
        return out
    keep = {k: v for k, v in out.items()
            if k in ("supporting_observations", "interpretation_template", "error",
                     "error_type", "reason", "assessable", "quality_limited",
                     "evidence_score", "evidence_strength", "signal_layers",
                     "min_mapq_applied", "position_provenance", "note",
                     # Phase 10: the attainable-ceiling echo. Kept for the same
                     # reason the observation prose is -- it is the part a model
                     # quotes -- and it is worthless if truncation eats it.
                     "attainable_note", "attainable_basis", "attainable_here",
                     "strong_band", "strong_band_reachable_here", "score_bands",
                     "max_all_layers", "max_with_flat_depth",
                     "attainable_ceiling_derivable", "attainable_ceiling_reason")}
    rest = {k: v for k, v in out.items() if k not in keep}
    keep["_truncated_fields"] = json.dumps(rest)[:budget // 2]
    return keep


SYSTEM_PROMPT = (
    "You inspect structural-variant breakpoints using the supplied tools.\n"
    "RULES:\n"
    "1. Every genomic fact you state must come from a tool call in THIS conversation. "
    "If no tool returned it, say you do not know and call the tool.\n"
    "2. Never state a number you did not receive from a tool.\n"
    "3. Do not use prior knowledge about genes, samples or variants.\n"
    "4. Datasets are referred to by LABEL. You do not have file paths and must not invent one.\n"
    "5. evidence_strength 'QUALITY-LIMITED' means the score was WITHHELD because too many reads "
    "are ambiguously mapped. It is NOT a low score. Report it as withheld.\n"
    "6. For a balanced translocation, flat read depth is EXPECTED and is not evidence against it.\n"
    "7. If a premise in the question is not supported by what the tools returned, say so plainly "
    "rather than accommodating it.\n"
    "8. If an image or panel could not be generated, say that. Never describe an image.\n"
)


# Ollama documents `think` as "a boolean or a thinking level ('low', 'medium',
# 'high', or 'max')". Coercing it with bool() silently turned every level into
# True, so a level could be requested and never reach the model. Validated here
# rather than coerced.
THINK_VALUES = (True, False, "low", "medium", "high", "max")


def normalise_think(v):
    if isinstance(v, str) and v.lower() in ("low", "medium", "high", "max"):
        return v.lower()
    if isinstance(v, str):
        return v.lower() not in ("0", "", "false", "no")
    return bool(v)


def run_turn(model, user_msg, tools, tool_names, exec_fn, num_ctx=8192,
             max_iters=8, system=SYSTEM_PROMPT, timeout=900, think=False):
    """One user turn. Returns the full trace: model messages, tool calls, and
    the verification pass over the final prose."""
    messages = [{"role": "system", "content": system},
                {"role": "user", "content": user_msg}]
    events, results_seen, call_ids = [], [], []
    text_call_failures, malformed = [], []
    think_supported = supports_thinking(model)
    think_requested = think
    if not think_supported:
        think = False
    produced_answer = False
    truncated = []
    done_reason = None
    last_prompt_tokens = 0
    t_start = time.time()
    stats = {"gen_tokens": 0, "gen_ns": 0, "prompt_tokens": 0}
    final_text = ""

    for it in range(max_iters):
        try:
            # `think` is sent explicitly. Ollama's docs state no default for it,
            # and measurement here showed it is effectively ON for qwen3.5: with
            # the 15 tool schemas costing ~5,110 prompt tokens, a thinking block
            # consumed the entire remaining generation budget and the model
            # returned empty content with done_reason='length'. Leaving it
            # unset silently starved the model.
            payload = {
                "model": model, "messages": messages, "tools": tools, "stream": False,
                "options": {"num_ctx": num_ctx}, "keep_alive": "600s"}
            # Sending `think` to a model without the capability is a hard 400 from
            # Ollama ("does not support thinking"), which killed 30 runs of
            # qwen2.5:7b outright and looked exactly like a model that answered
            # nothing. Ask what the model supports instead of assuming.
            if think_supported:
                payload["think"] = think
            resp = _post("/api/chat", payload, timeout=timeout)
        except Exception as e:
            events.append({"type": "error", "detail": f"ollama request failed: {e}"})
            break
        if isinstance(resp, dict) and resp.get("error"):
            events.append({"type": "error", "detail": str(resp["error"])})
            break
        msg = resp.get("message") or {}
        content = msg.get("content") or ""
        thinking = msg.get("thinking") or ""
        tcs = msg.get("tool_calls") or []
        stats["gen_tokens"] += resp.get("eval_count") or 0
        stats["gen_ns"] += resp.get("eval_duration") or 0
        stats["prompt_tokens"] += resp.get("prompt_eval_count") or 0
        done_reason = resp.get("done_reason")
        last_prompt_tokens = resp.get("prompt_eval_count") or 0
        if done_reason == "length":
            truncated.append({"iteration": it, "prompt_tokens": last_prompt_tokens,
                              "num_ctx": num_ctx,
                              "headroom": num_ctx - last_prompt_tokens})

        if not tcs:
            reason = detect_text_tool_call(content, tool_names)
            if reason:
                text_call_failures.append({"iteration": it, "reason": reason,
                                           "content": content[:600]})
        events.append({"type": "model", "iteration": it, "content": content,
                       "thinking": thinking[:1200], "n_tool_calls": len(tcs),
                       "done_reason": done_reason,
                       "prompt_tokens": last_prompt_tokens,
                       "headroom": num_ctx - last_prompt_tokens})
        messages.append({"role": "assistant", "content": content, "tool_calls": tcs})
        if not tcs:
            final_text = content
            produced_answer = True
            break

        for tc in tcs:
            fn = (tc or {}).get("function") or {}
            name = fn.get("name")
            args = fn.get("arguments")
            if isinstance(args, str):
                try:
                    args = json.loads(args)
                except Exception:
                    malformed.append({"name": name, "raw": args[:300],
                                      "reason": "arguments were a string that is not valid JSON"})
                    args = {}
            if not isinstance(args, dict):
                malformed.append({"name": name, "raw": str(args)[:300],
                                  "reason": f"arguments were {type(args).__name__}, not an object"})
                args = {}
            if name not in tool_names:
                malformed.append({"name": name, "raw": json.dumps(args)[:200],
                                  "reason": "named a tool that does not exist"})
                payload = {"error": f"no such tool {name!r}", "available": sorted(tool_names)}
                events.append({"type": "tool", "name": name, "params": args,
                               "call_id": None, "result": payload, "rejected": True})
                messages.append({"role": "tool", "tool_name": str(name),
                                 "content": json.dumps(payload)})
                continue
            rec, err = exec_fn(name, args)
            if err is not None:
                payload = {"error": err, "error_type": "bad_label"}
                events.append({"type": "tool", "name": name, "params": args,
                               "call_id": None, "result": payload, "rejected": True})
                messages.append({"role": "tool", "tool_name": name,
                                 "content": json.dumps(payload)})
                continue
            # FastMCP validates arguments against the generated schema. A
            # rejection here is an argument-FORMATTING failure by the model,
            # not a tool failure, and must be counted as such -- otherwise the
            # framework silently absorbs the model's mistakes and the
            # formatting-reliability number reads as perfect.
            res_err = ""
            if isinstance(rec.get("result"), dict):
                res_err = str(rec["result"].get("error", ""))
            if "validation error" in res_err.lower():
                malformed.append({"name": name, "raw": json.dumps(args)[:300],
                                  "reason": "arguments failed the tool's schema: "
                                            + res_err.split(chr(10))[0][:160]})
            call_ids.append(rec["id"])
            results_seen.append(rec["result"])
            events.append({"type": "tool", "name": name, "params": rec["params"],
                           "call_id": rec["id"], "result": rec["result"],
                           "is_error": rec["is_error"]})
            messages.append({"role": "tool", "tool_name": name,
                             "content": json.dumps(shrink_for_model(name, rec["result"]))})

    # A turn that never stopped requesting tools has no answer. Rendering the
    # empty string would look like a model that said nothing, which is a very
    # different claim from one that never finished.
    ended_without_answer = not produced_answer
    verification = verify_numbers(final_text, results_seen)
    gen_tps = (stats["gen_tokens"] / (stats["gen_ns"] / 1e9)) if stats["gen_ns"] else None
    return {
        "model": model, "num_ctx": num_ctx, "think": think,
        "think_requested": think_requested, "think_supported": think_supported,
        "max_iters": max_iters,
        "ended_without_answer": ended_without_answer,
        "context_truncated": truncated,
        "final_prompt_tokens": last_prompt_tokens,
        "final_headroom": num_ctx - last_prompt_tokens,
        "events": events, "final_text": final_text,
        "call_ids": call_ids, "n_tool_calls": len(call_ids),
        "text_tool_call_failures": text_call_failures,
        "malformed_tool_calls": malformed,
        "verification": verification,
        "wall_s": round(time.time() - t_start, 1),
        "gen_tokens": stats["gen_tokens"],
        "prompt_tokens": stats["prompt_tokens"],
        "gen_tokens_per_s": round(gen_tps, 1) if gen_tps else None,
        "iterations": len([e for e in events if e["type"] == "model"]),
    }


# ══════════════════════════════════════════════════════════════════════════════
# Anthropic API backend (Phase 9)
#
# A SECOND TRANSPORT, not a second path to the tools. Everything that decides
# what the model can see or say is shared with the Ollama path by construction:
# build_tools() produces the schemas, resolve_args() maps labels to paths inside
# the executor, the same exec_fn (ui._chat_exec) records and scrubs, and
# shrink_for_model / verify_numbers / detect_text_tool_call / SYSTEM_PROMPT are
# imported from above rather than restated. The only thing written twice is the
# request/response shape, because the two APIs disagree about it -- and that
# translation is one function, to_anthropic_tools(), which moves `parameters` to
# `input_schema` and unwraps the "function" envelope. Nothing else differs.
#
# This is the first phase in which data leaves the machine, so the leak surface
# is the WHOLE serialised request body, not just the tool returns: schemas,
# system prompt, user message, and every tool result. test_api_leak.py builds a
# real payload and asserts no registered path or basename occurs anywhere in it,
# with a negative control that creates the leak to prove the check can fail.
# ══════════════════════════════════════════════════════════════════════════════

# $ per 1M tokens, first-party Anthropic API rates. Sonnet 5's introductory
# $2/$10 ran through 2026-08-31 and has expired, so standard rates apply.
API_PRICES = {
    "claude-opus-5":   {"in": 5.00, "out": 25.00},
    "claude-sonnet-5": {"in": 3.00, "out": 15.00},
    "claude-opus-4-8": {"in": 5.00, "out": 25.00},
    "claude-haiku-4-5": {"in": 1.00, "out": 5.00},
}
CACHE_WRITE_MULT = 1.25    # 5-minute TTL
CACHE_READ_MULT = 0.10


def api_cost(model, usage):
    """USD for one turn from accumulated usage. Unknown model -> None, never a
    guessed price."""
    p = API_PRICES.get(model)
    if not p:
        return None
    return round(
        (usage["input_tokens"] * p["in"]
         + usage["cache_creation_input_tokens"] * p["in"] * CACHE_WRITE_MULT
         + usage["cache_read_input_tokens"] * p["in"] * CACHE_READ_MULT
         + usage["output_tokens"] * p["out"]) / 1e6, 6)


def to_anthropic_tools(tools, cache_last=False):
    """The ONE legitimate difference between the two paths: Ollama nests the
    schema under function.parameters, Anthropic takes it flat as input_schema.
    Same names, same descriptions, same schemas -- including the label enums,
    so the API model is offered exactly the parameter space the local one was."""
    out = []
    for t in tools:
        f = t["function"]
        out.append({"name": f["name"], "description": f["description"],
                    "input_schema": f["parameters"]})
    if cache_last and out:
        out[-1] = dict(out[-1], cache_control={"type": "ephemeral"})
    return out


def _blocks_to_plain(content):
    """SDK content blocks -> JSON-serialisable dicts, for the leak test and the
    stored trace. Never used to build a request: assistant turns are echoed back
    as the SDK objects so thinking blocks survive verbatim."""
    out = []
    for b in content:
        d = getattr(b, "model_dump", None)
        out.append(d(mode="json") if d else dict(b))
    return out


def _clear_cache_marks(messages):
    """Drop cache_control from every block we previously marked. Only plain
    dict blocks are touched: assistant turns hold SDK objects and are echoed
    back verbatim so thinking blocks survive."""
    for m in messages:
        c = m.get("content")
        if not isinstance(c, list):
            continue
        for i, b in enumerate(c):
            if isinstance(b, dict) and "cache_control" in b:
                c[i] = {k: v for k, v in b.items() if k != "cache_control"}


def run_turn_api(model, user_msg, tools, tool_names, exec_fn,
                 max_iters=40, system=SYSTEM_PROMPT, max_tokens=16000,
                 effort="high", thinking=True, api_key=None, timeout=1800,
                 client=None, dry_run=False):
    """One user turn against the Anthropic API. Returns the same trace shape as
    run_turn() so the Phase 8 signals() reads it unchanged, plus api_usage.

    dry_run=True builds and returns the first request body WITHOUT sending it.
    That is what the leak test inspects: the actual bytes, not a reconstruction.
    """
    import anthropic

    sys_blocks = [{"type": "text", "text": system,
                   "cache_control": {"type": "ephemeral"}}]
    api_tools = to_anthropic_tools(tools, cache_last=False)
    messages = [{"role": "user", "content": user_msg}]

    def build(msgs):
        req = {"model": model, "max_tokens": max_tokens,
               "system": sys_blocks, "tools": api_tools, "messages": msgs}
        if thinking:
            # Adaptive is the default on both control models; stated explicitly
            # so the record is unambiguous. display=summarized costs nothing --
            # thinking is billed identically under every display setting -- and
            # makes the reasoning readable for the ceiling question.
            req["thinking"] = {"type": "adaptive", "display": "summarized"}
        if effort:
            req["output_config"] = {"effort": effort}
        return req

    if dry_run:
        return build(messages)

    if client is None:
        client = anthropic.Anthropic(api_key=api_key, timeout=timeout, max_retries=2)

    events, results_seen, call_ids = [], [], []
    text_call_failures, malformed = [], []
    truncated, refusals = [], []
    usage = {"input_tokens": 0, "output_tokens": 0,
             "cache_creation_input_tokens": 0, "cache_read_input_tokens": 0}
    produced_answer = False
    final_text = ""
    stop_reason = None
    t_start = time.time()
    api_error = None

    for it in range(max_iters):
        try:
            resp = client.messages.create(**build(messages))
        except Exception as e:
            api_error = f"{type(e).__name__}: {e}"
            events.append({"type": "error", "detail": api_error})
            break

        u = resp.usage
        for k in usage:
            usage[k] += getattr(u, k, 0) or 0
        stop_reason = resp.stop_reason

        # A refusal is HTTP 200 with an empty or partial content array. Reading
        # content[0] unconditionally would crash here; treating it as an answer
        # would be worse. It is recorded as its own outcome. Server-side
        # fallbacks are deliberately NOT enabled: a silent switch to another
        # model would attribute that model's answer to this one.
        if stop_reason == "refusal":
            det = getattr(resp, "stop_details", None)
            refusals.append({"iteration": it,
                             "category": getattr(det, "category", None),
                             "explanation": getattr(det, "explanation", None)})
            events.append({"type": "refusal", "iteration": it,
                           "category": getattr(det, "category", None)})
            break
        if stop_reason == "max_tokens":
            truncated.append({"iteration": it, "max_tokens": max_tokens,
                              "output_tokens": u.output_tokens})

        content = resp.content
        text = "".join(b.text for b in content if b.type == "text")
        think_txt = "".join(getattr(b, "thinking", "") or ""
                            for b in content if b.type == "thinking")
        tcs = [b for b in content if b.type == "tool_use"]

        if not tcs:
            reason = detect_text_tool_call(text, tool_names)
            if reason:
                text_call_failures.append({"iteration": it, "reason": reason,
                                           "content": text[:600]})
        events.append({"type": "model", "iteration": it, "content": text,
                       "thinking": think_txt[:1200], "n_tool_calls": len(tcs),
                       "done_reason": stop_reason,
                       "prompt_tokens": (u.input_tokens or 0)
                                        + (u.cache_read_input_tokens or 0)
                                        + (u.cache_creation_input_tokens or 0),
                       "headroom": None})

        messages.append({"role": "assistant", "content": content})

        if not tcs:
            final_text = text
            produced_answer = True
            break

        # Every tool_use block must get exactly one tool_result, and they all go
        # back in ONE user message. Splitting them trains the model out of
        # parallel calls; omitting one is a hard 400.
        blocks = []
        for tc in tcs:
            name, args = tc.name, tc.input
            if not isinstance(args, dict):
                malformed.append({"name": name, "raw": str(args)[:300],
                                  "reason": f"arguments were {type(args).__name__}, not an object"})
                args = {}
            if name not in tool_names:
                malformed.append({"name": name, "raw": json.dumps(args)[:200],
                                  "reason": "named a tool that does not exist"})
                payload = {"error": f"no such tool {name!r}", "available": sorted(tool_names)}
                events.append({"type": "tool", "name": name, "params": args,
                               "call_id": None, "result": payload, "rejected": True})
                blocks.append({"type": "tool_result", "tool_use_id": tc.id,
                               "content": json.dumps(payload), "is_error": True})
                continue
            rec, err = exec_fn(name, args)
            if err is not None:
                payload = {"error": err, "error_type": "bad_label"}
                events.append({"type": "tool", "name": name, "params": args,
                               "call_id": None, "result": payload, "rejected": True})
                blocks.append({"type": "tool_result", "tool_use_id": tc.id,
                               "content": json.dumps(payload), "is_error": True})
                continue
            res_err = ""
            if isinstance(rec.get("result"), dict):
                res_err = str(rec["result"].get("error", ""))
            if "validation error" in res_err.lower():
                malformed.append({"name": name, "raw": json.dumps(args)[:300],
                                  "reason": "arguments failed the tool's schema: "
                                            + res_err.split(chr(10))[0][:160]})
            call_ids.append(rec["id"])
            results_seen.append(rec["result"])
            events.append({"type": "tool", "name": name, "params": rec["params"],
                           "call_id": rec["id"], "result": rec["result"],
                           "is_error": rec["is_error"]})
            blocks.append({"type": "tool_result", "tool_use_id": tc.id,
                           "content": json.dumps(shrink_for_model(name, rec["result"]))})
        # Incremental cache breakpoint on the growing history: the loop resends
        # everything every iteration, and Phase 8c saw turns of 22 calls.
        #
        # The breakpoint MOVES; it does not accumulate. Marking each new tool
        # result while leaving the previous marks in place hit the hard limit of
        # 4 cache_control blocks per request and 400'd the turn at iteration 4
        # -- which surfaced as a model that "ended without an answer", i.e. a
        # harness failure wearing a model failure's clothes. One permanent mark
        # on the system block plus one rolling mark here is two, always.
        _clear_cache_marks(messages)
        if blocks:
            blocks[-1] = dict(blocks[-1], cache_control={"type": "ephemeral"})
        messages.append({"role": "user", "content": blocks})

    verification = verify_numbers(final_text, results_seen)
    return {
        "model": model, "backend": "anthropic",
        "num_ctx": None, "think": effort if thinking else False,
        "think_requested": effort if thinking else False,
        "think_supported": True, "max_iters": max_iters,
        "ended_without_answer": not produced_answer,
        "context_truncated": truncated,
        "refusals": refusals,
        "api_error": api_error,
        "stop_reason": stop_reason,
        "final_prompt_tokens": None, "final_headroom": None,
        "events": events, "final_text": final_text,
        "call_ids": call_ids, "n_tool_calls": len(call_ids),
        "text_tool_call_failures": text_call_failures,
        "malformed_tool_calls": malformed,
        "verification": verification,
        "wall_s": round(time.time() - t_start, 1),
        "gen_tokens": usage["output_tokens"],
        "prompt_tokens": usage["input_tokens"] + usage["cache_read_input_tokens"]
                         + usage["cache_creation_input_tokens"],
        "gen_tokens_per_s": None,
        "iterations": len([e for e in events if e["type"] == "model"]),
        "api_usage": usage,
        "api_cost_usd": api_cost(model, usage),
    }
