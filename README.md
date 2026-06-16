# junction-ddpcr-designer

Design a **competitive ddPCR assay across a splice / trans-splice junction** — straight from two sequences.

You give it two junction sequences that differ only at a few discriminating bases:

- **WT** — the reference junction (e.g. native `exon-A | exon-B`)
- **TS** — the variant junction (e.g. `exon-A | recoded cargo`, or a SNP / edited allele)

…plus the index where the two exons meet. It **derives** (nothing hardcoded):

- a **competitive probe pair** (WT + TS), designed *independently* — different lengths and shifted/nested windows are allowed — so their **matched Tms land within 2 °C of each other** while each strongly rejects the opposite junction;
- a **shared forward + reverse primer pair** sitting in the regions identical to both junctions (one pair amplifies both products, and — placed across the intron — only spliced cDNA);
- an **upstream reference probe** inside the amplicon that fires on every molecule (total = WT + TS);
- a full **audit** of every oligo against ddPCR design rules (Bio-Rad Bulletin 6407 + IDT) and an HTML report.

It emits **bare sequences** — add FAM / HEX / Cy5 + ZEN/IBFQ at order time.

## Install

```bash
pip install -r requirements.txt    # biopython, primer3-py
```

## Use

```bash
# built-in worked example (human HEXA exon9 | exon10 vs trans-spliced cargo)
python junction_ddpcr_designer.py --example --out report.html

# your own junction
python junction_ddpcr_designer.py \
  --wt  GAAGTCC...GCTGCTGGAC...AAAG \
  --ts  GAAGTCC...GCTCCTAGAT...AAAG \
  --junction 87 \
  --intron 289 \
  --title "MYGENE exonA|exonB" \
  --out report.html

# or from a config file
python junction_ddpcr_designer.py --config myjunction.json --out report.html
```

`myjunction.json`:

```json
{ "wt": "....", "ts": "....", "junction": 87, "intron": 289, "title": "MYGENE exonA|exonB" }
```

Key options: `--probe-primer-gap` (default **6** — primers are designed this many °C *below*
the probes, per Bio-Rad's probe-over-primer rule; the run anneal is set near the primer Tm),
`--match-window` (2.0 °C), `--idt-offset` (7.5), `--intron` (enables the gDNA-by-length check).

> **Run the anneal near the primer Tm**, not at a generic 60 °C. With `--probe-primer-gap 6`
> the primers come out ~6 °C below the probes (e.g. ~55 °C); run them at a ~55 °C anneal (use the
> 55–65 °C gradient — the optimum sits at the low end). Low-Tm primers run at 60 °C go silent on cDNA.

## How it chooses (the principles)

**No temperature is preset — every Tm is derived from your sequences.** Only *constraints* are fixed
(the gap, the matched-window, a discrimination floor, a minimum anneal, the length cap). The actual
probe/primer/anneal temperatures emerge from what the junction allows: an AT-rich junction lands in
the high 50s, a GC-rich one higher.

**Constraints (fixed):** each probe crosses the junction and carries all discriminating bases; matched
Tms within `--match-window` (2 °C); each probe rejects the opposite junction by ≥ `--min-discrimination`
(10 °C); probe < 30 nt; 5′ ≠ G; more C than G; no runs > 3; no self-folding at the (derived) anneal.

**Objective (derives the temperature):** make the binding-limiting (colder) probe **as warm as the
sequences allow** for good occupancy → then a tight matched-Tm gap → more discrimination → concentric
nested windows (shorter probe nested in the longer, same centre → both compete at one site) → shorter
probes. The chosen probe Tm sets everything else: **primers = probes − `--probe-primer-gap` (6 °C),
and the run anneal = the primer Tm.**

## Tm model

- **Probes:** Biopython nearest-neighbor `Tm_NN`, **mismatch-aware** via `c_seq` — the only way to score a probe on the *opposite* junction (the discrimination ΔTm).
- **Primers:** primer3 `calc_tm`.
- Both reported on the **IDT OligoAnalyzer scale** by subtracting a calibration offset (default **7.5 °C**, fitted to the validated MECP2 RMD assay: probes 62.8 / 64.2, primers 61). Always confirm final Tm (and any LNA) in IDT OligoAnalyzer before ordering.

## What it does *not* do

- No dye/quencher chemistry is added — the output is bare DNA.
- No **BLAST specificity** check against paralogs / pseudogenes — run that separately before ordering.
- LNA balancing of an unusually cold probe is left to OligoAnalyzer (the offset estimate isn't precise enough to place LNA).

## Background

This grew out of designing a competitive RMD ddPCR assay to quantify 3′-exon-replacement trans-splicing — measuring trans-spliced (FAM) vs wild-type (HEX) junctions in one well, normalized to a total-transcript reference. The architecture mirrors Bio-Rad's Rare Mutation Detection (one shared primer pair, two competitive probes). See the built-in `--example`.

## License

MIT © 2026 Max Wraae
