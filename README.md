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

1. **Each probe crosses the junction and carries all discriminating bases** — so it reports the spliced product and discriminates.
2. **Matched Tms within 2 °C of each other, at or above the primers.** The cold (AT-rich) probe is allowed down to primer level; the anneal gradient sets the working temp.
3. **Rank by the *weaker* probe's discrimination first** (so neither probe is the weak link), then prefer **concentric windows** (the shorter probe nested in the longer, same centre → both compete at one site), then the smaller Tm gap, then shorter probes.

Probes shorter than 30 nt (quencher distance); 5′ ≠ G; more C than G; no runs > 3; hairpin Tm below the anneal temp.

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
