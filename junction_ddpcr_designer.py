#!/usr/bin/env python3
"""
junction-ddpcr-designer
=======================
Design a competitive ddPCR (RMD-style) assay across a splice / trans-splice junction.

You give it two junction sequences that differ only at a few discriminating bases:
  * WT  — the reference junction (e.g. native exon-A | exon-B)
  * TS  — the variant junction   (e.g. exon-A | recoded cargo, or a SNP/edit)
plus the index where the two exons meet. It DERIVES (nothing hardcoded):

  * a competitive probe PAIR (WT + TS) — designed independently (different lengths /
    shifted, nested windows allowed) so their MATCHED Tms land within 2 C of each other
    while each strongly rejects the opposite junction;
  * a SHARED forward + reverse primer pair sitting in the regions identical to both
    junctions (one pair amplifies both products);
  * an upstream REFERENCE probe inside the amplicon (counts every molecule = total);

and audits every oligo against ddPCR design rules (Bio-Rad Bulletin 6407 + IDT).

Tm model
--------
Probes: Biopython nearest-neighbor `Tm_NN` (mismatch-aware via `c_seq`) — this is the only
way to score a probe on the *opposite* junction. Primers: primer3 `calc_tm`. Both are reported
on the IDT-OligoAnalyzer scale by subtracting a calibration offset (default 7.5 C, fitted to the
validated MECP2 RMD assay: probes 62.8 / 64.2, primers 61).

No dye/quencher labels are emitted — add FAM / HEX / Cy5 + ZEN/IBFQ at order time, and confirm
final Tm/LNA in IDT OligoAnalyzer. BLAST specificity (paralogs, pseudogenes) is a separate step.

Requires: biopython, primer3-py  (`pip install -r requirements.txt`).
"""
import argparse, json, sys, warnings
warnings.filterwarnings("ignore")
from Bio.SeqUtils import MeltingTemp as mt
import primer3

# ----- HEXA exon9|exon10 vs trans-spliced-cargo: the built-in worked example -----
HEXA_WT = ("GAAGTCCAACCCAGAGATCCAGGACTTTATGAGGAAGAAAGGCTTCGGTGAGGACTTCAAGCAGCTGGAGTCCTTCTACATCCAGAC"
           "GCTGCTGGACATCGTCTCTTCTTATGGCAAGGGCTATGTGGTGTGGCAGGAGGTGTTTGATAATAAAGTAAAG")
HEXA_TS = ("GAAGTCCAACCCAGAGATCCAGGACTTTATGAGGAAGAAAGGCTTCGGTGAGGACTTCAAGCAGCTGGAGTCCTTCTACATCCAGAC"
           "GCTCCTAGATATCGTCTCTTCTTATGGCAAGGGCTATGTGGTGTGGCAGGAGGTGTTTGATAATAAAGTAAAG")
HEXA_JUNCTION = 87
HEXA_INTRON = 289   # intron between the two exons (bp); gDNA amplicon = cDNA amplicon + intron

_C = str.maketrans("ACGTacgt", "TGCAtgca")
comp = lambda s: s.translate(_C)
rc   = lambda s: s.translate(_C)[::-1]
gc   = lambda s: round(100*(s.count("G")+s.count("C"))/len(s), 1) if s else 0.0
def maxrun(s):
    b=c=1
    for i in range(1,len(s)): c=c+1 if s[i]==s[i-1] else 1; b=max(b,c)
    return b
def dinuc(s):
    best=0
    for i in range(len(s)-1):
        u=s[i:i+2]; n=1; j=i+2
        while s[j:j+2]==u: n+=1; j+=2
        best=max(best,n)
    return best


class Designer:
    def __init__(self, wt, ts, junction, **opt):
        self.WT, self.TS, self.J = wt.upper(), ts.upper(), junction
        if len(self.WT) != len(self.TS):
            sys.exit("ERROR: WT and TS must be the same length / coordinate frame.")
        self.MM = [i for i in range(len(self.WT)) if self.WT[i] != self.TS[i]]
        if not self.MM:
            sys.exit("ERROR: WT and TS are identical — no discriminating bases.")
        self.up_end, self.down_beg = min(self.MM), max(self.MM)+1
        # tunables
        self.PRIMER_TM   = opt.get("primer_tm", 61.0)
        self.ANNEAL      = opt.get("anneal", 60.0)
        self.IDT_OFF     = opt.get("idt_offset", 7.5)
        self.PROBE_MIN, self.PROBE_MAX   = opt.get("probe_len", (18, 29))
        self.PRIMER_MIN, self.PRIMER_MAX = opt.get("primer_len", (18, 28))
        self.MATCH_WINDOW = opt.get("match_window", 2.0)
        self.OVER_PRIMER  = opt.get("over_primer", (0.0, 11.0))
        self.AMPLICON     = opt.get("amplicon", (60, 200))
        self.INTRON       = opt.get("intron", None)
        self.SALT = dict(Na=50, Mg=3, dNTPs=0.8, dnac1=250, dnac2=0, saltcorr=7)

    # ---- Tm helpers (IDT scale) ----
    def ptm(self, probe, target): return round(mt.Tm_NN(probe, c_seq=comp(target), **self.SALT) - self.IDT_OFF, 1)
    def primer_tm(self, s):       return round(primer3.calc_tm(s, mv_conc=50, dv_conc=3.8, dntp_conc=0.8, dna_conc=250) - self.IDT_OFF, 1)
    def hp_tm(self, s):           return round(primer3.calc_hairpin(s).tm, 1)
    def hairpin(self, s):         return round(primer3.calc_hairpin(s).dg/1000, 1)
    def selfdimer(self, s):       return round(primer3.calc_homodimer(s).dg/1000, 1)
    def hetero(self, a, b):       return round(primer3.calc_heterodimer(a, b).dg/1000, 1)
    def end_stab(self, p):        return round(primer3.calc_heterodimer(p[-5:], rc(p[-5:])).dg/1000, 1)

    # ---- STEP 1: probe pair ----
    def _windows(self, seq):
        out=[]
        for L in range(self.PROBE_MIN, self.PROBE_MAX+1):
            for a in range(0, len(seq)-L+1):
                b=a+L
                if not (a < self.J < b): continue            # must cross the junction
                if not all(a<=m<b for m in self.MM): continue # carry all discriminating bases
                s=seq[a:b]
                if s[0]=="G": continue                       # no 5' G
                out.append((a,b,L,s))
        return out

    def design(self):
        lo, hi = self.PRIMER_TM+self.OVER_PRIMER[0], self.PRIMER_TM+self.OVER_PRIMER[1]
        WTw = [(a,b,L,s,self.ptm(s,s),self.ptm(s,self.TS[a:b])) for a,b,L,s in self._windows(self.WT)]
        TSw = [(a,b,L,s,self.ptm(s,s),self.ptm(s,self.WT[a:b])) for a,b,L,s in self._windows(self.TS)]
        WTw = [w for w in WTw if lo<=w[4]<=hi and self.hp_tm(w[3])<self.ANNEAL]
        TSw = [w for w in TSw if lo<=w[4]<=hi and self.hp_tm(w[3])<self.ANNEAL]
        if not WTw or not TSw:
            sys.exit("ERROR: no probe windows satisfy the Tm / hairpin constraints — loosen --over-primer or --probe-len.")
        pairs=[]
        for (aw,bw,Lw,wt,WTm,WTo) in WTw:
            for (at,bt,Lt,ts,TSm,TSo) in TSw:
                gap=round(abs(WTm-TSm),1)
                if gap>self.MATCH_WINDOW: continue
                pairs.append(dict(aw=aw,bw=bw,Lw=Lw,wt=wt,WTm=WTm,WTo=WTo,
                                  at=at,bt=bt,Lt=Lt,ts=ts,TSm=TSm,TSo=TSo,gap=gap,
                                  minD=round(min(WTm-WTo, TSm-TSo),1)))
        if not pairs:
            sys.exit(f"ERROR: no probe pair within {self.MATCH_WINDOW} C — raise --match-window.")
        # Ranking principles: (1) strongest weaker-probe discrimination; (2) windows concentric
        # (shorter nested in longer, same centre -> compete at one site); (3) smaller gap; (4) shorter.
        def nested(d):
            return (d["at"]<=d["aw"] and d["bw"]<=d["bt"]) if d["Lw"]<=d["Lt"] else (d["aw"]<=d["at"] and d["bt"]<=d["bw"])
        def concentric(d): return abs((d["aw"]+d["bw"])-(d["at"]+d["bt"]))/2
        pairs.sort(key=lambda d:(-d["minD"], 0 if nested(d) else 1, concentric(d), d["gap"], d["Lw"]+d["Lt"]))
        self.pairs=pairs; self.best=best=pairs[0]

        # ---- STEP 2: shared primers around the probes ----
        p5, p3 = min(best["aw"],best["at"]), max(best["bw"],best["bt"])
        self.FWD=(self._scan_primer(self.WT[:self.up_end],0,False,end_before=p5) or [None])[0]
        revs=self._scan_primer(self.WT[self.down_beg:],self.down_beg,True,start_after=p3)
        self.REV=None
        if self.FWD and revs:
            for d in revs:
                amp=d["ge"]-self.FWD["gs"]
                if self.AMPLICON[0]<=amp<=self.AMPLICON[1] and abs(d["tm"]-self.FWD["tm"])<=2.5:
                    self.REV=dict(d,amp=amp); break
            self.REV=self.REV or dict(revs[0],amp=revs[0]["ge"]-self.FWD["gs"])
        self.prim=round((self.FWD["tm"]+self.REV["tm"])/2,1) if (self.FWD and self.REV) else self.PRIMER_TM

        # ---- STEP 3: upstream reference probe (inside amplicon, before the junction probes) ----
        ref_lo=(self.FWD["ge"]+1) if self.FWD else 1
        self.REF=(self._scan_ref(self.WT[ref_lo:p5-1], ref_lo) or [None])[0]
        return best

    def _scan_primer(self, region, off, want_rc, end_before=None, start_after=None):
        out=[]
        for L in range(self.PRIMER_MIN, self.PRIMER_MAX+1):
            for i in range(0, len(region)-L+1):
                gs,ge=off+i,off+i+L
                if end_before is not None and ge>end_before: continue
                if start_after is not None and gs<start_after: continue
                p = rc(region[i:i+L]) if want_rc else region[i:i+L]
                if p[-1] not in "GC" or maxrun(p)>3: continue
                if sum(x in "GC" for x in p[-5:])>3: continue
                if not (40<=gc(p)<=60): continue
                out.append(dict(seq=p,gs=gs,ge=ge,L=L,tm=self.primer_tm(p),gc=gc(p)))
        out.sort(key=lambda d:(abs(d["tm"]-self.PRIMER_TM), abs(d["gc"]-53)))
        return out

    def _scan_ref(self, region, off, t_lo=60, t_hi=67):
        out=[]
        for L in range(20,27):
            for i in range(0,len(region)-L+1):
                s=region[i:i+L]
                s2 = rc(s) if (s[0]=="G" or s.count("C")<s.count("G")) else s
                if s2[0]=="G" or maxrun(s2)>3 or self.hp_tm(s2)>=self.ANNEAL: continue
                t=self.ptm(s2,s2)
                if t_lo<=t<=t_hi and 30<=gc(s2)<=80:
                    out.append((round(abs(t-63),2), s2, t, gc(s2), off+i, off+i+L))
        out.sort()
        return out

    # ---- audit ----
    def primer_checks(self, p):
        s=p["seq"]; last5=sum(x in "GC" for x in s[-5:])
        return [
            ("length 18-28", 18<=len(s)<=28, f"{len(s)} nt", False),
            ("GC 50-60%", 50<=p["gc"]<=60, f"{p['gc']}%", not(50<=p['gc']<=60) and 40<=p['gc']<=65),
            (f"Tm ~ anneal {self.ANNEAL:g}C", self.ANNEAL-2<=p["tm"]<=self.ANNEAL+3, f"{p['tm']}C", False),
            ("3' G/C clamp", s[-1] in "GC", s[-1], False),
            ("<=3 GC in last 5", last5<=3, f"{last5}/5", False),
            ("3' not over-stable", self.end_stab(s)>-9, f"{self.end_stab(s)} kcal", False),
            ("no run > 3", maxrun(s)<=3, f"run {maxrun(s)}", False),
            ("no dinuc repeat >=4", dinuc(s)<4, f"x{dinuc(s)}", False),
            ("hairpin Tm < anneal", self.hp_tm(s)<self.ANNEAL, f"{self.hp_tm(s)}C", self.ANNEAL-12<=self.hp_tm(s)<self.ANNEAL),
            ("self-dimer dG > -6", self.selfdimer(s)>-6, f"{self.selfdimer(s)} kcal", False),
        ]
    def probe_checks(self, s, tm, off=None):
        above=round(tm-self.prim,1)
        L=[
            ("length < 30", len(s)<30, f"{len(s)} nt", False),
            ("GC 30-80%", 30<=gc(s)<=80, f"{gc(s)}%", False),
            ("no 5' G", s[0]!="G", s[0], False),
            ("more C than G", s.count("C")>=s.count("G"), f"{s.count('C')}C/{s.count('G')}G", False),
            ("no run > 3", maxrun(s)<=3, f"run {maxrun(s)}", False),
            ("no dinuc repeat >=4", dinuc(s)<4, f"x{dinuc(s)}", False),
            ("hairpin Tm < anneal", self.hp_tm(s)<self.ANNEAL, f"{self.hp_tm(s)}C", self.ANNEAL-12<=self.hp_tm(s)<self.ANNEAL),
            ("Tm above primers (3-10 ideal)", above>=0, f"+{above}C", 0<=above<3),
        ]
        if off is not None:
            d=round(tm-off,1); L.append(("discrimination dTm >= 8", d>=8, f"{d}C", 8<=d<10))
        return L
    def audit(self):
        b=self.best; amp_seq=self.WT[self.FWD["gs"]:self.REV["ge"]] if (self.FWD and self.REV) else ""
        gdna=(self.REV["amp"]+self.INTRON) if (self.REV and self.INTRON) else None
        rows=[("Forward primer", self.FWD["seq"] if self.FWD else "-", self.primer_checks(self.FWD) if self.FWD else []),
              ("Reverse primer", self.REV["seq"] if self.REV else "-", self.primer_checks(self.REV) if self.REV else []),
              ("WT probe", b["wt"], self.probe_checks(b["wt"], b["WTm"], b["WTo"])),
              ("TS probe", b["ts"], self.probe_checks(b["ts"], b["TSm"], b["TSo"]))]
        if self.REF: rows.append(("Reference probe", self.REF[1], self.probe_checks(self.REF[1], self.REF[2])))
        assay=[
            ("amplicon 60-200 bp", bool(self.REV) and self.AMPLICON[0]<=self.REV["amp"]<=self.AMPLICON[1], f"{self.REV['amp'] if self.REV else '-'} bp", False),
            ("amplicon GC 40-60%", bool(amp_seq) and 40<=gc(amp_seq)<=60, f"{gc(amp_seq) if amp_seq else '-'}%", False),
            ("primer pair Tm within 2C", bool(self.FWD and self.REV) and abs(self.FWD["tm"]-self.REV["tm"])<=2, f"{round(abs(self.FWD['tm']-self.REV['tm']),1) if (self.FWD and self.REV) else '-'}C", False),
            ("FWD/REV hetero-dimer > -6", bool(self.FWD and self.REV) and self.hetero(self.FWD["seq"],self.REV["seq"])>-6, f"{self.hetero(self.FWD['seq'],self.REV['seq']) if (self.FWD and self.REV) else '-'} kcal", False),
            ("matched probe Tms within 2C", b["gap"]<=2, f"{b['gap']}C", False),
        ]
        if gdna is not None:
            assay.append(("gDNA-safe by length", gdna>1000, f"{gdna} bp (+intron) - DNase mandatory", gdna<=1000))
        nf=sum(1 for _,_,cs in rows for _,ok,_,w in cs if not ok and not w)+sum(1 for _,ok,_,w in assay if not ok and not w)
        nw=sum(1 for _,_,cs in rows for _,ok,_,w in cs if w)+sum(1 for _,ok,_,w in assay if w)
        return rows, assay, nf, nw


def render_text(d, rows, assay, nf, nw):
    b=d.best; SY=lambda ok,w:("OK " if not w else "~? ") if ok else ("~X " if w else "XX ")
    print("FINAL OLIGOS (bare; add FAM/HEX/Cy5 + ZEN/IBFQ at order)\n"+"-"*60)
    print(f"WT probe        {b['wt']}")
    print(f"TS probe        {b['ts']}")
    print(f"FORWARD         {d.FWD['seq'] if d.FWD else '-'}")
    print(f"REVERSE         {d.REV['seq'] if d.REV else '-'}")
    print(f"REFERENCE probe {d.REF[1] if d.REF else '- (none found)'}")
    print("-"*60)
    print(f"derived: WT {b['Lw']}nt Tm {b['WTm']} d{round(b['WTm']-b['WTo'],1)} | "
          f"TS {b['Lt']}nt Tm {b['TSm']} d{round(b['TSm']-b['TSo'],1)} | "
          f"matched-Tm gap {b['gap']} | primers {d.prim} | {len(d.pairs)} valid pairs")
    print("\nCHECKS  (OK pass / ~? flagged / ~X soft-fail / XX fail)\n"+"-"*60)
    for nm,seq,cs in rows:
        print(f"\n{nm}  {seq}")
        for lab,ok,val,w in cs: print(f"   {SY(ok,w)}{lab:<32}{val}")
    print("\nASSAY-LEVEL")
    for lab,ok,val,w in assay: print(f"   {SY(ok,w)}{lab:<32}{val}")
    print(f"\n{'ALL CHECKS PASS' if nf==0 else str(nf)+' FAILED'}{f'  ({nw} flagged)' if nw else ''}")


def render_html(d, rows, assay, nf, nw, title):
    b=d.best
    oligos=[("WT probe","HEX",b['wt'],b['WTm']),("TS probe","FAM",b['ts'],b['TSm']),
            ("Forward primer","",d.FWD['seq'] if d.FWD else "-",d.FWD['tm'] if d.FWD else "-"),
            ("Reverse primer","",d.REV['seq'] if d.REV else "-",d.REV['tm'] if d.REV else "-"),
            ("Reference probe","Cy5",d.REF[1] if d.REF else "-",d.REF[2] if d.REF else "-")]
    trows="".join(f"<tr><td>{n}</td><td>{dye}</td><td class=seq>{s}</td><td>{len(s) if s!='-' else '-'}</td><td>{tm}</td></tr>" for n,dye,s,tm in oligos)
    D0=min(b['aw'],b['at'],d.FWD['gs'] if d.FWD else 0)-1
    D1=max(b['bw'],b['bt'],d.REV['ge'] if d.REV else len(d.WT))+1
    W=D1-D0; pc=lambda x:round((x-D0)/W*100,2)
    bar=lambda a,bb,c,lab:f'<div class="bar {c}" style="left:{pc(a)}%;width:{pc(bb)-pc(a)}%">{lab}</div>'
    bars=(bar(d.FWD['gs'],d.FWD['ge'],'pr','FWD &rarr;') if d.FWD else '')+bar(b['at'],b['bt'],'ts','TS probe')+ \
         bar(b['aw'],b['bw'],'wt','WT probe')+(bar(d.REV['gs'],d.REV['ge'],'pr','&larr; REV') if d.REV else '')
    jx=pc(d.J); chip=lambda ok,w,lab,val:f'<span class="chip {"ok" if (ok and not w) else ("wn" if w else "no")}">{"&#10003;" if ok and not w else ("&#9888;" if w else "&#10007;")} {lab} <b>{val}</b></span>'
    cards="".join(f'<div class=oc><div class=on><b>{nm}</b><span class=sq>{sq}</span></div>'+"".join(chip(ok,w,lab,val) for lab,ok,val,w in cs)+'</div>' for nm,sq,cs in rows)
    cards+=f'<div class=oc><div class=on><b>Assay-level</b></div>{"".join(chip(ok,w,lab,val) for lab,ok,val,w in assay)}</div>'
    return f"""<!doctype html><meta charset=utf8><title>{title} — ddPCR oligos</title><style>
body{{font:15px/1.5 -apple-system,Segoe UI,sans-serif;margin:40px;color:#1a1a1a;max-width:980px}}
h1{{font-size:19px}}table{{border-collapse:collapse;margin-top:10px}}th,td{{border:1px solid #ccc;padding:8px 12px;text-align:left}}
th{{background:#2F5597;color:#fff}}.seq{{font-family:ui-monospace,Menlo,monospace}}td:nth-child(4),td:nth-child(5){{text-align:center}}
p{{color:#666;font-size:13px}}
.map{{position:relative;height:118px;margin:18px 0 8px;background:linear-gradient(90deg,#eef2fa {jx}%,#eafaf1 {jx}%);border:1px solid #ddd;border-radius:6px}}
.jl{{position:absolute;top:0;bottom:0;left:{jx}%;width:2px;background:#2F5597}}.jlab{{position:absolute;top:2px;left:{jx}%;transform:translateX(-50%);font-size:11px;color:#2F5597;font-weight:bold}}
.ex{{position:absolute;bottom:3px;font-size:11px;color:#88a}}
.bar{{position:absolute;height:20px;line-height:20px;font-size:11px;color:#fff;border-radius:3px;padding:0 5px;white-space:nowrap;overflow:hidden}}
.bar.wt{{background:#2F5597;top:30px}}.bar.ts{{background:#10a37f;top:56px}}.bar.pr{{background:#999;top:82px}}
.banner{{padding:9px 13px;border-radius:7px;font-weight:bold;margin:16px 0 6px}}.bpass{{background:#e3f5ec;color:#137a4b}}.bfail{{background:#fde6e6;color:#b3261e}}
.oc{{border:1px solid #e0e0e0;border-radius:7px;padding:9px 12px;margin:8px 0}}.on{{font-size:14px}}.on .sq{{color:#555;font-size:12.5px;margin-left:8px;font-family:ui-monospace,monospace}}
.chip{{display:inline-block;font-size:11.5px;padding:2px 7px;border-radius:11px;margin:4px 4px 0 0}}.chip.ok{{background:#e3f5ec;color:#137a4b}}.chip.no{{background:#fde6e6;color:#b3261e}}.chip.wn{{background:#fff4d6;color:#8a6100}}</style>
<h1>{title} — competitive ddPCR oligos</h1>
<div class=map><div class=jl></div><div class=jlab>junction</div><div class=ex style="left:6px">&larr; upstream exon</div><div class=ex style="right:6px">downstream exon &rarr;</div>{bars}</div>
<p style="margin-top:0">Both probes cross the junction — they only bind the spliced product. Probes designed independently (nested windows); matched Tms within 2&deg;C.</p>
<table><tr><th>Oligo</th><th>Dye</th><th>Sequence 5'&rarr;3'</th><th>nt</th><th>Tm</th></tr>{trows}</table>
<p>Bare sequences. Probes: add dye + ZEN/IBFQ double-quencher at order. Amplicon {d.REV['amp'] if d.REV else '-'} bp.</p>
<div class="banner {'bpass' if nf==0 else 'bfail'}">{'ALL CHECKS PASS' if nf==0 else str(nf)+' check(s) failed'}{f' &middot; {nw} flagged' if nw else ''}</div>
{cards}
<p style="color:#888">IDT-scale Tm (offset {d.IDT_OFF}C, MECP2-calibrated). Still pending: BLAST specificity vs paralogs/pseudogenes (separate step).</p>"""


def main():
    ap=argparse.ArgumentParser(description="Design a competitive ddPCR assay across a splice/trans-splice junction.")
    ap.add_argument("--wt", help="WT (reference) junction sequence")
    ap.add_argument("--ts", help="TS (variant) junction sequence, same length/frame as --wt")
    ap.add_argument("--junction", type=int, help="0-based index of the first base after the junction")
    ap.add_argument("--intron", type=int, default=None, help="intron length (bp) between the two exons, for the gDNA-safety check")
    ap.add_argument("--config", help="JSON file with {wt, ts, junction, intron, title}")
    ap.add_argument("--title", default="junction", help="label for the HTML report")
    ap.add_argument("--out", help="write an HTML report to this path")
    ap.add_argument("--primer-tm", type=float, default=61.0)
    ap.add_argument("--anneal", type=float, default=60.0)
    ap.add_argument("--idt-offset", type=float, default=7.5)
    ap.add_argument("--match-window", type=float, default=2.0, help="max Tm difference between the two matched probes")
    ap.add_argument("--example", action="store_true", help="run the built-in HEXA exon9|exon10 example")
    a=ap.parse_args()

    if a.config:
        cfg=json.load(open(a.config)); wt,ts,j=cfg["wt"],cfg["ts"],cfg["junction"]
        intron=cfg.get("intron"); title=cfg.get("title", a.title)
    elif a.wt and a.ts and a.junction is not None:
        wt,ts,j,intron,title=a.wt,a.ts,a.junction,a.intron,a.title
    else:
        if not a.example and (a.wt or a.ts):
            ap.error("provide --wt, --ts and --junction together (or --config, or --example).")
        wt,ts,j,intron,title=HEXA_WT,HEXA_TS,HEXA_JUNCTION,HEXA_INTRON,"HEXA exon9|exon10"

    d=Designer(wt,ts,j, primer_tm=a.primer_tm, anneal=a.anneal, idt_offset=a.idt_offset,
               match_window=a.match_window, intron=intron)
    d.design()
    rows,assay,nf,nw=d.audit()
    render_text(d,rows,assay,nf,nw)
    if a.out:
        open(a.out,"w").write(render_html(d,rows,assay,nf,nw,title))
        print("\nHTML:", a.out)


if __name__ == "__main__":
    main()
