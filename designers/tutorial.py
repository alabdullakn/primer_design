# designers/tutorial.py
import streamlit as st



_TOPICS = {
    "regular_pcr": "Regular PCR",
    "qpcr": "qPCR Primers",
    "splicing": "Splicing Primers",
}


def _nav_button(label: str, topic: str, active_topic: str):
    btn_type = "primary" if topic == active_topic else "secondary"
    if st.button(label, key=f"tut_nav_{topic}", type=btn_type, use_container_width=True):
        st.session_state["tutorial_topic"] = topic
        st.rerun()


def render():
    st.title("How to Use PrimerQ")

    # Topic selector (also settable from designer buttons)
    topic = st.session_state.get("tutorial_topic", "regular_pcr")

    # Style subtab buttons to match main nav tabs
    st.markdown("""
    <style>
    /* Tutorial subtab inactive */
    button[kind="secondary"], button[data-testid="baseButton-secondary"] {
        background: rgba(15, 23, 42, 0.82) !important;
        border: 1px solid rgba(148, 163, 184, 0.35) !important;
        border-radius: 10px !important;
        color: #e2e8f0 !important;
        font-weight: 700 !important;
        font-size: 1.02rem !important;
        width: 100% !important;
    }
    /* Tutorial subtab active */
    button[kind="primary"], button[data-testid="baseButton-primary"] {
        background: linear-gradient(180deg, rgba(37, 99, 235, 0.96), rgba(30, 64, 175, 0.96)) !important;
        border: 1px solid #7dd3fc !important;
        border-radius: 10px !important;
        color: #ffffff !important;
        font-weight: 700 !important;
        font-size: 1.02rem !important;
        box-shadow: 0 8px 20px rgba(37, 99, 235, 0.35) !important;
        width: 100% !important;
    }
    </style>
    """, unsafe_allow_html=True)

    cols = st.columns(3)
    with cols[0]:
        _nav_button("Regular PCR", "regular_pcr", topic)
    with cols[1]:
        _nav_button("Splicing Primers", "splicing", topic)
    with cols[2]:
        _nav_button("qPCR Primers", "qpcr", topic)

    st.divider()

    if topic == "regular_pcr":
        _tutorial_regular_pcr()
    elif topic == "qpcr":
        _tutorial_qpcr()
    elif topic == "splicing":
        _tutorial_splicing()


# ─────────────────────────────────────────────
# Regular PCR tutorial
# ─────────────────────────────────────────────

def _tutorial_regular_pcr():
    st.header("Regular PCR — Tutorial")
    st.write(
        "Designs a standard forward + reverse primer pair from a single template sequence "
        "for use in conventional (end-point) PCR."
    )

    st.subheader("Step 1 — Provide a template sequence")
    st.markdown("""
Choose one of three input methods:

| Method | When to use |
|--------|-------------|
| **Paste sequence** | You already have the DNA sequence as text (FASTA or raw A/C/G/T) |
| **Upload FASTA** | You have a `.fa` / `.fasta` file on your computer |
| **NCBI accession** | You know the GenBank accession (e.g. `NM_007294.4`). Click **Fetch from NCBI** to download it automatically. |

The tool accepts only unambiguous bases — anything that isn't A, C, G, or T is silently removed.
""")

    st.subheader("Step 2 — Set primer design parameters")
    st.markdown("""
Expand **Primer design parameters** to adjust:

**Primer length & Tm**
- `Min / Max primer length` — typical range is 18–25 bp. Longer primers have higher Tm and more specificity.
- `Target Tm (°C)` — the melting temperature you want. 60 °C is a good default for most polymerases.
- `Tm tolerance (±°C)` — how far from the target Tm is still acceptable. Widen this (e.g. to 7–10 °C) if no primers are found.

**Amplicon size**
- `Amplicon min / max (bp)` — the expected PCR product length. Keep realistic bounds. A 200–500 bp product is typical for gel-based PCR.

**GC content**
- `GC min / max (%)` — primers outside this range are rejected. The default 35–65% covers most cases.

**Search windows**
- `Forward search window` — how many bases from the **5' end** of the template to search for a forward primer.
- `Reverse search window` — how many bases from the **3' end** to search for a reverse primer.
- Narrow these to force primers into a specific region; widen if no pair is found.

**Reaction buffer** *(optional)*
- Enable **Edit reaction buffer** to enter your actual Na⁺, Mg²⁺, and dNTP concentrations. This corrects the nearest-neighbour Tm calculation for your specific conditions.
""")

    st.subheader("Step 3 — Click Design primers")
    st.markdown("""
The engine:
1. Scans every candidate within the search windows
2. Scores each primer (Tm deviation, GC content, homopolymer runs, hairpin/homodimer thermodynamics)
3. Picks the lowest-scoring *pair* that satisfies the amplicon size constraint and has no 3' complementarity (primer-dimer risk)
""")

    st.subheader("Step 4 — Interpret the results")
    st.markdown("""
**Primer table** — one row per primer with:
- `Tm (°C)` — nearest-neighbour melting temperature at your buffer conditions
- `GC (%)` — GC content percentage
- `Score` — lower is better (0 = perfect match to target Tm, no structural issues, GC clamp present)
- `Start (0-based)` — position in the template where the primer binds

**Amplicon length** — estimated size of the PCR product.

**Primer-BLAST link** — opens NCBI Primer-BLAST in your browser with both primers pre-filled. Use this to check genome-wide specificity.

**Dimer check** — a heuristic percentage of 3'-complementarity risk between the two primers. Values below ~15% are generally acceptable.

**In-app BLAST** — click *Check Specificity with BLAST* to query NCBI directly. Each primer is BLASTed against the nucleotide database filtered to your chosen organism. Results show identity %, query coverage, and an interpretation (Specific / Non-specific / No match).
""")

    st.subheader("Troubleshooting")
    st.markdown("""
| Error | Fix |
|-------|-----|
| *No primer pair found* | Widen Tm tolerance, widen amplicon range, or increase the search windows |
| *Forward/reverse window too small* | The window must be at least 2× the minimum primer length |
| *Amplicon min > max* | Fix the amplicon range inputs |
""")


# ─────────────────────────────────────────────
# qPCR tutorial
# ─────────────────────────────────────────────

def _tutorial_qpcr():
    st.header("qPCR Primers — Tutorial")
    st.write(
        "Designs primer pairs for quantitative PCR (qPCR) that amplify only spliced mRNA, "
        "not genomic DNA. Two strategies are available — choose the one that fits your gene and experiment."
    )

    # ── Strategy overview ──────────────────────────────────────────────────────
    st.subheader("Choosing a strategy")
    st.markdown("""
| Strategy | How it works | Best when |
|----------|-------------|-----------|
| **Junction-spanning** | One primer straddles the exon-exon junction — its 3' end sits on the downstream exon, so it cannot extend from unspliced pre-mRNA or genomic DNA | You know exactly which junction to target; the intron at that junction may be small |
| **Intron-flanking** | Both primers sit in different exons; the intron between them makes the genomic amplicon too large to amplify efficiently | The intron between your target exons is large (>500 bp); you want a simpler design without junction marking |

Both strategies produce a short cDNA amplicon (70–200 bp). Select the strategy using the radio button at the top of the qPCR tab.
""")

    st.divider()

    # ── STRATEGY A: Junction-spanning ─────────────────────────────────────────
    st.subheader("Strategy A — Junction-spanning")

    st.markdown("##### Step 1 — Prepare your sequence")
    st.markdown("""
Paste the **spliced mRNA / cDNA** sequence and mark the exon-exon junction with a caret `^`:

```
...AAGCTAGTCCAG^GTTAGCCTAAGG...
```

- Everything **before** `^` = end of the upstream exon
- Everything **after** `^` = start of the downstream exon
- Provide **at least 10 bases on each side** of `^` (50+ recommended for more candidates)
""")

    st.markdown("##### Step 2 — Chemistry")
    st.markdown("""
| Chemistry | What you get |
|-----------|-------------|
| **SYBR Green** | Primer pair only. Cheaper. Requires melt-curve analysis. |
| **TaqMan probe** | Primer pair + internal hydrolysis probe. More specific, no melt curve needed. |

TaqMan is only available in junction-spanning mode.
""")

    st.markdown("##### Step 3 — Junction primer direction")
    st.markdown("""
- `AUTO` — tries both FWD and REV spanning and picks the best-scoring pair
- `FWD` — forces the forward primer to span the junction (3' portion in the downstream exon)
- `REV` — forces the reverse primer to span the junction (3' portion in the upstream exon)

`AUTO` is recommended in most cases.
""")

    st.markdown("##### Step 4 — Junction overlap controls")
    st.markdown("""
- `Min overlap each side of junction` — minimum bases the junction primer must contribute from **each** exon.
  ≥6 bp is standard; less overlap risks priming on unspliced pre-mRNA.
- `Max distance of junction from 3' end` — limits how far the junction can sit from the primer's 3' end.
  The engine also **scores** junction position: primers where the junction is at the very last 1–2 bases
  are ranked higher, since polymerase is most sensitive to mismatches at the 3' end.
""")

    st.markdown("##### Step 5 — Results")
    st.markdown("""
- `Role` column shows **Junction** (the primer that straddles the splice site) vs **Partner**
- The junction primer should show **no perfect genomic DNA hit** in BLAST — this is expected and confirms mRNA specificity
- For TaqMan, a probe is designed between the two primers; it is also BLASTed in the in-app check
""")

    st.divider()

    # ── STRATEGY B: Intron-flanking ───────────────────────────────────────────
    st.subheader("Strategy B — Intron-flanking")

    st.markdown("##### How it works")
    st.markdown("""
```
Exon A                          Intron (large)                   Exon B
──────────────────────────────  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ──────────────────────
          [  FWD primer  ]──────────────────────────────────────────────[  REV primer  ]
              ↑ binds here (cDNA)                                         ↑ binds here (cDNA)

cDNA amplicon:  FWD ──────────────────── REV   (70–200 bp)
Genomic amplicon: FWD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ REV   (>>500 bp, does not amplify efficiently)
```

Because the intron is large, polymerase cannot complete the genomic amplicon under standard qPCR cycling conditions (short extension time). Only cDNA gives a product.
""")

    st.markdown("##### Step 1 — Paste exon sequences separately")
    st.markdown("""
Two text boxes appear when you select **Intron-flanking**:

| Box | What to paste |
|-----|---------------|
| **Exon A** | The exon where the FWD primer will bind — paste the last ~100 bp of this exon (the portion closest to the intron) |
| **Exon B** | The exon where the REV primer will bind — paste the first ~100 bp of this exon (the portion closest to the intron) |

The engine concatenates Exon A + Exon B internally to measure the cDNA amplicon size. The FWD primer is strictly constrained to stay within Exon A; REV is strictly within Exon B.
""")

    st.markdown("##### Step 2 — Parameters")
    st.markdown("""
All standard parameters apply (Tm, GC, length, amplicon size, homopolymer, reaction buffer).
There are no junction-specific controls in this mode.

- `Amplicon min / max` refers to the **cDNA** amplicon (Exon A tail + Exon B head). Keep it 70–200 bp.
- TaqMan probe is not available in this mode.
""")

    st.markdown("##### Step 3 — Results")
    st.markdown("""
- `Role` column shows **Exon A** (FWD) and **Exon B** (REV)
- The result banner shows the **cDNA amplicon size** and notes that the genomic amplicon will be larger by the intron size
- On a gel, you can confirm specificity: cDNA gives a band at the expected size; genomic DNA gives a much larger band or no band
- In BLAST, both primers should show a specific hit in your target gene. Unlike junction-spanning, you may see genomic hits — that is normal; the intron size is what provides discrimination, not the primer sequence itself
""")

    st.divider()

    # ── Shared parameters ─────────────────────────────────────────────────────
    st.subheader("Shared parameters (both strategies)")
    st.markdown("""
- `Amplicon min / max` — 70–200 bp is ideal for qPCR efficiency
- `Tm target / tolerance` — tight tolerance (±2 °C) ensures both primers amplify at the same efficiency
- `Max Tm difference pair` — both primers must be within this many °C of each other
- `GC 40–60%` — tighter than regular PCR; keeps melting behaviour predictable
- `Max homopolymer run` — runs of 4+ identical bases destabilise primers; 3 is a safe limit
- `Reaction buffer` — adjust if you know your master mix Na⁺/Mg²⁺/dNTP composition
""")

    st.subheader("Troubleshooting")
    st.markdown("""
| Error | Strategy | Fix |
|-------|----------|-----|
| *No qPCR junction primer pair found* | Junction-spanning | Widen Tm tolerance, increase amplicon max, reduce min overlap, or switch junction direction |
| *Contradictory settings* | Junction-spanning | `Max 3' distance` must be ≥ `Min overlap each side` |
| *No valid probe found* | Junction-spanning (TaqMan) | Widen probe Tm tolerance or increase probe length range |
| *No intron-flanking pair found* | Intron-flanking | Provide longer exon sequences (100+ bp each), widen Tm tolerance or GC range |
| *Exon A/B too short* | Intron-flanking | Paste more sequence — at least as many bases as the max primer length |
""")


# ─────────────────────────────────────────────
# Splicing tutorial
# ─────────────────────────────────────────────

def _tutorial_splicing():
    st.header("Splicing Primers — Tutorial")
    st.write(
        "Designs primers to detect alternative splicing events — exon skipping, intron retention, "
        "or the use of alternative splice sites."
    )

    st.subheader("The two conditions")
    st.markdown("""
The tool supports two experimental layouts, illustrated in the diagram on the Splicing page:

**Condition 1: FWD A + FWD B + REV A**

```
Exon 1   ──────────────────────────────────►
         [FWD A]
                    Exon 2
                    [FWD B]
                                     Exon 3
                              ◄─────────────[REV A]
```

Use this when you want to detect whether a **middle exon is included or skipped**:
- FWD A + REV A → amplifies the full-length isoform (all exons present)
- FWD B + REV A → amplifies the skipped isoform (Exon 2 excluded)

**Condition 2: FWD A + REV A + REV B**

Mirror image — use when the alternatively spliced region is on the reverse side.
""")

    st.subheader("Step 1 — Paste exon sequences")
    st.markdown("""
After choosing a condition, text boxes appear for the required exon sequences.

- Paste the **sense strand (5'→3')** for each exon
- The tool accepts raw A/C/G/T — FASTA headers are ignored, other characters are removed
- Provide enough sequence in each exon for the engine to find a good primer (at least 50 bp recommended)

**Which sequence goes where?**

| Input box | What to paste |
|-----------|---------------|
| Forward A | The exon where FWD primer A should bind (e.g. Exon 1) |
| Forward B | The exon where FWD primer B should bind (e.g. Exon 2) |
| Reverse A | The exon where REV primer A should bind (e.g. Exon 3) — paste the sense strand, the tool reverses it |
| Reverse B | The exon where REV primer B should bind |
""")

    st.subheader("Step 2 — Set parameters")
    st.markdown("""
Parameters work the same as Regular PCR:
- `Min / Max primer length`, `Target Tm`, `Tm tolerance`, `GC min/max`
- `3' dimer check window (k)` — the engine checks that the reverse primer's 3' end has no k-mer complementarity with either forward primer. k=4 is a reasonable default.

**Reaction buffer** — customize ionic conditions if needed.
""")

    st.subheader("Step 3 — Interpret the results")
    st.markdown("""
**Primer table** — one section per pair (e.g. "FWD A + REV A", "FWD B + REV A"). Each row shows primer sequence, Tm, GC%, and score.

**Expected band sizes** — you'll need to calculate these yourself based on the isoform exon sizes. Run on a gel and compare band sizes between isoforms to confirm splicing.

**Primer-BLAST** — check each pair. For splicing experiments, primers should not amplify other genes with similar sequences.

**Dimer check** — for each set (A or B), shows complementarity risk between the two forward primers and the shared reverse primer.

**In-app BLAST** — BLASTs each unique primer (FWD A, FWD B, REV A, REV B). All should show a single specific hit in your target gene.
""")

    st.subheader("Troubleshooting")
    st.markdown("""
| Error | Fix |
|-------|-----|
| *No primer found for Exon X* | Provide a longer exon sequence, widen Tm tolerance, or widen GC range |
| *No reverse primer avoids dimer* | Reduce dimer_k from 4 to 3, or widen Tm tolerance |
| *Missing: Forward A / Reverse A* | Fill in all required sequence boxes before clicking Design |
""")
