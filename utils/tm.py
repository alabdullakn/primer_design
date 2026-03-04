import primer3

# dg values from primer3 are in cal/mol (not kcal/mol).
# Thresholds are evaluated at annealing temperature (55°C), not 37°C,
# so only structures stable enough to compete at PCR annealing are rejected.
_HAIRPIN_DG_THRESHOLD = -2000.0   # -2 kcal/mol at 55°C
_HOMODIMER_DG_THRESHOLD = -6000.0  # -6 kcal/mol at 55°C
_STRUCTURE_TEMP_C = 55.0


def tm_nn(
    seq: str,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.2,
    dna_conc: float = 250.0,
) -> float:
    """
    Nearest-neighbor Tm (SantaLucia 1998) via primer3-py.

    Args:
        seq:       Primer sequence (5'->3'), DNA only.
        mv_conc:   Monovalent cation concentration (Na+) in mM. Default 50 mM.
        dv_conc:   Divalent cation concentration (Mg2+) in mM. Default 1.5 mM.
        dntp_conc: dNTP concentration in mM. Default 0.2 mM.
        dna_conc:  Oligo concentration in nM. Default 250 nM.

    Returns:
        Tm in degrees Celsius.
    """
    return primer3.calc_tm(
        seq.upper(),
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    )


def has_hairpin(
    seq: str,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.2,
    dna_conc: float = 250.0,
) -> bool:
    """
    Returns True if the primer is predicted to form a stable hairpin structure
    (dG < -2 kcal/mol at 55°C annealing temperature).
    """
    result = primer3.calc_hairpin(
        seq.upper(),
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=_STRUCTURE_TEMP_C,
    )
    return result.dg < _HAIRPIN_DG_THRESHOLD


def has_homodimer(
    seq: str,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.2,
    dna_conc: float = 250.0,
) -> bool:
    """
    Returns True if the primer is predicted to form a stable homodimer
    (dG < -6 kcal/mol at 55°C annealing temperature).
    """
    result = primer3.calc_homodimer(
        seq.upper(),
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=_STRUCTURE_TEMP_C,
    )
    return result.dg < _HOMODIMER_DG_THRESHOLD
