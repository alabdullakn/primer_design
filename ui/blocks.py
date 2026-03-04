import io
import streamlit as st


FOOTER_TEXT = (
    "This tool is free and open to everyone. "
    "It was designed and built by Khalid Alabdulla. "
    "If you find it useful, please consider sharing it. "
    "For feedback or bug reports, contact "
    "[alabdulla8932@gmail.com](mailto:alabdulla8932@gmail.com)."
)

def primer_diagram(
    template_len: int,
    fwd_start: int,
    fwd_len: int,
    rev_bind_start: int,
    rev_len: int,
    amp_len: int,
    junction_pos: int = None,
    display_width: int = 68,
) -> str:
    """Return an ASCII diagram string showing primer positions on the template."""
    W = display_width

    def scale(pos: int) -> int:
        return max(0, min(W - 1, round(pos * W / template_len)))

    fs = scale(fwd_start)
    fe = scale(fwd_start + fwd_len - 1)
    rs = scale(rev_bind_start)
    re = scale(rev_bind_start + rev_len - 1)

    # Template bar
    bar = ["─"] * W
    for k in range(fs, fe + 1):
        bar[k] = "█"
    for k in range(rs, re + 1):
        bar[k] = "█"
    if junction_pos is not None:
        j = scale(junction_pos)
        if 0 <= j < W:
            bar[j] = "│"

    template_line = "5'─" + "".join(bar) + "─3'"

    # Amplicon arrow
    a_s = scale(fwd_start)
    a_e = scale(rev_bind_start + rev_len - 1)
    arrow_body = "─" * max(0, a_e - a_s - 1)
    arrow_line = "   " + " " * a_s + "◄" + arrow_body + "►" + f"  amplicon: {amp_len} bp"

    # Position labels
    fwd_label = f"FWD: pos {fwd_start}–{fwd_start + fwd_len - 1} ({fwd_len} bp)"
    rev_label = f"REV: pos {rev_bind_start}–{rev_bind_start + rev_len - 1} ({rev_len} bp)"
    gap = max(1, (rs - fe) - len(fwd_label) + 3)
    label_line = "   " + " " * fs + fwd_label + " " * gap + rev_label

    lines = [template_line, arrow_line, label_line]
    if junction_pos is not None:
        lines.append(f"   Junction (│) at position {junction_pos}")

    return "\n".join(lines)


def render_diagram(
    template_len: int,
    fwd_start: int,
    fwd_len: int,
    rev_bind_start: int,
    rev_len: int,
    amp_len: int,
    junction_pos: int = None,
):
    """Render a primer position diagram in a Streamlit code block."""
    st.code(
        primer_diagram(template_len, fwd_start, fwd_len, rev_bind_start, rev_len, amp_len, junction_pos),
        language=None,
    )


def download_primers_csv(rows: list, filename: str = "primers.csv"):
    """Render a CSV download button for a list of row dicts."""
    import pandas as pd
    csv = pd.DataFrame(rows).to_csv(index=False)
    st.download_button(
        label="Download primers (CSV)",
        data=csv,
        file_name=filename,
        mime="text/csv",
    )


def add_footer():
    st.markdown("---")
    st.caption(FOOTER_TEXT)
