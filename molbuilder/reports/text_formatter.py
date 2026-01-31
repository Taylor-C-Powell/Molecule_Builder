"""ASCII text formatting utilities for report generation.

All output is cp1252-safe (pure ASCII printable characters only).
Used by molecule_report, synthesis_report, safety_report, and cost_report.
"""

from __future__ import annotations

import textwrap


# =====================================================================
#  Section headers
# =====================================================================

def section_header(title: str, width: int = 70, char: str = "=") -> str:
    """Generate a centered section header with border lines.

    Example::

        ======================================================================
                                 MOLECULE REPORT
        ======================================================================
    """
    border = char * width
    padded = title.upper().center(width)
    return "\n".join([border, padded, border])


def subsection_header(title: str, width: int = 70, char: str = "-") -> str:
    """Generate a subsection header.

    Example::

        --- Atom Composition --------------------------------------------------
    """
    prefix = char * 3 + " " + title + " "
    return prefix + char * max(0, width - len(prefix))


# =====================================================================
#  Tables
# =====================================================================

def ascii_table(headers: list[str], rows: list[list[str]],
                alignments: list[str] | None = None,
                min_widths: list[int] | None = None) -> str:
    """Generate a formatted ASCII table with column alignment.

    Parameters
    ----------
    headers : list[str]
        Column header labels.
    rows : list[list[str]]
        Table data (each row is a list of cell strings).
    alignments : list[str] | None
        Per-column alignment: ``'l'`` left, ``'r'`` right, ``'c'`` center.
        Defaults to left-aligned for every column.
    min_widths : list[int] | None
        Minimum column widths.  Actual widths expand to fit content.
    """
    n_cols = len(headers)
    if alignments is None:
        alignments = ["l"] * n_cols
    if min_widths is None:
        min_widths = [0] * n_cols

    # Compute column widths
    col_widths: list[int] = []
    for c in range(n_cols):
        w = max(len(headers[c]), min_widths[c])
        for row in rows:
            if c < len(row):
                w = max(w, len(str(row[c])))
        col_widths.append(w)

    def _fmt_cell(text: str, width: int, align: str) -> str:
        if align == "r":
            return text.rjust(width)
        if align == "c":
            return text.center(width)
        return text.ljust(width)

    sep = "  "
    header_line = sep.join(
        _fmt_cell(headers[c], col_widths[c], alignments[c])
        for c in range(n_cols)
    )
    divider = sep.join("-" * col_widths[c] for c in range(n_cols))

    lines = [header_line, divider]
    for row in rows:
        cells: list[str] = []
        for c in range(n_cols):
            val = str(row[c]) if c < len(row) else ""
            cells.append(_fmt_cell(val, col_widths[c], alignments[c]))
        lines.append(sep.join(cells))

    return "\n".join(lines)


# =====================================================================
#  Text utilities
# =====================================================================

def word_wrap(text: str, width: int = 70, indent: int = 0) -> str:
    """Word-wrap text to the given width with optional indent."""
    prefix = " " * indent
    wrapped = textwrap.fill(
        text, width=width, initial_indent=prefix,
        subsequent_indent=prefix,
    )
    return wrapped


def bullet_list(items: list[str], indent: int = 2, bullet: str = "-") -> str:
    """Format items as a bulleted list."""
    prefix = " " * indent + bullet + " "
    subsequent = " " * (indent + len(bullet) + 1)
    lines: list[str] = []
    for item in items:
        wrapped = textwrap.fill(
            item, width=70,
            initial_indent=prefix,
            subsequent_indent=subsequent,
        )
        lines.append(wrapped)
    return "\n".join(lines)


def key_value_block(pairs: list[tuple[str, str]], separator: str = ": ",
                    indent: int = 2) -> str:
    """Format key-value pairs aligned on the separator."""
    if not pairs:
        return ""
    max_key = max(len(k) for k, _ in pairs)
    prefix = " " * indent
    lines: list[str] = []
    for key, value in pairs:
        lines.append(f"{prefix}{key:<{max_key}}{separator}{value}")
    return "\n".join(lines)


# =====================================================================
#  Charts and number formatting
# =====================================================================

def horizontal_bar(value: float, max_value: float, width: int = 40,
                   char: str = "#") -> str:
    """Render a simple horizontal bar chart line.

    Returns a string of *char* characters proportional to *value / max_value*,
    padded to *width* with spaces.
    """
    if max_value <= 0:
        filled = 0
    else:
        ratio = max(0.0, min(1.0, value / max_value))
        filled = int(round(ratio * width))
    return char * filled + " " * (width - filled)


def format_currency(amount: float) -> str:
    """Format as USD with commas: ``$1,234.56``."""
    return "${:,.2f}".format(amount)


def format_percent(value: float, decimals: int = 1) -> str:
    """Format as percentage: ``85.0%``."""
    return "{:.{d}f}%".format(value, d=decimals)
