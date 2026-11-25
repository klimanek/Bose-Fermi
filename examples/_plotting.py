"""Shared plotting helper for the example scripts.

On a desktop session ``plt.show()`` opens a window.  On machines without a
GUI (a container, CI, or SSH without X forwarding) matplotlib silently falls
back to the non-interactive ``agg`` backend, where ``show()`` only emits
"FigureCanvasAgg is non-interactive".  In that case the figure is written to
``figures/`` instead, so the examples still produce something to look at.
"""

from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt

# Backends that can only write files, never open a window.
NON_INTERACTIVE_BACKENDS = {
    "agg",
    "cairo",
    "pdf",
    "pgf",
    "ps",
    "svg",
    "template",
}

FIGURE_DIR = Path(__file__).resolve().parent.parent / "figures"


def is_interactive() -> bool:
    """Return True when the current matplotlib backend can open a window."""
    return matplotlib.get_backend().lower() not in NON_INTERACTIVE_BACKENDS


def _target_path(filename: str) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    return FIGURE_DIR / filename


def _announce(path: Path) -> None:
    print(
        f"No interactive matplotlib backend ({matplotlib.get_backend()}); "
        f"figure saved to {path}"
    )


def show_or_save(
    fig=None, filename: str = "figure.png", dpi: int = 150
) -> Optional[Path]:
    """Show *fig* in a window, or save it to ``figures/<filename>`` if headless."""
    fig = fig if fig is not None else plt.gcf()
    if is_interactive():
        plt.show()
        return None
    path = _target_path(filename)
    fig.savefig(path, dpi=dpi)
    _announce(path)
    return path


def show_or_save_sympy_plot(
    plot, filename: str = "figure.png"
) -> Optional[Path]:
    """Same as :func:`show_or_save`, but for a ``sympy.plotting.Plot`` object."""
    if is_interactive():
        plot.show()
        return None
    path = _target_path(filename)
    plot.save(str(path))
    _announce(path)
    return path
