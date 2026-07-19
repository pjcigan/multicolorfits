#!/usr/bin/env python
"""Generate ``llms.txt`` and ``llms-full.txt`` from the package capability catalog.

Both files are rendered by :func:`multicolorfits._overview.render_llms` — the
SAME source of truth that ``mcf.overview()`` / ``mcf.recipes()`` read — so they
cannot drift from the code.  ``llms.txt`` carries the colorize-then-combine
model, the conventions, and the runnable recipe code; ``llms-full.txt`` adds a
function-signature reference.  ``tests/test_overview.py`` fails if either file
is out of sync (regenerate here to fix).

Run from the repo root:  python scripts/make_llms_txt.py
"""

from __future__ import annotations

import os

from multicolorfits._overview import render_llms

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> None:
    for name, full in (("llms.txt", False), ("llms-full.txt", True)):
        text = render_llms(full=full)
        path = os.path.join(ROOT, name)
        with open(path, "w") as fh:
            fh.write(text)
        print(f"wrote {path} ({len(text)} chars)")


if __name__ == "__main__":
    main()
