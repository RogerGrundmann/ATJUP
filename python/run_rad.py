#!/usr/bin/env python
# Ad-hoc driver: same compiled-default run as model.py, with a settable iteration cap and
# checkpoint cadence, so we can inspect the radiation / Q_rad diagnostics.
import sys
from pyatjup import Jupiter

nm         = int(sys.argv[1]) if len(sys.argv) > 1 else 100
checkpoint = int(sys.argv[2]) if len(sys.argv) > 2 else 10

j = Jupiter()
j.nm = nm
j.checkpoint = checkpoint
print("   ATJUP run: nm =", j.nm, " checkpoint =", j.checkpoint)
j.run()
print("   done")
