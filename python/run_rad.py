#!/usr/bin/env python
# Ad-hoc driver: same compiled-default run as model.py, with a settable iteration cap,
# checkpoint cadence and binary restart control.
#
#   run_rad.py [nm] [checkpoint] [checkpoint_save_iter] [restart_from_iter]
#
#   nm                    total iterations                         (default 100)
#   checkpoint            VTK/VTS + printMinMax cadence            (default 10; param.py uses 8)
#   checkpoint_save_iter  write output-Jupiter/jup_restart_<iter>.bin at this iteration,
#                         -1 = only the automatic every-100-iteration checkpoints (default -1)
#   restart_from_iter     resume from output-Jupiter/jup_restart_<iter>.bin instead of
#                         spinning up from scratch, -1 = from scratch (default -1)
#
# Examples:
#   python3 run_rad.py 500 8            # fresh 500-iteration run, restarts saved at 100..500
#   python3 run_rad.py 800 8 -1 500     # continue that run from iteration 500 to 800
import sys
from pyatjup import Jupiter

nm         = int(sys.argv[1]) if len(sys.argv) > 1 else 100
checkpoint = int(sys.argv[2]) if len(sys.argv) > 2 else 10
save_iter  = int(sys.argv[3]) if len(sys.argv) > 3 else -1
from_iter  = int(sys.argv[4]) if len(sys.argv) > 4 else -1

j = Jupiter()
j.nm = nm
j.checkpoint = checkpoint
j.checkpoint_save_iter = save_iter
j.restart_from_iter = from_iter
print("   ATJUP run: nm =", j.nm, " checkpoint =", j.checkpoint,
      " checkpoint_save_iter =", j.checkpoint_save_iter,
      " restart_from_iter =", j.restart_from_iter)
j.run()
print("   done")
