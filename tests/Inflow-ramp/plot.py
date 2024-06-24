import sys

try:
    import h5py as h5
    import numpy as np
    import matplotlib.figure as figure
except ImportError:
    print("Did not find h5py, matplotlib and/or numpy. Exiting without "
          "evaluating.")
    sys.exit(0)

fig = figure.Figure()
ax = fig.subplots(1, 2)
fig.suptitle("Inflow ramp example")

for idx, case in enumerate(["ramp", "ramp_inf"]):
    with h5.File(f"probes-{case}.h5", "r") as fh:
        time = fh["/PROBES/time"][...]
        U = fh["/PROBES/point/U"][...]

    ax[idx].plot(time["TIME"], U, label=case)
    ax[idx].set_xlabel("Time")
    ax[idx].set_title(case)

ax[0].set_ylabel("U")

fig.savefig("plot.pdf")
