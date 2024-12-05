import json
import subprocess
import sys


def postprocess(method):
    try:
        import h5py as h5
        import numpy as np
    except ImportError:
        print("Did not find h5py and/or numpy. Exiting without evaluating.")
        return

    error = np.zeros(len(nsteps))

    for idx, case in enumerate(nsteps):
        with h5.File(f"probes-{case}.h5", "r") as fh:
            time = fh["/PROBES/time"][...]
            phi = fh["/PROBES/point/PHI"][:, 0]

        # Analytical solution: y = exp(-t)
        phi_analytical = np.exp(-time["TIME"])

        # L2 norm of the error
        error[idx] = np.linalg.norm(phi - phi_analytical)

    return error


def plot(errors, labels):
    try:
        import matplotlib.figure as figure
    except ImportError:
        print("Did not find matplotlib and/or numpy. Exiting without "
              "evaluating.")
        return None

    fig = figure.Figure()
    axs = fig.subplots(1, 1)
    fig.suptitle("Scalar decay testcase")

    for idx, error in enumerate(errors):
        axs.loglog(nsteps, error, "o-", label=labels[idx])

    fig.suptitle("Scalar decay testcase")
    axs.set_xlabel("Time steps per time unit")
    axs.set_ylabel("L2 norm of error")
    axs.legend()
    axs.grid()

    fig.tight_layout()
    fig.savefig("plot.pdf")


def run(rkmethods, nsteps, tend, labels, mglet_bin):
    errors = []

    for method in rkmethods:
        # Run all timestep widths
        for case in nsteps:
            # Read parametrs-0.json and insert the right parameters
            with open("parameters-0.json", "r") as fh:
                params = json.load(fh)

            params["time"]["mtstep"] = int(tend*case)
            params["time"]["dt"] = 1.0/case
            params["time"]["tend"] = tend
            params["time"]["rkmethod"] = method
            params["probes"]["file"] = f"probes-{case}.h5"

            # Write the new parameters to a new file
            with open(f"parameters-{case}.json", "w") as fh:
                json.dump(params, fh, indent=4)

            # Run the simulation
            with open(f"mglet-{method}-{case}.OUT", "w") as outfile:
                subprocess.run([mglet_bin, f"parameters-{case}.json"],
                               stdout=outfile, stderr=outfile)

        # Run the post-processing
        error = postprocess(method)
        if error is not None:
            errors.append(error)

    # Plot the results
    if len(errors) > 0:
        plot(errors, labels)


mglet_bin = sys.argv[1]

rkmethods = ["williamson", "berland", "carpenter", "bernardini"]
labels = ["Williamson (3rd 3 step)", "Berland (4th 6 step)",
          "Carpenter (4th 5 step)", "Bernardini (2nd 5 step)"]
nsteps = [1, 10, 100, 1000, 10000]
tend = 10.0

run(rkmethods, nsteps, tend, labels, mglet_bin)
