import logging
import subprocess
import sys


def syscall(command, allow_fail=False):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code: {completed_process.returncode}")
    if (not allow_fail) and completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout and stderr:",
            completed_process.stdout,
            sep="\n",
            file=sys.stderr,
        )
        raise RuntimeError("Error in system call. Cannot continue")

    logging.info(f"stdout:\n{completed_process.stdout.rstrip()}")
    logging.info(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process
