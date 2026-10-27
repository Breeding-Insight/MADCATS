"""PyInstaller entry point for the standalone ``madcats`` binary.

``--jobs`` runs samples in a ProcessPoolExecutor.  In a frozen executable a
spawned worker re-launches the executable itself, so ``freeze_support`` must
run first: it turns that re-launch into a worker instead of a second CLI.
"""
import multiprocessing

from madcats import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    raise SystemExit(main())
