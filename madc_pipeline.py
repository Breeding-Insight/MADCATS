"""Compatibility shim: the tool is now MADCATS (madcats.py).

Everything that imported or invoked madc_pipeline continues to work;
new code should import madcats directly.
"""
from madcats import *          # noqa: F401,F403
from madcats import main, __version__  # noqa: F401

if __name__ == "__main__":
    import sys
    sys.exit(main())
