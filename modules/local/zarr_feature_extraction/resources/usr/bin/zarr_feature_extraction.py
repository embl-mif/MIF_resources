#!/usr/bin/env python3
import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from internal.cli import main

if __name__ == "__main__":
    sys.exit(main())
