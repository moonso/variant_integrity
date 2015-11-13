#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
variant_integrity.__main__
~~~~~~~~~~~~~~~~~~~~~

The main entry point for the command line interface.

Invoke as ``variant_integrity`` (if installed)
or ``python -m variant_integrity`` (no install required).
"""
import sys

from .cli import cli


if __name__ == '__main__':
    # exit using whatever exit code the CLI returned
    sys.exit(cli())
