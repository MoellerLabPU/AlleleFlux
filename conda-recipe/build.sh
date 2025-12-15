#!/bin/bash
set -ex

$PYTHON -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv
