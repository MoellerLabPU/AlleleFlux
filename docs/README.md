# AlleleFlux Documentation

This directory contains the Sphinx documentation for AlleleFlux.

## Directory Structure

```
docs/
├── README.md              # This file
├── Makefile               # Build shortcuts (make html, make clean)
├── environment.yml        # Conda environment for building docs
├── source/                # Documentation source files
│   ├── conf.py            # Sphinx configuration
│   ├── index.md           # Main entry point
│   ├── _static/           # Static assets (images, CSS)
│   ├── getting_started/   # Installation, overview, quickstart
│   ├── usage/             # Workflow guides, visualization, dN/dS
│   ├── examples/          # Tutorials and use cases
│   └── reference/         # CLI, configuration, input/output specs
└── build/                 # Generated output (git-ignored)
    └── html/              # HTML documentation
```

## Prerequisites

Create and activate the docs environment:

```bash
# From the docs/ directory
conda env create -f environment.yml
conda activate alleleflux-docs
```

Or install dependencies directly:

```bash
conda install -c conda-forge sphinx furo myst-parser sphinx-copybutton sphinx-autobuild
```

## Building Documentation

### Build HTML

```bash
cd docs
make html
```

Or using sphinx-build directly:

```bash
sphinx-build -b html source build/html
```

Output will be in `build/html/`.

### Clean Build

```bash
make clean html
```

### View Locally

Open `build/html/index.html` in your browser:

```bash
# macOS
open build/html/index.html

# Linux
xdg-open build/html/index.html
```

## Building on a Remote Server

When working on a remote server (e.g., via SSH), use one of these methods to view the docs:

### Option 1: Live Preview with Port Forwarding (Recommended)

Start the autobuild server on the remote:

```bash
cd docs
sphinx-autobuild source build/html --host 0.0.0.0 --port 8000
```

Then set up SSH port forwarding from your local machine:

```bash
ssh -L 8000:localhost:8000 user@remote-server
```

Open http://localhost:8000 in your local browser. Changes auto-reload.

### Option 2: Build and Download

Build on remote:

```bash
cd docs && make html
```

Download the built docs:

```bash
# From your local machine
scp -r user@remote-server:/path/to/AlleleFlux/docs/build/html ./alleleflux-docs
```

Open `alleleflux-docs/index.html` locally.

### Option 3: Simple HTTP Server

Build and serve with Python:

```bash
cd docs
make html
cd build/html
python -m http.server 8000
```

Access via SSH tunnel as in Option 1.

## Documentation Format

The documentation uses **MyST Markdown** (`.md` files) with Sphinx. Key syntax:

### Toctrees

```markdown
```{toctree}
:maxdepth: 1

getting_started/overview.md
getting_started/installation.md
```
```

### Admonitions

```markdown
:::{note}
This is a note.
:::

:::{warning}
This is a warning.
:::
```

### Cross-References

Use standard Markdown links:

```markdown
See [Installation](getting_started/installation.md) for details.
```

### Code Blocks

````markdown
```python
def example():
    return "Hello"
```
````

## Troubleshooting

**VS Code preview shows raw directives**: This is expected. VS Code's Markdown preview doesn't understand MyST/Sphinx syntax. Build with Sphinx to see the actual output.

**Build warnings about files not in toctree**: Legacy `*_old.md` files are intentionally excluded. These can be deleted if no longer needed.

**Missing dependencies**: Ensure you've activated the `alleleflux-docs` environment or installed all packages from `environment.yml`.
