# Documentation Migration Notes: RST to MyST Markdown

This document describes the migration of AlleleFlux documentation from reStructuredText (RST) to MyST Markdown.

## Summary of Changes

### What Changed

1. **Documentation Format**: All `.rst` files in `docs/source/` were converted to MyST Markdown (`.md`)
2. **Sphinx Configuration**: Updated `conf.py` to recognize both `.rst` and `.md` source files
3. **Dependencies**: `myst-parser` is now a required docs dependency (already was present)

### Files Converted

- `docs/source/index.md` (main entry point)
- `docs/source/getting_started/*.md` (3 files)
- `docs/source/usage/*.md` (5 files)
- `docs/source/examples/*.md` (2 files)
- `docs/source/reference/*.md` (7 files, including legacy `*_old.md` files)

### Legacy Files

The following `*_old.md` files are retained but not included in the toctree:
- `cli_reference_old.md`
- `inputs_old.md`
- `outputs_old.md`

These files contain `{eval-rst}` blocks for RST list-tables, which still work but could be converted to native Markdown tables in a future cleanup.

## Building Documentation Locally

### Prerequisites

```bash
# Install docs dependencies
pip install -r docs/requirements.txt
```

### Build HTML Documentation

```bash
cd docs
sphinx-build -b html source build/html
```

Or use the Makefile:

```bash
cd docs
make html
```

### Live Preview (with auto-reload)

```bash
cd docs
sphinx-autobuild source build/html
```

The docs will be served at `http://127.0.0.1:8000`.

## MyST Markdown Syntax Reference

### Directives (Fenced)

MyST uses fenced code blocks for directives:

```markdown
```{toctree}
:maxdepth: 2

getting_started/overview
getting_started/installation
```
```

### Admonitions

Use colon-fence syntax for admonitions:

```markdown
:::{note}
This is a note.
:::

:::{warning}
This is a warning.
:::
```

### Cross-References

- **Document links**: `{doc}`\`path/to/file\``
- **Section references**: `{ref}`\`label-name\``
- **Markdown links**: `[text](path/to/file.md)`

### Code Blocks

Standard Markdown fenced code blocks work:

````markdown
```python
def example():
    pass
```
````

## Known Limitations

1. **Legacy RST Tables**: Some files contain `{eval-rst}` blocks with RST list-tables. These render correctly but could be converted to native Markdown tables.

2. **Lexer Names**: Some code blocks use non-standard lexer names. Use `text` for TSV/CSV data instead of `tsv`.

3. **Build Warnings**: 5 warnings remain about legacy files not in toctree (expected, not errors).

## Conversion Tool Used

The conversion was performed using `rst-to-myst`:

```bash
pip install "rst-to-myst[sphinx]"
rst2myst convert docs/source/**/*.rst docs/source/index.rst
```

## Commands Run During Migration

```bash
# 1. Create migration branch
git checkout -b docs-myst-migration

# 2. Update conf.py with source_suffix mapping
# (Added source_suffix = {".rst": "restructuredtext", ".md": "markdown"})

# 3. Install conversion tool
pip install "rst-to-myst[sphinx]"

# 4. Convert all RST files
rst2myst convert docs/source/**/*.rst docs/source/index.rst

# 5. Remove old RST files
rm docs/source/index.rst docs/source/getting_started/*.rst \
   docs/source/usage/*.rst docs/source/examples/*.rst \
   docs/source/reference/*.rst

# 6. Fix post-conversion issues
# - Replace ```tsv and ```fasta with ```text
# - Update .rst cross-references to .md

# 7. Verify build
cd docs && sphinx-build -b html source build/html
```

## Future Improvements

1. Convert remaining `{eval-rst}` blocks to native MyST Markdown tables
2. Consider adding `myst_parser` extensions like `substitution`, `tasklist` if needed
3. Remove legacy `*_old.md` files after confirming they're not needed
