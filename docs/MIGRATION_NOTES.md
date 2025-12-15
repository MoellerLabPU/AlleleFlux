# Documentation Migration Notes: RST to MyST Markdown

This document describes the migration of AlleleFlux documentation from reStructuredText (RST) to MyST Markdown.

## Summary of Changes

### What Changed

1. **Documentation Format**: All `.rst` files in `docs/source/` were converted to MyST Markdown (`.md`)
2. **Sphinx Configuration**: Updated `conf.py` to recognize both `.rst` and `.md` source files
3. **Theme**: Migrated from `sphinx_rtd_theme` to **Furo** for a modern, clean appearance
4. **Extensions**: Added `sphinx-copybutton` for code block copy functionality
5. **Tables**: Converted all RST `list-table` directives to native Markdown tables

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

These contain legacy content and can be removed if no longer needed.

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

````markdown
```{toctree}
:maxdepth: 2

getting_started/overview
getting_started/installation
```
````

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

### Tables

Use standard Markdown tables:

```markdown
| Column 1 | Column 2 | Column 3 |
|----------|----------|----------|
| Value 1  | Value 2  | Value 3  |
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

1. **Build Warnings**: 5 warnings remain about legacy files not in toctree (expected, not errors)
2. **Definition Lists**: Require the `deflist` MyST extension (already enabled)

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

# 2. Update conf.py with source_suffix mapping and Furo theme
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
# - Convert {eval-rst} list-tables to Markdown tables

# 7. Install Furo theme and dependencies
pip install furo sphinx-copybutton

# 8. Verify build
cd docs && sphinx-build -b html source build/html
```

## Theme Configuration (Furo)

The Furo theme is configured in `conf.py`:

```python
html_theme = "furo"
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#2980B9",
        "color-brand-content": "#2980B9",
    },
    "dark_css_variables": {
        "color-brand-primary": "#56B4E9",
        "color-brand-content": "#56B4E9",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
}
```

## Dependencies

Current docs dependencies (`docs/requirements.txt`):

```
sphinx>=7.0.0
furo>=2024.1.0
sphinx-autodoc-typehints>=1.23.0
myst-parser>=2.0.0
sphinx-autobuild
sphinx-copybutton
```

## Future Improvements

1. Remove legacy `*_old.md` files after confirming they're not needed
2. Add more MyST extensions if needed (e.g., `substitution`, `tasklist`)
3. Consider adding custom CSS for additional styling
