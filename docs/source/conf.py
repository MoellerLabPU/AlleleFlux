import os
import sys

# Add project root to path for auto doc generation
sys.path.insert(0, os.path.abspath("../.."))

# Project information
project = "AlleleFlux"
copyright = "2024, Siddhartha Uppal, Andrew Moeller"
author = "Siddhartha Uppal, Andrew Moeller"
version = "0.0.1"
release = "0.0.1"

# Extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_copybutton",
]

# Source suffix configuration for both RST and MyST Markdown
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Templates
templates_path = ["_templates"]
exclude_patterns = []

# HTML output options
html_theme = "furo"
html_static_path = ["_static"]
html_title = "AlleleFlux Documentation"
html_favicon = None  # You can add a favicon later if needed

# Furo theme options
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

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_typehints_format = "short"

# Napoleon settings (for Google/NumPy style docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = True
napoleon_custom_sections = None

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# MyST parser settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3
