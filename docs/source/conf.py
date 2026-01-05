# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'completeness'
copyright = '2026, Quentin Groom'
author = 'Quentin Groom'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

import os
import sys
from datetime import datetime

# Make the repo root importable so autodoc can import completeness.py
sys.path.insert(0, os.path.abspath("../.."))

project = "completeness"
author = "Agentschap Plantentuin Meise"
copyright = f"{datetime.now().year}, {author}"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",      # parses Google/Numpy docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "myst_parser",              # allows Markdown pages
    "sphinx_autodoc_typehints",
]

autosummary_generate = True

html_theme = "sphinx_rtd_theme"
