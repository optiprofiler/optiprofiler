# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import re
from datetime import datetime

import OptiProfiler

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OptiProfiler'
author = 'Cunxin Huang, Tom M. Ragonneau, and Zaikun Zhang'
copyright = f'{datetime.now().year}, {author}'

# Short version (including .devX, rcX, b1 suffixes if present).
version = re.sub(r'(\d+\.\d+)\.\d+(.*)', r'\1\2', OptiProfiler.__version__)
version = re.sub(r'(\.dev\d+).*?$', r'\1', version)

# Full version, including alpha/beta/rc tags.
release = OptiProfiler.__version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'numpydoc',
    'sphinx_copybutton',
]

templates_path = ['_templates']

exclude_patterns = []

today_fmt = '%B %d, %Y'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'

html_static_path = ['_static']


# -- Generate autodoc summaries ----------------------------------------------

autosummary_generate = True
