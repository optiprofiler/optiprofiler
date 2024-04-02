# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import inspect
import re
import sys
from datetime import datetime
from pathlib import Path

import optiprofiler

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OptiProfiler'
author = 'Cunxin Huang, Tom M. Ragonneau, and Zaikun Zhang'
copyright = f'{2023}\u2013{datetime.now().year}, {author}'

# Short version (including .devX, rcX, b1 suffixes if present).
version = re.sub(r'(\d+\.\d+)\.\d+(.*)', r'\1\2', optiprofiler.__version__)
version = re.sub(r'(\.dev\d+).*?$', r'\1', version)

# Full version, including alpha/beta/rc tags.
release = optiprofiler.__version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.linkcode',
    'numpydoc',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
]

templates_path = ['_templates']

exclude_patterns = []

today_fmt = '%B %d, %Y'

default_role = 'autolink'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'

html_static_path = ['_static']

html_theme_options = {
    'repository_url': 'https://github.com/optiprofiler/optiprofiler',
    'repository_branch': 'main',
    'path_to_docs': 'doc/source',
    'use_repository_button': True,
    'use_source_button': True,
    'use_issues_button': True,
    'use_download_button': False,
    'max_navbar_depth': 2,
}

html_title = f'{project} v{version}'

htmlhelp_basename = project


# -- Options for LaTeX output ------------------------------------------------

latex_documents = [
    ('index', 'optiprofiler.tex', 'OptiProfiler Manual', author, 'manual'),
]

latex_elements = {
    'papersize': 'a4paper',
    'preamble': r'''
% Increase the default table of content depth.
\setcounter{tocdepth}{1}
'''
}


# -- Generate autodoc summaries ----------------------------------------------

autosummary_generate = True


# -- Link to other projects' documentation ------------------------------------

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    'python': ('https://docs.python.org/3/', None),
}


# -- BibTeX citations ---------------------------------------------------------

bibtex_bibfiles = ['_static/optiprofiler.bib']

bibtex_encoding = 'latin'

bibtex_default_style = 'plain'

bibtex_bibliography_header = '''.. only:: html or text

    .. rubric:: References
'''

bibtex_footbibliography_header = bibtex_bibliography_header


# -- Add external links to source code ----------------------------------------

def linkcode_resolve(domain, info):
    if domain != 'py':
        return None

    # Get the object indicated by the module name.
    obj = sys.modules.get(info['module'])
    if obj is None:
        return None
    for part in info['fullname'].split('.'):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None

    # Strip the decorators of the object.
    try:
        unwrap = inspect.unwrap
    except AttributeError:
        pass
    else:
        obj = unwrap(obj)

    # Get the relative path to the source of the object.
    try:
        fn = Path(inspect.getsourcefile(obj)).resolve(True)
    except TypeError:
        return None
    else:
        fn = fn.relative_to(Path(optiprofiler.__file__).resolve(True).parent)

    # Ignore re-exports as their source files are not within the repository.
    module = inspect.getmodule(obj)
    if module is not None and not module.__name__.startswith('optiprofiler'):
        return None

    # Get the line span of the object in the source file.
    try:
        source, lineno = inspect.getsourcelines(obj)
        lines = f'#L{lineno}-L{lineno + len(source) - 1}'
    except OSError:
        lines = ''

    repository = 'https://github.com/optiprofiler/optiprofiler'
    if 'dev' in release:
        return f'{repository}/blob/main/python/optiprofiler/{fn}{lines}'
    else:
        return f'{repository}/blob/v{release}/python/optiprofiler/{fn}{lines}'
