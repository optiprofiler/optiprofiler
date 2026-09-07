# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import inspect
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

# Import this checkout before any installed version, independently of the
# directory from which Sphinx was started.
_SOURCE_ROOT = Path(__file__).resolve().parents[2]
_PACKAGE_ROOT = _SOURCE_ROOT / 'python' / 'optiprofiler'
sys.path.insert(0, str(_SOURCE_ROOT / 'python'))
import optiprofiler

if Path(optiprofiler.__file__).resolve().parent != _PACKAGE_ROOT.resolve():
    raise RuntimeError(
        'The imported optiprofiler does not belong to this documentation checkout. '
        'Start Sphinx in a fresh Python process.')

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
    'sphinx_design',
    'sphinxcontrib.bibtex',
    'sphinxcontrib.matlab',
    'sphinxext.opengraph',
]

# -- Open Graph metadata (controls Google/social media preview images) --------
ogp_site_url = 'https://www.optprof.com/'
ogp_image = '_static/OP_logo.png'
ogp_image_alt = 'OptiProfiler Logo'
ogp_site_name = 'OptiProfiler'
ogp_description_length = 200

# primary_domain = 'mat'

# MATLAB domain configuration.
matlab_src_dir = str(_SOURCE_ROOT / 'matlab' / 'optiprofiler')
matlab_short_links = True

# Disable parallel reading
parallel_read_safe = False

templates_path = ['_templates']

exclude_patterns = ['dev']

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
    'use_issues_button': True,
    'use_download_button': False,
    'max_navbar_depth': 2,
    'home_page_in_toc': False,
    'announcement': (
        'Prefer running in the cloud? Try the hosted platform at '
        '<a href="https://app.optprof.com" target="_blank" '
        'rel="noopener noreferrer">app.optprof.com</a> &rarr;'
    ),
}

html_title = f'{project} v{version}'

html_favicon = '_static/favicon.ico'

html_css_files = [
    'custom.css',
]

html_js_files = [
    'system-theme.js',
]

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

autodoc_mock_imports = ['pycutest']


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

def _source_revision(root, override_name):
    """Resolve an exact source revision, or omit links when it is unknown."""
    override = os.environ.get(override_name)
    if override is not None:
        if not re.fullmatch(r'[0-9a-fA-F]{40}', override):
            raise ValueError(f'{override_name} must be a full 40-character hexadecimal commit SHA.')
        override = override.lower()

    revision = None
    # Git otherwise searches parents: a source archive or uninitialized provider
    # must never inherit an unrelated enclosing repository's revision.
    if (root / '.git').exists():
        env = os.environ.copy()
        for key in ('GIT_DIR', 'GIT_WORK_TREE', 'GIT_COMMON_DIR', 'GIT_INDEX_FILE'):
            env.pop(key, None)
        try:
            def git(*args):
                return subprocess.check_output(
                    ['git', '-C', str(root), 'rev-parse', *args],
                    env=env, text=True, stderr=subprocess.DEVNULL, timeout=5).strip()

            if Path(git('--show-toplevel')).resolve() == root.resolve():
                candidate = git('--verify', 'HEAD')
                if re.fullmatch(r'[0-9a-f]{40}', candidate):
                    revision = candidate
        except (OSError, subprocess.SubprocessError):
            pass  # Gitless builds may explicitly supply their archive's SHA.

    if revision and override and override != revision:
        raise ValueError(f'{override_name} does not match the checked-out revision.')
    return revision or override


# A maintenance branch can retain its release number while receiving fixes.
# Use the checkout's commit, not __version__, to identify the displayed source.
_core_revision = _source_revision(_SOURCE_ROOT, 'OPTIPROFILER_DOCS_SOURCE_REF')

# Bundled providers have their own repository identity, independent of core.
_submodule_sources = {
    'problem_libs/s2mpj/': ('https://github.com/optiprofiler/s2mpj_python', 'OPTIPROFILER_DOCS_S2MPJ_REF'),
    'problem_libs/pycutest/': ('https://github.com/optiprofiler/pycutest', 'OPTIPROFILER_DOCS_PYCUTEST_REF'),
    'problem_libs/solar/': ('https://github.com/optiprofiler/solar_python', 'OPTIPROFILER_DOCS_SOLAR_REF'),
}
_submodule_sources = {
    prefix: (url, _source_revision(_PACKAGE_ROOT / prefix, override))
    for prefix, (url, override) in _submodule_sources.items()
}


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
        fn = fn.relative_to(_PACKAGE_ROOT.resolve())
    except (TypeError, ValueError, OSError):
        return None

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

    fn_str = str(fn).replace('\\', '/')

    for prefix, (repo_url, revision) in _submodule_sources.items():
        if fn_str.startswith(prefix):
            sub_fn = quote(fn_str[len(prefix):], safe='/')
            return f'{repo_url}/blob/{revision}/{sub_fn}{lines}' if revision else None

    repository = 'https://github.com/optiprofiler/optiprofiler'
    if _core_revision:
        return f'{repository}/blob/{_core_revision}/python/optiprofiler/{quote(fn_str, safe="/")}{lines}'
    return None
