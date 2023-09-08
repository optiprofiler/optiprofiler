# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import inspect
import re
import sys
from datetime import datetime
from pathlib import Path

import OptiProfiler

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OptiProfiler'
author = 'Cunxin Huang, Tom M. Ragonneau, and Zaikun Zhang'
start_year = 2023
year = datetime.now().year
if year == start_year:
    copyright = f'{year}, {author}'
else:
    copyright = f'{start_year}\u2013{year}, {author}'

# Short version (including .devX, rcX, b1 suffixes if present).
version = re.sub(r'(\d+\.\d+)\.\d+(.*)', r'\1\2', OptiProfiler.__version__)
version = re.sub(r'(\.dev\d+).*?$', r'\1', version)

# Full version, including alpha/beta/rc tags.
release = OptiProfiler.__version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.doctest',
    'sphinx.ext.linkcode',
    'numpydoc',
    'sphinx_copybutton',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']

exclude_patterns = []

today_fmt = '%B %d, %Y'

default_role = 'autolink'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

# html_static_path = ['_static']

html_theme_options = {
    'navigation_depth': 2,
}

html_context = {
    'github_user': 'OptiProfiler',
    'github_repo': 'OptiProfiler',
    'github_version': 'main',
}

html_title = f'{project} v{version} Manual'

htmlhelp_basename = project


# -- Generate autodoc summaries ----------------------------------------------

autosummary_generate = True


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
        fn = fn.relative_to(Path(OptiProfiler.__file__).resolve(True).parent)

    # Ignore re-exports as their source files are not within the repository.
    module = inspect.getmodule(obj)
    if module is not None and not module.__name__.startswith('OptiProfiler'):
        return None

    # Get the line span of the object in the source file.
    try:
        source, lineno = inspect.getsourcelines(obj)
        lines = f'#L{lineno}-L{lineno + len(source) - 1}'
    except OSError:
        lines = ''

    repository = f'https://github.com/{html_context["github_user"]}/{html_context["github_repo"]}'
    if 'dev' in release:
        return f'{repository}/blob/{html_context["github_version"]}/python/OptiProfiler/{fn}{lines}'
    else:
        return f'{repository}/blob/v{release}/python/OptiProfiler/{fn}{lines}'
