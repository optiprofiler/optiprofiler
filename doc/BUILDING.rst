Building the documentation
==========================

From the repository root, use a Python environment with this checkout's
``doc`` dependencies (``python -m pip install -e ".[doc]"``) and initialize
the bundled problem-library submodules required by its API pages. Then run::

    python -m sphinx -b html -E -W doc/source /absolute/path/to/html-output

The configuration selects ``python/`` relative to its own location before
importing OptiProfiler; it does not depend on the shell's current directory.
If another OptiProfiler checkout was already imported in the Python process,
start Sphinx in a fresh process rather than mixing two source trees.

Source-link provenance
----------------------

Python API ``[source]`` links use the exact Git HEAD of this checkout, not the
package version number or a moving branch name. Bundled S2MPJ links use the
provider's own HEAD. The bundled PyCUTEst and SOLAR links follow the same rule.
The displayed package version is unchanged.

Build published documentation from clean, committed core and provider trees
whose commits are available in their corresponding GitHub repositories.
An uncommitted source edit is not represented by its HEAD link.

If Git or the relevant repository metadata is unavailable, the corresponding
source links are omitted; the build does not guess ``main``, a release tag, or
an enclosing repository's revision. For an archive build, a trusted build
operator who has verified the archive's exact source commits may explicitly
supply these environment variables:

* ``OPTIPROFILER_DOCS_SOURCE_REF``: core repository commit.
* ``OPTIPROFILER_DOCS_S2MPJ_REF``: bundled S2MPJ repository commit.
* ``OPTIPROFILER_DOCS_PYCUTEST_REF``: bundled PyCUTEst repository commit.
* ``OPTIPROFILER_DOCS_SOLAR_REF``: bundled SOLAR repository commit.

Each value must be a full 40-character hexadecimal commit SHA, not a branch,
tag, URL, or abbreviated SHA. If Git can identify the checkout, a disagreeing
override fails the build. In a Gitless archive the override is an explicit
provenance assertion by the build operator, not independent verification of
the archive; do not populate it from untrusted input.

Focused configuration tests
---------------------------

With Python and Git available, run from the repository root::

    python -m unittest discover -s doc/tests -v

These tests load the actual configuration against tiny disposable repositories.
They require neither Sphinx nor numerical solvers and do not modify the project
repository's commits or refs.
