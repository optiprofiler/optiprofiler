MATLAB output robustness
------------------------

Normal MATLAB sessions use native figures and PDF output. The ``-batch`` and
``-nodesktop`` flags do not disable the Java Virtual Machine (JVM).
On macOS and Linux, OptiProfiler also supports saving without a JVM:

* MATLAB releases with working native graphics (verified with R2026a)
  retain the native PDF route.
* If MATLAB cannot create a figure, as in R2023b with ``-nojvm``, OptiProfiler
  writes plain, browser-readable SVG charts and ``summary.html`` instead of
  PDF/FIG. These use the same computed profile curves and history
  preprocessing; they do not promise native styling, LaTeX interpretation,
  or interactive FIG files. Mean and errorband values remain the same; the
  portable history chart labels errorband boundaries as separate lines.
* A failing native PDF exporter can also use SVG as its last fallback. If
  only PDF merging fails, existing individual plots remain available through
  an HTML index and the complete backend diagnostics are saved.

The initial native-graphics probe is cached for the MATLAB process, avoiding
a new probe figure for every plot. A transient failure of that initial probe
therefore keeps subsequent plots on the SVG route until MATLAB is restarted.

On Linux with MATLAB R2026a, a separate safety check avoids a reproducible
native-export hang in the embedded Chromium renderer. Its observed temporary
socket name adds 46 bytes to ``TMPDIR``, while a Linux Unix-domain socket path
allows at most 107 bytes. An explicit temporary path longer than 61 UTF-8
bytes (ignoring trailing separators) therefore uses the existing SVG/HTML
route instead of native PDF/FIG.
Relative ``TMPDIR`` values also use this conservative fallback; use a short
absolute path for native output. This does not claim every relative path
fails in MATLAB. An unset ``TMPDIR`` uses the normal system default. An explicitly
empty ``TMPDIR`` also reproduced the hang, so it uses the portable route instead.
When MATLAB reports an empty value, a small read-only shell check distinguishes
unset from set-empty without requiring Java; an unsuccessful check also falls
back with an explicit diagnostic.

The message explains the fallback. To request native output, set a short
absolute temporary path, for example ``setenv('TMPDIR', '/tmp')``, before
calling ``benchmark``. OptiProfiler does not change that environment variable
or move experiment outputs itself. The path check is repeated on each graphics
capability query, so correcting the temporary path does not permanently retain
the unsafe-path fallback. The guard is specific to the tested Linux R2026a
backend; it is not a general timeout for arbitrary native graphics failures.
It also cannot help when an even longer temporary path prevents MATLAB itself
from starting; shorten that path before launching MATLAB.
When upgrading MATLAB, recheck the native exporter at the 61/62-byte boundary
and with unset versus set-empty ``TMPDIR`` before extending this release guard.

The SVG axes explicitly label any logarithmic coordinate transformation.
The display-only history limit and nonfinite replacements described above
are annotated; raw saved values and scores are never clipped for rendering.
``score_only=true`` neither allocates figures nor writes fallback artifacts,
even if summary options are enabled.

No-JVM file safety is implemented with exclusive POSIX directory creation and
a same-directory temporary file followed by atomic rename on macOS/Linux.
The JVM path retains Java's atomic file operations. Windows no-JVM output is
not supported: OptiProfiler reports an error asking you to restart MATLAB
without ``-nojvm``, rather than using an unsafe replacement operation. Normal
MATLAB startup enables the JVM; use ``usejava('jvm')`` to check your session.
The restriction applies to saving output, not ``score_only=true`` evaluation.
External PDF mergers (qpdf, pdfunite, or Ghostscript)
are optional fallbacks, not silently installed dependencies.

The full MATLAB regression suite additionally requires ``pdfunite`` and
``pdftotext`` (Ubuntu package ``poppler-utils``) to test actual PDF merging,
page order, and preservation of an existing summary when a merge fails.
Its CI job installs and checks these tools explicitly rather than relying
on MATLAB's bundled Java classes. This is a test prerequisite, not a new
requirement for users: the no-merger HTML fallback is tested separately.
On Linux CI, only these two system-tool subprocesses clear
``LD_LIBRARY_PATH`` so they do not load an incompatible MATLAB-bundled C++
library. MATLAB's own environment and user-installed tool configurations
are not changed by OptiProfiler.
In a normal Linux session, PDF merging first tries each tool with the caller's
environment unchanged. If a tool fails with a specific ``libstdc++.so``
``GLIBCXX`` or ``CXXABI`` required-version diagnostic, OptiProfiler retries
that tool once with ``LD_LIBRARY_PATH`` removed only in its child process.
The first successful isolated retry prints an informational message. Missing
tools, invalid PDFs, and other errors do not trigger this retry. If no merger
succeeds, the warning and saved diagnostics explain the failures; individual
plots and the HTML index remain available, and a previous summary PDF is not
replaced with a failed attempt. No tool is installed automatically, and no
experiment directory, raw data format, or PDF naming convention is changed.

The full-unit CI entry writes JUnit results, runner diagnostics, and coverage
outside disposable test fixtures. CI uploads these receipts even when a test
fails; an empty upload is an error, not evidence of a successful test run.

Completed experiment data is saved as a verified version-7.3 MAT file before
sequential history rendering. A timestamp becomes loadable only after this
save succeeds. Requested saving failures are errors, including in silent
mode; completed staging files are retained if final publication fails.
This is final-result preservation, not a checkpoint/resume mechanism for an
unfinished solver run or a guarantee against abrupt machine/disk failure.

Library-selection failures stop with their original diagnostic instead of
being treated as an empty or unfiltered selection. OptiProfiler restores the
caller's warning settings on success and exceptions, and borrows an existing
parallel pool without deleting or resizing it. A smaller existing pool may
use fewer workers than ``n_jobs``; this changes throughput, not the requested
problems or evaluation budgets.

If setup cannot save the default pathdef and no default startup location
exists, it keeps the installed paths available in the current session and
prints the commands needed for future sessions. It does not claim successful
persistence. This session-only fallback is not used when a persistence target
was explicitly specified. Uninstall removes only paths recorded as newly
added by setup, including the bundled S2MPJ runtime paths; pre-existing paths
are borrowed and preserved.
