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

The SVG axes explicitly label any logarithmic coordinate transformation.
The display-only history limit and nonfinite replacements described above
are annotated; raw saved values and scores are never clipped for rendering.
``score_only=true`` neither allocates figures nor writes fallback artifacts,
even if summary options are enabled.

No-JVM file safety is implemented with exclusive POSIX directory creation and
a same-directory temporary file followed by atomic rename on macOS/Linux.
The JVM path retains Java's atomic file operations. Windows no-JVM output is
not certified: enable the JVM there rather than relying on an unverified
replacement operation. External PDF mergers (qpdf, pdfunite, or Ghostscript)
are optional fallbacks, not silently installed dependencies.

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
