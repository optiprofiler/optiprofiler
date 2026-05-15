:html_theme.sidebar_secondary.remove:

.. toctree::
    :maxdepth: 2
    :hidden:

    User guide <user/index>
    Python API reference <ref/index>
    MATLAB API reference <matlab/index>


.. raw:: html

    <style>
        article.bd-article > section#optiprofiler > h1:first-child {display:none;}
    </style>

OptiProfiler
============

.. div:: hero-section

    .. raw:: html

        <div class="hero-grid">
            <div class="hero-copy">
                <div class="hero-logo-wrap">
                    <img src="_static/OP_logo.png" alt="OptiProfiler" class="hero-logo" />
                </div>
                <div class="hero-kicker">
                    Open-source benchmarking for optimization
                </div>
                <h1 class="hero-title">Benchmark your <span class="gradient-text">optimization solver</span></h1>
                <p class="hero-tagline">OptiProfiler standardizes solver benchmarking with curated problem libraries, shared interfaces, baseline comparisons, and reproducible performance profiles.</p>
                <div class="hero-actions">
                    <a class="hero-button hero-button-primary" href="user/installation.html#install">Get started</a>
                    <a class="hero-button hero-button-secondary" href="https://github.com/optiprofiler/optiprofiler" target="_blank" rel="noopener noreferrer">View on GitHub</a>
                </div>
            </div>
            <div class="hero-panel" aria-label="Quick start example">
                <div class="hero-panel-bar">
                    <span class="code-dot code-dot-red"></span>
                    <span class="code-dot code-dot-yellow"></span>
                    <span class="code-dot code-dot-green"></span>
                    <span class="hero-panel-title">quick_start.py</span>
                </div>
                <pre><code><span class="code-muted">from scipy.optimize import minimize</span>
        from optiprofiler import benchmark

        def bfgs(fun, x0):
            return minimize(
                fun, x0, method="BFGS").x

        def powell(fun, x0):
            return minimize(
                fun, x0, method="Powell").x

        scores = benchmark([bfgs, powell])
        print(scores)</code></pre>
            </div>
        </div>
        <div class="hero-metrics" aria-label="OptiProfiler capabilities">
            <div>
                <strong>Python + MATLAB</strong>
                <span>Benchmark solvers through either interface with the same workflow</span>
            </div>
            <div>
                <strong>Problem libraries</strong>
                <span>Use CUTEst, S2MPJ, custom problems, and diverse features</span>
            </div>
            <div>
                <strong>Reproducible outputs</strong>
                <span>Export profiles, scores, logs, and result files for reporting</span>
            </div>
        </div>

----

.. div:: section-eyebrow

    Two ways to use OptiProfiler

.. div:: section-heading

    Run it locally, or in the cloud

.. grid:: 1 1 2 2
    :gutter: 3

    .. grid-item-card::
        :class-card: use-card use-card-primary

        .. raw:: html

            <div class="use-eyebrow">Run locally</div>
            <h3>Use Python &amp; MATLAB packages</h3>
            <p>Install the open-source packages, run larger experiments, plug in your own problem libraries, and keep full control over private problems and compute.</p>

        .. raw:: html

            <div class="use-card-cta">

        :ref:`Read the installation guide → <install>`

        .. raw:: html

            </div>

    .. grid-item-card::
        :class-card: use-card

        .. raw:: html

            <div class="use-eyebrow">Run in the cloud</div>
            <h3>Try the hosted platform</h3>
            <p>Upload a Python solver, choose benchmark settings, compare against baseline solvers, and get publication-ready profiles with downloadable results in your browser.</p>

        .. raw:: html

            <div class="use-card-cta">
            <a href="https://app.optprof.com" target="_blank" rel="noopener noreferrer">Open app.optprof.com →</a>
            </div>

----

.. div:: section-heading

    Why OptiProfiler?

.. raw:: html

    <div class="feature-list">

        <div class="feature-item">
            <div class="feature-icon" aria-hidden="true">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><path d="M9.937 15.5A2 2 0 0 0 8.5 14.063l-6.135-1.582a.5.5 0 0 1 0-.962L8.5 9.937A2 2 0 0 0 9.937 8.5l1.582-6.135a.5.5 0 0 1 .963 0L14.063 8.5A2 2 0 0 0 15.5 9.937l6.135 1.582a.5.5 0 0 1 0 .962L15.5 14.063a2 2 0 0 0-1.437 1.437l-1.582 6.135a.5.5 0 0 1-.963 0z"/><path d="M20 3v4"/><path d="M22 5h-4"/><path d="M4 17v2"/><path d="M5 18H3"/></svg>
            </div>
            <div class="feature-content">
                <div class="feature-title"><span class="feature-keyword">Simple</span> usage for beginners</div>
                <div class="feature-desc">Easy <a href="user/installation.html">installation</a> and <a href="user/usage.html">quick start</a> with a few lines of code</div>
            </div>
        </div>

        <div class="feature-item">
            <div class="feature-icon" aria-hidden="true">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><line x1="21" y1="4" x2="14" y2="4"/><line x1="10" y1="4" x2="3" y2="4"/><line x1="21" y1="12" x2="12" y2="12"/><line x1="8" y1="12" x2="3" y2="12"/><line x1="21" y1="20" x2="16" y2="20"/><line x1="12" y1="20" x2="3" y2="20"/><line x1="14" y1="2" x2="14" y2="6"/><line x1="8" y1="10" x2="8" y2="14"/><line x1="16" y1="18" x2="16" y2="22"/></svg>
            </div>
            <div class="feature-content">
                <div class="feature-title"><span class="feature-keyword">Multiple</span> degrees of freedom for experts</div>
                <div class="feature-desc">Multiple built-in features and customization options for test suites</div>
            </div>
        </div>

        <div class="feature-item">
            <div class="feature-icon" aria-hidden="true">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><path d="M15 4V2"/><path d="M15 16v-2"/><path d="M8 9h2"/><path d="M20 9h2"/><path d="m17.8 11.8 1.4 1.4"/><path d="m17.8 6.2 1.4-1.4"/><path d="m3 21 9-9"/><path d="m12.2 6.2-1.4-1.4"/></svg>
            </div>
            <div class="feature-content">
                <div class="feature-title"><span class="feature-keyword">Automatic</span> generation of high-quality profiles</div>
                <div class="feature-desc">Publication-ready PDF visualizations with precise, clear, and aesthetically pleasing figures</div>
            </div>
        </div>

        <div class="feature-item">
            <div class="feature-icon" aria-hidden="true">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><path d="M20 13c0 5-3.5 7.5-7.66 8.95a1 1 0 0 1-.67-.01C7.5 20.5 4 18 4 13V6a1 1 0 0 1 1-1c2 0 4.5-1.2 6.24-2.72a1.17 1.17 0 0 1 1.52 0C14.51 3.81 17 5 19 5a1 1 0 0 1 1 1z"/><path d="m9 12 2 2 4-4"/></svg>
            </div>
            <div class="feature-content">
                <div class="feature-title"><span class="feature-keyword">Reliable</span> methodology for benchmarking</div>
                <div class="feature-desc">Based on the widely accepted performance profile and data profile</div>
            </div>
        </div>

        <div class="feature-item">
            <div class="feature-icon" aria-hidden="true">
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8"/><path d="M3 3v5h5"/><path d="M12 7v5l4 2"/></svg>
            </div>
            <div class="feature-content">
                <div class="feature-title"><span class="feature-keyword">Trackable</span> experimental results for reproducibility</div>
                <div class="feature-desc">Seed-controlled and time-stamped experimental results for easy reproducibility</div>
            </div>
        </div>

    </div>

----

.. div:: section-heading

    Multi-language support

.. grid:: 1 1 2 2
    :gutter: 3

    .. grid-item-card::
        :class-card: lang-card

        .. raw:: html

            <div class="lang-title-row">
                <img src="_static/python.svg" alt="Python" class="lang-inline-icon" />
                <strong>Python</strong>
            </div>

        The Python interface is ready to use.
        Check out the :ref:`installation guide <install>` to get started,
        explore :ref:`examples <use_python>` for common workflows,
        or dive into the :ref:`API reference <pythonapi>` for details.

    .. grid-item-card::
        :class-card: lang-card

        .. raw:: html

            <div class="lang-title-row">
                <img src="_static/matlab.png" alt="MATLAB" class="lang-inline-icon" />
                <strong>MATLAB</strong>
            </div>

        The full MATLAB interface is ready to use.
        Check out the :ref:`installation guide <install>` to get started,
        explore :ref:`examples <use>` for common workflows,
        or dive into the :ref:`API reference <matlabapi>` for details.
