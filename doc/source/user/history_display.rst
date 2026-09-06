History plots with extreme values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

History plots use a **display copy** of the recorded values. In both Python
and MATLAB, finite values outside ``[-1e100, 1e100]`` are clipped to that
interval **before** calculating the displayed means and error bands. This
large plotting limit leaves room for squared deviations, translations and
axis margins without floating-point overflow. It is not an assumption that
all valid objective values lie in that interval.

Each affected panel states the display limit and the number of clipped input
history entries. NaN and Inf entries are counted separately and displayed
above the finite range of the corresponding run (including its finite
initial value). If that run has no finite history or initial value, their
display placeholder is 1. A constant finite range gets a small positive gap
so its nonfinite placeholders remain distinguishable. These placeholders are
not successful function evaluations.

The clipping occurs before both raw-view and cumulative-minimum-view
statistics. The term "raw" in a plot filename denotes the non-cumulative
view, not an unmodified serialization of the experiment. Original saved
histories, oracle evaluations, output points, merit values and profile scores
are unchanged. For numerical analysis of extreme values, use the saved data
rather than recovering values from a capped plot.

Positive tiny values are not replaced by a fixed display floor. Nonpositive
histories are shifted when a logarithmic view is needed, and the shift is
shown in the axis label. Logarithmic limits and extreme-range ticks are kept
within representable floating-point values, so a finite history cannot be
silently replaced by a default ``[1, 10]`` viewing range.

For ``errorbar_type='meanstd'``, the two interfaces retain their established
standard-deviation conventions across the ``N`` runs: Python divides the
sum of squared deviations by ``N`` (population normalization), while MATLAB
divides by ``N-1`` when ``N > 1`` (sample normalization). With one run,
both use zero standard deviation. The same convention is used to select the
logarithmic shift and to draw the band within each interface. Therefore,
unequal-run ``meanstd`` bands, and sometimes their displayed shifts, are
not expected to match exactly between Python and MATLAB. This compatibility
choice does not change either interface's saved values or profile scores.
