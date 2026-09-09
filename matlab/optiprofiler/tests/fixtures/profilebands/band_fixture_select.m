function names = band_fixture_select(~)
%BAND_FIXTURE_SELECT Three-problem selector of the rendered-band test library (dimensions 2, 3 and 4).
%   The selector ignores its options on purpose: the framework's own problem_names filter selects, as for a real
%   library. Three problems of different dimension give the profile curves several jumps at pooled duplicate x,
%   the shape on which the legacy band polygon was rendered with spurious area.
    names = {'q2', 'q3', 'q4'};
end
