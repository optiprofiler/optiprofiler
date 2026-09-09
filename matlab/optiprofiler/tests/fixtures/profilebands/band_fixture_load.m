function problem = band_fixture_load(name)
%BAND_FIXTURE_LOAD Deterministic smooth-plus-nonsmooth quadratic of dimension 2, 3 or 4 for the rendered-band tests.
    dims = struct('q2', 2, 'q3', 3, 'q4', 4);
    n = dims.(name);
    shift = (1:n)';
    problem = Problem(struct('fun', @(x) sum((x(:) - shift).^2) + 0.1 * sum(abs(x(:))), 'x0', zeros(n, 1) + 0.5, 'name', name));
end
