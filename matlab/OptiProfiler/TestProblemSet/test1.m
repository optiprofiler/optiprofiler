function test1
    clc

    pb = loadTestProblem('pb2');
    % s = struct('fun', @(x) x.' * x, 'x0', 1);
    % pb = Problem(s);

    feature = Feature('noisy');
    
    ftpb = FeaturedProblem(pb, feature);

    ftpb.fun([1;1])
end