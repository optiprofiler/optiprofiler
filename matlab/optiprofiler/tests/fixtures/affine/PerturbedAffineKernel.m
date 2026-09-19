classdef PerturbedAffineKernel < optiprofiler_internal.FeatureKernel
%PERTURBEDAFFINEKERNEL Test fixture: a feature kernel whose affine pair is perturbed after it is built.
%   K = PerturbedAffineKernel(NAME, OPTIONS, PERTURB) behaves as
%   optiprofiler_internal.FeatureKernel(NAME, OPTIONS), except that
%   modifier_affine returns [A, b, inv] after [A, inv] = PERTURB(A, inv).
%
%   The framework builds both matrices of linearly_transformed itself, so
%   roundoff in its inverse cannot be requested through an option. This class
%   perturbs the pair where it is produced. The code under test, which decides
%   from the pair how bounds and linear constraints are represented
%   (modifier_bounds, modifier_linear_ub, modifier_linear_eq), is inherited
%   unchanged and reads the perturbed pair through ordinary method dispatch.
    properties (SetAccess = private)
        perturb
    end
    methods
        function obj = PerturbedAffineKernel(name, options, perturb)
            obj@optiprofiler_internal.FeatureKernel(name, options);
            obj.perturb = perturb;
        end
        function [A, b, inv] = modifier_affine(obj, seed, problem)
            [A, b, inv] = modifier_affine@optiprofiler_internal.FeatureKernel(obj, seed, problem);
            [A, inv] = obj.perturb(A, inv);
        end
    end
end
