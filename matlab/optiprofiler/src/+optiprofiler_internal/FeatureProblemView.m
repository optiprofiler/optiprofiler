classdef FeatureProblemView < Problem
%FEATUREPROBLEMVIEW Ordered oracle view; no solver budget or history.
% The Problem-facing methods expose observations to an immediate downstream
% callback. Reference reads have a separate path and never advance served
% counters. Effective dimensions are metadata, not hidden observed queries.
    properties (Access = private)
        predecessor
        stage
        kernel
        dimensions
        seeds
        affine_matrix
        affine_shift
        is_affine = false
        served = struct('fun',0,'cub',0,'ceq',0)
    end
    methods
        function obj = FeatureProblemView(predecessor,stage,run_seed)
            if nargin<2
                stage = [];
            end
            structure = struct('name',predecessor.name,'x0',predecessor.x0, ...
                'xl',predecessor.xl,'xu',predecessor.xu, ...
                'aub',predecessor.aub,'bub',predecessor.bub, ...
                'aeq',predecessor.aeq,'beq',predecessor.beq,'fun',@(x) NaN);
            dimensions = [predecessor.m_nonlinear_ub,predecessor.m_nonlinear_eq];
            kernel = [];
            seeds = [];
            A = eye(predecessor.n);
            b = zeros(predecessor.n,1);
            affine = false;
            if ~isempty(stage)
                options = stage.options;
                if strcmp(stage.name,'custom')
                    options = optiprofiler_internal.customRuntimeOptions(stage,predecessor);
                end
                % Composed views mix the per-query payload with the word fold
                % of matlab-stage-horner32-v2; construction streams are unchanged.
                kernel = optiprofiler_internal.FeatureKernel(stage.name,options,'horner32-words');
                channels = {'fun','cub','ceq','construction'};
                seeds = struct();
                for k = 1:4
                    seeds.(channels{k}) = optiprofiler_internal.deriveFeatureStageSeed( ...
                        run_seed,stage.code,stage.occurrence,k-1);
                end
                affine = ismember(stage.name,{'permuted','linearly_transformed','custom'});
                if affine
                    [A,b] = kernel.modifier_affine(seeds.construction,predecessor);
                    structure.x0 = kernel.modifier_x0(seeds.construction,predecessor);
                    [structure.xl,structure.xu] = kernel.modifier_bounds(seeds.construction,predecessor);
                    [structure.aub,structure.bub] = kernel.modifier_linear_ub(seeds.construction,predecessor);
                    [structure.aeq,structure.beq] = kernel.modifier_linear_eq(seeds.construction,predecessor);
                elseif strcmp(stage.name,'perturbed_x0')
                    structure.x0 = kernel.modifier_x0(seeds.construction,predecessor);
                end
            end
            obj@Problem(structure);
            obj.predecessor = predecessor;
            obj.stage = stage;
            obj.dimensions = dimensions;
            obj.kernel = kernel;
            obj.seeds = seeds;
            obj.affine_matrix = A;
            obj.affine_shift = b;
            obj.is_affine = affine;
        end
        function value = fun(obj,x)
            value = obj.observed('fun',x);
        end
        function value = cub(obj,x)
            value = obj.observed('cub',x);
        end
        function value = ceq(obj,x)
            value = obj.observed('ceq',x);
        end
        function value = reference(obj,channel,x)
            x = obj.point(x);
            if isempty(obj.stage)
                value = obj.predecessor.(channel)(x);
                if ~strcmp(channel,'fun')
                    value = value(:);
                end
            else
                x = obj.map(x);
                if strcmp(obj.stage.name,'quantized') && obj.stage.options.ground_truth
                    x = obj.quantize(x);
                end
                value = obj.predecessor.reference(channel,x);
            end
        end
        function varargout = maxcv(obj,x,detailed)
            if nargin<3
                detailed = false;
            end
            x = obj.point(x);
            [bounds,linear] = obj.observedStructural(x);
            nonlinear = obj.observedNonlinear(x);
            values = {max([bounds;linear;nonlinear],[],'includenan'),bounds,linear,nonlinear};
            if detailed
                varargout = values;
            else
                varargout = values(1);
            end
        end
        function varargout = referenceMaxcv(obj,x,detailed)
            if nargin<3
                detailed = false;
            end
            x = obj.point(x);
            [bounds,linear] = obj.referenceStructural(x);
            nonlinear = obj.referenceNonlinear(x);
            values = {max([bounds;linear;nonlinear],[],'includenan'),bounds,linear,nonlinear};
            if detailed
                varargout = values;
            else
                varargout = values(1);
            end
        end
        function [bounds,linear] = observedStructural(obj,x)
            bounds = 0;
            linear = 0;
            if any(isfinite(obj.xl))
                bounds = max([obj.xl-x;0],[],'includenan');
            end
            if any(isfinite(obj.xu))
                bounds = max([x-obj.xu;bounds],[],'includenan');
            end
            if ~isempty(obj.aub)
                linear = max([obj.aub*x-obj.bub;0],[],'includenan');
            end
            if ~isempty(obj.aeq)
                linear = max([abs(obj.aeq*x-obj.beq);linear],[],'includenan');
            end
        end
        function value = observedNonlinear(obj,x)
            value = 0;
            if obj.dimensions(1)>0
                value = max([obj.cub(x);0],[],'includenan');
            end
            if obj.dimensions(2)>0
                value = max([abs(obj.ceq(x));value],[],'includenan');
            end
        end
        function [bounds,linear] = referenceStructural(obj,x)
            if isempty(obj.stage)
                [bounds,linear] = obj.observedStructural(x);
            else
                [bounds,linear] = obj.predecessor.referenceStructural(obj.map(x));
            end
        end
        function value = referenceNonlinear(obj,x)
            if isempty(obj.stage)
                value = 0;
                if obj.dimensions(1)>0
                    value = max([obj.reference('cub',x);0],[],'includenan');
                end
                if obj.dimensions(2)>0
                    value = max([abs(obj.reference('ceq',x));value],[],'includenan');
                end
            else
                x = obj.map(x);
                if strcmp(obj.stage.name,'quantized') && obj.stage.options.ground_truth
                    x = obj.quantize(x);
                end
                value = obj.predecessor.referenceNonlinear(x);
            end
        end
        function x = toOriginalCoordinates(obj,x)
            x = obj.point(x);
            if ~isempty(obj.stage)
                x = obj.predecessor.toOriginalCoordinates(obj.map(x));
            end
        end
        function records = runtimeStages(obj)
            if isempty(obj.stage)
                records = {};
                return;
            end
            records = obj.predecessor.runtimeStages();
            records{end+1} = struct('identity',obj.stage.identity,'name',obj.stage.name, ...
                'code',obj.stage.code,'occurrence',obj.stage.occurrence, ...
                'seeds',obj.seeds,'served',obj.served);
        end
        function value = grad(obj,x)
            value = obj.noDerivatives(x);
        end
        function value = hess(obj,x)
            value = obj.noDerivatives(x);
        end
        function value = jcub(obj,x)
            value = obj.noDerivatives(x);
        end
        function value = jceq(obj,x)
            value = obj.noDerivatives(x);
        end
        function value = hcub(obj,x)
            value = obj.noDerivatives(x);
        end
        function value = hceq(obj,x)
            value = obj.noDerivatives(x);
        end
    end
    methods (Access = protected)
        function value = constraintDimension(obj,channel)
            if strcmp(channel,'cub')
                value = obj.dimensions(1);
            else
                value = obj.dimensions(2);
            end
        end
    end
    methods (Access = private)
        function value = observed(obj,channel,x)
            x = obj.point(x);
            if isempty(obj.stage)
                value = obj.predecessor.(channel)(x);
                if ~strcmp(channel,'fun')
                    value = value(:);
                end
                return
            end
            index = obj.served.(channel);
            obj.served.(channel) = index+1;
            x = obj.map(x);
            if strcmp(obj.stage.name,'quantized')
                % A local mesh reads its predecessor once. It is not the
                % legacy kernel's eager unsnapped pre-read plus snapped read.
                value = obj.predecessor.(channel)(obj.quantize(x));
                return
            end
            if strcmp(obj.stage.name,'unrelaxable_constraints') && strcmp(channel,'fun')
                value = obj.predecessor.fun(x);
                [bounds,linear] = obj.predecessor.observedStructural(x);
                options = obj.stage.options;
                if (options.unrelaxable_bounds && bounds>0) || ...
                        (options.unrelaxable_linear_constraints && linear>0)
                    value = Inf;
                elseif options.unrelaxable_nonlinear_constraints && obj.predecessor.observedNonlinear(x)>0
                    value = Inf;
                end
                return
            end
            method = ['modifier_',channel];
            value = obj.kernel.(method)(x,obj.seeds.(channel),obj.predecessor,index);
        end
        function x = quantize(obj,x)
            % Same MATLAB mesh formula as the retained one-stage kernel.
            mesh = obj.stage.options.mesh_size;
            if strcmp(obj.stage.options.mesh_type,'relative')
                mesh = mesh.*max(1,abs(x));
            end
            x = mesh.*round(x./mesh);
        end
        function x = map(obj,x)
            if obj.is_affine
                x = obj.affine_matrix*x+obj.affine_shift;
            end
        end
        function x = point(obj,x)
            if ~(isnumeric(x) && isreal(x) && isvector(x) && numel(x)==obj.n)
                error('MATLAB:FeaturedProblem:InvalidPoint','A real vector of length %d is required.',obj.n);
            end
            x = x(:);
        end
        function value = noDerivatives(~,~) %#ok<STOUT>
            error('MATLAB:FeaturedProblem:UnsupportedCompositeDerivative', ...
                'Derivatives of composite feature views are not supported.');
        end
    end
end
