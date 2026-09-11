classdef TestFeatureReportBounds < matlab.unittest.TestCase
% Public bounded metadata must not constrain the reusable native Feature.
    methods (Test)
        function identityDeclarationBoundary(testCase)
            for count = [256, 257]
                feature = Feature(strjoin(repmat({'plain'}, 1, count), '+'));
                report = runReport(feature, count);
                value = report.configuration.effective.feature.declared;
                if count == 256
                    testCase.verifyNumElements(value, 256);
                else
                    testCase.verifyEqual(value.total_items, 257);
                    testCase.verifyEqual(value.reason, 'metadata_item_limit');
                    testCase.verifyNumElements(value.values, 256);
                end
                testCase.verifyEmpty(report.configuration.effective.feature.stages);
                testCase.verifyNumElements(feature.declared.entries, count);
                testCase.verifyEmpty(feature.stages);
            end
        end

        function effectiveChainAndRetainedProjectionAreBounded(testCase)
            feature = Feature(strjoin(repmat({'truncated'}, 1, 257), '+'));
            report = runReport(feature, 257);
            block = report.configuration.effective.feature;
            testCase.verifyEqual(block.declared.total_items, 257);
            testCase.verifyEqual(block.stages.total_items, 257);
            testCase.verifyNumElements(block.stages.values, 256);
            testCase.verifyEqual(block.stages.values(end).position, 255);
            testCase.verifyNumElements(feature.stages, 257);
            % The same safe feature projection is used for retained metadata;
            % numerical/problem lists outside feature blocks are not capped.
            nested = struct('feature', struct('effective_name', feature.name, ...
                'declared', {feature.declared.entries}, 'stages', {feature.stages}), ...
                'ordinary_list', {num2cell(1:300)});
            projected = jsondecode(optiprofiler_internal.EvalReport.encodeMetadata(nested));
            testCase.verifyEqual(projected.feature.stages.total_items, 257);
            testCase.verifyNumElements(projected.ordinary_list, 300);
        end
    end
end

function report = runReport(feature, count)
    output = tempname(getenv('OP_ARTIFACTS')); mkdir(output);
    problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', 1, 'name', 'BOUNDED_FEATURE'));
    options = struct('problem', problem, 'feature', feature, 'n_runs', 1, ...
        'n_jobs', 1, 'max_eval_factor', 1, 'max_tol_order', 1, 'seed', 17, ...
        'solver_names', {{'stay', 'zero'}}, 'score_only', true, 'silent', true, ...
        'report_path', fullfile(output, sprintf('feature-%d.json', count)));
    benchmark({@stay, @zero}, options);
    report = jsondecode(fileread(options.report_path));
end

function x = stay(fun, x0)
    fun(x0); x = x0;
end

function x = zero(fun, x0)
    x = zeros(size(x0)); fun(x);
end
