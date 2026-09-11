classdef TestFeatureNativeV2 < matlab.unittest.TestCase
% Trusted native transport and explicit import of the frozen ac4 layout.
    methods (Test)
        function legacyImportAndCurrentTransportKeepOwnership(testCase)
            here = fileparts(mfilename('fullpath'));
            fixture = fullfile(here,'..','fixtures','feature-v2','native-legacy-enum-feature.mat');
            before = readBytes(fixture);
            restored = load(fixture);
            testCase.verifyClass(restored.enum_n_runs,'FeatureOptionKey');
            testCase.verifyEqual(restored.enum_n_runs.value,'n_runs');
            visible = enumeration('FeatureOptionKey');
            testCase.verifyFalse(any(arrayfun(@(v) strcmp(v.value,'n_runs'),visible)));
            testCase.verifyClass(restored.legacy_feature,'optiprofiler_internal.LegacyFeatureEnvelope');
            [feature, request] = optiprofiler_internal.importLegacyFeature(restored.legacy_feature);
            testCase.verifyClass(feature,'Feature');
            testCase.verifyEqual(feature.name,'noisy');
            testCase.verifyEqual(feature.stages{1}.options.noise_level,.001);
            testCase.verifyFalse(isfield(feature.stages{1}.options,'n_runs'));
            testCase.verifyEqual(request.n_runs,3);
            testCase.verifyEqual(request.origin,'retained_resolved_value');
            testCase.verifyEqual(request.requested_origin,'unknown');
            testCase.verifyEmpty(feature.declared);
            testCase.verifyEmpty(feature.declared_name);
            testCase.verifyEqual(restored.native_callback,@sin);
            testCase.verifyEqual(readBytes(fixture),before);

            feature = Feature({'plain',struct('name','noisy','options', ...
                struct('noise_mode','deterministic','noise_map',@sum)),'noisy'});
            native = feature.saveobj();
            testCase.verifyEqual(native.schema,'matlab-feature-native-v2');
            testCase.verifyFalse(isfield(native.specification,'n_runs'));
            file = [tempname,'.mat']; cleanup = onCleanup(@() deleteIfPresent(file));
            save(file,'feature','-v7');
            roundtrip = load(file);
            testCase.verifyClass(roundtrip.feature,'Feature');
            testCase.verifyEqual(roundtrip.feature.stages,feature.stages);
            testCase.verifyEqual(roundtrip.feature.declared,feature.declared);
            testCase.verifyEqual(roundtrip.feature.stages{1}.options.noise_map,@sum);
            testCase.verifyError(@() Feature.loadobj(struct('schema','unsupported')), ...
                'MATLAB:Feature:InvalidNativeState');
        end

        function nestedEnumsAndIdentityUnknownDeclarationRoundtrip(testCase)
            here = fileparts(mfilename('fullpath'));
            fixture = fullfile(here,'..','fixtures','feature-v2','native-legacy-enum-feature.mat');
            restored = load(fixture);
            testCase.verifyClass(restored.enum_as_key_record.key,'FeatureOptionKey');
            testCase.verifyEqual(restored.enum_as_key_record.key.value,'n_runs');
            testCase.verifyEqual(restored.enum_as_key_record.value,3);
            testCase.verifyClass(restored.enum_as_value_record.value,'FeatureOptionKey');
            testCase.verifyEqual(restored.enum_as_value_record.value.value,'n_runs');
            identity = Feature({'plain','plain'});
            [imported, request] = optiprofiler_internal.importLegacyFeature(restored.legacy_feature);
            file = [tempname,'.mat']; cleanup = onCleanup(@() deleteIfPresent(file));
            save(file,'identity','imported','-v7');
            roundtrip = load(file);
            testCase.verifyTrue(roundtrip.identity.is_identity);
            testCase.verifyEmpty(roundtrip.identity.stages);
            testCase.verifyEqual(roundtrip.identity.declared,identity.declared);
            testCase.verifyEqual(roundtrip.identity.declared_name,'plain+plain');
            testCase.verifyEqual(roundtrip.imported.stages,imported.stages);
            testCase.verifyEmpty(roundtrip.imported.declared);
            testCase.verifyEmpty(roundtrip.imported.declared_name);
            testCase.verifyFalse(isfield(roundtrip.imported.stages{1}.options,'n_runs'));
            testCase.verifyEqual(request.n_runs,3);
        end
    end
end

function value = readBytes(file)
    fid = fopen(file,'rb'); assert(fid>=0,'Fixture cannot be read.');
    cleanup = onCleanup(@() fclose(fid));
    value = fread(fid,'*uint8');
end
function deleteIfPresent(file)
    if isfile(file), delete(file); end
end
