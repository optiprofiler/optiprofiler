function seed = deriveFeatureStageSeed(run_seed, stage_code, occurrence, channel_tag)
%DERIVEFEATURESTAGESEED matlab-stage-horner32-v1, for true multistage only.
% This exact finite identity mixer is not an independence or unrestricted
% collision-free guarantee. Keep the legacy MATLAB numerical RNG kernel.
    words = {run_seed,stage_code,occurrence,channel_tag};
    for k=1:numel(words)
        w=words{k};
        if ~(isnumeric(w) && isreal(w) && isscalar(w) && isfinite(w) ...
                && w>=0 && w<2^32 && w==floor(w))
            error('MATLAB:Feature:InvalidSeedWord','Composite seed words must be finite uint32-range integers.');
        end
    end
    if ~ismember(stage_code,[1,2,3,4,5,6,7,8,9,10]) || ~ismember(channel_tag,[0,1,2,3])
        error('MATLAB:Feature:InvalidSeedIdentity','Unknown frozen stage code or channel tag.');
    end
    seed=double(run_seed);
    for word=[1,double(stage_code),double(occurrence),double(channel_tag)]
        seed=mod(65599*seed+word,2^32);
    end
end
