classdef RawReferenceProblem < Problem
%RAWREFERENCEPROBLEM Test fixture: a Problem whose stored reference skips validation.
%   P = RawReferenceProblem(S, RAW) builds Problem(S) and then stores RAW as
%   the reference record through the protected setter, exactly as loading a
%   file does. The Problem constructor itself never stores an unvalidated
%   record; this class exists only so that tests can put one into an object
%   (an unknown mapping token written by another version, a non-finite merit,
%   the superseded five-field layout) and pin that it reads as unknown.
    methods
        function obj = RawReferenceProblem(s, raw)
            obj@Problem(s);
            obj.reference = raw;
        end
    end
end
