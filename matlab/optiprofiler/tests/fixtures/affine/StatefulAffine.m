classdef StatefulAffine < handle
%STATEFULAFFINE Test fixture: a mod_affine callback with a state.
%   C = StatefulAffine(KIND) counts in C.calls how often C.transform is asked
%   and answers according to KIND:
%
%   'alternating'  a valid diagonal transform and a valid dense one in turn;
%   'drifting'     another valid rotation on every call, as when user code
%                  draws from the global stream instead of the one it is handed;
%   'counting'     always the same dense transform.
%
%   C = StatefulAffine(KIND, CALLS) starts the count at CALLS. Use it as
%   @(stream, problem) C.transform(stream, problem). The framework must ask
%   once per problem and seed: a second answer may be a different map, and the
%   bounds, the linear constraints and the evaluations have to share one.
    properties
        calls = 0
    end
    properties (SetAccess = private)
        kind
    end
    methods
        function obj = StatefulAffine(kind, calls)
            obj.kind = kind;
            if nargin > 1
                obj.calls = calls;
            end
        end
        function [A, b, inverse] = transform(obj, ~, ~)
            obj.calls = obj.calls + 1;
            switch obj.kind
                case 'alternating'
                    if mod(obj.calls, 2) == 1
                        [A, b, inverse] = TestAffineStructure.exactDiagonal();
                    else
                        [A, b, inverse] = TestAffineStructure.dense();
                    end
                case 'drifting'
                    [A, b, inverse] = StatefulAffine.rotation(obj.calls);
                otherwise
                    [A, b, inverse] = TestAffineStructure.dense();
            end
        end
    end
    methods (Static)
        function [A, b, inverse] = rotation(call)
            c = cos(0.3 * call); s = sin(0.3 * call);
            A = [c, -s, 0; s, c, 0; 0, 0, 1];
            b = TestAffineStructure.B;
            inverse = A';
        end
    end
end
