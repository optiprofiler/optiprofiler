function varargout = NUFFIELD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NUFFIELD
%    *********
% 
%    A problem from economics.
%    Maximize a 2-D integral representing consumer surplus subject to 
%    linear and quadratic constraints representing incentive compatibility
% 
%    Let v( . , . ) : R^2 -> R, Omega = [a,a+1] x [a,a+1], and
%    the corners A, B, C, D be as follows:
% 
%            (a+1,a+1)
%        A *-----* B
%          |     |
%          |     |
%        D *-----* C
%        (a,a)  
% 
%    The problem is to maximize
% 
%       (a+1) line integral_{AB U BC} v(w)dw 
%        - a line integral_{CD U DA} v(w)dw
%        - 3 volume integral_{Omega} v(w)dw
% 
%    subject to v being symmetric (i.e., v(x,y) = v(y,x))
%               v(a,a) = 0
%               nabla_w v(w) >= 0
%               < e, nabla_w v(w) > <= 1
%         and   nabla_ww v(w) positive definite
% 
%    this last constraint is guaranteed by ensuring that
% 
%               d^2 v/dx^2 >= 0
%               d^2 v/dy^2 >= 0
%               ( d^2 v/dx^2 )( d^2 v/dy^2 ) >= ( d^2 v/dxdy )^2
% 
%    Symmetry is ensured by only considering v(x,y) for x <= y
% 
%    Here v(x,y) is the consumer surplus. that is if the consumer values good 
%    1 at x pounds and good 2 at y pounds then they will have a utility 
%    equivalent to v(x,y) pounds after being faced with the optimal monopoly 
%    pricing strategy. (Apparently, from this we can infer what the optimal 
%    pricing strategy was... ).
% 
%    More background is available from
% 
%    "Optimal Selling Strategies: When to haggle, when to hold firm",
%      Riley and Zeckhauser. The Quarterly Journal of Economics, 1983, and
% 
%    "Multidimensional Incentive Compatibility and Mechanism Design", 
%      McAfee and McMillan. The Journal of Economic Theory, 1988.
% 
%    Source: John Thanassoulis <john.thanassoulis@nuffield.oxford.ac.uk>
% 
%    Standard finite-differences are used to ap[proximate derivatives, and 
%    1- and 2-D trapezoidal rules to approximate integrals
% 
%    SIF input: Nick Gould, February 2001
% 
%    classification = 'C-CLQR2-AN-V-V'
% 
%    The parameter a
% 
%       Alternative values for the SIF file parameters:
% RE A                   5.0            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NUFFIELD';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        if(nargs<1)
            v_('A') = 5.0;  %  SIF file default value
        else
            v_('A') = varargin{1};
        end
% IE N                   10            $-PARAMETER
% IE N                   20            $-PARAMETER
% IE N                   30            $-PARAMETER
% IE N                   40            $-PARAMETER
        if(nargs<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{2};
        end
% IE N                   100           $-PARAMETER
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('1/H') = v_('RN');
        v_('-1/H') = -1.0*v_('1/H');
        v_('H**2') = v_('H')*v_('H');
        v_('1/H**2') = v_('1/H')*v_('1/H');
        v_('-2/H**2') = -2.0*v_('1/H**2');
        v_('1/H**4') = v_('1/H**2')*v_('1/H**2');
        v_('A+1') = 1.0+v_('A');
        v_('-A-1') = -1.0*v_('A+1');
        v_('C2') = 3.0*v_('H');
        v_('C3') = 0.5*v_('C2');
        v_('C4') = v_('C3')+v_('A');
        v_('C1') = v_('C3')+v_('-A-1');
        v_('C5') = -1.0+v_('C3');
        v_('C5') = 0.5*v_('C5');
        v_('C6') = 0.5*v_('C3');
        v_('C6') = v_('C6')+v_('-A-1');
        v_('C6') = 0.5*v_('C6');
        v_('C1') = v_('C1')*v_('H');
        v_('C2') = v_('C2')*v_('H');
        v_('C3') = v_('C3')*v_('H');
        v_('C4') = v_('C4')*v_('H');
        v_('C5') = v_('C5')*v_('H');
        v_('C6') = v_('C6')*v_('H');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('N')
            for J=v_('0'):I
                [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('N'))),',',int2str(J)]);
            valA(end+1) = v_('C1');
        end
        for I=v_('2'):v_('N-1')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('C2');
            end
        end
        for I=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(I)]);
            valA(end+1) = v_('C3');
        end
        for I=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('0')))]);
            valA(end+1) = v_('C4');
        end
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('0')))]);
        valA(end+1) = v_('C5');
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        valA(end+1) = v_('C6');
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            for J=v_('0'):I
                [ig,ig_] = s2mpjlib('ii',['VX',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['VX',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(J)]);
                valA(end+1) = v_('1/H');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-1/H');
                [ig,ig_] = s2mpjlib('ii',['VV',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['VV',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(J)]);
                valA(end+1) = v_('1/H');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-1/H');
            end
        end
        for J=v_('0'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['VX',int2str(round(v_('N'))),',',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['VX',int2str(round(v_('N'))),',',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('N'))),',',int2str(J)]);
            valA(end+1) = v_('1/H');
            [ig,ig_] = s2mpjlib('ii',['VX',int2str(round(v_('N'))),',',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['VX',int2str(round(v_('N'))),',',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('N-1'))),',',int2str(J)]);
            valA(end+1) = v_('-1/H');
            [ig,ig_] = s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(J)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('N'))),',',int2str(J)]);
            valA(end+1) = v_('1/H');
            [ig,ig_] = s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(J)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('N-1'))),',',int2str(J)]);
            valA(end+1) = v_('-1/H');
        end
        [ig,ig_] =...
              s2mpjlib('ii',['VX',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['VX',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        valA(end+1) = v_('1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VX',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['VX',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))]);
        valA(end+1) = v_('-1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        valA(end+1) = v_('1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))]);
        valA(end+1) = v_('-1/H');
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('0'):v_('I-1')
                v_('J+1') = 1+J;
                [ig,ig_] = s2mpjlib('ii',['VY',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['VY',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('J+1')))]);
                valA(end+1) = v_('1/H');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-1/H');
                [ig,ig_] = s2mpjlib('ii',['VV',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['VV',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('J+1')))]);
                valA(end+1) = v_('1/H');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-1/H');
            end
        end
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['VY',int2str(I),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['VY',int2str(I),',',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(I)]);
            valA(end+1) = v_('1/H');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(I)]);
            valA(end+1) = v_('-1/H');
            [ig,ig_] = s2mpjlib('ii',['VV',int2str(I),',',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['VV',int2str(I),',',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(I)]);
            valA(end+1) = v_('1/H');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(I)]);
            valA(end+1) = v_('-1/H');
        end
        [ig,ig_] =...
              s2mpjlib('ii',['VY',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['VY',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        valA(end+1) = v_('1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VY',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['VY',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))]);
        valA(end+1) = v_('-1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        valA(end+1) = v_('1/H');
        [ig,ig_] =...
              s2mpjlib('ii',['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['VV',int2str(round(v_('N'))),',',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  =...
              ix_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))]);
        valA(end+1) = v_('-1/H');
        for I=v_('1'):v_('N-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('0'):v_('I-1')
                [ig,ig_] = s2mpjlib('ii',['VXX',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['VXX',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(J)]);
                valA(end+1) = v_('1/H**2');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-2/H**2');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I-1'))),',',int2str(J)]);
                valA(end+1) = v_('1/H**2');
            end
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                [ig,ig_] = s2mpjlib('ii',['VYY',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['VYY',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('J+1')))]);
                valA(end+1) = v_('1/H**2');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('-2/H**2');
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('J-1')))]);
                valA(end+1) = v_('1/H**2');
            end
        end
        for I=v_('1'):v_('N-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['VXX',int2str(I),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['VXX',int2str(I),',',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(I)]);
            valA(end+1) = v_('1/H**2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(I)]);
            valA(end+1) = v_('-2/H**2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('I-1')))]);
            valA(end+1) = v_('1/H**2');
            [ig,ig_] = s2mpjlib('ii',['VYY',int2str(I),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['VYY',int2str(I),',',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(I)]);
            valA(end+1) = v_('1/H**2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(I)]);
            valA(end+1) = v_('-2/H**2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('I-1')))]);
            valA(end+1) = v_('1/H**2');
        end
        for I=v_('1'):v_('N-1')
            for J=v_('1'):I
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('0'):v_('N')
            for J=v_('0'):I
                pbm.gconst(ig_(['VV',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['V',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['V',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCONVEX',iet_);
        elftv{it}{1} = 'VIP1J';
        elftv{it}{2} = 'VIJP1';
        elftv{it}{3} = 'VIJ';
        elftv{it}{4} = 'VIM1J';
        elftv{it}{5} = 'VIJM1';
        elftv{it}{6} = 'VIPJP';
        elftv{it}{7} = 'VIPJM';
        elftv{it}{8} = 'VIMJM';
        elftv{it}{9} = 'VIMJP';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            for J=v_('1'):I
                v_('J+1') = 1+J;
                v_('J-1') = -1+J;
                ename = ['C',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eCONVEX';
                ielftype(ie) = iet_('eCONVEX');
                vname = ['V',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIP1J',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIM1J',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIJP1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIJM1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I+1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIPJP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I-1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIMJM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I+1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIPJM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I-1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('VIMJP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N-1')
            for J=v_('1'):I
                ig = ig_(['C',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('1/H**4');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solutions (may be local!)
% LO SOLTN               -2.512312500   $ (n=10)
% LO SOLTN               -2.512359371   $ (n=20)
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLQR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eCONVEX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,9);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-2;
        U_(1,4) = U_(1,4)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)-2;
        U_(2,5) = U_(2,5)+1;
        U_(3,6) = U_(3,6)+2.500000e-01;
        U_(3,8) = U_(3,8)+2.500000e-01;
        U_(3,9) = U_(3,9)-2.500000e-01;
        U_(3,7) = U_(3,7)-2.500000e-01;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)-IV_(3)*IV_(3);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            g_(3,1) = -2.0*IV_(3);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                H_(3,3) = -2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

