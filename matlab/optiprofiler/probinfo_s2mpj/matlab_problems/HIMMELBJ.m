function varargout = HIMMELBJ(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELBJ
%    *********
% 
%    An chemical equilibrium problem by A.P. Jones.
%    It has a nonlinear objective and linear constraints
% 
%    Source: problem 6 in
%    D.H. Himmelblau,
%    "Applied nonlinear programming",
%    McGraw-Hill, New-York, 1972.
% 
%    SIF input: Ph. Toint, March 1991.
% 
%    classification = 'C-COLR2-MY-45-14'
% 
%    Number of variable sets
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELBJ';

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
        v_('NSETS') = 7;
        v_('NS1') = 4.0;
        v_('NS2') = 13.0;
        v_('NS3') = 18.0;
        v_('NS4') = 3.0;
        v_('NS5') = 3.0;
        v_('NS6') = 2.0;
        v_('NS7') = 2.0;
        v_('NEQ') = 16;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('13') = 13;
        v_('18') = 18;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for K=v_('1'):v_('NSETS')
            v_('RNSK') = v_(['NS',int2str(K)]);
            v_('NSK') = fix(v_('RNSK'));
            for J=v_('1'):v_('NSK')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(J),',',int2str(K)],ix_);
                pb.xnames{iv} = ['X',int2str(J),',',int2str(K)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,1');
        valA(end+1) = -7.69;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,1');
        valA(end+1) = -11.52;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,1');
        valA(end+1) = -36.60;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,2');
        valA(end+1) = -10.94;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8,2');
        valA(end+1) = 2.5966;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,2');
        valA(end+1) = -39.39;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,2');
        valA(end+1) = -21.35;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,2');
        valA(end+1) = -32.84;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,2');
        valA(end+1) = 6.26;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,3');
        valA(end+1) = 10.45;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,3');
        valA(end+1) = -0.5;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7,3');
        valA(end+1) = 2.2435;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,3');
        valA(end+1) = -39.39;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,3');
        valA(end+1) = -21.49;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,3');
        valA(end+1) = -32.84;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,3');
        valA(end+1) = 6.12;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15,3');
        valA(end+1) = -1.9028;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16,3');
        valA(end+1) = -2.8889;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17,3');
        valA(end+1) = -3.3622;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18,3');
        valA(end+1) = -7.4854;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,4');
        valA(end+1) = -15.639;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,4');
        valA(end+1) = 21.81;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,5');
        valA(end+1) = -16.79;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,5');
        valA(end+1) = 18.9779;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,6');
        valA(end+1) = 11.959;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,7');
        valA(end+1) = 12.899;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16,3');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17,3');
        valA(end+1) = 3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18,3');
        valA(end+1) = 4.0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,6');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,7');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,4');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,5');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,6');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,7');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17,3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18,3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5,2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6,2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10,2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12,2');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13,2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14,3');
        valA(end+1) = -4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15,3');
        valA(end+1) = -3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16,3');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17,3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C13';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15,3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16,3');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17,3');
        valA(end+1) = -3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18,3');
        valA(end+1) = -4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,4');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C14';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,5');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C15';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,4');
        valA(end+1) = -4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,6');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,6');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C16';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3,5');
        valA(end+1) = -4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1,7');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2,7');
        valA(end+1) = 1.0;
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
        pbm.gconst(ig_('C1')) = 0.652981;
        pbm.gconst(ig_('C2')) = 0.281941;
        pbm.gconst(ig_('C3')) = 3.705233;
        pbm.gconst(ig_('C4')) = 47.00022;
        pbm.gconst(ig_('C5')) = 47.02972;
        pbm.gconst(ig_('C6')) = 0.08005;
        pbm.gconst(ig_('C7')) = 0.08813;
        pbm.gconst(ig_('C8')) = 0.04829;
        pbm.gconst(ig_('C11')) = 0.0022725;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 1.e-12*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_('X13,2'),1) = 0.0155;
        pb.xupper(ix_('X13,2'),1) = 0.0155;
        pb.xlower(ix_('X13,3'),1) = 0.0211275;
        pb.xupper(ix_('X13,3'),1) = 0.0211275;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.1*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX2',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX3',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        elftv{it}{4} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX4',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        elftv{it}{4} = 'Y4';
        elftv{it}{5} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX13',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        elftv{it}{4} = 'Y4';
        elftv{it}{5} = 'Y5';
        elftv{it}{6} = 'Y6';
        elftv{it}{7} = 'Y7';
        elftv{it}{8} = 'Y8';
        elftv{it}{9} = 'Y9';
        elftv{it}{10} = 'Y10';
        elftv{it}{11} = 'Y11';
        elftv{it}{12} = 'Y12';
        elftv{it}{13} = 'Y13';
        elftv{it}{14} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX18',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        elftv{it}{4} = 'Y4';
        elftv{it}{5} = 'Y5';
        elftv{it}{6} = 'Y6';
        elftv{it}{7} = 'Y7';
        elftv{it}{8} = 'Y8';
        elftv{it}{9} = 'Y9';
        elftv{it}{10} = 'Y10';
        elftv{it}{11} = 'Y11';
        elftv{it}{12} = 'Y12';
        elftv{it}{13} = 'Y13';
        elftv{it}{14} = 'Y14';
        elftv{it}{15} = 'Y15';
        elftv{it}{16} = 'Y16';
        elftv{it}{17} = 'Y17';
        elftv{it}{18} = 'Y18';
        elftv{it}{19} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for K=v_('1'):v_('NSETS')
            v_('RNSK') = v_(['NS',int2str(K)]);
            v_('NSK') = fix(v_('RNSK'));
            for J=v_('1'):v_('NSK')
                ename = ['A',int2str(J),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXLOGX';
                ielftype(ie) = iet_('eXLOGX');
                vname = ['X',int2str(J),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        v_('K') = 1;
        for J=v_('1'):v_('4')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX4';
            ielftype(ie) = iet_('eXLOGX4');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3,1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4,1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 2;
        for J=v_('1'):v_('13')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX13';
            ielftype(ie) = iet_('eXLOGX13');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X5,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X6,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X7,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X8,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X9,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y9',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X10,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y10',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X11,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y11',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X12,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y12',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X13,2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y13',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 3;
        for J=v_('1'):v_('18')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX18';
            ielftype(ie) = iet_('eXLOGX18');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X5,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X6,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X7,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X8,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X9,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y9',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X10,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y10',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X11,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y11',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X12,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y12',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X13,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y13',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X14,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y14',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X15,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y15',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X16,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y16',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X17,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y17',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X18,3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y18',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 4;
        for J=v_('1'):v_('3')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX3';
            ielftype(ie) = iet_('eXLOGX3');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3,4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 5;
        for J=v_('1'):v_('3')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX3';
            ielftype(ie) = iet_('eXLOGX3');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3,5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 6;
        for J=v_('1'):v_('2')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX2';
            ielftype(ie) = iet_('eXLOGX2');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        v_('K') = 7;
        for J=v_('1'):v_('2')
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX2';
            ielftype(ie) = iet_('eXLOGX2');
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(J),',',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1,7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(J),',',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2,7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],0.1);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for K=v_('1'):v_('NSETS')
            v_('RNSK') = v_(['NS',int2str(K)]);
            v_('NSK') = fix(v_('RNSK'));
            for J=v_('1'):v_('NSK')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(J),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(J),',',int2str(K)]);
                pbm.grelw{ig}(posel) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -1910.344724
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-MY-45-14';
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

    case 'eXLOGX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LOGX = log(EV_(1));
        varargout{1} = EV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX+1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0/EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eXLOGX2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(1,3) = U_(1,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXLOGX3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(1,4) = U_(1,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXLOGX4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,5);
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(1,5) = U_(1,5)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXLOGX13'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,14);
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)+1;
        U_(2,7) = U_(2,7)+1;
        U_(2,8) = U_(2,8)+1;
        U_(2,9) = U_(2,9)+1;
        U_(2,10) = U_(2,10)+1;
        U_(2,11) = U_(2,11)+1;
        U_(2,12) = U_(2,12)+1;
        U_(2,13) = U_(2,13)+1;
        U_(1,14) = U_(1,14)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXLOGX18'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,19);
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)+1;
        U_(2,7) = U_(2,7)+1;
        U_(2,8) = U_(2,8)+1;
        U_(2,9) = U_(2,9)+1;
        U_(2,10) = U_(2,10)+1;
        U_(2,11) = U_(2,11)+1;
        U_(2,12) = U_(2,12)+1;
        U_(2,13) = U_(2,13)+1;
        U_(2,14) = U_(2,14)+1;
        U_(2,15) = U_(2,15)+1;
        U_(2,16) = U_(2,16)+1;
        U_(2,17) = U_(2,17)+1;
        U_(2,18) = U_(2,18)+1;
        U_(1,19) = U_(1,19)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
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

