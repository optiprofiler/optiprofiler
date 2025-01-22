function varargout = HYDCAR20(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HYDCAR20
%    *********
% 
%    The hydrocarbon-20 problem by Fletcher.
% 
%    Source: Problem 2b in
%    J.J. More',"A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer Seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input : N. Gould and Ph. Toint, Feb 1991.
% 
%    classification = 'C-CNOR2-AN-99-99'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HYDCAR20';

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
        v_('N') = 20;
        v_('M') = 3;
        v_('K') = 9;
        v_('A1') = 9.647;
        v_('B1') = -2998.00;
        v_('C1') = 230.66;
        v_('A2') = 9.953;
        v_('B2') = -3448.10;
        v_('C2') = 235.88;
        v_('A3') = 9.466;
        v_('B3') = -3347.25;
        v_('C3') = 215.31;
        v_('AL1') = 0.0;
        v_('ALp1') = 37.6;
        v_('ALpp1') = 0.0;
        v_('AL2') = 0.0;
        v_('ALp2') = 48.2;
        v_('ALpp2') = 0.0;
        v_('AL3') = 0.0;
        v_('ALp3') = 45.4;
        v_('ALpp3') = 0.0;
        v_('BE1') = 8425.0;
        v_('BEp1') = 24.2;
        v_('BEpp1') = 0.0;
        v_('BE2') = 9395.0;
        v_('BEp2') = 35.6;
        v_('BEpp2') = 0.0;
        v_('BE3') = 10466.0;
        v_('BEp3') = 31.9;
        v_('BEpp3') = 0.0;
        v_('FL1') = 30.0;
        v_('FL2') = 30.0;
        v_('FL3') = 40.0;
        v_('FV1') = 0.0;
        v_('FV2') = 0.0;
        v_('FV3') = 0.0;
        v_('TF') = 100.0;
        v_('B') = 40.0;
        v_('D') = 60.0;
        v_('Q') = 2500000.0;
        v_('0') = 0;
        v_('1') = 1;
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('K-') = -1+v_('K');
        v_('K+') = 1+v_('K');
        for I=v_('0'):v_('N-1')
            v_(['PI',int2str(I)]) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
            v_(['INVPI',int2str(I)]) = 1.0/v_(['PI',int2str(I)]);
            for J=v_('1'):v_('M')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('0'):v_('N-2')
            [iv,ix_] = s2mpjlib('ii',['V',int2str(I)],ix_);
            pb.xnames{iv} = ['V',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['2.1-',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['2.1-',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('0'))),',',int2str(J)]);
            valA(end+1) = v_('B');
            pbm.gscale(ig,1) = 1.0e+2;
            [ig,ig_] = s2mpjlib('ii',['2.3-',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['2.3-',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('N-1'))),',',int2str(J)]);
            valA(end+1) = -1.0;
            for I=v_('1'):v_('N-2')
                [ig,ig_] = s2mpjlib('ii',['2.2-',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['2.2-',int2str(I),',',int2str(J)];
                pbm.gscale(ig,1) = 1.0e+2;
            end
        end
        for I=v_('0'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['2.7-',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['2.7-',int2str(I)];
        end
        [ig,ig_] = s2mpjlib('ii','2.8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '2.8';
        pbm.gscale(ig,1) = 1.0e+5;
        for I=v_('1'):v_('N-2')
            [ig,ig_] = s2mpjlib('ii',['2.9-',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['2.9-',int2str(I)];
            pbm.gscale(ig,1) = 1.0e+5;
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
        v_('SMALLHF') = 0.0e+0;
        v_('BIGHF') = 0.0e+0;
        for J=v_('1'):v_('M')
            pbm.gconst(ig_(['2.2-',int2str(round(v_('K'))),',',int2str(J)])) =...
                  v_(['FL',int2str(J)]);
            pbm.gconst(ig_(['2.2-',int2str(round(v_('K+'))),',',int2str(J)])) =...
                  v_(['FV',int2str(J)]);
            v_('TFTF') = v_('TF')*v_('TF');
            v_('TEMP1') = v_('TFTF')*v_(['ALpp',int2str(J)]);
            v_('TEMP2') = v_('TF')*v_(['ALp',int2str(J)]);
            v_('TEMP1') = v_('TEMP1')+v_('TEMP2');
            v_('TEMP1') = v_('TEMP1')+v_(['AL',int2str(J)]);
            v_('TEMP1') = v_('TEMP1')*v_(['FL',int2str(J)]);
            v_('SMALLHF') = v_('SMALLHF')+v_('TEMP1');
            v_('TEMP1') = v_('TFTF')*v_(['BEpp',int2str(J)]);
            v_('TEMP2') = v_('TF')*v_(['BEp',int2str(J)]);
            v_('TEMP1') = v_('TEMP1')+v_('TEMP2');
            v_('TEMP1') = v_('TEMP1')+v_(['BE',int2str(J)]);
            v_('TEMP1') = v_('TEMP1')*v_(['FV',int2str(J)]);
            v_('BIGHF') = v_('BIGHF')+v_('TEMP1');
        end
        for I=v_('0'):v_('N-1')
            pbm.gconst(ig_(['2.7-',int2str(I)])) = 1.0;
        end
        pbm.gconst(ig_('2.8')) = v_('Q');
        pbm.gconst(ig_(['2.9-',int2str(round(v_('K')))])) = v_('SMALLHF');
        pbm.gconst(ig_(['2.9-',int2str(round(v_('K+')))])) = v_('BIGHF');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X0,1'))
            pb.x0(ix_('X0,1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X0,1')),1) = 0.0;
        end
        if(isKey(ix_,'X0,2'))
            pb.x0(ix_('X0,2'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X0,2')),1) = 0.3;
        end
        if(isKey(ix_,'X0,3'))
            pb.x0(ix_('X0,3'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X0,3')),1) = 0.1;
        end
        if(isKey(ix_,'X1,1'))
            pb.x0(ix_('X1,1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1,1')),1) = 0.0;
        end
        if(isKey(ix_,'X1,2'))
            pb.x0(ix_('X1,2'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X1,2')),1) = 0.3;
        end
        if(isKey(ix_,'X1,3'))
            pb.x0(ix_('X1,3'),1) = 0.9;
        else
            pb.y0(find(pbm.congrps==ig_('X1,3')),1) = 0.9;
        end
        if(isKey(ix_,'X2,1'))
            pb.x0(ix_('X2,1'),1) = 0.01;
        else
            pb.y0(find(pbm.congrps==ig_('X2,1')),1) = 0.01;
        end
        if(isKey(ix_,'X2,2'))
            pb.x0(ix_('X2,2'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X2,2')),1) = 0.3;
        end
        if(isKey(ix_,'X2,3'))
            pb.x0(ix_('X2,3'),1) = 0.9;
        else
            pb.y0(find(pbm.congrps==ig_('X2,3')),1) = 0.9;
        end
        if(isKey(ix_,'X3,1'))
            pb.x0(ix_('X3,1'),1) = 0.02;
        else
            pb.y0(find(pbm.congrps==ig_('X3,1')),1) = 0.02;
        end
        if(isKey(ix_,'X3,2'))
            pb.x0(ix_('X3,2'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X3,2')),1) = 0.4;
        end
        if(isKey(ix_,'X3,3'))
            pb.x0(ix_('X3,3'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X3,3')),1) = 0.8;
        end
        if(isKey(ix_,'X4,1'))
            pb.x0(ix_('X4,1'),1) = 0.05;
        else
            pb.y0(find(pbm.congrps==ig_('X4,1')),1) = 0.05;
        end
        if(isKey(ix_,'X4,2'))
            pb.x0(ix_('X4,2'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X4,2')),1) = 0.4;
        end
        if(isKey(ix_,'X4,3'))
            pb.x0(ix_('X4,3'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X4,3')),1) = 0.8;
        end
        if(isKey(ix_,'X5,1'))
            pb.x0(ix_('X5,1'),1) = 0.07;
        else
            pb.y0(find(pbm.congrps==ig_('X5,1')),1) = 0.07;
        end
        if(isKey(ix_,'X5,2'))
            pb.x0(ix_('X5,2'),1) = 0.45;
        else
            pb.y0(find(pbm.congrps==ig_('X5,2')),1) = 0.45;
        end
        if(isKey(ix_,'X5,3'))
            pb.x0(ix_('X5,3'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X5,3')),1) = 0.8;
        end
        if(isKey(ix_,'X6,1'))
            pb.x0(ix_('X6,1'),1) = 0.09;
        else
            pb.y0(find(pbm.congrps==ig_('X6,1')),1) = 0.09;
        end
        if(isKey(ix_,'X6,2'))
            pb.x0(ix_('X6,2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X6,2')),1) = 0.5;
        end
        if(isKey(ix_,'X6,3'))
            pb.x0(ix_('X6,3'),1) = 0.7;
        else
            pb.y0(find(pbm.congrps==ig_('X6,3')),1) = 0.7;
        end
        if(isKey(ix_,'X7,1'))
            pb.x0(ix_('X7,1'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X7,1')),1) = 0.1;
        end
        if(isKey(ix_,'X7,2'))
            pb.x0(ix_('X7,2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X7,2')),1) = 0.5;
        end
        if(isKey(ix_,'X7,3'))
            pb.x0(ix_('X7,3'),1) = 0.7;
        else
            pb.y0(find(pbm.congrps==ig_('X7,3')),1) = 0.7;
        end
        if(isKey(ix_,'X8,1'))
            pb.x0(ix_('X8,1'),1) = 0.15;
        else
            pb.y0(find(pbm.congrps==ig_('X8,1')),1) = 0.15;
        end
        if(isKey(ix_,'X8,2'))
            pb.x0(ix_('X8,2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X8,2')),1) = 0.5;
        end
        if(isKey(ix_,'X8,3'))
            pb.x0(ix_('X8,3'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X8,3')),1) = 0.6;
        end
        if(isKey(ix_,'X9,1'))
            pb.x0(ix_('X9,1'),1) = 0.2;
        else
            pb.y0(find(pbm.congrps==ig_('X9,1')),1) = 0.2;
        end
        if(isKey(ix_,'X9,2'))
            pb.x0(ix_('X9,2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X9,2')),1) = 0.5;
        end
        if(isKey(ix_,'X9,3'))
            pb.x0(ix_('X9,3'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X9,3')),1) = 0.6;
        end
        if(isKey(ix_,'X10,1'))
            pb.x0(ix_('X10,1'),1) = 0.25;
        else
            pb.y0(find(pbm.congrps==ig_('X10,1')),1) = 0.25;
        end
        if(isKey(ix_,'X10,2'))
            pb.x0(ix_('X10,2'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X10,2')),1) = 0.6;
        end
        if(isKey(ix_,'X10,3'))
            pb.x0(ix_('X10,3'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X10,3')),1) = 0.5;
        end
        if(isKey(ix_,'X11,1'))
            pb.x0(ix_('X11,1'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X11,1')),1) = 0.3;
        end
        if(isKey(ix_,'X11,2'))
            pb.x0(ix_('X11,2'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X11,2')),1) = 0.6;
        end
        if(isKey(ix_,'X11,3'))
            pb.x0(ix_('X11,3'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X11,3')),1) = 0.5;
        end
        if(isKey(ix_,'X12,1'))
            pb.x0(ix_('X12,1'),1) = 0.35;
        else
            pb.y0(find(pbm.congrps==ig_('X12,1')),1) = 0.35;
        end
        if(isKey(ix_,'X12,2'))
            pb.x0(ix_('X12,2'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X12,2')),1) = 0.6;
        end
        if(isKey(ix_,'X12,3'))
            pb.x0(ix_('X12,3'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X12,3')),1) = 0.5;
        end
        if(isKey(ix_,'X13,1'))
            pb.x0(ix_('X13,1'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X13,1')),1) = 0.4;
        end
        if(isKey(ix_,'X13,2'))
            pb.x0(ix_('X13,2'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('X13,2')),1) = 0.6;
        end
        if(isKey(ix_,'X13,3'))
            pb.x0(ix_('X13,3'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X13,3')),1) = 0.4;
        end
        if(isKey(ix_,'X14,1'))
            pb.x0(ix_('X14,1'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X14,1')),1) = 0.4;
        end
        if(isKey(ix_,'X14,2'))
            pb.x0(ix_('X14,2'),1) = 0.7;
        else
            pb.y0(find(pbm.congrps==ig_('X14,2')),1) = 0.7;
        end
        if(isKey(ix_,'X14,3'))
            pb.x0(ix_('X14,3'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X14,3')),1) = 0.4;
        end
        if(isKey(ix_,'X15,1'))
            pb.x0(ix_('X15,1'),1) = 0.42;
        else
            pb.y0(find(pbm.congrps==ig_('X15,1')),1) = 0.42;
        end
        if(isKey(ix_,'X15,2'))
            pb.x0(ix_('X15,2'),1) = 0.7;
        else
            pb.y0(find(pbm.congrps==ig_('X15,2')),1) = 0.7;
        end
        if(isKey(ix_,'X15,3'))
            pb.x0(ix_('X15,3'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X15,3')),1) = 0.3;
        end
        if(isKey(ix_,'X16,1'))
            pb.x0(ix_('X16,1'),1) = 0.45;
        else
            pb.y0(find(pbm.congrps==ig_('X16,1')),1) = 0.45;
        end
        if(isKey(ix_,'X16,2'))
            pb.x0(ix_('X16,2'),1) = 0.75;
        else
            pb.y0(find(pbm.congrps==ig_('X16,2')),1) = 0.75;
        end
        if(isKey(ix_,'X16,3'))
            pb.x0(ix_('X16,3'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('X16,3')),1) = 0.3;
        end
        if(isKey(ix_,'X17,1'))
            pb.x0(ix_('X17,1'),1) = 0.45;
        else
            pb.y0(find(pbm.congrps==ig_('X17,1')),1) = 0.45;
        end
        if(isKey(ix_,'X17,2'))
            pb.x0(ix_('X17,2'),1) = 0.75;
        else
            pb.y0(find(pbm.congrps==ig_('X17,2')),1) = 0.75;
        end
        if(isKey(ix_,'X17,3'))
            pb.x0(ix_('X17,3'),1) = 0.2;
        else
            pb.y0(find(pbm.congrps==ig_('X17,3')),1) = 0.2;
        end
        if(isKey(ix_,'X18,1'))
            pb.x0(ix_('X18,1'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X18,1')),1) = 0.5;
        end
        if(isKey(ix_,'X18,2'))
            pb.x0(ix_('X18,2'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X18,2')),1) = 0.8;
        end
        if(isKey(ix_,'X18,3'))
            pb.x0(ix_('X18,3'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X18,3')),1) = 0.1;
        end
        if(isKey(ix_,'X19,1'))
            pb.x0(ix_('X19,1'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('X19,1')),1) = 0.5;
        end
        if(isKey(ix_,'X19,2'))
            pb.x0(ix_('X19,2'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X19,2')),1) = 0.8;
        end
        if(isKey(ix_,'X19,3'))
            pb.x0(ix_('X19,3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X19,3')),1) = 0.0;
        end
        for I=v_('0'):v_('N-1')
            if(isKey(ix_,['T',int2str(I)]))
                pb.x0(ix_(['T',int2str(I)]),1) = 100.0;
            else
                pb.y0(find(pbm.congrps==ig_(['T',int2str(I)])),1) = 100.0;
            end
        end
        for I=v_('0'):v_('N-2')
            if(isKey(ix_,['V',int2str(I)]))
                pb.x0(ix_(['V',int2str(I)]),1) = 300.0;
            else
                pb.y0(find(pbm.congrps==ig_(['V',int2str(I)])),1) = 300.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        [it,iet_] = s2mpjlib( 'ii', 'ePOLY1PRD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P6';
        elftp{it}{3} = 'P7';
        elftp{it}{4} = 'P8';
        [it,iet_] = s2mpjlib( 'ii', 'ePOLY2PRD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P6';
        elftp{it}{4} = 'P7';
        elftp{it}{5} = 'P8';
        [it,iet_] = s2mpjlib( 'ii', 'eEXP2PROD',iet_);
        elftv{it}{1} = 'V2';
        elftv{it}{2} = 'V3';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        elftp{it}{5} = 'P5';
        [it,iet_] = s2mpjlib( 'ii', 'eEXP3PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        elftp{it}{5} = 'P5';
        [it,iet_] = s2mpjlib( 'ii', 'eEXP4PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        elftp{it}{5} = 'P5';
        elftp{it}{6} = 'P6';
        elftp{it}{7} = 'P7';
        elftp{it}{8} = 'P8';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        v_('-D') = -1.0*v_('D');
        for J=v_('1'):v_('M')
            ename = ['E11-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PROD';
            ielftype(ie) = iet_('en2PROD');
            vname = ['X',int2str(round(v_('1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['V',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
            ename = ['E12-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXP3PROD';
            ielftype(ie) = iet_('eEXP3PROD');
            vname = ['V',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('0'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['INVPI',int2str(round(v_('0')))]);
            [~,posep] = ismember('P3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
            [~,posep] = ismember('P4',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
            [~,posep] = ismember('P5',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
            for I=v_('1'):v_('N-2')
                v_('I-1') = -1+I;
                v_('I+1') = 1+I;
                ename = ['E21-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PROD';
                ielftype(ie) = iet_('en2PROD');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
                ename = ['E22-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXP3PROD';
                ielftype(ie) = iet_('eEXP3PROD');
                vname = ['V',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['INVPI',int2str(round(v_('I-1')))]);
                [~,posep] = ismember('P3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
                [~,posep] = ismember('P4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
                [~,posep] = ismember('P5',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
                ename = ['E23-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PROD';
                ielftype(ie) = iet_('en2PROD');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                ename = ['E24-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXP3PROD';
                ielftype(ie) = iet_('eEXP3PROD');
                vname = ['V',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['INVPI',int2str(I)]);
                [~,posep] = ismember('P3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
                [~,posep] = ismember('P4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
                [~,posep] = ismember('P5',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
            end
            for I=v_('1'):v_('K-')
                ename = ['E21-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('B');
                ename = ['E23-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('B');
            end
            ename = ['E21-',int2str(round(v_('K'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-D');
            ename = ['E23-',int2str(round(v_('K'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
            for I=v_('K+'):v_('N-2')
                ename = ['E21-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('-D');
                ename = ['E23-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('-D');
            end
            ename = ['E31-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXP2PROD';
            ielftype(ie) = iet_('eEXP2PROD');
            vname = ['X',int2str(round(v_('N-2'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(round(v_('N-2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['INVPI',int2str(round(v_('N-2')))]);
            [~,posep] = ismember('P3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
            [~,posep] = ismember('P4',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
            [~,posep] = ismember('P5',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
        end
        for J=v_('1'):v_('M')
            for I=v_('0'):v_('N-1')
                ename = ['E71-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXP2PROD';
                ielftype(ie) = iet_('eEXP2PROD');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['INVPI',int2str(I)]);
                [~,posep] = ismember('P3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
                [~,posep] = ismember('P4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
                [~,posep] = ismember('P5',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
            end
        end
        for J=v_('1'):v_('M')
            ename = ['E81-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXP4PROD';
            ielftype(ie) = iet_('eEXP4PROD');
            vname = ['V',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('0'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['INVPI',int2str(round(v_('0')))]);
            [~,posep] = ismember('P3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
            [~,posep] = ismember('P4',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
            [~,posep] = ismember('P5',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
            [~,posep] = ismember('P6',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['BE',int2str(J)]);
            [~,posep] = ismember('P7',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['BEp',int2str(J)]);
            [~,posep] = ismember('P8',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['BEpp',int2str(J)]);
            ename = ['E82-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePOLY1PRD';
            ielftype(ie) = iet_('ePOLY1PRD');
            vname = ['X',int2str(round(v_('0'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
            [~,posep] = ismember('P6',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['AL',int2str(J)]);
            [~,posep] = ismember('P7',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['ALp',int2str(J)]);
            [~,posep] = ismember('P8',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['ALpp',int2str(J)]);
            ename = ['E83-',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePOLY2PRD';
            ielftype(ie) = iet_('ePOLY2PRD');
            vname = ['X',int2str(round(v_('1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['V',int2str(round(v_('0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
            [~,posep] = ismember('P6',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['AL',int2str(J)]);
            [~,posep] = ismember('P7',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['ALp',int2str(J)]);
            [~,posep] = ismember('P8',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['ALpp',int2str(J)]);
            for I=v_('1'):v_('N-2')
                v_('I-1') = -1+I;
                v_('I+1') = 1+I;
                ename = ['E91-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXP4PROD';
                ielftype(ie) = iet_('eEXP4PROD');
                vname = ['V',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['INVPI',int2str(I)]);
                [~,posep] = ismember('P3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
                [~,posep] = ismember('P4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
                [~,posep] = ismember('P5',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
                [~,posep] = ismember('P6',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BE',int2str(J)]);
                [~,posep] = ismember('P7',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BEp',int2str(J)]);
                [~,posep] = ismember('P8',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BEpp',int2str(J)]);
                ename = ['E92-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePOLY2PRD';
                ielftype(ie) = iet_('ePOLY2PRD');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                [~,posep] = ismember('P6',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['AL',int2str(J)]);
                [~,posep] = ismember('P7',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['ALp',int2str(J)]);
                [~,posep] = ismember('P8',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['ALpp',int2str(J)]);
                ename = ['E93-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXP4PROD';
                ielftype(ie) = iet_('eEXP4PROD');
                vname = ['V',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['INVPI',int2str(round(v_('I-1')))]);
                [~,posep] = ismember('P3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(J)]);
                [~,posep] = ismember('P4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['B',int2str(J)]);
                [~,posep] = ismember('P5',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['C',int2str(J)]);
                [~,posep] = ismember('P6',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BE',int2str(J)]);
                [~,posep] = ismember('P7',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BEp',int2str(J)]);
                [~,posep] = ismember('P8',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['BEpp',int2str(J)]);
                ename = ['E94-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePOLY2PRD';
                ielftype(ie) = iet_('ePOLY2PRD');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['T',int2str(round(v_('I+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('P1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
                [~,posep] = ismember('P6',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['AL',int2str(J)]);
                [~,posep] = ismember('P7',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['ALp',int2str(J)]);
                [~,posep] = ismember('P8',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['ALpp',int2str(J)]);
            end
            for I=v_('1'):v_('K-')
                ename = ['E92-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('B');
                ename = ['E94-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('B');
            end
            ename = ['E92-',int2str(round(v_('K'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
            ename = ['E94-',int2str(round(v_('K'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-D');
            for I=v_('K+'):v_('N-2')
                ename = ['E92-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('-D');
                ename = ['E94-',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('P2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('-D');
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('1'):v_('M')
            ig = ig_(['2.1-',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E11-',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['E12-',int2str(J)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['2.3-',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E31-',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_('2.8');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E81-',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['E82-',int2str(J)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E83-',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            for I=v_('1'):v_('N-2')
                ig = ig_(['2.2-',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E21-',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['E22-',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E23-',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['E24-',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['2.9-',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E91-',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['E92-',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E93-',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['E94-',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
            for I=v_('0'):v_('N-1')
                ig = ig_(['2.7-',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E71-',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-AN-99-99';
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

    case 'en2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*(EV_(2)+pbm.elpar{iel_}(2));
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*(EV_(2)+pbm.elpar{iel_}(2));
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'ePOLY1PRD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        POLY = pbm.elpar{iel_}(2)+pbm.elpar{iel_}(3)*EV_(2)+pbm.elpar{iel_}(4)*...
             EV_(2)*EV_(2);
        DPOLY = pbm.elpar{iel_}(3)+2.0*pbm.elpar{iel_}(4)*EV_(2);
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*POLY;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*POLY;
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1)*DPOLY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1)*DPOLY;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*EV_(1)*2.0e+0*pbm.elpar{iel_}(4);
                varargout{3} = H_;
            end
        end

    case 'ePOLY2PRD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        POLY = pbm.elpar{iel_}(3)+pbm.elpar{iel_}(4)*EV_(3)+pbm.elpar{iel_}(5)*...
             EV_(3)*EV_(3);
        DPOLY = pbm.elpar{iel_}(4)+2.0*pbm.elpar{iel_}(5)*EV_(3);
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*(pbm.elpar{iel_}(2)+EV_(2))*POLY;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(2)+EV_(2))*POLY;
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1)*POLY;
            g_(3,1) = pbm.elpar{iel_}(1)*EV_(1)*(pbm.elpar{iel_}(2)+EV_(2))*DPOLY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = pbm.elpar{iel_}(1)*POLY;
                H_(2,1) = H_(1,2);
                H_(1,3) = pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(2)+EV_(2))*DPOLY;
                H_(3,1) = H_(1,3);
                H_(2,3) = pbm.elpar{iel_}(1)*EV_(1)*DPOLY;
                H_(3,2) = H_(2,3);
                H_(3,3) = pbm.elpar{iel_}(1)*EV_(1)*(pbm.elpar{iel_}(2)+EV_(2))*2.0e+0*...
                     pbm.elpar{iel_}(5);
                varargout{3} = H_;
            end
        end

    case 'eEXP2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPROD =...
              pbm.elpar{iel_}(1)*pbm.elpar{iel_}(2)*exp(pbm.elpar{iel_}(3)+(pbm.elpar{iel_}(4)/(EV_(2)+pbm.elpar{iel_}(5))));
        F = EV_(1)*EXPROD;
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = EXPROD;
            g_(2,1) = -EV_(1)*EXPROD*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(2))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -EXPROD*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(2))^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = F*(pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(2))^2)^2+...
                     2.0e+0*F*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(2))^3;
                varargout{3} = H_;
            end
        end

    case 'eEXP3PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPROD =...
              pbm.elpar{iel_}(1)*pbm.elpar{iel_}(2)*exp(pbm.elpar{iel_}(3)+(pbm.elpar{iel_}(4)/(EV_(3)+pbm.elpar{iel_}(5))));
        F = EV_(1)*EV_(2)*EXPROD;
        TERM = -pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^2;
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = EV_(2)*EXPROD;
            g_(2,1) = EV_(1)*EXPROD;
            g_(3,1) = F*TERM;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EXPROD;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EXPROD*TERM;
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1)*EXPROD*TERM;
                H_(3,2) = H_(2,3);
                H_(3,3) =...
                      F*(TERM*TERM+2.0e+0*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^3);
                varargout{3} = H_;
            end
        end

    case 'eEXP4PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPROD =...
              pbm.elpar{iel_}(1)*pbm.elpar{iel_}(2)*exp(pbm.elpar{iel_}(3)+(pbm.elpar{iel_}(4)/(EV_(3)+pbm.elpar{iel_}(5))));
        F = EV_(1)*EV_(2)*EXPROD;
        POLY = pbm.elpar{iel_}(6)+pbm.elpar{iel_}(7)*EV_(3)+pbm.elpar{iel_}(8)*...
             EV_(3)*EV_(3);
        DPOLY = pbm.elpar{iel_}(7)+2.0*pbm.elpar{iel_}(8)*EV_(3);
        TERM = DPOLY-POLY*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^2;
        varargout{1} = F*POLY;
        if(nargout>1)
            g_(1,1) = EV_(2)*EXPROD*POLY;
            g_(2,1) = EV_(1)*EXPROD*POLY;
            g_(3,1) = F*TERM;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EXPROD*POLY;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EXPROD*TERM;
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1)*EXPROD*TERM;
                H_(3,2) = H_(2,3);
                H_(3,3) =...
                      F*(-(pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^2)*TERM+2.0*pbm.elpar{iel_}(8)-DPOLY*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^2+2.0e+0*POLY*pbm.elpar{iel_}(4)/(pbm.elpar{iel_}(5)+EV_(3))^3);
                varargout{3} = H_;
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

