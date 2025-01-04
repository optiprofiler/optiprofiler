function varargout = METHANB8(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    The methanol-8 problem by Fletcher, with ti values from Table 2.7
% 
%    Source: Problem 2d in
%    J.J. More',"A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer Seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: N. Gould and Ph. Toint, Feb 1991.
% 
%    classification = 'C-CNOR2-MN-31-31'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'METHANB8';

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
        v_('N') = 8;
        v_('M') = 2;
        v_('K') = 2;
        v_('0') = 0;
        v_('1') = 1;
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('K-') = -1+v_('K');
        v_('K+') = 1+v_('K');
        v_('A1') = 18.5751;
        v_('B1') = -3632.649;
        v_('C1') = 239.2;
        v_('A2') = 18.3443;
        v_('B2') = -3841.2203;
        v_('C2') = 228.0;
        v_('AL1') = 0.0;
        v_('ALp1') = 15.97;
        v_('ALpp1') = 0.0422;
        v_('AL2') = 0.0;
        v_('ALp2') = 18.1;
        v_('ALpp2') = 0.0;
        v_('BE1') = 9566.67;
        v_('BEp1') = -1.59;
        v_('BEpp1') = 0.0422;
        v_('BE2') = 10834.67;
        v_('BEp2') = 8.74;
        v_('BEpp2') = 0.0;
        v_('FL1') = 451.25;
        v_('FL2') = 684.25;
        v_('FV1') = 0.0;
        v_('FV2') = 0.0;
        v_('TF') = 89.0;
        v_('B') = 693.37;
        v_('D') = 442.13;
        v_('Q') = 8386200.0;
        v_('PI0') = 1210.0;
        v_('PI1') = 1200.0;
        v_('PI2') = 1190.0;
        v_('PI3') = 1180.0;
        v_('PI4') = 1170.0;
        v_('PI5') = 1160.0;
        v_('PI6') = 1150.0;
        v_('PI7') = 1140.0;
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
            pb.x0(ix_('X0,1'),1) = 0.09203;
        else
            pb.y0(find(pbm.congrps==ig_('X0,1')),1) = 0.09203;
        end
        if(isKey(ix_,'X0,2'))
            pb.x0(ix_('X0,2'),1) = 0.908;
        else
            pb.y0(find(pbm.congrps==ig_('X0,2')),1) = 0.908;
        end
        if(isKey(ix_,'X1,1'))
            pb.x0(ix_('X1,1'),1) = 0.1819;
        else
            pb.y0(find(pbm.congrps==ig_('X1,1')),1) = 0.1819;
        end
        if(isKey(ix_,'X1,2'))
            pb.x0(ix_('X1,2'),1) = 0.8181;
        else
            pb.y0(find(pbm.congrps==ig_('X1,2')),1) = 0.8181;
        end
        if(isKey(ix_,'X2,1'))
            pb.x0(ix_('X2,1'),1) = 0.284;
        else
            pb.y0(find(pbm.congrps==ig_('X2,1')),1) = 0.284;
        end
        if(isKey(ix_,'X2,2'))
            pb.x0(ix_('X2,2'),1) = 0.716;
        else
            pb.y0(find(pbm.congrps==ig_('X2,2')),1) = 0.716;
        end
        if(isKey(ix_,'X3,1'))
            pb.x0(ix_('X3,1'),1) = 0.3051;
        else
            pb.y0(find(pbm.congrps==ig_('X3,1')),1) = 0.3051;
        end
        if(isKey(ix_,'X3,2'))
            pb.x0(ix_('X3,2'),1) = 0.6949;
        else
            pb.y0(find(pbm.congrps==ig_('X3,2')),1) = 0.6949;
        end
        if(isKey(ix_,'X4,1'))
            pb.x0(ix_('X4,1'),1) = 0.3566;
        else
            pb.y0(find(pbm.congrps==ig_('X4,1')),1) = 0.3566;
        end
        if(isKey(ix_,'X4,2'))
            pb.x0(ix_('X4,2'),1) = 0.6434;
        else
            pb.y0(find(pbm.congrps==ig_('X4,2')),1) = 0.6434;
        end
        if(isKey(ix_,'X5,1'))
            pb.x0(ix_('X5,1'),1) = 0.468;
        else
            pb.y0(find(pbm.congrps==ig_('X5,1')),1) = 0.468;
        end
        if(isKey(ix_,'X5,2'))
            pb.x0(ix_('X5,2'),1) = 0.532;
        else
            pb.y0(find(pbm.congrps==ig_('X5,2')),1) = 0.532;
        end
        if(isKey(ix_,'X6,1'))
            pb.x0(ix_('X6,1'),1) = 0.6579;
        else
            pb.y0(find(pbm.congrps==ig_('X6,1')),1) = 0.6579;
        end
        if(isKey(ix_,'X6,2'))
            pb.x0(ix_('X6,2'),1) = 0.3421;
        else
            pb.y0(find(pbm.congrps==ig_('X6,2')),1) = 0.3421;
        end
        if(isKey(ix_,'X7,1'))
            pb.x0(ix_('X7,1'),1) = 0.8763;
        else
            pb.y0(find(pbm.congrps==ig_('X7,1')),1) = 0.8763;
        end
        if(isKey(ix_,'X7,2'))
            pb.x0(ix_('X7,2'),1) = 0.1237;
        else
            pb.y0(find(pbm.congrps==ig_('X7,2')),1) = 0.1237;
        end
        if(isKey(ix_,'T0'))
            pb.x0(ix_('T0'),1) = 107.47;
        else
            pb.y0(find(pbm.congrps==ig_('T0')),1) = 107.47;
        end
        if(isKey(ix_,'T1'))
            pb.x0(ix_('T1'),1) = 102.4;
        else
            pb.y0(find(pbm.congrps==ig_('T1')),1) = 102.4;
        end
        if(isKey(ix_,'T2'))
            pb.x0(ix_('T2'),1) = 97.44;
        else
            pb.y0(find(pbm.congrps==ig_('T2')),1) = 97.44;
        end
        if(isKey(ix_,'T3'))
            pb.x0(ix_('T3'),1) = 96.3;
        else
            pb.y0(find(pbm.congrps==ig_('T3')),1) = 96.3;
        end
        if(isKey(ix_,'T4'))
            pb.x0(ix_('T4'),1) = 93.99;
        else
            pb.y0(find(pbm.congrps==ig_('T4')),1) = 93.99;
        end
        if(isKey(ix_,'T5'))
            pb.x0(ix_('T5'),1) = 89.72;
        else
            pb.y0(find(pbm.congrps==ig_('T5')),1) = 89.72;
        end
        if(isKey(ix_,'T6'))
            pb.x0(ix_('T6'),1) = 83.71;
        else
            pb.y0(find(pbm.congrps==ig_('T6')),1) = 83.71;
        end
        if(isKey(ix_,'T7'))
            pb.x0(ix_('T7'),1) = 78.31;
        else
            pb.y0(find(pbm.congrps==ig_('T7')),1) = 78.31;
        end
        if(isKey(ix_,'V0'))
            pb.x0(ix_('V0'),1) = 886.37;
        else
            pb.y0(find(pbm.congrps==ig_('V0')),1) = 886.37;
        end
        if(isKey(ix_,'V1'))
            pb.x0(ix_('V1'),1) = 910.01;
        else
            pb.y0(find(pbm.congrps==ig_('V1')),1) = 910.01;
        end
        if(isKey(ix_,'V2'))
            pb.x0(ix_('V2'),1) = 922.52;
        else
            pb.y0(find(pbm.congrps==ig_('V2')),1) = 922.52;
        end
        if(isKey(ix_,'V3'))
            pb.x0(ix_('V3'),1) = 926.46;
        else
            pb.y0(find(pbm.congrps==ig_('V3')),1) = 926.46;
        end
        if(isKey(ix_,'V4'))
            pb.x0(ix_('V4'),1) = 935.56;
        else
            pb.y0(find(pbm.congrps==ig_('V4')),1) = 935.56;
        end
        if(isKey(ix_,'V5'))
            pb.x0(ix_('V5'),1) = 952.83;
        else
            pb.y0(find(pbm.congrps==ig_('V5')),1) = 952.83;
        end
        if(isKey(ix_,'V6'))
            pb.x0(ix_('V6'),1) = 975.73;
        else
            pb.y0(find(pbm.congrps==ig_('V6')),1) = 975.73;
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
        pb.pbclass = 'C-CNOR2-MN-31-31';
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

