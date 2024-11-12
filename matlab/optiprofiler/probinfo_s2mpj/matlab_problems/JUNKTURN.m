function varargout = JUNKTURN(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : JUNKTURN
%    *********
% 
%    The spacecraft orientation problem by Junkins and Turner. This is a
%    nonlinear optimal control problem.
% 
%    The problem is not convex.
% 
%    Source:
%    A.I Tyatushkin, A.I. Zholudev and N. M. Erinchek,
%    "The gradient method for solving optimal control problems with phase
%    constraints", 
%    in "System Modelling and Optimization", P. Kall, ed., pp. 456--464,
%    Springer Verlag, Lecture Notes in Control and Information Sciences 180, 1992.
%    This reference itself refers to:
%    I.L. Junkins and I.D. Turner,
%    "Optimal continuous torque attitude maneuvers",
%    AIAA/AAS Astrodynamics Conference, Palo Alto, 1978.
% 
%    SIF input: Ph. Toint, February 1994.
% 
%    classification = 'C-CQQR2-MN-V-V'
% 
%    Number of discretized points in [0,100] - 1.
%    The number of variables is    10 * ( N + 1 )
%    The number of constraints is  7 * N
%    N should be large enough to ensure feasibility.
% 
%       Alternative values for the SIF file parameters:
% IE N                   50             $-PARAMETER n =     510, m =    350
% IE N                   100            $-PARAMETER n =    1010, m =    700
% IE N                   500            $-PARAMETER n =    5010, m =   3500
% IE N                   1000           $-PARAMETER n =   10010, m =   7000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'JUNKTURN';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
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
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   10000          $-PARAMETER n =  100010, m =  70000
% IE N                   20000          $-PARAMETER n =  200010, m = 140000
% IE N                   100000         $-PARAMETER n = 1000010, m = 700000
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 100.0/v_('RN');
        v_('H/2') = 0.5*v_('H');
        v_('6H/5') = 1.2*v_('H');
        v_('SH') = 1.0909*v_('H');
        v_('S1H') = -0.08333*v_('H');
        v_('S2H') = 0.18182*v_('H');
        v_('-H/2') = -0.5*v_('H');
        v_('-H') = -1.0*v_('H');
        v_('-H/10') = -0.1*v_('H');
        v_('H/4') = 0.25*v_('H');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('7')
            for T=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(T)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(T)];
            end
        end
        for I=v_('1'):v_('3')
            for T=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(T)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(T)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for T=v_('1'):v_('N')
            v_('T-1') = -1+T;
            for I=v_('1'):v_('7')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(T)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['C',int2str(I),',',int2str(T)];
                iv = ix_(['X',int2str(I),',',int2str(T)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('T-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('5'))),',',int2str(T)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(round(v_('5'))),',',int2str(T)];
            iv = ix_(['U',int2str(round(v_('1'))),',',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('H');
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('6'))),',',int2str(T)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(round(v_('6'))),',',int2str(T)];
            iv = ix_(['U',int2str(round(v_('2'))),',',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('6H/5')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('6H/5');
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('7'))),',',int2str(T)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(round(v_('7'))),',',int2str(T)];
            iv = ix_(['U',int2str(round(v_('3'))),',',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SH');
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 1.0;
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 1.0;
        pb.xlower(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('3'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('3'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('4'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('4'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('5'))),',',int2str(round(v_('0')))]),1) = 0.01;
        pb.xupper(ix_(['X',int2str(round(v_('5'))),',',int2str(round(v_('0')))]),1) = 0.01;
        pb.xlower(ix_(['X',int2str(round(v_('6'))),',',int2str(round(v_('0')))]),1) = 0.005;
        pb.xupper(ix_(['X',int2str(round(v_('6'))),',',int2str(round(v_('0')))]),1) = 0.005;
        pb.xlower(ix_(['X',int2str(round(v_('7'))),',',int2str(round(v_('0')))]),1) = 0.001;
        pb.xupper(ix_(['X',int2str(round(v_('7'))),',',int2str(round(v_('0')))]),1) = 0.001;
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('N')))]),1) = 0.43047;
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('N')))]),1) = 0.43047;
        pb.xlower(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('N')))]),1) = 0.70106;
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('N')))]),1) = 0.70106;
        pb.xlower(ix_(['X',int2str(round(v_('3'))),',',int2str(round(v_('N')))]),1) = 0.0923;
        pb.xupper(ix_(['X',int2str(round(v_('3'))),',',int2str(round(v_('N')))]),1) = 0.0923;
        pb.xlower(ix_(['X',int2str(round(v_('4'))),',',int2str(round(v_('N')))]),1) = 0.56098;
        pb.xupper(ix_(['X',int2str(round(v_('4'))),',',int2str(round(v_('N')))]),1) = 0.56098;
        pb.xlower(ix_(['X',int2str(round(v_('5'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('5'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('6'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('6'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('7'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('7'))),',',int2str(round(v_('N')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) =...
              1.0;
        pb.x0(ix_(['X',int2str(round(v_('5'))),',',int2str(round(v_('0')))]),1) =...
              0.01;
        pb.x0(ix_(['X',int2str(round(v_('6'))),',',int2str(round(v_('0')))]),1) =...
              0.005;
        pb.x0(ix_(['X',int2str(round(v_('7'))),',',int2str(round(v_('0')))]),1) =...
              0.001;
        pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('N')))]),1) =...
              0.43047;
        pb.x0(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('N')))]),1) =...
              0.70106;
        pb.x0(ix_(['X',int2str(round(v_('3'))),',',int2str(round(v_('N')))]),1) =...
              0.0923;
        pb.x0(ix_(['X',int2str(round(v_('4'))),',',int2str(round(v_('N')))]),1) =...
              0.56098;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for T=v_('0'):v_('N')
            ename = ['U1S',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(round(v_('1'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['U2S',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(round(v_('2'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['U3S',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(round(v_('3'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for T=v_('1'):v_('N')
            ename = ['P15',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('1'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P16',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('1'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P17',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('1'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P25',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('2'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P26',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('2'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P27',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('2'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P35',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('3'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P36',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('3'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P37',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('3'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P45',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('4'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P46',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('4'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P47',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('4'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P56',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P57',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('5'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P67',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('6'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('7'))),',',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['U1S',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('H/4');
        for T=v_('1'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U1S',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U2S',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U3S',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
        end
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['U1S',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('H/4');
        for T=v_('1'):v_('N')
            ig = ig_(['C',int2str(round(v_('1'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P25',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P36',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P47',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            ig = ig_(['C',int2str(round(v_('2'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P15',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P37',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P46',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            ig = ig_(['C',int2str(round(v_('3'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P16',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P27',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P45',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            ig = ig_(['C',int2str(round(v_('4'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P17',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P26',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P35',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            ig = ig_(['C',int2str(round(v_('5'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P67',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('S1H');
            ig = ig_(['C',int2str(round(v_('6'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P57',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/10');
            ig = ig_(['C',int2str(round(v_('7'))),',',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P56',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('S2H');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(500)          7.417771100D-5
% LO SOLTN(1000)         1.224842784D-5
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQQR2-MN-V-V';
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

% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

