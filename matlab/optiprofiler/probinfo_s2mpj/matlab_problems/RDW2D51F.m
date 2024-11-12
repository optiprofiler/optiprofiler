function varargout = RDW2D51F(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : RDW2D51F
%    *********
% 
%    A finite-element approximation to the distributed optimal control problem
% 
%       min 1/2||u-v||_L2^2 + beta ||f||_L2^2
% 
%    subject to - nabla^2 u = f
% 
%    where v is given on and within the boundary of a unit [0,1] box in 
%    2 dimensions, and u = v on its boundary. The discretization uses 
%    quadrilateral elememts. There are simple bounds on the controls f
% 
%    The problem is stated as a quadratic program
% 
%    Source:  example 5.1 in 
%     T. Rees, H. S. Dollar and A. J. Wathen
%     "Optimal solvers for PDE-constrained optimization"
%     SIAM J. Sci. Comp. (to appear) 2009
% 
%    with the control bounds as specified in 
% 
%     M. Stoll and A. J. Wathen
%     "Preconditioning for PDE constrained optimization with 
%      control constraints"
%     OUCL Technical Report 2009
% 
%    SIF input: Nick Gould, May 2009
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'C-CQOR2-AN-V-V'
% 
%    Number of nodes in each direction (a power of 2)
% 
%       Alternative values for the SIF file parameters:
% IE N                   2             $-PARAMETER
% IE N                   4             $-PARAMETER
% IE N                   8             $-PARAMETER
% IE N                   16            $-PARAMETER
% IE N                   32            $-PARAMETER
% IE N                   64            $-PARAMETER
% IE N                   128           $-PARAMETER
% IE N                   256           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'RDW2D51F';

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
            v_('N') = 4;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   512           $-PARAMETER
% IE N                   1024          $-PARAMETER
% IE N                   2048          $-PARAMETER
% IE N                   4096          $-PARAMETER
% IE N                   8192          $-PARAMETER
% IE N                   16384         $-PARAMETER
        if(nargs<2)
            v_('BETA') = 0.01;  %  SIF file default value
        else
            v_('BETA') = varargin{2};
        end
        v_('ZERO') = 0.0;
        v_('ONE') = 1.0;
        v_('TWO') = 2.0;
        v_('SIX') = 6.0;
        v_('THIRTYSIX') = 36.0;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N-') = -1+v_('N');
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N/2+1') = v_('N/2')+v_('1');
        v_('RN') = v_('N');
        v_('H') = v_('ONE')/v_('RN');
        v_('H**2') = v_('H')*v_('H');
        v_('H**2/36') = v_('H**2')/v_('THIRTYSIX');
        v_('-H**2/36') = -1.0*v_('H**2/36');
        v_('2BETA') = 2.0*v_('BETA');
        v_('2BH**2/36') = v_('2BETA')*v_('H**2/36');
        v_('1/6') = v_('ONE')/v_('SIX');
        for I=v_('0'):v_('N/2')
            v_('RI') = I;
            v_('2RI') = 2.0*v_('RI');
            v_('2RIH') = v_('2RI')*v_('H');
            v_('2RIH-1') = v_('2RIH')-v_('ONE');
            v_('2RIH-1S') = v_('2RIH-1')*v_('2RIH-1');
            for J=v_('0'):v_('N/2')
                v_('RJ') = J;
                v_('2RJ') = 2.0*v_('RJ');
                v_('2RJH') = v_('2RJ')*v_('H');
                v_('2RJH-1') = v_('2RJH')-v_('ONE');
                v_('2RJH-1S') = v_('2RJH-1')*v_('2RJH-1');
                v_('V') = v_('2RIH-1S')*v_('2RJH-1S');
                v_(['V',int2str(I),',',int2str(J)]) = v_('V');
            end
        end
        for I=v_('N/2+1'):v_('N')
            for J=v_('N/2+1'):v_('N')
                v_(['V',int2str(I),',',int2str(J)]) = v_('ZERO');
            end
        end
        for I=v_('0'):v_('N/2')
            for J=v_('N/2+1'):v_('N')
                v_(['V',int2str(I),',',int2str(J)]) = v_('ZERO');
            end
        end
        for I=v_('N/2+1'):v_('N')
            for J=v_('0'):v_('N/2')
                v_(['V',int2str(I),',',int2str(J)]) = v_('ZERO');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['F',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['F',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N-1')
            for J=v_('1'):v_('N-1')
                [ig,ig_] = s2mpjlib('ii',['L',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['L',int2str(I),',',int2str(J)];
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
        pb.xlower(ix_(['U',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = v_(['V',int2str(round(v_('0'))),',',int2str(round(v_('0')))]);
        pb.xupper(ix_(['U',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = v_(['V',int2str(round(v_('0'))),',',int2str(round(v_('0')))]);
        pb.xlower(ix_(['U',int2str(round(v_('N'))),',',int2str(round(v_('0')))]),1) = v_(['V',int2str(round(v_('N'))),',',int2str(round(v_('0')))]);
        pb.xupper(ix_(['U',int2str(round(v_('N'))),',',int2str(round(v_('0')))]),1) = v_(['V',int2str(round(v_('N'))),',',int2str(round(v_('0')))]);
        pb.xlower(ix_(['U',int2str(round(v_('0'))),',',int2str(round(v_('N')))]),1) = v_(['V',int2str(round(v_('0'))),',',int2str(round(v_('N')))]);
        pb.xupper(ix_(['U',int2str(round(v_('0'))),',',int2str(round(v_('N')))]),1) = v_(['V',int2str(round(v_('0'))),',',int2str(round(v_('N')))]);
        pb.xlower(ix_(['U',int2str(round(v_('N'))),',',int2str(round(v_('N')))]),1) = v_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        pb.xupper(ix_(['U',int2str(round(v_('N'))),',',int2str(round(v_('N')))]),1) = v_(['V',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        for I=v_('1'):v_('N-1')
            pb.xlower(ix_(['U',int2str(round(v_('0'))),',',int2str(I)]),1) =...
                  v_(['V',int2str(round(v_('0'))),',',int2str(I)]);
            pb.xupper(ix_(['U',int2str(round(v_('0'))),',',int2str(I)]),1) =...
                  v_(['V',int2str(round(v_('0'))),',',int2str(I)]);
            pb.xlower(ix_(['U',int2str(round(v_('N'))),',',int2str(I)]),1) =...
                  v_(['V',int2str(round(v_('N'))),',',int2str(I)]);
            pb.xupper(ix_(['U',int2str(round(v_('N'))),',',int2str(I)]),1) =...
                  v_(['V',int2str(round(v_('N'))),',',int2str(I)]);
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('0')))]),1) =...
                  v_(['V',int2str(I),',',int2str(round(v_('0')))]);
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('0')))]),1) =...
                  v_(['V',int2str(I),',',int2str(round(v_('0')))]);
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('N')))]),1) =...
                  v_(['V',int2str(I),',',int2str(round(v_('N')))]);
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('N')))]),1) =...
                  v_(['V',int2str(I),',',int2str(round(v_('N')))]);
        end
        pb.xlower(ix_(['F',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['F',int2str(round(v_('0'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['F',int2str(round(v_('N'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['F',int2str(round(v_('N'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['F',int2str(round(v_('0'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['F',int2str(round(v_('0'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xlower(ix_(['F',int2str(round(v_('N'))),',',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['F',int2str(round(v_('N'))),',',int2str(round(v_('N')))]),1) = 0.0;
        for I=v_('1'):v_('N-1')
            pb.xlower(ix_(['F',int2str(round(v_('0'))),',',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['F',int2str(round(v_('0'))),',',int2str(I)]),1) = 0.0;
            pb.xlower(ix_(['F',int2str(round(v_('N'))),',',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['F',int2str(round(v_('N'))),',',int2str(I)]),1) = 0.0;
            pb.xlower(ix_(['F',int2str(I),',',int2str(round(v_('0')))]),1) = 0.0;
            pb.xupper(ix_(['F',int2str(I),',',int2str(round(v_('0')))]),1) = 0.0;
            pb.xlower(ix_(['F',int2str(I),',',int2str(round(v_('N')))]),1) = 0.0;
            pb.xupper(ix_(['F',int2str(I),',',int2str(round(v_('N')))]),1) = 0.0;
        end
        for I=v_('1'):v_('N-1')
            v_('RI') = I;
            v_('X1') = v_('RI')*v_('H');
            v_('X1**2') = v_('X1')*v_('X1');
            v_('-X1**2') = -1.0*v_('X1**2');
            v_('2-X1') = v_('TWO')-v_('X1');
            v_('0.1(2-X1)') = 0.1*v_('2-X1');
            for J=v_('1'):v_('N-1')
                v_('RJ') = J;
                v_('X2') = v_('RJ')*v_('H');
                v_('X2**2') = v_('X2')*v_('X2');
                v_('ARG') = v_('-X1**2')-v_('X2**2');
                v_('EARG') = exp(v_('ARG'));
                v_('UA') = v_('0.1(2-X1)')*v_('EARG');
                pb.xlower(ix_(['F',int2str(I),',',int2str(J)]),1) = v_('UA');
            end
            for J=v_('1'):v_('N/2')
                pb.xupper(ix_(['F',int2str(I),',',int2str(J)])) = 0.6;
            end
            for J=v_('N/2+1'):v_('N-1')
                pb.xupper(ix_(['F',int2str(I),',',int2str(J)])) = 0.9;
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                v_('I+J') = I+J;
                v_('RI+J') = v_('I+J');
                v_('RI+J/N') = v_('RI+J')/v_('RN');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eM',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftv{it}{4} = 'U4';
        elftp{it}{1} = 'V1';
        elftp{it}{2} = 'V2';
        elftp{it}{3} = 'V3';
        elftp{it}{4} = 'V4';
        [it,iet_] = s2mpjlib( 'ii', 'eM0',iet_);
        elftv{it}{1} = 'F1';
        elftv{it}{2} = 'F2';
        elftv{it}{3} = 'F3';
        elftv{it}{4} = 'F4';
        [it,iet_] = s2mpjlib( 'ii', 'eA',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftv{it}{4} = 'U4';
        [it,iet_] = s2mpjlib( 'ii', 'eB',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftv{it}{4} = 'U4';
        [it,iet_] = s2mpjlib( 'ii', 'eC',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftv{it}{4} = 'U4';
        [it,iet_] = s2mpjlib( 'ii', 'eD',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftv{it}{4} = 'U4';
        [it,iet_] = s2mpjlib( 'ii', 'eP',iet_);
        elftv{it}{1} = 'F1';
        elftv{it}{2} = 'F2';
        elftv{it}{3} = 'F3';
        elftv{it}{4} = 'F4';
        [it,iet_] = s2mpjlib( 'ii', 'eQ',iet_);
        elftv{it}{1} = 'F1';
        elftv{it}{2} = 'F2';
        elftv{it}{3} = 'F3';
        elftv{it}{4} = 'F4';
        [it,iet_] = s2mpjlib( 'ii', 'eR',iet_);
        elftv{it}{1} = 'F1';
        elftv{it}{2} = 'F2';
        elftv{it}{3} = 'F3';
        elftv{it}{4} = 'F4';
        [it,iet_] = s2mpjlib( 'ii', 'eS',iet_);
        elftv{it}{1} = 'F1';
        elftv{it}{2} = 'F2';
        elftv{it}{3} = 'F3';
        elftv{it}{4} = 'F4';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('0'):v_('N-1')
            v_('I+') = I+v_('1');
            for J=v_('0'):v_('N-1')
                v_('J+') = J+v_('1');
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eM';
                ielftype(ie) = iet_('eM');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('V1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['V',int2str(I),',',int2str(J)]);
                [~,posep] = ismember('V2',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['V',int2str(I),',',int2str(round(v_('J+')))]);
                [~,posep] = ismember('V3',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['V',int2str(round(v_('I+'))),',',int2str(J)]);
                [~,posep] = ismember('V4',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) =...
                      v_(['V',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))]);
            end
        end
        for I=v_('0'):v_('N-1')
            v_('I+') = I+v_('1');
            for J=v_('0'):v_('N-1')
                v_('J+') = J+v_('1');
                ename = ['F',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eM0';
                ielftype(ie) = iet_('eM0');
                vname = ['F',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('0'):v_('N-1')
            v_('I+') = I+v_('1');
            for J=v_('0'):v_('N-1')
                v_('J+') = J+v_('1');
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA';
                ielftype(ie) = iet_('eA');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eB';
                ielftype(ie) = iet_('eB');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['C',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eC';
                ielftype(ie) = iet_('eC');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['D',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eD';
                ielftype(ie) = iet_('eD');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('0'):v_('N-1')
            v_('I+') = I+v_('1');
            for J=v_('0'):v_('N-1')
                v_('J+') = J+v_('1');
                ename = ['P',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eP';
                ielftype(ie) = iet_('eP');
                vname = ['F',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Q',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eQ';
                ielftype(ie) = iet_('eQ');
                vname = ['F',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['R',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eR';
                ielftype(ie) = iet_('eR');
                vname = ['F',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['S',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eS';
                ielftype(ie) = iet_('eS');
                vname = ['F',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['F',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('F4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('0'):v_('N-1')
            for J=v_('0'):v_('N-1')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('H**2/36');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['F',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('2BH**2/36');
            end
        end
        for I=v_('1'):v_('N-2')
            v_('I+') = I+v_('1');
            for J=v_('1'):v_('N-2')
                v_('J+') = J+v_('1');
                ig = ig_(['L',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('1/6');
                ig = ig_(['L',int2str(I),',',int2str(round(v_('J+')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('1/6');
                ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('1/6');
                ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['D',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('1/6');
                ig = ig_(['L',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-H**2/36');
                ig = ig_(['L',int2str(I),',',int2str(round(v_('J+')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-H**2/36');
                ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['R',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-H**2/36');
                ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('J+')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-H**2/36');
            end
        end
        for I=v_('1'):v_('N-2')
            v_('I+') = I+v_('1');
            ig = ig_(['L',int2str(I),',',int2str(round(v_('N-')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(round(v_('N-')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(I),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(round(v_('0')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('N-')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(round(v_('N-')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(I),',',int2str(round(v_('0')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(I),',',int2str(round(v_('N-')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(round(v_('N-')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(I),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(round(v_('0')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('N-')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['R',int2str(I),',',int2str(round(v_('N-')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(round(v_('I+'))),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(round(v_('0')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
        end
        for J=v_('1'):v_('N-2')
            v_('J+') = J+v_('1');
            ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('N-'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('J+')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(round(v_('N-'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('1'))),',',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('0'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('J+')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('0'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/6');
            ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('J+')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(round(v_('1'))),',',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['R',int2str(round(v_('0'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
            ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('J+')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('0'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H**2/36');
        end
        ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('N-')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['A',int2str(round(v_('N-'))),',',int2str(round(v_('N-')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/6');
        ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['B',int2str(round(v_('N-'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/6');
        ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('N-')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['C',int2str(round(v_('0'))),',',int2str(round(v_('N-')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/6');
        ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['D',int2str(round(v_('0'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/6');
        ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('N-')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['P',int2str(round(v_('N-'))),',',int2str(round(v_('N-')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-H**2/36');
        ig = ig_(['L',int2str(round(v_('N-'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['Q',int2str(round(v_('N-'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-H**2/36');
        ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('N-')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['R',int2str(round(v_('0'))),',',int2str(round(v_('N-')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-H**2/36');
        ig = ig_(['L',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('0'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-H**2/36');
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-AN-V-V';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eM'

        EV_  = varargin{1};
        iel_ = varargin{2};
        UV1 = EV_(1)-pbm.elpar{iel_}(1);
        UV2 = EV_(2)-pbm.elpar{iel_}(2);
        UV3 = EV_(3)-pbm.elpar{iel_}(3);
        UV4 = EV_(4)-pbm.elpar{iel_}(4);
        varargout{1} = 2.0*UV1^2+2.0*UV2^2+2.0*UV3^2+2.0*UV4^2+2.0*UV1*UV2+...
             2.0*UV1*UV3+UV1*UV4+UV2*UV3+2.0*UV2*UV4+2.0*UV3*UV4;
        if(nargout>1)
            g_(1,1) = 4.0*UV1+2.0*UV2+2.0*UV3+UV4;
            g_(2,1) = 2.0*UV1+4.0*UV2+UV3+2.0*UV4;
            g_(3,1) = 2.0*UV1+UV2+4.0*UV3+2.0*UV4;
            g_(4,1) = UV1+2.0*UV2+2.0*UV3+4.0*UV4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 4.0;
                H_(1,2) = 2.0;
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0;
                H_(3,1) = H_(1,3);
                H_(1,4) = 1.0;
                H_(4,1) = H_(1,4);
                H_(2,2) = 4.0;
                H_(2,3) = 1.0;
                H_(3,2) = H_(2,3);
                H_(2,4) = 2.0;
                H_(4,2) = H_(2,4);
                H_(3,3) = 4.0;
                H_(3,4) = 2.0;
                H_(4,3) = H_(3,4);
                H_(4,4) = 4.0;
                varargout{3} = H_;
            end
        end

    case 'eM0'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 2.0*EV_(1)^2+2.0*EV_(2)^2+2.0*EV_(3)^2+2.0*EV_(4)^2+...
             2.0*EV_(1)*EV_(2)+2.0*EV_(1)*EV_(3)+EV_(1)*EV_(4)+EV_(2)*EV_(3)+2.0*EV_(2)*EV_(4)+2.0*EV_(3)*EV_(4);
        if(nargout>1)
            g_(1,1) = 4.0*EV_(1)+2.0*EV_(2)+2.0*EV_(3)+EV_(4);
            g_(2,1) = 2.0*EV_(1)+4.0*EV_(2)+EV_(3)+2.0*EV_(4);
            g_(3,1) = 2.0*EV_(1)+EV_(2)+4.0*EV_(3)+2.0*EV_(4);
            g_(4,1) = EV_(1)+2.0*EV_(2)+2.0*EV_(3)+4.0*EV_(4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 4.0;
                H_(1,2) = 2.0;
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0;
                H_(3,1) = H_(1,3);
                H_(1,4) = 1.0;
                H_(4,1) = H_(1,4);
                H_(2,2) = 4.0;
                H_(2,3) = 1.0;
                H_(3,2) = H_(2,3);
                H_(2,4) = 2.0;
                H_(4,2) = H_(2,4);
                H_(3,3) = 4.0;
                H_(3,4) = 2.0;
                H_(4,3) = H_(3,4);
                H_(4,4) = 4.0;
                varargout{3} = H_;
            end
        end

    case 'eA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = 4.0;
        C2 = -1.0;
        C3 = -1.0;
        C4 = -2.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = -1.0;
        C2 = 4.0;
        C3 = -2.0;
        C4 = -1.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = -1.0;
        C2 = -2.0;
        C3 = 4.0;
        C4 = -1.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = -2.0;
        C2 = -1.0;
        C3 = -1.0;
        C4 = 4.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = 4.0;
        C2 = 2.0;
        C3 = 2.0;
        C4 = 1.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = 2.0;
        C2 = 4.0;
        C3 = 1.0;
        C4 = 2.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = 2.0;
        C2 = 1.0;
        C3 = 4.0;
        C4 = 2.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
                varargout{3} = H_;
            end
        end

    case 'eS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C1 = 1.0;
        C2 = 2.0;
        C3 = 2.0;
        C4 = 4.0;
        varargout{1} = C1*EV_(1)+C2*EV_(2)+C3*EV_(3)+C4*EV_(4);
        if(nargout>1)
            g_(1,1) = C1;
            g_(2,1) = C2;
            g_(3,1) = C3;
            g_(4,1) = C4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 0.0;
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

