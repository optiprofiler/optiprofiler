function varargout = WATSON(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : WATSON
%    *********
% 
%    Watson problem in 12 variables.
% 
%    This function  is a nonlinear least squares with 31 groups.  Each
%    group has 1 nonlinear and 1 linear elements.
% 
%    Source:  problem 20 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#128 (p. 100).
% 
%    SIF input: Ph. Toint, Dec 1989.
%    (bug fix July 2007)
% 
%    classification = 'SUR2-AN-V-0'
% 
%    The number of variables can be varied, but should be smaller than
%    31
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   12             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'WATSON';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargs<1)
            v_('N') = 12;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   31             $-PARAMETER
        v_('M') = 31;
        v_('1') = 1;
        v_('2') = 2;
        v_('29') = 29;
        v_('30') = 30;
        v_('29') = 29.0;
        v_('1/29') = 1.0/v_('29');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('29')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('1/29');
            v_('LNTI') = log(v_('TI'));
            for J=v_('2'):v_('N')
                v_('RJ') = J;
                v_('RJ-1') = -1.0+v_('RJ');
                v_('RJ-2') = -2.0+v_('RJ');
                v_('AE') = v_('RJ-2')*v_('LNTI');
                v_('C0') = exp(v_('AE'));
                v_('C') = v_('C0')*v_('RJ-1');
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('C');
                end
            end
        end
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('30')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('M')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0*ones(ngrp,1);
        pbm.gconst(ig_(['G',int2str(round(v_('30')))])) = 0.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eMWSQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        elftv{it}{9} = 'V9';
        elftv{it}{10} = 'V10';
        elftv{it}{11} = 'V11';
        elftv{it}{12} = 'V12';
        elftp{it}{1} = 'T1';
        elftp{it}{2} = 'T2';
        elftp{it}{3} = 'T3';
        elftp{it}{4} = 'T4';
        elftp{it}{5} = 'T5';
        elftp{it}{6} = 'T6';
        elftp{it}{7} = 'T7';
        elftp{it}{8} = 'T8';
        elftp{it}{9} = 'T9';
        elftp{it}{10} = 'T10';
        elftp{it}{11} = 'T11';
        elftp{it}{12} = 'T12';
        [it,iet_] = s2mpjlib( 'ii', 'eMSQ',iet_);
        elftv{it}{1} = 'V1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('29')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('1/29');
            v_('LNTI') = log(v_('TI'));
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eMWSQ';
            ielftype(ie) = iet_('eMWSQ');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X8';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X9';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V9',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X10';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V10',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X11';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V11',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X12';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V12',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            for J=v_('1'):v_('N')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('CE0') = v_('RJ-1')*v_('LNTI');
                v_(['CE',int2str(J)]) = exp(v_('CE0'));
            end
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('T1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE1');
            [~,posep] = ismember('T2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE2');
            [~,posep] = ismember('T3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE3');
            [~,posep] = ismember('T4',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE4');
            [~,posep] = ismember('T5',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE5');
            [~,posep] = ismember('T6',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE6');
            [~,posep] = ismember('T7',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE7');
            [~,posep] = ismember('T8',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE8');
            [~,posep] = ismember('T9',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE9');
            [~,posep] = ismember('T10',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE10');
            [~,posep] = ismember('T11',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE11');
            [~,posep] = ismember('T12',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('CE12');
        end
        ename = ['E',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eMSQ';
        ielftype(ie) = iet_('eMSQ');
        ename = ['E',int2str(round(v_('M')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('29')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_(['G',int2str(round(v_('M')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('M')))]);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(12)           2.27559922D-9
% LO SOLTN(31)           1.53795068D-9
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-0';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eMSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = -EV_(1)-EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -2.0;
                varargout{3} = H_;
            end
        end

    case 'eMWSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U = pbm.elpar{iel_}(1)*EV_(1)+pbm.elpar{iel_}(2)*EV_(2)+pbm.elpar{iel_}(3)*...
             EV_(3)+pbm.elpar{iel_}(4)*EV_(4)+pbm.elpar{iel_}(5)*EV_(5)+pbm.elpar{iel_}(6)*EV_(6)+pbm.elpar{iel_}(7)*EV_(7)+pbm.elpar{iel_}(8)*EV_(8)+pbm.elpar{iel_}(9)*EV_(9)+pbm.elpar{iel_}(10)*EV_(10)+pbm.elpar{iel_}(11)*EV_(11)+pbm.elpar{iel_}(12)*EV_(12);
        TWOT1 = pbm.elpar{iel_}(1)+pbm.elpar{iel_}(1);
        TWOT2 = pbm.elpar{iel_}(2)+pbm.elpar{iel_}(2);
        TWOT3 = pbm.elpar{iel_}(3)+pbm.elpar{iel_}(3);
        TWOT4 = pbm.elpar{iel_}(4)+pbm.elpar{iel_}(4);
        TWOT5 = pbm.elpar{iel_}(5)+pbm.elpar{iel_}(5);
        TWOT6 = pbm.elpar{iel_}(6)+pbm.elpar{iel_}(6);
        TWOT7 = pbm.elpar{iel_}(7)+pbm.elpar{iel_}(7);
        TWOT8 = pbm.elpar{iel_}(8)+pbm.elpar{iel_}(8);
        TWOT9 = pbm.elpar{iel_}(9)+pbm.elpar{iel_}(9);
        TWOT10 = pbm.elpar{iel_}(10)+pbm.elpar{iel_}(10);
        TWOT11 = pbm.elpar{iel_}(11)+pbm.elpar{iel_}(11);
        TWOT12 = pbm.elpar{iel_}(12)+pbm.elpar{iel_}(12);
        varargout{1} = -U*U;
        if(nargout>1)
            g_(1,1) = -TWOT1*U;
            g_(2,1) = -TWOT2*U;
            g_(3,1) = -TWOT3*U;
            g_(4,1) = -TWOT4*U;
            g_(5,1) = -TWOT5*U;
            g_(6,1) = -TWOT6*U;
            g_(7,1) = -TWOT7*U;
            g_(8,1) = -TWOT8*U;
            g_(9,1) = -TWOT9*U;
            g_(10,1) = -TWOT10*U;
            g_(11,1) = -TWOT11*U;
            g_(12,1) = -TWOT12*U;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(12,12);
                H_(1,1) = -TWOT1*pbm.elpar{iel_}(1);
                H_(1,2) = -TWOT1*pbm.elpar{iel_}(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = -TWOT1*pbm.elpar{iel_}(3);
                H_(3,1) = H_(1,3);
                H_(1,4) = -TWOT1*pbm.elpar{iel_}(4);
                H_(4,1) = H_(1,4);
                H_(1,5) = -TWOT1*pbm.elpar{iel_}(5);
                H_(5,1) = H_(1,5);
                H_(1,6) = -TWOT1*pbm.elpar{iel_}(6);
                H_(6,1) = H_(1,6);
                H_(1,7) = -TWOT1*pbm.elpar{iel_}(7);
                H_(7,1) = H_(1,7);
                H_(1,8) = -TWOT1*pbm.elpar{iel_}(8);
                H_(8,1) = H_(1,8);
                H_(1,9) = -TWOT1*pbm.elpar{iel_}(9);
                H_(9,1) = H_(1,9);
                H_(1,10) = -TWOT1*pbm.elpar{iel_}(10);
                H_(10,1) = H_(1,10);
                H_(1,11) = -TWOT1*pbm.elpar{iel_}(11);
                H_(11,1) = H_(1,11);
                H_(1,12) = -TWOT1*pbm.elpar{iel_}(12);
                H_(12,1) = H_(1,12);
                H_(2,2) = -TWOT2*pbm.elpar{iel_}(2);
                H_(2,3) = -TWOT2*pbm.elpar{iel_}(3);
                H_(3,2) = H_(2,3);
                H_(2,4) = -TWOT2*pbm.elpar{iel_}(4);
                H_(4,2) = H_(2,4);
                H_(2,5) = -TWOT2*pbm.elpar{iel_}(5);
                H_(5,2) = H_(2,5);
                H_(2,6) = -TWOT2*pbm.elpar{iel_}(6);
                H_(6,2) = H_(2,6);
                H_(2,7) = -TWOT2*pbm.elpar{iel_}(7);
                H_(7,2) = H_(2,7);
                H_(2,8) = -TWOT2*pbm.elpar{iel_}(8);
                H_(8,2) = H_(2,8);
                H_(2,9) = -TWOT2*pbm.elpar{iel_}(8);
                H_(9,2) = H_(2,9);
                H_(2,10) = -TWOT2*pbm.elpar{iel_}(10);
                H_(10,2) = H_(2,10);
                H_(2,11) = -TWOT2*pbm.elpar{iel_}(11);
                H_(11,2) = H_(2,11);
                H_(2,12) = -TWOT2*pbm.elpar{iel_}(12);
                H_(12,2) = H_(2,12);
                H_(3,3) = -TWOT3*pbm.elpar{iel_}(3);
                H_(3,4) = -TWOT3*pbm.elpar{iel_}(4);
                H_(4,3) = H_(3,4);
                H_(3,5) = -TWOT3*pbm.elpar{iel_}(5);
                H_(5,3) = H_(3,5);
                H_(3,6) = -TWOT3*pbm.elpar{iel_}(6);
                H_(6,3) = H_(3,6);
                H_(3,7) = -TWOT3*pbm.elpar{iel_}(7);
                H_(7,3) = H_(3,7);
                H_(3,8) = -TWOT3*pbm.elpar{iel_}(8);
                H_(8,3) = H_(3,8);
                H_(3,9) = -TWOT3*pbm.elpar{iel_}(8);
                H_(9,3) = H_(3,9);
                H_(3,10) = -TWOT3*pbm.elpar{iel_}(10);
                H_(10,3) = H_(3,10);
                H_(3,11) = -TWOT3*pbm.elpar{iel_}(11);
                H_(11,3) = H_(3,11);
                H_(3,12) = -TWOT3*pbm.elpar{iel_}(12);
                H_(12,3) = H_(3,12);
                H_(4,4) = -TWOT4*pbm.elpar{iel_}(4);
                H_(4,5) = -TWOT4*pbm.elpar{iel_}(5);
                H_(5,4) = H_(4,5);
                H_(4,6) = -TWOT4*pbm.elpar{iel_}(6);
                H_(6,4) = H_(4,6);
                H_(4,7) = -TWOT4*pbm.elpar{iel_}(7);
                H_(7,4) = H_(4,7);
                H_(4,8) = -TWOT4*pbm.elpar{iel_}(8);
                H_(8,4) = H_(4,8);
                H_(4,9) = -TWOT4*pbm.elpar{iel_}(8);
                H_(9,4) = H_(4,9);
                H_(4,10) = -TWOT4*pbm.elpar{iel_}(10);
                H_(10,4) = H_(4,10);
                H_(4,11) = -TWOT4*pbm.elpar{iel_}(11);
                H_(11,4) = H_(4,11);
                H_(4,12) = -TWOT4*pbm.elpar{iel_}(12);
                H_(12,4) = H_(4,12);
                H_(5,5) = -TWOT5*pbm.elpar{iel_}(5);
                H_(5,6) = -TWOT5*pbm.elpar{iel_}(6);
                H_(6,5) = H_(5,6);
                H_(5,7) = -TWOT5*pbm.elpar{iel_}(7);
                H_(7,5) = H_(5,7);
                H_(5,8) = -TWOT5*pbm.elpar{iel_}(8);
                H_(8,5) = H_(5,8);
                H_(5,9) = -TWOT5*pbm.elpar{iel_}(8);
                H_(9,5) = H_(5,9);
                H_(5,10) = -TWOT5*pbm.elpar{iel_}(10);
                H_(10,5) = H_(5,10);
                H_(5,11) = -TWOT5*pbm.elpar{iel_}(11);
                H_(11,5) = H_(5,11);
                H_(5,12) = -TWOT5*pbm.elpar{iel_}(12);
                H_(12,5) = H_(5,12);
                H_(6,6) = -TWOT6*pbm.elpar{iel_}(6);
                H_(6,7) = -TWOT6*pbm.elpar{iel_}(7);
                H_(7,6) = H_(6,7);
                H_(6,8) = -TWOT6*pbm.elpar{iel_}(8);
                H_(8,6) = H_(6,8);
                H_(6,9) = -TWOT6*pbm.elpar{iel_}(8);
                H_(9,6) = H_(6,9);
                H_(6,10) = -TWOT6*pbm.elpar{iel_}(10);
                H_(10,6) = H_(6,10);
                H_(6,11) = -TWOT6*pbm.elpar{iel_}(11);
                H_(11,6) = H_(6,11);
                H_(6,12) = -TWOT6*pbm.elpar{iel_}(12);
                H_(12,6) = H_(6,12);
                H_(7,7) = -TWOT7*pbm.elpar{iel_}(7);
                H_(7,8) = -TWOT7*pbm.elpar{iel_}(8);
                H_(8,7) = H_(7,8);
                H_(7,9) = -TWOT7*pbm.elpar{iel_}(8);
                H_(9,7) = H_(7,9);
                H_(7,10) = -TWOT7*pbm.elpar{iel_}(10);
                H_(10,7) = H_(7,10);
                H_(7,11) = -TWOT7*pbm.elpar{iel_}(11);
                H_(11,7) = H_(7,11);
                H_(7,12) = -TWOT7*pbm.elpar{iel_}(12);
                H_(12,7) = H_(7,12);
                H_(8,8) = -TWOT8*pbm.elpar{iel_}(8);
                H_(8,9) = -TWOT8*pbm.elpar{iel_}(8);
                H_(9,8) = H_(8,9);
                H_(8,10) = -TWOT8*pbm.elpar{iel_}(10);
                H_(10,8) = H_(8,10);
                H_(8,11) = -TWOT8*pbm.elpar{iel_}(11);
                H_(11,8) = H_(8,11);
                H_(8,12) = -TWOT8*pbm.elpar{iel_}(12);
                H_(12,8) = H_(8,12);
                H_(9,9) = -TWOT9*pbm.elpar{iel_}(9);
                H_(9,10) = -TWOT9*pbm.elpar{iel_}(10);
                H_(10,9) = H_(9,10);
                H_(9,11) = -TWOT9*pbm.elpar{iel_}(11);
                H_(11,9) = H_(9,11);
                H_(9,12) = -TWOT9*pbm.elpar{iel_}(12);
                H_(12,9) = H_(9,12);
                H_(10,10) = -TWOT10*pbm.elpar{iel_}(10);
                H_(10,11) = -TWOT10*pbm.elpar{iel_}(11);
                H_(11,10) = H_(10,11);
                H_(10,12) = -TWOT10*pbm.elpar{iel_}(12);
                H_(12,10) = H_(10,12);
                H_(11,11) = -TWOT11*pbm.elpar{iel_}(11);
                H_(11,12) = -TWOT11*pbm.elpar{iel_}(12);
                H_(12,11) = H_(11,12);
                H_(12,12) = -TWOT12*pbm.elpar{iel_}(12);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

