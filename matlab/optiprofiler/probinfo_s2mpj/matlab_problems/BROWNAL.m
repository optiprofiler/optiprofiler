function varargout = BROWNAL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BROWNAL
%    *********
%    Brown almost linear least squares problem.
%    This problem is a sum of n least-squares groups, the last one of
%    which has a nonlinear element.
%    It Hessian matrix is dense.
% 
%    Source: Problem 27 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#79
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SUR2-AN-V-0'
% 
%    N is the number of free variables (variable).
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   100            $-PARAMETER
% IE N                   200            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BROWNAL';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   1000           $-PARAMETER
        v_('1') = 1;
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        v_('RN+1') = v_('N+1');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('1'):v_('I-1')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 2.0;
            end
            for J=v_('I+1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['G',int2str(I)])) = v_('RN+1');
        end
        pbm.gconst(ig_(['G',int2str(round(v_('N')))])) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
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
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
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
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_(['G',int2str(round(v_('N')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-0';
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

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V12 = EV_(1)*EV_(2);
        V34 = EV_(3)*EV_(4);
        V56 = EV_(5)*EV_(6);
        V78 = EV_(7)*EV_(8);
        V910 = EV_(9)*EV_(10);
        varargout{1} = V12*V34*V56*V78*V910;
        if(nargout>1)
            g_(1,1) = EV_(2)*V34*V56*V78*V910;
            g_(2,1) = EV_(1)*V34*V56*V78*V910;
            g_(3,1) = V12*EV_(4)*V56*V78*V910;
            g_(4,1) = V12*EV_(3)*V56*V78*V910;
            g_(5,1) = V12*V34*EV_(6)*V78*V910;
            g_(6,1) = V12*V34*EV_(5)*V78*V910;
            g_(7,1) = V12*V34*V56*EV_(8)*V910;
            g_(8,1) = V12*V34*V56*EV_(7)*V910;
            g_(9,1) = V12*V34*V56*V78*EV_(10);
            g_(10,1) = V12*V34*V56*V78*EV_(9);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(10,10);
                H_(1,2) = V34*V56*V78*V910;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4)*V56*V78*V910;
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3)*V56*V78*V910;
                H_(4,1) = H_(1,4);
                H_(1,5) = EV_(2)*V34*EV_(6)*V78*V910;
                H_(5,1) = H_(1,5);
                H_(1,6) = EV_(2)*V34*EV_(5)*V78*V910;
                H_(6,1) = H_(1,6);
                H_(1,7) = EV_(2)*V34*V56*EV_(8)*V910;
                H_(7,1) = H_(1,7);
                H_(1,8) = EV_(2)*V34*V56*EV_(7)*V910;
                H_(8,1) = H_(1,8);
                H_(1,9) = EV_(2)*V34*V56*V78*EV_(10);
                H_(9,1) = H_(1,9);
                H_(1,10) = EV_(2)*V34*V56*V78*EV_(9);
                H_(10,1) = H_(1,10);
                H_(2,3) = EV_(1)*EV_(4)*V56*V78*V910;
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3)*V56*V78*V910;
                H_(4,2) = H_(2,4);
                H_(2,5) = EV_(1)*V34*EV_(6)*V78*V910;
                H_(5,2) = H_(2,5);
                H_(2,6) = EV_(1)*V34*EV_(5)*V78*V910;
                H_(6,2) = H_(2,6);
                H_(2,7) = EV_(1)*V34*V56*EV_(8)*V910;
                H_(7,2) = H_(2,7);
                H_(2,8) = EV_(1)*V34*V56*EV_(7)*V910;
                H_(8,2) = H_(2,8);
                H_(2,9) = EV_(1)*V34*V56*V78*EV_(10);
                H_(9,2) = H_(2,9);
                H_(2,10) = EV_(1)*V34*V56*V78*EV_(9);
                H_(10,2) = H_(2,10);
                H_(3,4) = V12*V56*V78*V910;
                H_(4,3) = H_(3,4);
                H_(3,5) = V12*EV_(4)*EV_(6)*V78*V910;
                H_(5,3) = H_(3,5);
                H_(3,6) = V12*EV_(4)*EV_(5)*V78*V910;
                H_(6,3) = H_(3,6);
                H_(3,7) = V12*EV_(4)*V56*EV_(8)*V910;
                H_(7,3) = H_(3,7);
                H_(3,8) = V12*EV_(4)*V56*EV_(7)*V910;
                H_(8,3) = H_(3,8);
                H_(3,9) = V12*EV_(4)*V56*V78*EV_(10);
                H_(9,3) = H_(3,9);
                H_(3,10) = V12*EV_(4)*V56*V78*EV_(9);
                H_(10,3) = H_(3,10);
                H_(4,5) = V12*EV_(3)*EV_(6)*V78*V910;
                H_(5,4) = H_(4,5);
                H_(4,6) = V12*EV_(3)*EV_(5)*V78*V910;
                H_(6,4) = H_(4,6);
                H_(4,7) = V12*EV_(3)*V56*EV_(8)*V910;
                H_(7,4) = H_(4,7);
                H_(4,8) = V12*EV_(3)*V56*EV_(7)*V910;
                H_(8,4) = H_(4,8);
                H_(4,9) = V12*EV_(3)*V56*V78*EV_(10);
                H_(9,4) = H_(4,9);
                H_(4,10) = V12*EV_(3)*V56*V78*EV_(9);
                H_(10,4) = H_(4,10);
                H_(5,6) = V12*V34*V78*V910;
                H_(6,5) = H_(5,6);
                H_(5,7) = V12*V34*EV_(6)*EV_(8)*V910;
                H_(7,5) = H_(5,7);
                H_(5,8) = V12*V34*EV_(6)*EV_(7)*V910;
                H_(8,5) = H_(5,8);
                H_(5,9) = V12*V34*EV_(6)*V78*EV_(10);
                H_(9,5) = H_(5,9);
                H_(5,10) = V12*V34*EV_(6)*V78*EV_(9);
                H_(10,5) = H_(5,10);
                H_(6,7) = V12*V34*EV_(5)*EV_(8)*V910;
                H_(7,6) = H_(6,7);
                H_(6,8) = V12*V34*EV_(5)*EV_(7)*V910;
                H_(8,6) = H_(6,8);
                H_(6,9) = V12*V34*EV_(5)*V78*EV_(10);
                H_(9,6) = H_(6,9);
                H_(6,10) = V12*V34*EV_(5)*V78*EV_(9);
                H_(10,6) = H_(6,10);
                H_(7,8) = V12*V34*V56*V910;
                H_(8,7) = H_(7,8);
                H_(7,9) = V12*V34*V56*EV_(8)*EV_(10);
                H_(9,7) = H_(7,9);
                H_(7,10) = V12*V34*V56*EV_(8)*EV_(9);
                H_(10,7) = H_(7,10);
                H_(8,9) = V12*V34*V56*EV_(7)*EV_(10);
                H_(9,8) = H_(8,9);
                H_(8,10) = V12*V34*V56*EV_(7)*EV_(9);
                H_(10,8) = H_(8,10);
                H_(9,10) = V12*V34*V56*V78;
                H_(10,9) = H_(9,10);
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

