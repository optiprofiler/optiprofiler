function varargout = DJTL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DJTL
%    *********
% 
%    Source: modified version of problem 19 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
%    that is meant to simulate the Lagrangian barrier objective function
%    for particular values of the shifts and multipliers
% 
%    SIF input: A.R. Conn August 1993
% 
%    classification = 'OUR2-AN-2-0'
% 
%    Define multipliers and shifts
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DJTL';

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
        v_('LL1') = 1.0;
        v_('LL2') = 1.0;
        v_('LL3') = 1.0;
        v_('LL4') = 1.0;
        v_('SL1') = 1.0;
        v_('SL2') = 1.0;
        v_('SL3') = 1.0;
        v_('SL4') = 1.0;
        v_('LU1') = 1.0;
        v_('LU2') = 1.0;
        v_('LU3') = 1.0;
        v_('LU4') = 1.0;
        v_('SU1') = 1.0;
        v_('SU2') = 1.0;
        v_('SU3') = 1.0;
        v_('SU4') = 1.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CONU1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CONL1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CONU2',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CONL2',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','BNDU1',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','BNDL1',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','BNDU2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','BNDL2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X2');
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('CONU1')) = -200.0;
        pbm.gconst(ig_('CONL1')) = 100.0;
        pbm.gconst(ig_('CONU2')) = 0.0;
        pbm.gconst(ig_('CONL2')) = -82.81;
        pbm.gconst(ig_('BNDU1')) = -100.0;
        pbm.gconst(ig_('BNDL1')) = 13.0;
        pbm.gconst(ig_('BNDU2')) = -100.0;
        pbm.gconst(ig_('BNDL2')) = 0.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 15.0;
        pb.x0(ix_('X2'),1) = 6.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCBm10',iet_);
        elftv{it}{1} = 'V1';
        [it,iet_] = s2mpjlib( 'ii', 'eCBm20',iet_);
        elftv{it}{1} = 'V1';
        [it,iet_] = s2mpjlib( 'ii', 'eSQm5',iet_);
        elftv{it}{1} = 'V1';
        [it,iet_] = s2mpjlib( 'ii', 'eSQm6',iet_);
        elftv{it}{1} = 'V1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCBm10';
        ielftype(ie) = iet_('eCBm10');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCBm20';
        ielftype(ie) = iet_('eCBm20');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQm5';
        ielftype(ie) = iet_('eSQm5');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQm5';
        ielftype(ie) = iet_('eSQm5');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQm6';
        ielftype(ie) = iet_('eSQm6');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        grftp{it}{1} = 'P1';
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        grftp{it}{2} = 'P2';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('CONL1');
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SL1');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LL1');
        ig = ig_('CONU1');
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = -1.0;
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SU1');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LU1');
        ig = ig_('CONL2');
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        pbm.grelw{ig}(posel) = -1.0;
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SL2');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LL2');
        ig = ig_('CONU2');
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        pbm.grelw{ig}(posel) = 1.;
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SU2');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LU2');
        ig = ig_('BNDL1');
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SL3');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LL3');
        ig = ig_('BNDU1');
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SU3');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LU3');
        ig = ig_('BNDL2');
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SL4');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LL4');
        ig = ig_('BNDU2');
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P1',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('SU4');
        [~,posgp] = ismember('P2',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('LU4');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -8951.54472
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-2-0';
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

    case 'eCBm10'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DIF = EV_(1)-10.0;
        varargout{1} = DIF^3;
        if(nargout>1)
            g_(1,1) = 3.0*DIF*DIF;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*DIF;
                varargout{3} = H_;
            end
        end

    case 'eCBm20'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DIF = EV_(1)-20.0;
        varargout{1} = DIF^3;
        if(nargout>1)
            g_(1,1) = 3.0*DIF*DIF;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*DIF;
                varargout{3} = H_;
            end
        end

    case 'eSQm5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DIF = EV_(1)-5.0;
        varargout{1} = DIF^2;
        if(nargout>1)
            g_(1,1) = 2.0*DIF;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eSQm6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DIF = EV_(1)-6.0;
        varargout{1} = DIF^2;
        if(nargout>1)
            g_(1,1) = 2.0*DIF;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gLOG'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        APP1 = GVAR_+pbm.grpar{igr_}(1);
        P1P2 = pbm.grpar{igr_}(1)*pbm.grpar{igr_}(2);
        ARG0 = APP1<=0.0;
        BIG = 1.0000e+10;
        if(ARG0)
            FF = BIG*GVAR_^2;
        end
        if(~ARG0)
            FF = -P1P2*log(APP1);
        end
        if(ARG0)
            GG = 2.0*BIG*GVAR_;
        end
        if(~ARG0)
            GG = -P1P2/APP1;
        end
        if(ARG0)
            HH = 2.0*BIG;
        end
        if(~ARG0)
            HH = P1P2/APP1^2;
        end
        varargout{1} = FF;
        if(nargout>1)
            g_ = GG;
            varargout{2} = g_;
            if(nargout>2)
                H_ = HH;
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

