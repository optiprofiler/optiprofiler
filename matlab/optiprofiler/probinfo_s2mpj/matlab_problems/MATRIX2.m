function varargout = MATRIX2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MATRIX2
%    *********
%    Find the pair of 2 by 2 symmetric matrices that are closest in
%    Frobenius norm, subject to the first (X) being positive semidefinite and
%    the second (Y) being negative semidefinite.
% 
%    Source : a problem made up for the SIF/LANCELOT user' manual.
% 
%    SIF input: Ph. Toint, Jan 91.
% 
%    classification = 'QOR2-AY-6-2'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MATRIX2';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X11',ix_);
        pb.xnames{iv} = 'X11';
        [iv,ix_] = s2mpjlib('ii','X12',ix_);
        pb.xnames{iv} = 'X12';
        [iv,ix_] = s2mpjlib('ii','X22',ix_);
        pb.xnames{iv} = 'X22';
        [iv,ix_] = s2mpjlib('ii','Y11',ix_);
        pb.xnames{iv} = 'Y11';
        [iv,ix_] = s2mpjlib('ii','Y12',ix_);
        pb.xnames{iv} = 'Y12';
        [iv,ix_] = s2mpjlib('ii','Y22',ix_);
        pb.xnames{iv} = 'Y22';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','FROBDIST',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','XPOSDEF',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'XPOSDEF';
        [ig,ig_] = s2mpjlib('ii','YNEGDEF',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'YNEGDEF';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X12')) = -Inf;
        pb.xupper(ix_('X12'),1) = +Inf;
        pb.xlower(ix_('Y12')) = -Inf;
        pb.xupper(ix_('Y12'),1) = +Inf;
        pb.xlower(ix_('Y11'),1) = -Inf;
        pb.xlower(ix_('Y22'),1) = -Inf;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'XY11SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XY12SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XY22SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'X22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X1122';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X12SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Y1122';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Y11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Y12SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'Y12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('FROBDIST');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XY11SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XY22SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XY12SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        ig = ig_('XPOSDEF');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X1122');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X12SQ');
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('YNEGDEF');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Y1122');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Y12SQ');
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QOR2-AY-6-2';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE ELEMENTS *
%  ROUTINE             *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eISQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
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

