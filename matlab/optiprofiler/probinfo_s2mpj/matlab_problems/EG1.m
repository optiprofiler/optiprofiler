function varargout = EG1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A simple nonlinear problem given as an example in Section 1.2.3 of
%    the LANCELOT Manual.
% 
%    Source:
%    A.R. Conn, N. Gould and Ph.L. Toint,
%    "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
%    (Release A)"
%    Springer Verlag, 1992.
% 
%    SIF input: N. Gould and Ph. Toint, June 1994.
% 
%    classification = 'OBR2-AY-3-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'EG1';

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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','GROUP1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','GROUP2',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','GROUP3',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = numEntries(ig_);
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        ig = ig_('GROUP1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        ig = ig_('GROUP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1')) = -Inf;
        pb.xupper(ix_('X1'),1) = +Inf;
        pb.xlower(ix_('X2'),1) = -1.0;
        pb.xlower(ix_('X3'),1) = 1.0;
        pb.xupper(ix_('X2')) = 1.0;
        pb.xupper(ix_('X3')) = 2.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eETYPE1',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'eETYPE2',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'G2E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eETYPE1';
        ielftype(ie) = iet_('eETYPE1');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'G3E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eETYPE2';
        ielftype(ie) = iet_('eETYPE2');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'G3E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eETYPE1';
        ielftype(ie) = iet_('eETYPE1');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gGTYPE1',igt_);
        [it,igt_] = s2mpjlib('ii','gGTYPE2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('GROUP1');
        pbm.grftype{ig} = 'gGTYPE1';
        ig = ig_('GROUP2');
        pbm.grftype{ig} = 'gGTYPE2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G2E1');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('GROUP3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G3E1');
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('G3E2');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-AY-3-0';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eETYPE1'

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

    case 'eETYPE2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        CS = cos(IV_(2));
        SN = sin(IV_(2));
        varargout{1} = IV_(1)*SN;
        if(nargout>1)
            g_(1,1) = SN;
            g_(2,1) = IV_(1)*CS;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = CS;
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)*SN;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gGTYPE1'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        TWO = 2.0;
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = TWO*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = TWO;
                varargout{3} = H_;
            end
        end

    case 'gGTYPE2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        ALPHA2 = GVAR_*GVAR_;
        varargout{1} = ALPHA2*ALPHA2;
        if(nargout>1)
            g_ = 4.0*ALPHA2*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*ALPHA2;
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

