function varargout = LOGHAIRY(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LOGHAIRY
%    *********
% 
%    A more difficult variant of the HAIRY problem in two variables.  
%    It is defined by a logarithmic transformation of the HAIRY surface,
%    which is defined by this function has a large number of relatively
%    sharp hills between  which a valley leads to the minimizer.
%    This problem contains a large number of saddle points.
% 
%    The problem is nonconvex.
% 
%    Source:
%    Ph. Toint, private communication,
% 
%    SIF input: Ph. Toint, April 1997.
% 
%    classification = 'OUR2-AN-2-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LOGHAIRY';

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
        v_('HLENGTH') = 30.0;
        v_('CSLOPE') = 100.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','FURCUP',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = -500.0;
        pb.x0(ix_('X2'),1) = -700.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eFUR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'DENS';
        [it,iet_] = s2mpjlib( 'ii', 'eDCUP',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'SMOOTH';
        [it,iet_] = s2mpjlib( 'ii', 'en1CUP',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'SMOOTH';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'HAIR';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eFUR';
        ielftype(ie) = iet_('eFUR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('DENS',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.0;
        ename = 'DBOWL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDCUP';
        ielftype(ie) = iet_('eDCUP');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('SMOOTH',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.01;
        ename = '1BOWL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en1CUP';
        ielftype(ie) = iet_('en1CUP');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('SMOOTH',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.01;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('FURCUP');
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('HAIR');
        pbm.grelw{ig}(posel) = v_('HLENGTH');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('DBOWL');
        pbm.grelw{ig}(posel) = v_('CSLOPE');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('1BOWL');
        pbm.grelw{ig}(posel) = v_('CSLOPE');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.1823216
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

    case 'eFUR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DV1 = pbm.elpar{iel_}(1)*EV_(1);
        DV2 = pbm.elpar{iel_}(1)*EV_(2);
        TDV1 = DV1+DV1;
        TDV2 = DV2+DV2;
        TDL2 = 2.0*pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        S1SQ = sin(DV1)^2;
        C2SQ = cos(DV2)^2;
        STDV1 = sin(TDV1);
        STDV2 = sin(TDV2);
        varargout{1} = S1SQ*C2SQ;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*STDV1*C2SQ;
            g_(2,1) = -pbm.elpar{iel_}(1)*S1SQ*STDV2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = TDL2*cos(TDV1)*C2SQ;
                H_(1,2) = -pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*STDV1*STDV2;
                H_(2,1) = H_(1,2);
                H_(2,2) = -TDL2*S1SQ*cos(TDV2);
                varargout{3} = H_;
            end
        end

    case 'eDCUP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        VSQ = IV_(1)*IV_(1);
        ARG = pbm.elpar{iel_}(1)+VSQ;
        SQARG = sqrt(ARG);
        DEN = 1.0/SQARG;
        varargout{1} = SQARG;
        if(nargout>1)
            g_(1,1) = IV_(1)*DEN;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = (1.0-VSQ/ARG)*DEN;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'en1CUP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        VSQ = EV_(1)*EV_(1);
        ARG = pbm.elpar{iel_}(1)+VSQ;
        SQARG = sqrt(ARG);
        DEN = 1.0/SQARG;
        varargout{1} = SQARG;
        if(nargout>1)
            g_(1,1) = EV_(1)*DEN;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = (1.0-VSQ/ARG)*DEN;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gLOG'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        S = 1.0e2;
        varargout{1} = log((S+GVAR_)/S);
        if(nargout>1)
            g_ = 1.0/(S+GVAR_);
            varargout{2} = g_;
            if(nargout>2)
                H_ = -1.0/(S+GVAR_)^2;
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

