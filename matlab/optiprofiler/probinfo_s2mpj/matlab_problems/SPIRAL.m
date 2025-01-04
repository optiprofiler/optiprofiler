function varargout = SPIRAL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPIRAL
%    *********
% 
%    A nonlinear minmax problem.
% 
%    Source:
%    E. Polak, J.E. Higgins and D. Mayne,
%    "A barrier function for minmax problems",
%    Mathematical Programming, vol.54(2), pp. 155-176, 1992.
% 
%    SIF input: Ph. Toint, April 1992.
% 
%    classification = 'C-CLOR2-AN-3-2'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SPIRAL';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U');
        valA(end+1) = 1.0;
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
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('X1'),1) = 1.41831;
        pb.x0(ix_('X2'),1) = -4.79462;
        pb.x0(ix_('U'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eBADCOS',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eBADSIN',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'X1SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X2SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'BC';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBADCOS';
        ielftype(ie) = iet_('eBADCOS');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'BS';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBADSIN';
        ielftype(ie) = iet_('eBADSIN');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('BC');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X1SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.005;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X2SQ');
        pbm.grelw{ig}(posel) = -0.005;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('BS');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X1SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.005;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X2SQ');
        pbm.grelw{ig}(posel) = -0.005;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN              0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-AN-3-2';
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

    case 'eBADCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XX = EV_(1)*EV_(1);
        YY = EV_(2)*EV_(2);
        R = sqrt(XX+YY);
        DRDX = EV_(1)/R;
        DRDY = EV_(2)/R;
        R3 = R^3;
        D2RDXX = 1.0/R-XX/R3;
        D2RDYY = 1.0/R-YY/R3;
        D2RDXY = -EV_(1)*EV_(2)/R3;
        C = cos(R);
        S = sin(R);
        DCDX = -S*DRDX;
        DCDY = -S*DRDY;
        D2CDXX = -C*DRDX*DRDX-S*D2RDXX;
        D2CDYY = -C*DRDY*DRDY-S*D2RDYY;
        D2CDXY = -C*DRDX*DRDY-S*D2RDXY;
        DSDX = C*DRDX;
        DSDY = C*DRDY;
        D2SDXX = -S*DRDX*DRDX+C*D2RDXX;
        D2SDYY = -S*DRDY*DRDY+C*D2RDYY;
        D2SDXY = -S*DRDX*DRDY+C*D2RDXY;
        Z = EV_(1)-R*C;
        DZDX = 1.0-DRDX*C-R*DCDX;
        DZDY = -DRDY*C-R*DCDY;
        D2ZDXX = -D2RDXX*C-2.0*DRDX*DCDX-R*D2CDXX;
        D2ZDYY = -D2RDYY*C-2.0*DRDY*DCDY-R*D2CDYY;
        D2ZDXY = -D2RDXY*C-DRDX*DCDY-DRDY*DCDX-R*D2CDXY;
        varargout{1} = Z*Z;
        if(nargout>1)
            g_(1,1) = 2.0*DZDX*Z;
            g_(2,1) = 2.0*DZDY*Z;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*(D2ZDXX*Z+DZDX*DZDX);
                H_(1,2) = 2.0*(D2ZDXY*Z+DZDX*DZDY);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*(D2ZDYY*Z+DZDY*DZDY);
                varargout{3} = H_;
            end
        end

    case 'eBADSIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XX = EV_(1)*EV_(1);
        YY = EV_(2)*EV_(2);
        R = sqrt(XX+YY);
        DRDX = EV_(1)/R;
        DRDY = EV_(2)/R;
        R3 = R^3;
        D2RDXX = 1.0/R-XX/R3;
        D2RDYY = 1.0/R-YY/R3;
        D2RDXY = -EV_(1)*EV_(2)/R3;
        C = cos(R);
        S = sin(R);
        DCDX = -S*DRDX;
        DCDY = -S*DRDY;
        D2CDXX = -C*DRDX*DRDX-S*D2RDXX;
        D2CDYY = -C*DRDY*DRDY-S*D2RDYY;
        D2CDXY = -C*DRDX*DRDY-S*D2RDXY;
        DSDX = C*DRDX;
        DSDY = C*DRDY;
        D2SDXX = -S*DRDX*DRDX+C*D2RDXX;
        D2SDYY = -S*DRDY*DRDY+C*D2RDYY;
        D2SDXY = -S*DRDX*DRDY+C*D2RDXY;
        Z = EV_(2)-R*S;
        DZDX = -DRDX*S-R*DSDX;
        DZDY = 1.0-DRDY*S-R*DSDY;
        D2ZDXX = -D2RDXX*S-2.0*DRDX*DSDX-R*D2SDXX;
        D2ZDYY = -D2RDYY*S-2.0*DRDY*DSDY-R*D2SDYY;
        D2ZDXY = -D2RDXY*S-DRDX*DSDY-DRDY*DSDX-R*D2SDXY;
        varargout{1} = Z*Z;
        if(nargout>1)
            g_(1,1) = 2.0*DZDX*Z;
            g_(2,1) = 2.0*DZDY*Z;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*(D2ZDXX*Z+DZDX*DZDX);
                H_(1,2) = 2.0*(D2ZDXY*Z+DZDX*DZDY);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*(D2ZDYY*Z+DZDY*DZDY);
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

