function varargout = TWOBARS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Structureal analysis of the simplest two bar scheme.  The structure has
%    the following simple symmetric shape
% 
%                                 *
%                                / \
%                               /   \
%                              /     \
%                            """     """
% 
%    and a force is applied at the top node.  The unknown are the distance
%    of the left and right feet wrt to the projection of the top node and the
%    weight of the bars.
% 
%    Source:
%    an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.
% 
%    SIF input: Ph. Toint, November 1994
% 
%    classification = 'C-COOR2-MN-2-2'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TWOBARS';

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
        [ig,ig_] = s2mpjlib('ii','CONS1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CONS1';
        [ig,ig_] = s2mpjlib('ii','CONS2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CONS2';
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
        pbm.gconst(ig_('CONS1')) = 1.0;
        pbm.gconst(ig_('CONS2')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.2;
        pb.xupper(ix_('X1')) = 4.0;
        pb.xlower(ix_('X2'),1) = 0.1;
        pb.xupper(ix_('X2')) = 1.6;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOE',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        [it,iet_] = s2mpjlib( 'ii', 'eCE1',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        [it,iet_] = s2mpjlib( 'ii', 'eCE2',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'OBEL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE';
        ielftype(ie) = iet_('eOE');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'COEL1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCE1';
        ielftype(ie) = iet_('eCE1');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'COEL2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCE2';
        ielftype(ie) = iet_('eCE2');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OBEL');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('CONS1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('COEL1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.124;
        ig = ig_('CONS2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('COEL2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.124;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               1.5086379655
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-MN-2-2';
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

    case 'eOE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = 1.0+EV_(2)*EV_(2);
        RA = sqrt(A);
        varargout{1} = EV_(1)*RA;
        if(nargout>1)
            g_(1,1) = RA;
            g_(2,1) = EV_(1)*EV_(2)/RA;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = EV_(2)/RA;
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)/(A*RA);
                varargout{3} = H_;
            end
        end

    case 'eCE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = 1.0+EV_(2)*EV_(2);
        RA = sqrt(A);
        B = 8.0/EV_(1);
        DB = -8.0/EV_(1)^2;
        D2B = 16.0/EV_(1)^3;
        C = 1.0/(EV_(1)*EV_(2));
        DCDX = -1.0/(EV_(1)^2*EV_(2));
        DCDY = -1.0/(EV_(2)^2*EV_(1));
        D2CDXX = 2.0/(EV_(1)^3*EV_(2));
        D2CDXY = 1.0/(EV_(1)*EV_(2))^2;
        D2CDYY = 2.0/(EV_(1)*EV_(2)^3);
        BC = B+C;
        varargout{1} = RA*BC;
        if(nargout>1)
            g_(1,1) = RA*(DB+DCDX);
            g_(2,1) = EV_(2)*BC/RA+RA*DCDY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = RA*(D2B+D2CDXX);
                H_(1,2) = RA*D2CDXY+EV_(2)*(DB+DCDX)/RA;
                H_(2,1) = H_(1,2);
                H_(2,2) = (BC+2.0*EV_(2)*DCDY-EV_(2)*EV_(2)*BC/A)/RA+RA*D2CDYY;
                varargout{3} = H_;
            end
        end

    case 'eCE2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = 1.0+EV_(2)*EV_(2);
        RA = sqrt(A);
        B = 8.0/EV_(1);
        DB = -8.0/EV_(1)^2;
        D2B = 16.0/EV_(1)^3;
        C = 1.0/(EV_(1)*EV_(2));
        DCDX = -1.0/(EV_(1)^2*EV_(2));
        DCDY = -1.0/(EV_(2)^2*EV_(1));
        D2CDXX = 2.0/(EV_(1)^3*EV_(2));
        D2CDXY = 1.0/(EV_(1)*EV_(2))^2;
        D2CDYY = 2.0/(EV_(1)*EV_(2)^3);
        BC = B-C;
        varargout{1} = RA*BC;
        if(nargout>1)
            g_(1,1) = RA*(DB-DCDX);
            g_(2,1) = EV_(2)*BC/RA-RA*DCDY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = RA*(D2B-D2CDXX);
                H_(1,2) = -RA*D2CDXY+EV_(2)*(DB-DCDX)/RA;
                H_(2,1) = H_(1,2);
                H_(2,2) = (BC-2.0*EV_(2)*DCDY-EV_(2)*EV_(2)*BC/A)/RA-RA*D2CDYY;
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

