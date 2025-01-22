function varargout = DIPIGRI(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DIPIGRI
%    *********
% 
%    A problem proposed by Di Pillo and Grippo.
% 
%    Source:
%    G. Di Pillo and L. Grippo,
%    "An new augmented Lagrangian function for inequality constraints
%    in nonlinear programming problems",
%    JOTA, vol. 36, pp. 495-519, 1982.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'C-COOR2-AN-7-4'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DIPIGRI';

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
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        [iv,ix_] = s2mpjlib('ii','X5',ix_);
        pb.xnames{iv} = 'X5';
        [iv,ix_] = s2mpjlib('ii','X6',ix_);
        pb.xnames{iv} = 'X6';
        [iv,ix_] = s2mpjlib('ii','X7',ix_);
        pb.xnames{iv} = 'X7';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -10.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -8.0;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 5.0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 7.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 23.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -8.0;
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -11.0;
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
        pbm.gconst(ig_('C1')) = 127.0;
        pbm.gconst(ig_('C2')) = 282.0;
        pbm.gconst(ig_('C3')) = 196.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 1.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 2.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 4.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 4.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 0.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 1.0;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSSQ',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'S';
        [it,iet_] = s2mpjlib( 'ii', 'eP4',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eP6',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eOE',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eCE',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'O1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSSQ';
        ielftype(ie) = iet_('eSSQ');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('S',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 10.0;
        ename = 'O2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSSQ';
        ielftype(ie) = iet_('eSSQ');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('S',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 12.0;
        ename = 'O3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eP4';
        ielftype(ie) = iet_('eP4');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSSQ';
        ielftype(ie) = iet_('eSSQ');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('S',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 11.0;
        ename = 'O5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eP6';
        ielftype(ie) = iet_('eP6');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE';
        ielftype(ie) = iet_('eOE');
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eP4';
        ielftype(ie) = iet_('eP4');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C21';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C31';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C32';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C41';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCE';
        ielftype(ie) = iet_('eCE');
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
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('O2');
        pbm.grelw{ig}(posel) = 5.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('O4');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 10.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('O6');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C12');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C21');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 10.0;
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C31');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C32');
        pbm.grelw{ig}(posel) = 6.0;
        ig = ig_('C4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C41');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C21');
        pbm.grelw{ig}(posel) = 2.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               680.630
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-7-4';
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

    case 'eSSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-pbm.elpar{iel_}(1))^2;
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-pbm.elpar{iel_}(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

    case 'eP4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^4;
        if(nargout>1)
            g_(1,1) = 4.0*EV_(1)^3;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 12.0*EV_(1)^2;
                varargout{3} = H_;
            end
        end

    case 'eP6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^6;
        if(nargout>1)
            g_(1,1) = 6.0*EV_(1)^5;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 30.0*EV_(1)^4;
                varargout{3} = H_;
            end
        end

    case 'eOE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 7.0*EV_(1)^2+EV_(2)^4-4.0*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 14.0*EV_(1)-4.0*EV_(2);
            g_(2,1) = 4.0*EV_(2)^3-4.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 14.0;
                H_(1,2) = -4.0;
                H_(2,1) = H_(1,2);
                H_(2,2) = 12.0*EV_(2)^2;
                varargout{3} = H_;
            end
        end

    case 'eCE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 4.0*EV_(1)^2+EV_(2)^2-3.0*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 8.0*EV_(1)-3.0*EV_(2);
            g_(2,1) = 2.0*EV_(2)-3.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 8.0;
                H_(1,2) = -3.0;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0;
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

