function varargout = HS69(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS69
%    *********
% 
%    This is a cost optimal inspection plan.
% 
%    Source: problem 69 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'C-COOR2-MN-4-2'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS69';

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
        v_('N') = 4;
        v_('A') = 0.1;
        v_('B') = 1000.0;
        v_('D') = 1.0;
        v_('NN') = 4.0;
        v_('1') = 1;
        v_('AN') = v_('A')*v_('NN');
        v_('ROOTN') = sqrt(v_('NN'));
        v_('DROOTN') = v_('D')*v_('ROOTN');
        v_('-DROOTN') = -v_('DROOTN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0e+0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0e+0;
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.0001;
        pb.xupper(ix_('X1')) = 100.0;
        pb.xupper(ix_('X2')) = 100.0;
        pb.xupper(ix_('X3')) = 2.0;
        pb.xupper(ix_('X4')) = 2.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.y0 = 1.0*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eRECIP',iet_);
        elftv{it}{1} = 'X1';
        [it,iet_] = s2mpjlib( 'ii', 'eNASTYEXP',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X3';
        elftv{it}{3} = 'X4';
        elftp{it}{1} = 'B';
        [it,iet_] = s2mpjlib( 'ii', 'ePHI',iet_);
        elftv{it}{1} = 'X2';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OE1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eRECIP';
        ielftype(ie) = iet_('eRECIP');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OE2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eNASTYEXP';
        ielftype(ie) = iet_('eNASTYEXP');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('B',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('B');
        ename = 'C1E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePHI';
        ielftype(ie) = iet_('ePHI');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0e+0;
        ename = 'C2E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePHI';
        ielftype(ie) = iet_('ePHI');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('DROOTN');
        ename = 'C2E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePHI';
        ielftype(ie) = iet_('ePHI');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('-DROOTN');
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('AN');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0e+0;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.0e+0;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C2E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C2E2');
        pbm.grelw{ig}(posel) = -1.0e+0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION            -956.71288
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-MN-4-2';
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 3.9894228040143270e-01;
        varargout{1} = pbm;

    case 'eRECIP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 1.0e+0/EV_(1);
        if(nargout>1)
            g_(1,1) = -1.0e+0/EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0/EV_(1)^3;
                varargout{3} = H_;
            end
        end

    case 'eNASTYEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E = exp(EV_(1));
        R = pbm.elpar{iel_}(1)*(E-1.0e+0)-EV_(2);
        S = E-1.0e+0+EV_(3);
        F1 = -(-(R*EV_(3)/S))/(EV_(1)*EV_(1));
        F2 = -R/(S*EV_(1));
        F3 = -EV_(3)/(S*EV_(1));
        F4 = (EV_(3)*R)/(EV_(1)*S*S);
        F11 = 2.0e+0*(-(R*EV_(3)/S))/(EV_(1)*EV_(1)*EV_(1));
        F12 = E*EV_(3)*((R/S)-pbm.elpar{iel_}(1))/(EV_(1)*S);
        F11 = F11+F12;
        F12 = R/(S*EV_(1)*EV_(1));
        F13 = EV_(3)/(S*EV_(1)*EV_(1));
        F14 = R/(S*S*EV_(1));
        F15 = -1.0e+0/(S*EV_(1));
        F16 = -EV_(3)*R/(S*S*EV_(1)*EV_(1));
        F17 = EV_(3)/(S*S*EV_(1));
        F18 = -2.0e+0*EV_(3)*R/(EV_(1)*S*S*S);
        varargout{1} = EV_(3)*R/(S*EV_(1));
        if(nargout>1)
            g_(1,1) = -F1-(pbm.elpar{iel_}(1)*E*F3)-(E*F4);
            g_(2,1) = F3;
            g_(3,1) = -F2-F4;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = -F11-(2.0e+0*pbm.elpar{iel_}(1)*E*F13)-(2.0e+0*E*F16)...
                          -(2.0e+0*pbm.elpar{iel_}(1)*E*E*F17)-(E*E*F18);
                H_(2,1) = F13+(E*F17);
                H_(1,2) = H_(2,1);
                H_(3,1) = -(E*F14)-(pbm.elpar{iel_}(1)*E*F15)...
                          -F16-(pbm.elpar{iel_}(1)*E*F17)-(E*F18)-F12;
                H_(1,3) = H_(3,1);
                H_(3,2) = F15+F17;
                H_(2,3) = H_(3,2);
                H_(3,3) = -(2.0e+0*F14)-F18;
                varargout{3} = H_;
            end
        end

    case 'ePHI'

        x2    = varargin{1};
        iel_  = varargin{2};
        mx2pp = -x2+pbm.elpar{iel_}(1);
        E     = exp(-0.5*mx2pp^2);
        arg   =  7.071067811865475e-1 * abs( mx2pp );
        if ( mx2pp >= 0 )
           phi = 0.5 + 0.5 * erf( arg );
        else
           phi = 0.5 - 0.5 * erf( arg );
        end
        varargout{1} = phi;
        if(nargout>1)
            g_(1,1) = -pbm.efpar(1)*E;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -pbm.efpar(1)*E*mx2pp;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
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

