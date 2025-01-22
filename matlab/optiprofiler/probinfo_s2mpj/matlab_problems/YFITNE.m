function varargout = YFITNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    This problem arises in measuring angles and distances to a vibrating beam
%    using a laser-Doppler velocimeter.
%    This is a nonlinear equation variant of the bounded constrained
%    problem YFIT.
% 
%    Source:
%    an exercize for L. Watson course on LANCELOT in the Spring 1993.
% 
%    SIF input: B. E. Lindholm, Virginia Tech., Spring 1993,
%               modified by Ph. Toint, March 1994.
%               derivatives corrected by Nick Gould, June 2019.
% 
%    classification = 'C-CNOR2-MN-3-17'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'YFITNE';

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
        v_('zero') = 0;
        v_('p') = 16;
        v_('realp') = 16.0;
        v_('y0') = 21.158931;
        v_('y1') = 17.591719;
        v_('y2') = 14.046854;
        v_('y3') = 10.519732;
        v_('y4') = 7.0058392;
        v_('y5') = 3.5007293;
        v_('y6') = 0.0000000;
        v_('y7') = -3.5007293;
        v_('y8') = -7.0058392;
        v_('y9') = -10.519732;
        v_('y10') = -14.046854;
        v_('y11') = -17.591719;
        v_('y12') = -21.158931;
        v_('y13') = -24.753206;
        v_('y14') = -28.379405;
        v_('y15') = -32.042552;
        v_('y16') = -35.747869;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','alpha',ix_);
        pb.xnames{iv} = 'alpha';
        [iv,ix_] = s2mpjlib('ii','beta',ix_);
        pb.xnames{iv} = 'beta';
        [iv,ix_] = s2mpjlib('ii','dist',ix_);
        pb.xnames{iv} = 'dist';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for i=v_('zero'):v_('p')
            [ig,ig_] = s2mpjlib('ii',['diff',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['diff',int2str(i)];
        end
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
        for i=v_('zero'):v_('p')
            pbm.gconst(ig_(['diff',int2str(i)])) = v_(['y',int2str(i)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('alpha'),1) = 0.60;
        pb.x0(ix_('beta'),1) = -0.60;
        pb.x0(ix_('dist'),1) = 20.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'etanab',iet_);
        elftv{it}{1} = 'a1';
        elftv{it}{2} = 'b1';
        elftv{it}{3} = 'd1';
        elftp{it}{1} = 'point';
        elftp{it}{2} = 'count';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for i=v_('zero'):v_('p')
            v_('index') = i;
            ename = ['est',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'etanab';
            ielftype(ie) = iet_('etanab');
            vname = 'alpha';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('a1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'beta';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('b1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'dist';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('d1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('point',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('index');
            [~,posep] = ismember('count',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('realp');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for i=v_('zero'):v_('p')
            ig = ig_(['diff',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['est',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-3-17';
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

    case 'etanab'

        EV_  = varargin{1};
        iel_ = varargin{2};
        frac = pbm.elpar{iel_}(1)/pbm.elpar{iel_}(2);
        ttan = tan(EV_(1)*(1.0-frac)+EV_(2)*frac);
        tsec = 1.0/cos(EV_(1)*(1.0-frac)+EV_(2)*frac);
        tsec2 = tsec*tsec;
        varargout{1} = EV_(3)*ttan;
        if(nargout>1)
            g_(1,1) = EV_(3)*(1.0-frac)*tsec2;
            g_(2,1) = EV_(3)*frac*tsec2;
            g_(3,1) = ttan;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0*EV_(3)*((1.0-frac)^2)*tsec2*ttan;
                H_(2,2) = 2.0*EV_(3)*(frac^2)*tsec2*ttan;
                H_(1,2) = 2.0*EV_(3)*(1.0-frac)*frac*tsec2*ttan;
                H_(2,1) = H_(1,2);
                H_(1,3) = (1.0-frac)*tsec2;
                H_(3,1) = H_(1,3);
                H_(2,3) = frac*tsec2;
                H_(3,2) = H_(2,3);
                H_(3,3) = 0.0;
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

