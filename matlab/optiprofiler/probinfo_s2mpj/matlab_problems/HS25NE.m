function varargout = HS25NE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS25NE
%    *********
% 
%    A nonlinear least squares problem with bounds.
% 
%    Source: problem 25 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: J-M Collin, Mar 1990.
%    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
% 
%    classification = 'NOR2-AN-3-99'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS25NE';

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
        v_('N') = 3;
        v_('1') = 1;
        v_('99') = 99;
        v_('2/3') = 0.6666666666;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('99')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['O',int2str(I)];
        end
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('99')
            v_('IR') = I;
            v_('I/100') = 0.01*v_('IR');
            pbm.gconst(ig_(['O',int2str(I)])) = v_('I/100');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.1;
        pb.xupper(ix_('X1')) = 100.0;
        pb.xlower(ix_('X2'),1) = 0.0;
        pb.xupper(ix_('X2')) = 25.6;
        pb.xlower(ix_('X3'),1) = 0.0;
        pb.xupper(ix_('X3')) = 5.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 100.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 100.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 12.5;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 12.5;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 3.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 3.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eWFI',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        elftp{it}{1} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('99')
            v_('IR') = I;
            v_('I/100') = 0.01*v_('IR');
            v_('LOG01I') = log(v_('I/100'));
            v_('M50LOG') = -50.0*v_('LOG01I');
            v_('EXPLOG') = log(v_('M50LOG'));
            v_('EXPL2/3') = v_('EXPLOG')*v_('2/3');
            v_('EXP2/3') = exp(v_('EXPL2/3'));
            v_('UI') = 25.0+v_('EXP2/3');
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWFI';
            ielftype(ie) = iet_('eWFI');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('UI');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('99')
            ig = ig_(['O',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO HS25                0.0
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-3-99';
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

    case 'eWFI'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XI = 1.0/EV_(1);
        X2I = XI*XI;
        X3I = X2I*XI;
        WMY = pbm.elpar{iel_}(1)-EV_(2);
        WMYEZ = WMY^EV_(3);
        LWMY = log(WMY);
        EXPO = exp(-XI*WMYEZ);
        varargout{1} = EXPO;
        if(nargout>1)
            g_(1,1) = X2I*WMYEZ*EXPO;
            g_(2,1) = XI*EV_(3)*WMY^(EV_(3)-1.0)*EXPO;
            g_(3,1) = -XI*LWMY*WMYEZ*EXPO;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = EXPO*WMYEZ*X3I*(-2.0+XI*WMY^EV_(3));
                H_(1,2) = EXPO*EV_(3)*X2I*WMY^(EV_(3)-1.0)*(-1.0+XI*WMYEZ);
                H_(2,1) = H_(1,2);
                H_(1,3) = EXPO*X2I*WMYEZ*LWMY*(1.0-XI*WMYEZ);
                H_(3,1) = H_(1,3);
                H_(2,2) = EXPO*XI*WMY^(EV_(3)-2.0)*EV_(3)*(-EV_(3)+1.0+XI*EV_(3)*WMYEZ);
                H_(2,3) = EXPO*XI*WMY^(EV_(3)-1.0)*(1.0+EV_(3)*LWMY*(1.0-XI*WMYEZ));
                H_(3,2) = H_(2,3);
                H_(3,3) = EXPO*WMYEZ*XI*LWMY^2*(-1.0+XI*WMYEZ);
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

