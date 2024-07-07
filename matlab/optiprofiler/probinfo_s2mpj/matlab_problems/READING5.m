function varargout = READING5(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : READING5
%    *********
% 
%    A nonlinear optimal control problem from Nancy Nichols
%    with a given initial condition.
%    This problem arises in tide modelling.
% 
%    Source: a variant upon a problem in
%    S. Lyle and N.K. Nichols,
%    "Numerical Methods for Optimal Control Problems with State Constraints",
%    Numerical Analysis Report 8/91, Dept of Mathematics, 
%    University of Reading, UK.
% 
%    SIF input: Ph. Toint, Aug 1992
% 
%    classification = 'OOR2-MN-V-V'
% 
%    Number of discretized points in [0,1]
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER  n =    3, m =    2
% IE N                   50             $-PARAMETER  n =   51, m =   50
% IE N                   100            $-PARAMETER  n =  101, m =  100
% IE N                   500            $-PARAMETER  n =  501, m =  500
% IE N                   1000           $-PARAMETER  n = 1001, m = 1000
% IE N                   5000           $-PARAMETER  n = 5001, m = 5000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'READING5';

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
        if(nargs<1)
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('2/H') = 2.0*v_('RN');
        v_('H/2') = 0.5*v_('H');
        v_('1/H') = 1.0*v_('RN');
        v_('-1/H') = -1.0*v_('RN');
        v_('PI') = 3.1415926535;
        v_('2PI') = 2.0*v_('PI');
        v_('A') = 0.07716;
        v_('1/A') = 1.0/v_('A');
        v_('1/2A') = 0.5*v_('1/A');
        v_('2A') = 2.0*v_('A');
        v_('H/2A') = v_('H')*v_('1/2A');
        v_('2A/H') = 1.0/v_('H/2A');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii',['X',int2str(round(v_('0')))],ix_);
        pb.xnames{iv} = ['X',int2str(round(v_('0')))];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','J',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('1/A');
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['U',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['U',int2str(I)];
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 0.25;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 0.25;
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = -0.5;
            pb.xupper(ix_(['X',int2str(I)])) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eUC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'XP';
        elftp{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'eENERGY',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'XP';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            v_('I-1') = -1+I;
            ename = ['I',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eENERGY';
            ielftype(ie) = iet_('eENERGY');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XP',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            ename = ['UC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eUC';
            ielftype(ie) = iet_('eUC');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XP',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('J');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['I',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        for I=v_('2'):v_('N-1')
            ig = ig_('J');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['I',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -2.0;
        end
        ig = ig_('J');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['I',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        for I=v_('1'):v_('N')
            ig = ig_(['U',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['UC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('2A/H');
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-V-V';
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

    case 'eENERGY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        F = cos(2.0*3.141592653589*pbm.elpar{iel_}(1));
        varargout{1} = (F-EV_(1))*(EV_(1)-EV_(2));
        if(nargout>1)
            g_(1,1) = -2.0*EV_(1)+EV_(2)+F;
            g_(2,1) = -(F-EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -2.0;
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eUC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        F = cos(2.0*3.141592653589*pbm.elpar{iel_}(1));
        C = (EV_(1)-EV_(2))/(F-EV_(1));
        D = (1.0+C)/(F-EV_(1));
        varargout{1} = C;
        if(nargout>1)
            g_(1,1) = D;
            g_(2,1) = -1.0/(F-EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*D/(F-EV_(1));
                H_(1,2) = -1.0/(F-EV_(1))^2;
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

