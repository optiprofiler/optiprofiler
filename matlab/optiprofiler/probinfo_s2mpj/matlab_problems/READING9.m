function varargout = READING9(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : READING9
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
%    SIF input: Nick Gould and Ph. Toint, March 1995
% 
%    classification = 'OOR2-MN-V-V'
% 
%    Number of discretized points in [0,1] - 1
% 
%       Alternative values for the SIF file parameters:
% IE N+1                 3              $-PARAMETER n=6, m=2
% IE N+1                 51             $-PARAMETER n=102, m=50
% IE N+1                 101            $-PARAMETER n=202, m=100
% IE N+1                 201            $-PARAMETER n=402, m=200
% IE N+1                 501            $-PARAMETER n=1002, m=500  original value
% IE N+1                 1001           $-PARAMETER n=2002, m=1000
% IE N+1                 5001           $-PARAMETER n=10002, m= 5000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'READING9';

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
            v_('N+1') = 3;  %  SIF file default value
        else
            v_('N+1') = varargin{1};
        end
        v_('N') = -1+v_('N+1');
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('-H/2') = -0.5*v_('H');
        v_('K1') = 0.07716;
        v_('K2') = 0.07716;
        v_('K1H') = v_('K1')*v_('H');
        v_('K1H+1') = 1.0+v_('K1H');
        v_('-K1H-1') = -1.0*v_('K1H+1');
        v_('K2H') = v_('K2')*v_('H');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['S',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['S',int2str(I)];
            iv = ix_(['P',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['P',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-K1H-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-K1H-1');
            end
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
        for I=v_('0'):v_('N-1')
            v_('RI') = I;
            v_('T') = v_('RI')*v_('H');
            v_('SINT') = sin(v_('T'));
            v_('HSINT') = v_('H')*v_('SINT');
            pbm.gconst(ig_(['S',int2str(I)])) = v_('HSINT');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('0'):v_('N')
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
            pb.xlower(ix_(['P',int2str(I)])) = -Inf;
            pb.xupper(ix_(['P',int2str(I)]),1) = +Inf;
        end
        pb.xlower(ix_(['P',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['P',int2str(round(v_('0')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.2*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('N')
            ename = ['OE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD2';
            ielftype(ie) = iet_('ePROD2');
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.2);
            posev = find(strcmp('P',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.2);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('0'):v_('N-1')
            ename = ['CE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.2);
            posev = find(strcmp('P',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.2);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['OE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['OE',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            ig = ig_(['S',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('K2H');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION            -4.41677D-02   $ (n=500)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

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

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)*EV_(2))^2;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*EV_(2)^2;
            g_(2,1) = 2.0*EV_(2)*EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*EV_(2)^2;
                H_(1,2) = 4.0*EV_(1)*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1)^2;
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

