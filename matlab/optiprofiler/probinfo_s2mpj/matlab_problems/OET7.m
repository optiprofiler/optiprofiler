function varargout = OET7(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OET7
%    *********
% 
%    A nonlinear programming formulation of a discretization of
%    a nonlinear Chebychev problem.
% 
%    The problem is
% 
%        min  max || phi(x,w) ||, for all w in the interval I.
%         x    w
% 
%    I is discretized, and the problem solved over the
%    discrete points.
% 
%    Nonlinear programming formulation
%        min   u     s.t.  u - phi >= 0, u + phi >= 0
%        x,u
% 
%    Specific problem: I = [-0.5,0.5]
%    phi(x,w) = 1/1+w - x1 exp(w x4) - x2 exp(w x5) - x3 exp(w x6)
% 
%    Source: K. Oettershagen
%    "Ein superlinear knonvergenter algorithmus zur losung 
%     semi-infiniter optimierungsproblem",
%     Ph.D thesis, Bonn University, 1982
% 
%    SIF input: Nick Gould, February, 1994.
% 
%    classification = 'LOR2-AN-7-V'
% 
%    Discretization
% 
% IE M                   2
% IE M                   100
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OET7';

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
        v_('M') = 500;
        v_('LOWER') = -0.5;
        v_('UPPER') = 0.5;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('DIFF') = v_('UPPER')-v_('LOWER');
        v_('RM') = v_('M');
        v_('H') = v_('DIFF')/v_('RM');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('0'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['LO',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['LO',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['UP',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['UP',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
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
        for I=v_('0'):v_('M')
            v_('RI') = I;
            v_('W') = v_('RI')*v_('H');
            v_('W') = v_('W')+v_('LOWER');
            v_('1+W') = 1.0+v_('W');
            v_('1/1+W') = 1.0/v_('1+W');
            v_('-1/1+W') = -1.0*v_('1/1+W');
            pbm.gconst(ig_(['LO',int2str(I)])) = v_('-1/1+W');
            pbm.gconst(ig_(['UP',int2str(I)])) = v_('1/1+W');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXEYW',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('0'):v_('M')
            v_('RI') = I;
            v_('W') = v_('RI')*v_('H');
            v_('W') = v_('W')+v_('LOWER');
            ename = ['E',int2str(round(v_('1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXEYW';
            ielftype(ie) = iet_('eXEYW');
            ename = ['E',int2str(round(v_('1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('W');
            ename = ['E',int2str(round(v_('2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXEYW';
            ielftype(ie) = iet_('eXEYW');
            ename = ['E',int2str(round(v_('2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('W');
            ename = ['E',int2str(round(v_('3'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXEYW';
            ielftype(ie) = iet_('eXEYW');
            ename = ['E',int2str(round(v_('3'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('3'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('3'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('W');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('M')
            ig = ig_(['LO',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['UP',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-AN-7-V';
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

    case 'eXEYW'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EYW = exp(EV_(2)*pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*EYW;
        if(nargout>1)
            g_(1,1) = EYW;
            g_(2,1) = EV_(1)*pbm.elpar{iel_}(1)*EYW;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1)*EYW;
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)*pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*EYW;
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

