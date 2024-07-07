function varargout = SCW1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SCW1
%    *********
% 
%    Source: a discretization of an infinite-demsional problem proposed 
%    by Simon Chandler-Wilde (U. Reading):
% 
%    Given a function u in C[0,2 pi] with ||u||_infty <= 1, find the 
%    supremum of c^2(u) + s^2(u), where
%      c(u) = int_0^2 pi cos(t)u(t) dt and
%      s(u) = int_0^2 pi sin(t)u(t) dt      
% 
%    The discretized version ignores the required continuity, and 
%    posits a piecewise constant solution that oscilates between
%    plus and minus one. The anticipated solution is -16.
% 
%    SIF input: Nick Gould, July 2020
% 
%    classification = 'SLR2-MN-V-V'
% 
%    Number of internal knots
% 
%       Alternative values for the SIF file parameters:
% IE K                   1              $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SCW1';

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
            v_('K') = 7;  %  SIF file default value
        else
            v_('K') = varargin{1};
        end
% IE K                   3              $-PARAMETER
% IE K                   7              $-PARAMETER     original value
% IE K                   10             $-PARAMETER
% IE K                   100            $-PARAMETER
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('K+1') = 1+v_('K');
        v_('RK+1') = v_('K+1');
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('2PI') = 2.0*v_('PI');
        v_('2PI/K+1') = v_('2PI')/v_('RK+1');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('K+1')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','S',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','C',ig_);
        gtype{ig} = '<>';
        for I=v_('0'):v_('K')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['CON',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CON',int2str(I)];
            iv = ix_(['T',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['T',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('K')
            pb.xlower(ix_(['T',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['T',int2str(I)])) = v_('2PI');
        end
        pb.xlower(ix_(['T',int2str(round(v_('K+1')))]),1) = v_('2PI');
        pb.xupper(ix_(['T',int2str(round(v_('K+1')))]),1) = v_('2PI');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,['T',int2str(round(v_('1')))]))
            pb.x0(ix_(['T',int2str(round(v_('1')))]),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_(['T',int2str(round(v_('1')))])),1) = 1.0;
        end
        for I=v_('2'):v_('K')
            v_('RI') = I;
            v_('START') = v_('RI')*v_('2PI/K+1');
            if(isKey(ix_,['T',int2str(I)]))
                pb.x0(ix_(['T',int2str(I)]),1) = 0.0;
            else
                pb.y0(find(pbm.congrps==ig_(['T',int2str(I)])),1) = 0.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSINT',iet_);
        elftv{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'eCOST',iet_);
        elftv{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('K+1')
            ename = ['S',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSINT';
            ielftype(ie) = iet_('eSINT');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOST';
            ielftype(ie) = iet_('eCOST');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gMAXSQ',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('C');
        pbm.grftype{ig} = 'gMAXSQ';
        for I=v_('0'):v_('2'):v_('K')
            v_('I+1') = 1+I;
            ig = ig_('C');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        for I=v_('1'):v_('2'):v_('K')
            v_('I+1') = 1+I;
            ig = ig_('C');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
        end
        ig = ig_('S');
        pbm.grftype{ig} = 'gMAXSQ';
        for I=v_('0'):v_('2'):v_('K')
            v_('I+1') = 1+I;
            ig = ig_('S');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
        end
        for I=v_('1'):v_('2'):v_('K')
            v_('I+1') = 1+I;
            ig = ig_('S');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SCW                 0.0
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'SLR2-MN-V-V';
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

    case 'eSINT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S = sin(EV_(1));
        varargout{1} = S;
        if(nargout>1)
            g_(1,1) = cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -S;
                varargout{3} = H_;
            end
        end

    case 'eCOST'

        EV_  = varargin{1};
        iel_ = varargin{2};
        C = cos(EV_(1));
        varargout{1} = C;
        if(nargout>1)
            g_(1,1) = -sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -C;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gMAXSQ'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = -GVAR_*GVAR_;
        if(nargout>1)
            g_ = -GVAR_-GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -2.0e+0;
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

