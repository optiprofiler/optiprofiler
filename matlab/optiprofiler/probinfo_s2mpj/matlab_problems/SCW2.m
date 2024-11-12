function varargout = SCW2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SCW2
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
%    posits a piecewise constant solution that varies anywhere between
%    plus and minus one. The anticipated solution is -16.
% 
%    SIF input: Nick Gould, July 2020
% 
%    classification = 'C-CSLR2-MN-V-V'
% 
%    Number of internal knots
% 
%       Alternative values for the SIF file parameters:
% IE K                   1              $-PARAMETER
% IE K                   10             $-PARAMETER
% IE K                   100            $-PARAMETER
% IE K                   1000           $-PARAMETER     original value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SCW2';

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
        if(nargs<1)
            v_('K') = 10;  %  SIF file default value
        else
            v_('K') = varargin{1};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('ONE') = 1.0;
        v_('K+1') = 1+v_('K');
        v_('RK') = v_('K');
        v_('RK+1') = v_('K+1');
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('2PI') = 2.0*v_('PI');
        v_('2PI/K+1') = v_('2PI')/v_('RK+1');
        v_('1/K') = v_('ONE')/v_('RK');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('K+1')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
        end
        for I=v_('0'):v_('K')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
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
        pb.xlower(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('K')
            pb.xlower(ix_(['T',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['T',int2str(I)])) = v_('2PI');
        end
        pb.xlower(ix_(['T',int2str(round(v_('K+1')))]),1) = v_('2PI');
        pb.xupper(ix_(['T',int2str(round(v_('K+1')))]),1) = v_('2PI');
        for I=v_('0'):v_('K')
            pb.xlower(ix_(['U',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
        end
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
        for I=v_('0'):v_('K')
            v_('RI') = I;
            v_('START') = v_('RI')*v_('1/K');
            if(isKey(ix_,['U',int2str(I)]))
                pb.x0(ix_(['U',int2str(I)]),1) = v_('START');
            else
                pb.y0(find(pbm.congrps==ig_(['U',int2str(I)])),1) = v_('START');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eUSINT',iet_);
        elftv{it}{1} = 'T';
        elftv{it}{2} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'eUCOST',iet_);
        elftv{it}{1} = 'T';
        elftv{it}{2} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('K')
            v_('I+1') = 1+I;
            ename = ['US',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eUSINT';
            ielftype(ie) = iet_('eUSINT');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['USP',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eUSINT';
            ielftype(ie) = iet_('eUSINT');
            vname = ['T',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['UC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eUCOST';
            ielftype(ie) = iet_('eUCOST');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['UCP',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eUCOST';
            ielftype(ie) = iet_('eUCOST');
            vname = ['T',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gMAXSQ',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('C');
        pbm.grftype{ig} = 'gMAXSQ';
        for I=v_('0'):v_('K')
            ig = ig_('C');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['USP',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['US',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        ig = ig_('S');
        pbm.grftype{ig} = 'gMAXSQ';
        for I=v_('0'):v_('K')
            ig = ig_('S');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['UCP',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['UC',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
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
        pb.pbclass = 'C-CSLR2-MN-V-V';
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

    case 'eUSINT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S = sin(EV_(1));
        C = cos(EV_(1));
        varargout{1} = EV_(2)*S;
        if(nargout>1)
            g_(1,1) = EV_(2)*C;
            g_(2,1) = S;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -EV_(2)*S;
                H_(1,2) = C;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eUCOST'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S = sin(EV_(1));
        C = cos(EV_(1));
        varargout{1} = EV_(2)*C;
        if(nargout>1)
            g_(1,1) = -EV_(2)*S;
            g_(2,1) = C;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -EV_(2)*C;
                H_(1,2) = -S;
                H_(2,1) = H_(1,2);
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

