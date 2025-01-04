function varargout = MINSURF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MINSURF
%    *********
%    Variable dimension full rank linear problem
%    A version of the minimum surface problem
%    on the unit square with simple boundary conditions.
% 
%    SIF input: Ph. Toint, Jan 1991.
% 
%    classification = 'C-COXR2-MY-64-0'
% 
%    Discretization parameter
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MINSURF';

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
        v_('P') = 7;
        v_('1') = 1;
        v_('P+1') = 1+v_('P');
        v_('RP') = v_('P');
        v_('RPSQ') = v_('RP')*v_('RP');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for i=v_('1'):v_('P+1')
            for j=v_('1'):v_('P+1')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(i),',',int2str(j)],ix_);
                pb.xnames{iv} = ['X',int2str(i),',',int2str(j)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for i=v_('1'):v_('P')
            for j=v_('1'):v_('P')
                [ig,ig_] = s2mpjlib('ii',['S',int2str(i),',',int2str(j)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('RPSQ');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for i=v_('1'):v_('P')
            for j=v_('1'):v_('P')
                pbm.gconst(ig_(['S',int2str(i),',',int2str(j)])) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        v_('2') = 2;
        for i=v_('2'):v_('P')
            for j=v_('2'):v_('P')
                pb.xlower(ix_(['X',int2str(i),',',int2str(j)])) = -Inf;
                pb.xupper(ix_(['X',int2str(i),',',int2str(j)]),1) = +Inf;
            end
        end
        for i=v_('1'):v_('P+1')
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(i)]),1) = 1.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(i)]),1) = 1.0;
            pb.xlower(ix_(['X',int2str(round(v_('P+1'))),',',int2str(i)]),1) = 1.0;
            pb.xupper(ix_(['X',int2str(round(v_('P+1'))),',',int2str(i)]),1) = 1.0;
            pb.xlower(ix_(['X',int2str(i),',',int2str(round(v_('1')))]),1) = 1.0;
            pb.xupper(ix_(['X',int2str(i),',',int2str(round(v_('1')))]),1) = 1.0;
            pb.xlower(ix_(['X',int2str(i),',',int2str(round(v_('P+1')))]),1) = 1.0;
            pb.xupper(ix_(['X',int2str(i),',',int2str(round(v_('P+1')))]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for i=v_('1'):v_('P')
            v_('i+1') = 1+i;
            for j=v_('1'):v_('P')
                v_('j+1') = 1+j;
                ename = ['A',int2str(i),',',int2str(j)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(i),',',int2str(j)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('i+1'))),',',int2str(round(v_('j+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('W',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(i),',',int2str(j)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(i),',',int2str(round(v_('j+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('i+1'))),',',int2str(j)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('W',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSQROOT',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        v_('WEIGHT') = 0.5*v_('RPSQ');
        for i=v_('1'):v_('P')
            for j=v_('1'):v_('P')
                ig = ig_(['S',int2str(i),',',int2str(j)]);
                pbm.grftype{ig} = 'gSQROOT';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(i),',',int2str(j)]);
                pbm.grelw{ig}(posel) = v_('WEIGHT');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(i),',',int2str(j)]);
                pbm.grelw{ig}(posel) = v_('WEIGHT');
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COXR2-MY-64-0';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eISQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQROOT'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        SQRAL = sqrt(GVAR_);
        varargout{1} = SQRAL;
        if(nargout>1)
            g_ = 0.5/SQRAL;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -0.25/(GVAR_*SQRAL);
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

