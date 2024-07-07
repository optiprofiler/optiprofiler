function varargout = OSCIGRAD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OSCIGRAD
%    *********
% 
%    The roots of the gradient of Yurii Nesterov's "oscillating path" problem
% 
%    SIF input: Nick Gould, June 2011.
% 
%    classification = 'SUR2-AN-V-0'
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER
% IE N                   5              $-PARAMETER
% IE N                   10             $-PARAMETER
% IE N                   15             $-PARAMETER
% IE N                   25             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   10000          $-PARAMETER
% IE N                   100000         $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OSCIGRAD';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% RE RHO                 1.0            $-PARAMETER    Nesterov's original value
        if(nargs<2)
            v_('RHO') = 500.0;  %  SIF file default value
        else
            v_('RHO') = varargin{2};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('2RHO') = 2.0*v_('RHO');
        v_('-4RHO') = -4.0*v_('RHO');
        v_('N-1') = -1+v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','G1',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.5;
        end
        for I=v_('2'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('G1')) = 0.5;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_(['X',int2str(round(v_('1')))]),1) = -2.0;
        for I=v_('2'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eA',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eB',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'B1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eB';
        ielftype(ie) = iet_('eB');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('-4RHO');
        for I=v_('2'):v_('N-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA';
            ielftype(ie) = iet_('eA');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('2RHO');
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eB';
            ielftype(ie) = iet_('eB');
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-4RHO');
        end
        ename = ['A',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA';
        ielftype(ie) = iet_('eA');
        ename = ['A',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['A',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['A',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('2RHO');
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_('G1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('B1');
        pbm.grelw{ig}(posel) = 1.0;
        for I=v_('2'):v_('N-1')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
        end
        ig = ig_(['G',int2str(round(v_('N')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-0';
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

    case 'eA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*(EV_(2)-2.0*EV_(1)^2+1.0);
        if(nargout>1)
            g_(2,1) = pbm.elpar{iel_}(1);
            g_(1,1) = -4.0*pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -4.0*pbm.elpar{iel_}(1);
                varargout{3} = H_;
            end
        end

    case 'eB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*(EV_(2)-2.0*EV_(1)^2+1.0)*EV_(1);
        if(nargout>1)
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            g_(1,1) = pbm.elpar{iel_}(1)*(EV_(2)-6.0*EV_(1)^2+1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = pbm.elpar{iel_}(1);
                H_(1,2) = H_(2,1);
                H_(1,1) = -12.0*pbm.elpar{iel_}(1)*EV_(1);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
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

