function varargout = FLETCHBV(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : FLETCHBV
%    *********
% 
%    Another Boundary Value problem.
% 
%    Source:  The first problem given by
%    R. Fletcher,
%    "An optimal positive definite update for sparse Hessian matrices"
%    Numerical Analysis report NA/145, University of Dundee, 1992.
% 
%    N.B. This formulation is incorrect. See FLETCBV2.SIF for
%         the correct version.
% 
%    SDIF input: Nick Gould, Oct 1992.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    The number of variables is N.
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'FLETCHBV';

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
% IE N                   10000          $-PARAMETER
        if(nargs<2)
            v_('KAPPA') = 1.0;  %  SIF file default value
        else
            v_('KAPPA') = varargin{2};
        end
% RE KAPPA               0.0            $-PARAMETER
        v_('OBJSCALE') = 1.0e+0;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('1.0') = 1.0;
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        v_('RN+1') = v_('N+1');
        v_('H') = v_('1.0')/v_('RN+1');
        v_('H2') = v_('H')*v_('H');
        v_('1/H2') = v_('RN+1')*v_('RN+1');
        v_('KAPPA/H2') = v_('1/H2')*v_('KAPPA');
        v_('-KAPPA/H2') = -1.0*v_('KAPPA/H2');
        v_('-2/H2') = -2.0*v_('1/H2');
        v_('-1-2/H2') = -1.0*v_('-2/H2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('0')))],ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('OBJSCALE');
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('0')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('OBJSCALE');
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('OBJSCALE');
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['L',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('OBJSCALE');
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2/H2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2/H2');
            end
        end
        [ig,ig_] = s2mpjlib('ii',['L',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('OBJSCALE');
        [ig,ig_] = s2mpjlib('ii',['L',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-1-2/H2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-1-2/H2');
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('OBJSCALE');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('IH') = v_('RI')*v_('H');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('IH');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCOS',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['C',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eCOS';
                ielftype(ie) = iet_('eCOS');
            end
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gHALFL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('N')
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gHALFL2';
        end
        for I=v_('1'):v_('N')
            ig = ig_(['C',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('-KAPPA/H2');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                ??
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-V-0';
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

    case 'eCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = cos(EV_(1));
        if(nargout>1)
            g_(1,1) = -sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -cos(EV_(1));
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gHALFL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 5.0e-1*GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 1.0e+0;
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

