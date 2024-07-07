function varargout = POWELLSG(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : POWELLSG
%    *********
% 
%    The extended Powell singular problem.
%    This problem is a sum of n/4 sets of four terms, each of which is
%    assigned its own group.
% 
%    Source:  Problem 13 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Toint#19, Buckley#34 (p.85)
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    N is the number of free variables, and should be a multiple of 4
% 
%       Alternative values for the SIF file parameters:
% IE N                   4              $-PARAMETER     original value
% IE N                   8              $-PARAMETER
% IE N                   16             $-PARAMETER
% IE N                   20             $-PARAMETER
% IE N                   36             $-PARAMETER
% IE N                   40             $-PARAMETER
% IE N                   60             $-PARAMETER
% IE N                   80             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'POWELLSG';

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
            v_('N') = 12;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   10000          $-PARAMETER
        v_('1') = 1;
        v_('4') = 4;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('4'):v_('N')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 10.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 10.0;
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('I+1')))],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('I+1')))],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.2;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('I+2')))],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0;
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('I+3')))],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('I+3')))],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.1;
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
        for I=v_('1'):v_('4'):v_('N')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            pb.x0(ix_(['X',int2str(I)]),1) = 3.0;
            pb.x0(ix_(['X',int2str(round(v_('I+1')))]),1) = -1.0;
            pb.x0(ix_(['X',int2str(round(v_('I+2')))]),1) = 0.0;
            pb.x0(ix_(['X',int2str(round(v_('I+3')))]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        [it,igt_] = s2mpjlib('ii','gL4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('4'):v_('N')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['G',int2str(round(v_('I+1')))]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['G',int2str(round(v_('I+2')))]);
            pbm.grftype{ig} = 'gL4';
            ig = ig_(['G',int2str(round(v_('I+3')))]);
            pbm.grftype{ig} = 'gL4';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
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
% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

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

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
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

