function varargout = NCVXBQP2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NCVXBQP2
%    *********
% 
%    A non-convex bound constrained quadratic program.
% 
%    SIF input: Nick Gould, July 1995
% 
%    classification = 'QBR2-AN-V-0'
% 
%    The number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER    original value
% IE N                   10000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NCVXBQP2';

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
% IE N                   100000         $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('M') = fix(v_('N')/v_('2'));
        v_('NPLUS') = fix(v_('N')/v_('2'));
        v_('NPLUS+1') = 1+v_('NPLUS');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            v_('J') = 2*I;
            v_('J') = -1+v_('J');
            v_('K') = fix(v_('J')/v_('N'));
            v_('K') = v_('K')*v_('N');
            v_('J') = v_('J')-v_('K');
            v_('J') = 1+v_('J');
            iv = ix_(['X',int2str(round(v_('J')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            v_('J') = 3*I;
            v_('J') = -1+v_('J');
            v_('K') = fix(v_('J')/v_('N'));
            v_('K') = v_('K')*v_('N');
            v_('J') = v_('J')-v_('K');
            v_('J') = 1+v_('J');
            iv = ix_(['X',int2str(round(v_('J')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.1;
            pb.xupper(ix_(['X',int2str(I)])) = 10.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gSQR',igt_);
        [it,igt_] = s2mpjlib('ii','gSQR',igt_);
        grftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NPLUS')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gSQR';
            v_('RI') = I;
            [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = v_('RI');
        end
        for I=v_('NPLUS+1'):v_('N')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gSQR';
            v_('RI') = I;
            v_('RI') = -1.0*v_('RI');
            [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = v_('RI');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -1.33305D+06   $ (n=100)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'QBR2-AN-V-0';
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

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQR'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 0.5*pbm.grpar{igr_}(1)*GVAR_*GVAR_;
        if(nargout>1)
            g_ = pbm.grpar{igr_}(1)*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = pbm.grpar{igr_}(1);
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

