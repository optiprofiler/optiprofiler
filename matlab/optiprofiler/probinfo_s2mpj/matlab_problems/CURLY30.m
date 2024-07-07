function varargout = CURLY30(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CURLY30
%    --------
% 
%    A banded function with semi-bandwidth 30 and
%    negative curvature near the starting point
% 
%    Source: Nick Gould
% 
%    SIF input: Nick Gould, September 1997.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   10000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CURLY30';

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
            v_('N') = 35;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('K') = 30;
        v_('1') = 1;
        v_('N-K') = v_('N')-v_('K');
        v_('N-K+1') = 1+v_('N-K');
        v_('RN') = v_('N');
        v_('RN+1') = 1+v_('RN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-K')
            v_('I+K') = I+v_('K');
            for J=I:v_('I+K')
                [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for I=v_('N-K+1'):v_('N')
            for J=I:v_('N')
                [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
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
            v_('T') = v_('RI')/v_('RN+1');
            v_('T') = 0.0001*v_('T');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('T');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gP4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gP4';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% ZL SOLTN               -1.003163D+5   $ (n=1000)
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

    case 'gP4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        APB = 2.0e+1;
        varargout{1} = GVAR_*(GVAR_*(GVAR_^2-APB)-1.0e-1);
        if(nargout>1)
            g_ = 2.0e+0*GVAR_*(2.0e+0*GVAR_^2-APB)-1.0e-1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 1.2e+1*GVAR_^2-2.0e+0*APB;
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

