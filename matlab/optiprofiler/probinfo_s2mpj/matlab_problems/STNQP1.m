function varargout = STNQP1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : STNQP1
%    *********
% 
%    A non-convex quadratic program with some structure.
% 
%    The objective function is of the form
%       sum (i=0,n) x_i^2 - 0.5 sum (l=1,n/p) sum(i=1,p) sum(k;i) x_{k+l}^2,
%    where n = 2^p and (k;i) means k takes the values of the first i powers of 2
%    eg, (k:3) = {k = {1,2,4}} and (k:7) = {k = {1,2,4,8,16,32}}.
%    There are equality constraints of the form
%    
%       sum(k;i) x_{k+l-1} = i, where l=1,n/p,2 and i=1,p.
%    Finally, there are simple bounds
%          2 <= x_i, y_i <= 2    (i=0,n).
% 
%    SIF input: Nick Gould, May 1996
% 
%    classification = 'QLR2-AN-V-V'
% 
%    There will be 2**p + 1 variables
% 
%       Alternative values for the SIF file parameters:
% IE P                   2              $-PARAMETER n = 5
% IE P                   4              $-PARAMETER n = 17
% IE P                   6              $-PARAMETER n = 65
% IE P                   8              $-PARAMETER n = 257
% IE P                   10             $-PARAMETER n = 1025
% IE P                   12             $-PARAMETER n = 4097     original value
% IE P                   13             $-PARAMETER n = 8193
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'STNQP1';

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
            v_('P') = 4;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
% IE P                   14             $-PARAMETER n = 16395
% IE P                   15             $-PARAMETER n = 32769
% IE P                   16             $-PARAMETER n = 65537
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N') = 1;
        for I=v_('1'):v_('P')
            v_('N') = v_('N')*v_('2');
        end
        v_('N/P') = fix(v_('N')/v_('P'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('0'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for L=v_('1'):v_('N/P')
            for I=v_('1'):v_('P')
                v_('K') = v_('1');
                for J=v_('1'):I
                    v_('K+L') = v_('K')+L;
                    [ig,ig_] = s2mpjlib('ii',['N',int2str(I),',',int2str(L)],ig_);
                    gtype{ig} = '<>';
                    iv = ix_(['X',int2str(round(v_('K+L')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1.0;
                    end
                    v_('K') = v_('K')*v_('2');
                end
            end
        end
        for L=v_('1'):v_('2'):v_('N/P')
            for I=v_('1'):v_('P')
                v_('K') = v_('1');
                for J=v_('1'):I
                    v_('K-1') = v_('K')-v_('1');
                    v_('K+L-1') = v_('K-1')+L;
                    [ig,ig_] = s2mpjlib('ii',['E',int2str(I),',',int2str(L)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['E',int2str(I),',',int2str(L)];
                    iv = ix_(['X',int2str(round(v_('K+L-1')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1.0;
                    end
                    v_('K') = v_('K')*v_('2');
                end
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
        for L=v_('1'):v_('2'):v_('N/P')
            for I=v_('1'):v_('P')
                v_('RI') = I;
                pbm.gconst(ig_(['E',int2str(I),',',int2str(L)])) = v_('RI');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('0'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = -2.0;
            pb.xupper(ix_(['X',int2str(I)])) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gPSQR',igt_);
        [it,igt_] = s2mpjlib('ii','gPSQR',igt_);
        grftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('N')
            ig = ig_(['O',int2str(I)]);
            pbm.grftype{ig} = 'gPSQR';
            [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = 1.0;
        end
        for L=v_('1'):v_('N/P')
            for I=v_('1'):v_('P')
                ig = ig_(['N',int2str(I),',',int2str(L)]);
                pbm.grftype{ig} = 'gPSQR';
                [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
                pbm.grpar{ig}(posgp) = -0.5;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLUTION            -1.361565E+5   $ (P=12)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'QLR2-AN-V-V';
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

    case 'gPSQR'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = pbm.grpar{igr_}(1)*GVAR_*GVAR_;
        if(nargout>1)
            g_ = 2.0*pbm.grpar{igr_}(1)*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0*pbm.grpar{igr_}(1);
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

