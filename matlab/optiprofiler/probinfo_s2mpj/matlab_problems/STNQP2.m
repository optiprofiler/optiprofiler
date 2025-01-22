function varargout = STNQP2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : STNQP2
%    *********
% 
%    Another non-convex quadratic program with some structure.
% 
%    The objective function is of the form
%       sum (i=0,n) x_i^2 - 0.5 sum (l=1,n/p) sum(i=1,p) sum(k;i) x_{k+l}^2,
%    where n = 2^p and (k;i) means k takes the values of the first i powers of 2
%    eg, (k:3) = {k = {1,2,4}} and (k:7) = {k = {1,2,4,8,16,32}}.
%    There are equality constraints of the form
%    
%       sum(j=1,i) x_{(l-1)p+i} = i, where l=1,n/p,2 and i=1,p.
%    Finally, there are simple bounds
%          2 <= x_i, y_i <= 2    (i=0,n).
% 
%    SIF input: Nick Gould, May 1996
% 
%    classification = 'C-CQLR2-AN-V-V'
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
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'STNQP2';

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
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('0'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
        end
        for L=v_('1'):v_('N/P')
            for I=v_('1'):v_('P')
                v_('K') = v_('1');
                for J=v_('1'):I
                    v_('K+L') = v_('K')+L;
                    [ig,ig_] = s2mpjlib('ii',['N',int2str(I),',',int2str(L)],ig_);
                    gtype{ig} = '<>';
                    irA(end+1)  = ig;
                    icA(end+1)  = ix_(['X',int2str(round(v_('K+L')))]);
                    valA(end+1) = 1.0;
                    v_('K') = v_('K')*v_('2');
                end
            end
        end
        for L=v_('1'):v_('2'):v_('N/P')
            v_('LL') = L*v_('P');
            v_('LL') = v_('LL')-v_('P');
            for I=v_('1'):v_('P')
                for J=v_('1'):I
                    v_('LL+J') = v_('LL')+J;
                    [ig,ig_] = s2mpjlib('ii',['E',int2str(I),',',int2str(L)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['E',int2str(I),',',int2str(L)];
                    irA(end+1)  = ig;
                    icA(end+1)  = ix_(['X',int2str(round(v_('LL+J')))]);
                    valA(end+1) = 1.0;
                end
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
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gPSQR',igt_);
        [it,igt_] = s2mpjlib('ii','gPSQR',igt_);
        grftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
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
% LO SOLUTION            -2.476395E+5   $ (P=12)
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CQLR2-AN-V-V';
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

