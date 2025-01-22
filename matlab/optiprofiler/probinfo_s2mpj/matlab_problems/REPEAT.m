function varargout = REPEAT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : REPEAT
%    *********
% 
%    This problem is to find the nearest feasible point to 2n+1 inconsistent
%    linear equations subject to bounds
% 
%    Source: blue-cheese delerium
% 
%    SIF input: Nick Gould, December 2020.
% 
%    classification = 'C-CNLR2-AN-V-V'
% 
%    N is the number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'REPEAT';

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
            v_('N') = 100;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   10000          $-PARAMETER
% IE N                   100000         $-PARAMETER     original value
% IE N                   1000000        $-PARAMETER
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('100') = 100;
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N/100') = fix(v_('N')/v_('100'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
        valA(end+1) = 1.0;
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['R',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['R',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['R',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
        valA(end+1) = 1.0;
        for I=v_('1'):v_('N/100'):v_('N')
            v_('RI') = I;
            [ig,ig_] = s2mpjlib('ii','E',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'E';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('RI');
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
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['C',int2str(I)])) = 2.0;
        end
        pbm.gconst(ig_(['C',int2str(round(v_('N')))])) = 1.0;
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['R',int2str(I)])) = 4.0;
        end
        pbm.gconst(ig_(['R',int2str(round(v_('N')))])) = 3.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = -1.0*ones(pb.n,1);
        pb.xupper(ix_(['X',int2str(round(v_('1')))])) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('2')))]),1) = 3.0;
        pb.xupper(ix_(['X',int2str(round(v_('N/2')))])) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('N-1')))]),1) = 3.0;
        pb.xupper(ix_(['X',int2str(round(v_('N')))])) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CNLR2-AN-V-V';
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

