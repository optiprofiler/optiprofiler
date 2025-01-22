function varargout = QPBAND(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : QPBAND
%    *********
% 
%    A banded QP
%    SIF input: Nick Gould, December 1999.
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'C-CQLR2-AN-V-V'
% 
%       Alternative values for the SIF file parameters:
% IE N                   10000          $-PARAMETER
% IE N                   50000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'QPBAND';

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
% IE N                   100000         $-PARAMETER
% IE N                   200000         $-PARAMETER
% IE N                   400000         $-PARAMETER
% IE N                   500000         $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        v_('ONE') = 1.0;
        v_('M') = fix(v_('N')/v_('2'));
        v_('N-1') = v_('N')-v_('1');
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
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('RI/RN') = v_('RI')/v_('RN');
            v_('-RI/RN') = -1.0*v_('RI/RN');
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('-RI/RN');
        end
        for I=v_('1'):v_('M')
            v_('M+I') = v_('M')+I;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('M+I')))]);
            valA(end+1) = 1.0;
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['C',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 2.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        irH  = [];
        icH  = [];
        valH = [];
        for I=v_('1'):v_('N-1')
            v_('I+1') = I+v_('1');
            irH(end+1)  =  ix_(['X',int2str(I)]);
            icH(end+1)  =  ix_(['X',int2str(I)]);
            valH(end+1) =  2.0;
            irH(end+1)  =  ix_(['X',int2str(I)]);
            icH(end+1)  =  ix_(['X',int2str(round(v_('I+1')))]);
            valH(end+1) =  -1.0;
            irH(end+1)  =  ix_(['X',int2str(round(v_('I+1')))]);
            icH(end+1)  =  ix_(['X',int2str(I)]);
            valH(end+1) =  -1.0;
        end
        irH(end+1)  =  ix_(['X',int2str(round(v_('N')))]);
        icH(end+1)  =  ix_(['X',int2str(round(v_('N')))]);
        valH(end+1) =  2.0;
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        pbm.H = sparse(irH,icH,valH,pb.n,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CQLR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
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

