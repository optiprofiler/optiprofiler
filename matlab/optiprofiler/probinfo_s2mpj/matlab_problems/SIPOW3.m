function varargout = SIPOW3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SIPOW3
%    *********
% 
%    This is a discretization of a one sided approximation problem of
%    approximating the function xi * xi * eta by a linear polynomial
%    on the boundary of the unit square [0,1]x[0,1].
% 
%    Source: problem 3 in
%    M. J. D. Powell,
%    "Log barrier methods for semi-infinite programming calculations"
%    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
% 
%    SIF input: A. R. Conn and Nick Gould, August 1993
% 
%    classification = 'C-CLLR2-AN-4-V'
% 
%    Problem variants: they are identified by the values of M (even)
% 
% IE M                   20 
% IE M                   100 
% IE M                   500 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SIPOW3';

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
        v_('M') = 2000;
        v_('1') = 1;
        v_('2') = 2;
        v_('RM') = v_('M');
        v_('RM/8') = 0.125*v_('RM');
        v_('RM/4') = 0.25*v_('RM');
        v_('3RM/8') = 0.375*v_('RM');
        v_('M/8') = fix(v_('RM/8'));
        v_('M/4') = fix(v_('RM/4'));
        v_('3M/8') = fix(v_('3RM/8'));
        v_('M/8+1') = 1+v_('M/8');
        v_('M/4+1') = 1+v_('M/4');
        v_('3M/8+1') = 1+v_('3M/8');
        v_('M/2') = fix(v_('M')/v_('2'));
        v_('M/2+1') = 1+v_('M/2');
        v_('RM') = v_('M');
        v_('STEP') = 8.0/v_('RM');
        for J=v_('1'):v_('M/2')
            v_('I') = -1+J;
            v_('RI') = v_('I');
            v_(['XI',int2str(J)]) = v_('RI')*v_('STEP');
        end
        for J=v_('1'):v_('M/8')
            v_('RJ') = J;
            v_(['ETA',int2str(J)]) = v_(['XI',int2str(J)]);
            v_(['XI',int2str(J)]) = 0.0;
        end
        for J=v_('M/8+1'):v_('M/4')
            v_('RJ') = J;
            v_(['XI',int2str(J)]) = -1.0+v_(['XI',int2str(J)]);
            v_(['ETA',int2str(J)]) = 1.0;
        end
        for J=v_('M/4+1'):v_('3M/8')
            v_('RJ') = J;
            v_(['ETA',int2str(J)]) = -2.0+v_(['XI',int2str(J)]);
            v_(['XI',int2str(J)]) = 1.0;
        end
        for J=v_('3M/8+1'):v_('M/2')
            v_('RJ') = J;
            v_(['XI',int2str(J)]) = -3.0+v_(['XI',int2str(J)]);
            v_(['ETA',int2str(J)]) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        for J=v_('1'):v_('M/2')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X4');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X2');
            valA(end+1) = v_(['XI',int2str(J)]);
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X3');
            valA(end+1) = v_(['ETA',int2str(J)]);
        end
        for J=v_('1'):v_('M/2')
            v_('J+') = v_('M/2')+J;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X1');
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X2');
            valA(end+1) = v_(['XI',int2str(J)]);
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('J+')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('J+')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('X3');
            valA(end+1) = v_(['ETA',int2str(J)]);
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
        for J=v_('1'):v_('M/2')
            v_('J+') = v_('M/2')+J;
            v_('XIXI') = v_(['XI',int2str(J)])*v_(['XI',int2str(J)]);
            v_('XIXIETA') = v_('XIXI')*v_(['ETA',int2str(J)]);
            pbm.gconst(ig_(['C',int2str(J)])) = v_('XIXIETA');
            pbm.gconst(ig_(['C',int2str(round(v_('J+')))])) = v_('XIXIETA');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = -0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = -0.1;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 0.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 1.2;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 1.2;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            3.0315716D-1 ! m = 20
% LO SOLUTION            5.0397238D-1 ! m = 100
% LO SOLUTION            5.3016386D-1 ! m = 500
% LO SOLUTION            5.3465470D-1 ! m = 2000
% LO SOLUTION            5.3564207D-1 ! m = 10000
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CLLR2-AN-4-V';
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

