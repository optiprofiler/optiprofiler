function varargout = LOTSCHD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A simple quadratic program inspired by the economic lot scheduling
%    problem.
% 
%    Source:
%    an exercize for L. Watson course on LANCELOT in the Spring 1993.
% 
%    SIF input: T. Kuan, Virginia Tech., Spring 1993.
% 
%    classification = 'C-CQLR2-AN-12-7'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LOTSCHD';

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
        v_('X1') = 1.502;
        v_('X2') = 1.126;
        v_('X3') = 0.815;
        v_('X4') = 1.268;
        v_('X5') = 1.502;
        v_('X6') = 0.740;
        v_('A1') = 1.8;
        v_('A2') = 3.2;
        v_('A3') = 6.1;
        v_('A4') = 3.2;
        v_('A5') = 1.8;
        v_('A6') = 7.4;
        v_('C1') = 11.0;
        v_('C2') = 3.0;
        v_('C3') = 20.0;
        v_('C4') = 17.0;
        v_('C5') = 9.0;
        v_('C6') = 20.0;
        v_('C7') = 126.1;
        v_('1') = 1;
        v_('2') = 2;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('6')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('6')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = v_(['X',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii','CONS7',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CONS7';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['CONS',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CONS',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = v_(['A',int2str(I)]);
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for I=v_('2'):v_('4')
            [ig,ig_] = s2mpjlib('ii',['CONS',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CONS',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CONS2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CONS2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('T3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U3');
        valA(end+1) = -1.0;
        for I=v_('1'):v_('2')
            [ig,ig_] = s2mpjlib('ii','CONS3',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CONS3';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for I=v_('4'):v_('6')
            [ig,ig_] = s2mpjlib('ii','CONS3',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CONS3';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CONS4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CONS4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('T1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U1');
        valA(end+1) = -1.0;
        for I=v_('5'):v_('6')
            [ig,ig_] = s2mpjlib('ii','CONS4',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CONS4';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CONS5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CONS5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('T1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U1');
        valA(end+1) = -1.0;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii','CONS6',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CONS6';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -1.0;
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
        for I=v_('1'):v_('7')
            pbm.gconst(ig_(['CONS',int2str(I)])) = v_(['C',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSQUARE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('6')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gSQUARE';
        end
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CQLR2-AN-12-7';
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


    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQUARE'

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

