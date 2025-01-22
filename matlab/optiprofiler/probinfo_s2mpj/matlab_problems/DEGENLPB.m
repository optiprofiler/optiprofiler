function varargout = DEGENLPB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DEGENLPB
%    *********
% 
%    A linear program with some degeneracy.
% 
%    Source:
%    T.C.T. Kotiah and D.I. Steinberg,
%    "Occurences of cycling and other phenomena arising in a class of
%    linear programming models",
%    Communications of the ACM, vol. 20, pp. 107-112, 1977.
% 
%    SIF input: Ph. Toint, Aug 1990.
% 
%    classification = 'C-CLLR2-AN-20-15'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DEGENLPB';

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
        v_('N') = 20;
        v_('M') = 15;
        v_('1') = 1;
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
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -0.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -33.333;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -0.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -33.343;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -100.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -33.333;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = -133.33;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = -100.0;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 300.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.09;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.03;
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 0.326;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -101.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 200.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 0.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 0.02;
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 0.0066667;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -1.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 200.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 0.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 0.02;
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 6.6667e-4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -1.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 200.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 0.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 0.02;
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.978;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -201.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.489;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -101.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.001;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.489;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -101.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C9';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.001;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.01;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = -1.04;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C10';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -1.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.002;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = -1.02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 100.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 0.03;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = 0.01;
        [ig,ig_] = s2mpjlib('ii','C12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = -2.5742e-6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 0.00252;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = -0.61975;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = 1.03;
        [ig,ig_] = s2mpjlib('ii','C13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C13';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = -0.00257;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 0.25221;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = -6.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = 1.09;
        [ig,ig_] = s2mpjlib('ii','C14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C14';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 0.00629;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = -0.20555;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = -4.1106;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 101.04;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 505.1;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = -256.72;
        [ig,ig_] = s2mpjlib('ii','C15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C15';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 0.00841;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = -0.08406;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = -0.20667;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 20.658;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = 1.07;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = -10.5;
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
        pbm.gconst(ig_('C1')) = 0.70785;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 1.0*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               3.06435
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CLLR2-AN-20-15';
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

