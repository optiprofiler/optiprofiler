function varargout = DALLASS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DALLASS
%    *********
% 
%    The small Dallas water distribution problem
%    The problem is also named "W30" in some references.
%    This is a nonlinear network problem with conditioning of
%    the order of 10**4.
% 
%    Source:
%    R. Dembo,
%    private communication, 1986.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'C-CONR2-MN-46-31'
% 
%    Number of arcs
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DALLASS';

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
        v_('N') = 46;
        v_('NODES') = 31;
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
        icA(end+1)  = ix_('X42');
        valA(end+1) = -6.38400e+02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X43');
        valA(end+1) = -6.33000e+02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X44');
        valA(end+1) = -5.54500e+02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X45');
        valA(end+1) = -5.05000e+02;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X46');
        valA(end+1) = -4.36900e+02;
        [ig,ig_] = s2mpjlib('ii','N1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X46');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X41');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X45');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X44');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N9';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N10';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N13';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X42');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N14';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X21');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N15';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X43');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X21');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N16';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X23');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X22');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N17',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N17';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X23');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X25');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X24');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N18';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X31');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X25');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X22');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X26');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N19',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N19';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X26');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X28');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X27');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N20';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X28');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','N21',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N21';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X31');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X30');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X29');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N22',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N22';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X30');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X27');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','N23',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N23';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X24');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X32');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N24',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N24';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X38');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X29');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X34');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X33');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N25',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N25';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X32');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X35');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N26',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N26';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X35');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X37');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X36');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N27',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N27';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X37');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X34');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','N28',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N28';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X36');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X40');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X39');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X38');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N29',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N29';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X39');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X33');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','N30',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N30';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X40');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X41');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','N31',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'N31';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X46');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X45');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X44');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X43');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X42');
        valA(end+1) = -1.0;
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
        pbm.gconst(ig_('N5')) = 2.80000;
        pbm.gconst(ig_('N7')) = 4.03000e-01;
        pbm.gconst(ig_('N8')) = 5.92000e-01;
        pbm.gconst(ig_('N9')) = 1.15600;
        pbm.gconst(ig_('N10')) = 2.00000e-01;
        pbm.gconst(ig_('N11')) = 4.95000e-01;
        pbm.gconst(ig_('N16')) = 3.13000e-01;
        pbm.gconst(ig_('N17')) = 8.44000e-01;
        pbm.gconst(ig_('N18')) = 3.31000e-01;
        pbm.gconst(ig_('N19')) = 5.30000e-02;
        pbm.gconst(ig_('N21')) = 2.72000e-01;
        pbm.gconst(ig_('N22')) = 8.83000e-01;
        pbm.gconst(ig_('N23')) = 5.71000e-01;
        pbm.gconst(ig_('N24')) = 7.55000e-01;
        pbm.gconst(ig_('N26')) = 5.27000e-01;
        pbm.gconst(ig_('N29')) = 1.00000e-03;
        pbm.gconst(ig_('N31')) = -1.01960e+01;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -2.00000e+02*ones(pb.n,1);
        pb.xupper = 2.00000e+02*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.00000;
        pb.xupper(ix_('X1')) = 2.11673e+01;
        pb.xlower(ix_('X2'),1) = 0.00000;
        pb.xupper(ix_('X2')) = 4.37635e+01;
        pb.xlower(ix_('X3'),1) = 0.00000;
        pb.xupper(ix_('X3')) = 3.28255e+01;
        pb.xlower(ix_('X19'),1) = 0.00000;
        pb.xupper(ix_('X19')) = 2.20120e+01;
        pb.xlower(ix_('X21'),1) = 0.00000;
        pb.xupper(ix_('X21')) = 1.36703e+01;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = -2.00000e+02*ones(pb.n,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 2.11673e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 2.11673e+01;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 4.37635e+01;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 4.37635e+01;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 3.28255e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 3.28255e+01;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 1.42109e-14;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 1.42109e-14;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 1.68826e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 1.68826e+02;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 2.81745e+01;
        else
            pb.y0(find(pbm.congrps==ig('X7')),1) = 2.81745e+01;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 8.75603e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X8')),1) = 8.75603e+01;
        end
        if(isKey(ix_,'X9'))
            pb.x0(ix_('X9'),1) = -5.93858e+01;
        else
            pb.y0(find(pbm.congrps==ig('X9')),1) = -5.93858e+01;
        end
        if(isKey(ix_,'X10'))
            pb.x0(ix_('X10'),1) = -5.97888e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X10')),1) = -5.97888e+01;
        end
        if(isKey(ix_,'X11'))
            pb.x0(ix_('X11'),1) = 1.83383e+02;
        else
            pb.y0(find(pbm.congrps==ig('X11')),1) = 1.83383e+02;
        end
        if(isKey(ix_,'X13'))
            pb.x0(ix_('X13'),1) = -1.68331e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X13')),1) = -1.68331e+02;
        end
        if(isKey(ix_,'X15'))
            pb.x0(ix_('X15'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X15')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X16'))
            pb.x0(ix_('X16'),1) = 2.00000e-01;
        else
            pb.y0(find(pbm.congrps==ig_('X16')),1) = 2.00000e-01;
        end
        if(isKey(ix_,'X17'))
            pb.x0(ix_('X17'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X17')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X18'))
            pb.x0(ix_('X18'),1) = -7.67574e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X18')),1) = -7.67574e+01;
        end
        if(isKey(ix_,'X19'))
            pb.x0(ix_('X19'),1) = 2.20120e+01;
        else
            pb.y0(find(pbm.congrps==ig('X19')),1) = 2.20120e+01;
        end
        if(isKey(ix_,'X20'))
            pb.x0(ix_('X20'),1) = 1.36703e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X20')),1) = 1.36703e+01;
        end
        if(isKey(ix_,'X21'))
            pb.x0(ix_('X21'),1) = 1.36703e+01;
        else
            pb.y0(find(pbm.congrps==ig('X21')),1) = 1.36703e+01;
        end
        if(isKey(ix_,'X22'))
            pb.x0(ix_('X22'),1) = -1.98461e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X22')),1) = -1.98461e+02;
        end
        if(isKey(ix_,'X23'))
            pb.x0(ix_('X23'),1) = 1.81531e+02;
        else
            pb.y0(find(pbm.congrps==ig('X23')),1) = 1.81531e+02;
        end
        if(isKey(ix_,'X24'))
            pb.x0(ix_('X24'),1) = -1.93133e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X24')),1) = -1.93133e+01;
        end
        if(isKey(ix_,'X25'))
            pb.x0(ix_('X25'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X25')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X26'))
            pb.x0(ix_('X26'),1) = -1.98792e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X26')),1) = -1.98792e+02;
        end
        if(isKey(ix_,'X27'))
            pb.x0(ix_('X27'),1) = 1.15500;
        else
            pb.y0(find(pbm.congrps==ig('X27')),1) = 1.15500;
        end
        if(isKey(ix_,'X28'))
            pb.x0(ix_('X28'),1) = 0.00000;
        else
            pb.y0(find(pbm.congrps==ig_('X28')),1) = 0.00000;
        end
        if(isKey(ix_,'X29'))
            pb.x0(ix_('X29'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X29')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X30'))
            pb.x0(ix_('X30'),1) = 2.72000e-01;
        else
            pb.y0(find(pbm.congrps==ig_('X30')),1) = 2.72000e-01;
        end
        if(isKey(ix_,'X32'))
            pb.x0(ix_('X32'),1) = -1.98843e+01;
        else
            pb.y0(find(pbm.congrps==ig('X32')),1) = -1.98843e+01;
        end
        if(isKey(ix_,'X33'))
            pb.x0(ix_('X33'),1) = 1.78834e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X33')),1) = 1.78834e+02;
        end
        if(isKey(ix_,'X34'))
            pb.x0(ix_('X34'),1) = -1.79589e+02;
        else
            pb.y0(find(pbm.congrps==ig('X34')),1) = -1.79589e+02;
        end
        if(isKey(ix_,'X35'))
            pb.x0(ix_('X35'),1) = -1.98843e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X35')),1) = -1.98843e+01;
        end
        if(isKey(ix_,'X37'))
            pb.x0(ix_('X37'),1) = 1.79589e+02;
        else
            pb.y0(find(pbm.congrps==ig('X37')),1) = 1.79589e+02;
        end
        if(isKey(ix_,'X40'))
            pb.x0(ix_('X40'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X40')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X41'))
            pb.x0(ix_('X41'),1) = 2.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X41')),1) = 2.00000e+02;
        end
        if(isKey(ix_,'X42'))
            pb.x0(ix_('X42'),1) = 9.87694e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X42')),1) = 9.87694e+01;
        end
        if(isKey(ix_,'X43'))
            pb.x0(ix_('X43'),1) = 1.36703e+01;
        else
            pb.y0(find(pbm.congrps==ig('X43')),1) = 1.36703e+01;
        end
        if(isKey(ix_,'X44'))
            pb.x0(ix_('X44'),1) = 3.28255e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X44')),1) = 3.28255e+01;
        end
        if(isKey(ix_,'X45'))
            pb.x0(ix_('X45'),1) = 4.37635e+01;
        else
            pb.y0(find(pbm.congrps==ig('X45')),1) = 4.37635e+01;
        end
        if(isKey(ix_,'X46'))
            pb.x0(ix_('X46'),1) = -1.78833e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X46')),1) = -1.78833e+02;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eT1',iet_);
        elftv{it}{1} = 'ARC';
        elftp{it}{1} = 'C1';
        elftp{it}{2} = 'C2';
        elftp{it}{3} = 'C3';
        [it,iet_] = s2mpjlib( 'ii', 'eT4',iet_);
        elftv{it}{1} = 'ARC';
        elftp{it}{1} = 'C1';
        elftp{it}{2} = 'C2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eT4';
        ielftype(ie) = iet_('eT4');
        vname = 'X1';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.48060e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.51200e+02;
        ename = 'E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eT4';
        ielftype(ie) = iet_('eT4');
        vname = 'X2';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.91526e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.46300e+01;
        ename = 'E3';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eT4';
        ielftype(ie) = iet_('eT4');
        vname = 'X3';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.07752e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.81400e+01;
        ename = 'E4';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X4';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.90000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.22000e+02;
        ename = 'E5';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X5';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        ename = 'E6';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X6';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.63000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+02;
        ename = 'E7';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X7';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.10000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.22000e+02;
        ename = 'E8';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X8';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.45000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.00000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+02;
        ename = 'E9';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X9';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.40000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.22000e+02;
        ename = 'E10';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X10';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.00000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.50000e+01;
        ename = 'E11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X11';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.00000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.07000e+02;
        ename = 'E12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X12';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.20000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.80000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X13';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.00000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.80000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X14';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.00000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        ename = 'E15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X15';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.00000e+01;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.12200e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.30000e+02;
        ename = 'E16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X16';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.50000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.22000e+02;
        ename = 'E17';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X17';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.10000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E18';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X18';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.00000e+01;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.80000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.18000e+02;
        ename = 'E19';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eT4';
        ielftype(ie) = iet_('eT4');
        vname = 'X19';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.84530e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.12970e+02;
        ename = 'E20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X20';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.60000e+04;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.80000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E21';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eT4';
        ielftype(ie) = iet_('eT4');
        vname = 'X21';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.86880e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.60610e+02;
        ename = 'E22';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X22';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.20000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.36100e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.30000e+02;
        ename = 'E23';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X23';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.60000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.50000e+01;
        ename = 'E24';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X24';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.50000e+01;
        ename = 'E25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X25';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.60000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E26';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X26';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.30000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E27';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X27';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.20000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.24000e+02;
        ename = 'E28';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X28';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X29';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.90000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.13000e+02;
        ename = 'E30';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X30';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.80000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.13000e+02;
        ename = 'E31';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X31';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.70000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E32';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X32';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.10000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.50000e+01;
        ename = 'E33';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X33';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        ename = 'E34';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X34';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.30000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.13000e+02;
        ename = 'E35';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X35';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.20000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.50000e+01;
        ename = 'E36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X36';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.80000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E37';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X37';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.40000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+02;
        ename = 'E38';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X38';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.31000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.00000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        ename = 'E39';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X39';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.65000e+02;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+02;
        ename = 'E40';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X40';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.10000e+03;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.60000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.20000e+02;
        ename = 'E41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eT1';
            ielftype(ie) = iet_('eT1');
        end
        vname = 'X41';
        [iv,ix_,pb] =...
              s2mpjlib('nlx',vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02);
        posev = find(strcmp('ARC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('C1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.23000e+01;
        [~,posep] = ismember('C2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+01;
        [~,posep] = ismember('C3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.00000e+02;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E12');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E14');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E16');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E17');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E18');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E19');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E20');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E21');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E22');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E23');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E24');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E26');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E27');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E28');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E29');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E30');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E31');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E32');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E33');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E34');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E35');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E36');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E37');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E38');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E39');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E40');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E41');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -3.2393D+04
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CONR2-MN-46-31';
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

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eT1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP = 850559.0e0/2.85*pbm.elpar{iel_}(1);
        TMP = TMP/(pbm.elpar{iel_}(3)^1.85);
        TMP = TMP/(pbm.elpar{iel_}(2)^4.87);
        X = abs(EV_(1));
        XEXP = X^0.85;
        varargout{1} = TMP*X^2*XEXP;
        if(nargout>1)
            g_(1,1) = 2.85*TMP*EV_(1)*XEXP;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 5.2725*TMP*XEXP;
                varargout{3} = H_;
            end
        end

    case 'eT4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EPS2 = 1.0e-14;
        SQC1 = sqrt(pbm.elpar{iel_}(1));
        X = min(EV_(1),SQC1);
        TMP = pbm.elpar{iel_}(2)*(pbm.elpar{iel_}(1)-X*X);
        TMP = sqrt(max(TMP,EPS2));
        TMP2 = sqrt(pbm.elpar{iel_}(2))*asin(X/SQC1);
        varargout{1} = 0.5*(-X*TMP-pbm.elpar{iel_}(1)*TMP2);
        if(nargout>1)
            g_(1,1) = -TMP;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = pbm.elpar{iel_}(2)*X/TMP;
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

