function varargout = AVION2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : AVION2
%    *********
% 
%    Dassault France avion (airplane design) problem
% 
%    SIF input:  A. R. Conn, June 1993.
% 
%    classification = 'C-COLR2-RN-49-15'
% 
%    Define useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'AVION2';

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
        v_('790') = 790.0;
        v_('1/790') = 1.0/v_('790');
        v_('24000') = 24000.0;
        v_('1/24000') = 1.0/v_('24000');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','SR',ix_);
        pb.xnames{iv} = 'SR';
        [iv,ix_] = s2mpjlib('ii','LR',ix_);
        pb.xnames{iv} = 'LR';
        [iv,ix_] = s2mpjlib('ii','PK',ix_);
        pb.xnames{iv} = 'PK';
        [iv,ix_] = s2mpjlib('ii','EF',ix_);
        pb.xnames{iv} = 'EF';
        [iv,ix_] = s2mpjlib('ii','SX',ix_);
        pb.xnames{iv} = 'SX';
        [iv,ix_] = s2mpjlib('ii','LX',ix_);
        pb.xnames{iv} = 'LX';
        [iv,ix_] = s2mpjlib('ii','SD',ix_);
        pb.xnames{iv} = 'SD';
        [iv,ix_] = s2mpjlib('ii','SK',ix_);
        pb.xnames{iv} = 'SK';
        [iv,ix_] = s2mpjlib('ii','ST',ix_);
        pb.xnames{iv} = 'ST';
        [iv,ix_] = s2mpjlib('ii','SF',ix_);
        pb.xnames{iv} = 'SF';
        [iv,ix_] = s2mpjlib('ii','LF',ix_);
        pb.xnames{iv} = 'LF';
        [iv,ix_] = s2mpjlib('ii','AM',ix_);
        pb.xnames{iv} = 'AM';
        [iv,ix_] = s2mpjlib('ii','CA',ix_);
        pb.xnames{iv} = 'CA';
        [iv,ix_] = s2mpjlib('ii','CB',ix_);
        pb.xnames{iv} = 'CB';
        [iv,ix_] = s2mpjlib('ii','SO',ix_);
        pb.xnames{iv} = 'SO';
        [iv,ix_] = s2mpjlib('ii','SS',ix_);
        pb.xnames{iv} = 'SS';
        [iv,ix_] = s2mpjlib('ii','IMPDER',ix_);
        pb.xnames{iv} = 'IMPDER';
        [iv,ix_] = s2mpjlib('ii','IMPK',ix_);
        pb.xnames{iv} = 'IMPK';
        [iv,ix_] = s2mpjlib('ii','IMPFUS',ix_);
        pb.xnames{iv} = 'IMPFUS';
        [iv,ix_] = s2mpjlib('ii','QI',ix_);
        pb.xnames{iv} = 'QI';
        [iv,ix_] = s2mpjlib('ii','PT',ix_);
        pb.xnames{iv} = 'PT';
        [iv,ix_] = s2mpjlib('ii','MV',ix_);
        pb.xnames{iv} = 'MV';
        [iv,ix_] = s2mpjlib('ii','MC',ix_);
        pb.xnames{iv} = 'MC';
        [iv,ix_] = s2mpjlib('ii','MD',ix_);
        pb.xnames{iv} = 'MD';
        [iv,ix_] = s2mpjlib('ii','PD',ix_);
        pb.xnames{iv} = 'PD';
        [iv,ix_] = s2mpjlib('ii','NS',ix_);
        pb.xnames{iv} = 'NS';
        [iv,ix_] = s2mpjlib('ii','VS',ix_);
        pb.xnames{iv} = 'VS';
        [iv,ix_] = s2mpjlib('ii','CR',ix_);
        pb.xnames{iv} = 'CR';
        [iv,ix_] = s2mpjlib('ii','PM',ix_);
        pb.xnames{iv} = 'PM';
        [iv,ix_] = s2mpjlib('ii','DV',ix_);
        pb.xnames{iv} = 'DV';
        [iv,ix_] = s2mpjlib('ii','MZ',ix_);
        pb.xnames{iv} = 'MZ';
        [iv,ix_] = s2mpjlib('ii','VN',ix_);
        pb.xnames{iv} = 'VN';
        [iv,ix_] = s2mpjlib('ii','QV',ix_);
        pb.xnames{iv} = 'QV';
        [iv,ix_] = s2mpjlib('ii','QF',ix_);
        pb.xnames{iv} = 'QF';
        [iv,ix_] = s2mpjlib('ii','IMPTRAIN',ix_);
        pb.xnames{iv} = 'IMPTRAIN';
        [iv,ix_] = s2mpjlib('ii','IMPMOT',ix_);
        pb.xnames{iv} = 'IMPMOT';
        [iv,ix_] = s2mpjlib('ii','IMPNMOT',ix_);
        pb.xnames{iv} = 'IMPNMOT';
        [iv,ix_] = s2mpjlib('ii','IMPPET',ix_);
        pb.xnames{iv} = 'IMPPET';
        [iv,ix_] = s2mpjlib('ii','IMPPIL',ix_);
        pb.xnames{iv} = 'IMPPIL';
        [iv,ix_] = s2mpjlib('ii','IMPCAN',ix_);
        pb.xnames{iv} = 'IMPCAN';
        [iv,ix_] = s2mpjlib('ii','IMPSNA',ix_);
        pb.xnames{iv} = 'IMPSNA';
        [iv,ix_] = s2mpjlib('ii','MS',ix_);
        pb.xnames{iv} = 'MS';
        [iv,ix_] = s2mpjlib('ii','EL',ix_);
        pb.xnames{iv} = 'EL';
        [iv,ix_] = s2mpjlib('ii','DE',ix_);
        pb.xnames{iv} = 'DE';
        [iv,ix_] = s2mpjlib('ii','DS',ix_);
        pb.xnames{iv} = 'DS';
        [iv,ix_] = s2mpjlib('ii','IMPVOIL',ix_);
        pb.xnames{iv} = 'IMPVOIL';
        [iv,ix_] = s2mpjlib('ii','NM',ix_);
        pb.xnames{iv} = 'NM';
        [iv,ix_] = s2mpjlib('ii','NP',ix_);
        pb.xnames{iv} = 'NP';
        [iv,ix_] = s2mpjlib('ii','NG',ix_);
        pb.xnames{iv} = 'NG';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','E1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SD');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SR');
        valA(end+1) = -0.13;
        [ig,ig_] = s2mpjlib('ii','E2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SX');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SR');
        valA(end+1) = -0.7;
        [ig,ig_] = s2mpjlib('ii','E3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LX');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LR');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E4',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SK');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SF');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ST');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SD');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SX');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SK');
        valA(end+1) = -2.0;
        [ig,ig_] = s2mpjlib('ii','E6',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('CA');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E7',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('AM');
        valA(end+1) = -2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SO');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SS');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E8',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('AM');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E9',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPDER');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SD');
        valA(end+1) = -27.5;
        [ig,ig_] = s2mpjlib('ii','E10',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPK');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SK');
        valA(end+1) = -70.0;
        [ig,ig_] = s2mpjlib('ii','E11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPFUS');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SF');
        valA(end+1) = -20.0;
        [ig,ig_] = s2mpjlib('ii','E12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MD');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MV');
        valA(end+1) = -2.0;
        [ig,ig_] = s2mpjlib('ii','E13',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QI');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E14',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PT');
        valA(end+1) = 1000.0;
        [ig,ig_] = s2mpjlib('ii','E15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E15';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QF');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QI');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QV');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E16',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('VN');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('VS');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QF');
        valA(end+1) = v_('1/790');
        [ig,ig_] = s2mpjlib('ii','E17',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E17';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPTRAIN');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MV');
        valA(end+1) = -0.137;
        [ig,ig_] = s2mpjlib('ii','E18',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPMOT');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E19',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E19';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPNMOT');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NM');
        valA(end+1) = -35.0;
        [ig,ig_] = s2mpjlib('ii','E20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E20';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPPET');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QI');
        valA(end+1) = -0.043;
        [ig,ig_] = s2mpjlib('ii','E21',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E21';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPPIL');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NP');
        valA(end+1) = -200.0;
        [ig,ig_] = s2mpjlib('ii','E22',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E22';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPCAN');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NG');
        valA(end+1) = -120.0;
        [ig,ig_] = s2mpjlib('ii','E23',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E23';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPSNA');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NS');
        valA(end+1) = -300.0;
        [ig,ig_] = s2mpjlib('ii','E24',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E24';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MC');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MV');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NP');
        valA(end+1) = 95.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NG');
        valA(end+1) = 70.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NM');
        valA(end+1) = 660.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QI');
        valA(end+1) = 0.5;
        [ig,ig_] = s2mpjlib('ii','E25',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E25';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MZ');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPTRAIN');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPNMOT');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPPET');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPPIL');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPCAN');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPSNA');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E26',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ST');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E27',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SR');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E28',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('QV');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E29',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SO');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E30',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('SS');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E31',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('CB');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E32',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IMPVOIL');
        valA(end+1) = 1.0;
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
        pbm.gconst(ig_('E13')) = 1000.0;
        pbm.gconst(ig_('E16')) = -2.0;
        pbm.gconst(ig_('E23')) = 400.0;
        pbm.gconst(ig_('E24')) = 380.0;
        pbm.gconst(ig_('E25')) = -290.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('SR'),1) = 10.0;
        pb.xlower(ix_('LR'),1) = 0.0;
        pb.xlower(ix_('PK'),1) = 0.0;
        pb.xlower(ix_('EF'),1) = 0.0;
        pb.xlower(ix_('SX'),1) = 7.0;
        pb.xlower(ix_('LX'),1) = 1.5;
        pb.xlower(ix_('SD'),1) = 2.0;
        pb.xlower(ix_('SK'),1) = 2.0;
        pb.xlower(ix_('ST'),1) = 30.0;
        pb.xlower(ix_('SF'),1) = 20.0;
        pb.xlower(ix_('LF'),1) = 0.001;
        pb.xlower(ix_('AM'),1) = 0.0;
        pb.xlower(ix_('CA'),1) = -0.2;
        pb.xlower(ix_('CB'),1) = 0.1;
        pb.xlower(ix_('SO'),1) = 0.0;
        pb.xlower(ix_('SS'),1) = 0.0;
        pb.xlower(ix_('IMPDER'),1) = 100.0;
        pb.xlower(ix_('IMPK'),1) = 500.0;
        pb.xlower(ix_('IMPFUS'),1) = 500.0;
        pb.xlower(ix_('QI'),1) = 1000.0;
        pb.xlower(ix_('PT'),1) = 2.0;
        pb.xlower(ix_('MV'),1) = 2000.0;
        pb.xlower(ix_('MC'),1) = 3000.0;
        pb.xlower(ix_('MD'),1) = 5000.0;
        pb.xlower(ix_('PD'),1) = 0.2;
        pb.xlower(ix_('NS'),1) = 1.0;
        pb.xlower(ix_('VS'),1) = 0.0;
        pb.xlower(ix_('CR'),1) = 100.0;
        pb.xlower(ix_('PM'),1) = 4.0;
        pb.xlower(ix_('DV'),1) = 0.0;
        pb.xlower(ix_('MZ'),1) = 500.0;
        pb.xlower(ix_('VN'),1) = 10.0;
        pb.xlower(ix_('QV'),1) = 250.0;
        pb.xlower(ix_('QF'),1) = 750.0;
        pb.xlower(ix_('IMPTRAIN'),1) = 250.0;
        pb.xlower(ix_('IMPMOT'),1) = 10.0;
        pb.xlower(ix_('IMPNMOT'),1) = 35.0;
        pb.xlower(ix_('IMPPET'),1) = 100.0;
        pb.xlower(ix_('IMPPIL'),1) = 200.0;
        pb.xlower(ix_('IMPCAN'),1) = 120.0;
        pb.xlower(ix_('IMPSNA'),1) = 700.0;
        pb.xlower(ix_('MS'),1) = 100.0;
        pb.xlower(ix_('EL'),1) = 2.0;
        pb.xlower(ix_('DE'),1) = 0.0;
        pb.xlower(ix_('DS'),1) = 0.0;
        pb.xlower(ix_('IMPVOIL'),1) = 500.0;
        pb.xupper(ix_('SR')) = 150.0;
        pb.xupper(ix_('LR')) = 10.0;
        pb.xupper(ix_('PK')) = 10.0;
        pb.xupper(ix_('EF')) = 5.0;
        pb.xupper(ix_('SX')) = 120.0;
        pb.xupper(ix_('LX')) = 8.0;
        pb.xupper(ix_('SD')) = 20.0;
        pb.xupper(ix_('SK')) = 30.0;
        pb.xupper(ix_('ST')) = 500.0;
        pb.xupper(ix_('SF')) = 200.0;
        pb.xupper(ix_('LF')) = 20.0;
        pb.xupper(ix_('AM')) = 10.0;
        pb.xupper(ix_('CA')) = -0.001;
        pb.xupper(ix_('CB')) = 2.0;
        pb.xupper(ix_('SO')) = 1.0;
        pb.xupper(ix_('SS')) = 2.0;
        pb.xupper(ix_('IMPDER')) = 1000.0;
        pb.xupper(ix_('IMPK')) = 5000.0;
        pb.xupper(ix_('IMPFUS')) = 5000.0;
        pb.xupper(ix_('QI')) = 20000.0;
        pb.xupper(ix_('PT')) = 30.0;
        pb.xupper(ix_('MV')) = 20000.0;
        pb.xupper(ix_('MC')) = 30000.0;
        pb.xupper(ix_('MD')) = 50000.0;
        pb.xupper(ix_('PD')) = 0.8;
        pb.xupper(ix_('NS')) = 5.0;
        pb.xupper(ix_('VS')) = 20.0;
        pb.xupper(ix_('CR')) = 400.0;
        pb.xupper(ix_('PM')) = 15.0;
        pb.xupper(ix_('DV')) = 10.0;
        pb.xupper(ix_('MZ')) = 10000.0;
        pb.xupper(ix_('VN')) = 50.0;
        pb.xupper(ix_('QV')) = 5000.0;
        pb.xupper(ix_('QF')) = 15000.0;
        pb.xupper(ix_('IMPTRAIN')) = 3000.0;
        pb.xupper(ix_('IMPMOT')) = 5000.0;
        pb.xupper(ix_('IMPNMOT')) = 70.0;
        pb.xupper(ix_('IMPPET')) = 3000.0;
        pb.xupper(ix_('IMPPIL')) = 400.0;
        pb.xupper(ix_('IMPCAN')) = 240.0;
        pb.xupper(ix_('IMPSNA')) = 1900.0;
        pb.xupper(ix_('MS')) = 1000.0;
        pb.xupper(ix_('EL')) = 20.0;
        pb.xupper(ix_('DE')) = 1.0;
        pb.xupper(ix_('DS')) = 2.0;
        pb.xupper(ix_('IMPVOIL')) = 5000.0;
        pb.xupper(ix_('NM')) = 2.0;
        pb.xlower(ix_('NM'),1) = 1.0;
        pb.xupper(ix_('NP')) = 2.0;
        pb.xlower(ix_('NP'),1) = 1.0;
        pb.xupper(ix_('NG')) = 2.0;
        pb.xlower(ix_('NG'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('SR'),1) = 2.7452e+01;
        pb.x0(ix_('LR'),1) = 1.5000;
        pb.x0(ix_('PK'),1) = 1.0000e+01;
        pb.x0(ix_('EF'),1) = 0.0000;
        pb.x0(ix_('SX'),1) = 1.9217e+01;
        pb.x0(ix_('LX'),1) = 1.5000;
        pb.x0(ix_('SD'),1) = 3.5688;
        pb.x0(ix_('SK'),1) = 4.0696;
        pb.x0(ix_('ST'),1) = 3.4315e+01;
        pb.x0(ix_('SF'),1) = 8.8025e+01;
        pb.x0(ix_('LF'),1) = 5.1306;
        pb.x0(ix_('AM'),1) = 0.0000;
        pb.x0(ix_('CA'),1) = -1.4809e-01;
        pb.x0(ix_('CB'),1) = 7.5980e-01;
        pb.x0(ix_('SO'),1) = 0.0000;
        pb.x0(ix_('SS'),1) = 0.0000;
        pb.x0(ix_('IMPDER'),1) = 1.1470e+02;
        pb.x0(ix_('IMPK'),1) = 5.0000e+02;
        pb.x0(ix_('IMPFUS'),1) = 1.7605e+03;
        pb.x0(ix_('QI'),1) = 2.3256e+03;
        pb.x0(ix_('PT'),1) = 5.6788;
        pb.x0(ix_('MV'),1) = 1.4197e+04;
        pb.x0(ix_('MC'),1) = 1.2589e+04;
        pb.x0(ix_('MD'),1) = 2.8394e+04;
        pb.x0(ix_('PD'),1) = 2.0000e-01;
        pb.x0(ix_('NS'),1) = 1.0000;
        pb.x0(ix_('VS'),1) = 0.0000;
        pb.x0(ix_('CR'),1) = 1.0000e+02;
        pb.x0(ix_('PM'),1) = 1.5000e+01;
        pb.x0(ix_('DV'),1) = 0.0000;
        pb.x0(ix_('MZ'),1) = 5.0000e+02;
        pb.x0(ix_('VN'),1) = 1.0000e+01;
        pb.x0(ix_('QV'),1) = 8.1490e+02;
        pb.x0(ix_('QF'),1) = 3.1405e+03;
        pb.x0(ix_('IMPTRAIN'),1) = 1.9450e+03;
        pb.x0(ix_('IMPMOT'),1) = 1.9085e+02;
        pb.x0(ix_('IMPNMOT'),1) = 3.5000e+01;
        pb.x0(ix_('IMPPET'),1) = 1.0000e+02;
        pb.x0(ix_('IMPPIL'),1) = 2.0000e+02;
        pb.x0(ix_('IMPCAN'),1) = 1.2000e+02;
        pb.x0(ix_('IMPSNA'),1) = 7.0000e+02;
        pb.x0(ix_('MS'),1) = 1.0000e+03;
        pb.x0(ix_('EL'),1) = 4.9367;
        pb.x0(ix_('DE'),1) = 0.0000;
        pb.x0(ix_('DS'),1) = 0.0000;
        pb.x0(ix_('IMPVOIL'),1) = 5.0000e+03;
        pb.x0(ix_('NM'),1) = 1.0000;
        pb.x0(ix_('NP'),1) = 1.0000;
        pb.x0(ix_('NG'),1) = 1.0000;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eQDdSQ',iet_);
        elftv{it}{1} = 'W';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        elftv{it}{4} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en12',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en12d1',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eQT',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en1dLIN',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSQRT',iet_);
        elftv{it}{1} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSURD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSQPRD',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eCBdSQQD',iet_);
        elftv{it}{1} = 'W';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        elftv{it}{4} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSREL',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        elftv{it}{3} = 'X';
        elftv{it}{4} = 'Y';
        elftv{it}{5} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'EL1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'PK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SR';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQDdSQ';
        ielftype(ie) = iet_('eQDdSQ');
        vname = 'SS';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SO';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en12';
        ielftype(ie) = iet_('en12');
        vname = 'EF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en12d1';
        ielftype(ie) = iet_('en12d1');
        vname = 'SO';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CA';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'SD';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'SK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'MV';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'MD';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PD';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'MZ';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CR';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'DV';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en1dLIN';
        ielftype(ie) = iet_('en1dLIN');
        vname = 'PT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRT';
        ielftype(ie) = iet_('eSQRT');
        vname = 'PT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'SR';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL14';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'MD';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MS';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL15';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSURD';
        ielftype(ie) = iet_('eSURD');
        vname = 'SX';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'EL';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LX';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL16';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQPRD';
        ielftype(ie) = iet_('eSQPRD');
        vname = 'DE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQPRD';
        ielftype(ie) = iet_('eSQPRD');
        vname = 'DS';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCBdSQQD';
        ielftype(ie) = iet_('eCBdSQQD');
        vname = 'VN';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CA';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SO';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSREL';
        ielftype(ie) = iet_('eSREL');
        vname = 'SX';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MC';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LX';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SR';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'EL';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSQUARE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('E4');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.01;
        ig = ig_('E6');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('E7');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.01;
        ig = ig_('E8');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.25;
        ig = ig_('E9');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.3;
        ig = ig_('E10');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.6;
        ig = ig_('E13');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/24000');
        ig = ig_('E14');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('E16');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EL10');
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('E18');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1000.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EL12');
        pbm.grelw{ig}(posel) = -12.0;
        ig = ig_('E26');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.25;
        ig = ig_('E27');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('E28');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.4;
        ig = ig_('E29');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL16');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.785;
        ig = ig_('E30');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL17');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.785;
        ig = ig_('E31');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.0;
        ig = ig_('E32');
        pbm.grftype{ig} = 'gSQUARE';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EL19');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.15;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               9.46801297093018D+07
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-RN-49-15';
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

% ***********************
%  SET UP THE FUNCTIONS *
% ***********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eQDdSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        QD = EV_(1)-EV_(2)-EV_(3)*EV_(4);
        SQ = EV_(4)^2;
        RSQ = 1.0e0/SQ;
        QDOSQ = QD/SQ;
        varargout{1} = QDOSQ;
        if(nargout>1)
            g_(1,1) = RSQ;
            g_(2,1) = -RSQ;
            g_(3,1) = -1.0e0/EV_(4);
            g_(4,1) = -EV_(3)/SQ-2.0e0*QDOSQ/EV_(4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,4) = -2.0e0/(SQ*EV_(4));
                H_(4,1) = H_(1,4);
                H_(2,4) = 2.0e0/(SQ*EV_(4));
                H_(4,2) = H_(2,4);
                H_(3,4) = RSQ;
                H_(4,3) = H_(3,4);
                H_(4,4) = (4.0e0*EV_(3))/(SQ*EV_(4))+6.0e0*QDOSQ/SQ;
                varargout{3} = H_;
            end
        end

    case 'en12'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)/EV_(2);
        if(nargout>1)
            g_(1,1) = 1.0/EV_(2);
            g_(2,1) = -EV_(1)/EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0e0/EV_(2)^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = (2.0e0*EV_(1))/EV_(2)^3;
                varargout{3} = H_;
            end
        end

    case 'en12d1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ZSQ = EV_(3)^2;
        YSQ = EV_(2)^2;
        varargout{1} = (EV_(1)*YSQ)/EV_(3);
        if(nargout>1)
            g_(1,1) = YSQ/EV_(3);
            g_(2,1) = (2.0e0*EV_(1)*EV_(2))/EV_(3);
            g_(3,1) = -(EV_(1)*YSQ)/ZSQ;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = (2.0e0*EV_(2))/EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = -YSQ/ZSQ;
                H_(3,1) = H_(1,3);
                H_(2,2) = (2.0e0*EV_(1))/EV_(3);
                H_(2,3) = -(2.0e0*EV_(1)*EV_(2))/ZSQ;
                H_(3,2) = H_(2,3);
                H_(3,3) = (2.0e0*EV_(1)*YSQ)/(ZSQ*EV_(3));
                varargout{3} = H_;
            end
        end

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0e0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e0;
                varargout{3} = H_;
            end
        end

    case 'eQT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)/EV_(2);
        if(nargout>1)
            g_(1,1) = 1.0/EV_(2);
            g_(2,1) = -EV_(1)/EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0e0/EV_(2)^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = (2.0e0*EV_(1))/EV_(2)^3;
                varargout{3} = H_;
            end
        end

    case 'en1dLIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LIN = EV_(2)+20.0e0;
        varargout{1} = EV_(1)/LIN;
        if(nargout>1)
            g_(1,1) = 1.0/LIN;
            g_(2,1) = -EV_(1)/LIN^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0e0/LIN^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = (2.0e0*EV_(1))/LIN^3;
                varargout{3} = H_;
            end
        end

    case 'eSQRT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        RTZ = sqrt(EV_(1));
        varargout{1} = RTZ;
        if(nargout>1)
            g_(1,1) = 0.5e0/RTZ;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -0.25e0/(EV_(1)*RTZ);
                varargout{3} = H_;
            end
        end

    case 'eSURD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        RTX = sqrt(EV_(1));
        RTZ = sqrt(EV_(3));
        XRTX = EV_(1)*sqrt(EV_(1));
        ZRTZ = EV_(3)*sqrt(EV_(3));
        varargout{1} = XRTX*EV_(2)/RTZ;
        if(nargout>1)
            g_(1,1) = 1.5e0*RTX*EV_(2)/RTZ;
            g_(2,1) = XRTX/RTZ;
            g_(3,1) = -(0.5e0*XRTX*EV_(2))/ZRTZ;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 0.75e0*EV_(2)/(RTZ*RTX);
                H_(1,2) = 1.5e0*RTX/RTZ;
                H_(2,1) = H_(1,2);
                H_(1,3) = -(0.75e0*RTX*EV_(2))/ZRTZ;
                H_(3,1) = H_(1,3);
                H_(2,3) = -(0.5e0*XRTX)/ZRTZ;
                H_(3,2) = H_(2,3);
                H_(3,3) = (0.75e0*XRTX*EV_(2))/(ZRTZ*EV_(3));
                varargout{3} = H_;
            end
        end

    case 'eSQPRD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^2*EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0e0*EV_(1)*EV_(2);
            g_(2,1) = EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0e0*EV_(2);
                H_(1,2) = 2.0e0*EV_(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eCBdSQQD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        YCB = EV_(3)^3;
        CB = EV_(1)-EV_(2)*YCB;
        TMZY = 3.0e0-EV_(4)*EV_(3);
        SQQD = EV_(3)^2*TMZY;
        DYSQQD = 3.0e0*EV_(3)*(2.0e0-EV_(4)*EV_(3));
        varargout{1} = CB/SQQD;
        if(nargout>1)
            g_(1,1) = 1.0e0/SQQD;
            g_(2,1) = -YCB/SQQD;
            g_(3,1) = -3.0e0*EV_(2)*EV_(3)^2/SQQD-(CB*DYSQQD)/(SQQD*SQQD);
            g_(4,1) = (YCB*CB)/SQQD^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,3) = -DYSQQD/SQQD^2;
                H_(3,1) = H_(1,3);
                H_(1,4) = YCB/SQQD^2;
                H_(4,1) = H_(1,4);
                H_(2,3) = -3.0e0*(1.0e0-(TMZY-1.0e0)/TMZY)/TMZY;
                H_(3,2) = H_(2,3);
                H_(2,4) = -(YCB^2*CB)/SQQD^2;
                H_(4,2) = H_(2,4);
                H_(3,3) = -6.0e0*EV_(2)*EV_(3)/SQQD+(3.0e0*EV_(2)*EV_(3)^2*DYSQQD)/...
                     SQQD^2+(2.0e0*CB*DYSQQD^2)/SQQD^3-(6.0e0*CB*(1.0-EV_(4)*EV_(3)))/SQQD^2;
                H_(3,4) = -(3.0e0*EV_(2)*EV_(3))/TMZY^2-(6.0e0*EV_(3)^4*(TMZY-1.0e0)*CB)/...
                     SQQD^3+3.0e0*CB*EV_(3)^2/SQQD^2;
                H_(4,3) = H_(3,4);
                H_(4,4) = 2.0e0*(YCB^2*CB)/SQQD^3;
                varargout{3} = H_;
            end
        end

    case 'eSREL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SRLIN = 15.0e0+0.15e0*EV_(1);
        SRPD = EV_(2)*EV_(3);
        SRQD = 50.0e0*EV_(4)*EV_(5);
        SRQT = SRPD/SRQD;
        SRRT = sqrt(SRQT);
        SRSUM = 15.0e0+0.3e0*EV_(1);
        SRVLN = EV_(1)*SRLIN;
        varargout{1} = EV_(1)*SRLIN*(SRQT*SRRT+8.0e0);
        if(nargout>1)
            g_(1,1) = SRSUM*(SRQT*SRRT+8.0e0);
            g_(2,1) = 1.5e0*(SRVLN*SRRT*EV_(3)/SRQD);
            g_(3,1) = 1.5e0*(SRVLN*SRRT*EV_(2)/SRQD);
            g_(4,1) = -1.5e0*SRVLN*SRRT*SRQT/EV_(4);
            g_(5,1) = -1.5e0*SRVLN*SRRT*SRQT/EV_(5);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = 0.3e0*(SRQT*SRRT+8.0);
                H_(1,2) = 1.5e0*(SRSUM*SRRT*EV_(3)/SRQD);
                H_(2,1) = H_(1,2);
                H_(1,3) = 1.5e0*(SRSUM*SRRT*EV_(2)/SRQD);
                H_(3,1) = H_(1,3);
                H_(1,4) = -1.5e0*SRSUM*SRRT*SRQT/EV_(4);
                H_(4,1) = H_(1,4);
                H_(1,5) = -1.5e0*SRSUM*SRRT*SRQT/EV_(5);
                H_(5,1) = H_(1,5);
                H_(2,2) = (0.75e0*SRVLN*EV_(3)^2)/(SRQD^2*SRRT);
                H_(2,3) = SRVLN*((0.75e0*SRPD)/(SRQD^2*SRRT)+(1.5e0*SRRT)/SRQD);
                H_(3,2) = H_(2,3);
                H_(2,4) = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_(2)*EV_(4));
                H_(4,2) = H_(2,4);
                H_(2,5) = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_(2)*EV_(5));
                H_(5,2) = H_(2,5);
                H_(3,3) = (SRVLN*0.75e0*EV_(2)*EV_(2))/(SRRT*SRQD^2);
                H_(3,4) = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_(3)*EV_(4));
                H_(4,3) = H_(3,4);
                H_(3,5) = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_(3)*EV_(5));
                H_(5,3) = H_(3,5);
                H_(4,4) = (SRVLN*3.75e0*SRRT*SRQT)/EV_(4)^2;
                H_(4,5) = (SRVLN*2.25e0*SRRT*SRQT)/(EV_(4)*EV_(5));
                H_(5,4) = H_(4,5);
                H_(5,5) = (SRVLN*3.75e0*SRRT*SRQT)/EV_(5)^2;
                varargout{3} = H_;
            end
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
                H_ = 2.0e0;
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

