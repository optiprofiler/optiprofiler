function varargout = LAUNCH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    The objective function to be minimized represents the total cost of
%    the development and launching of a 3 stages space launching vehicle.
%    Constraints are imposed on physical interrelations between the variables
%    and performance.
% 
%    The problem is highly non-convex. 
% 
%    Source:
%    B. Rush, J. Bracken and G. McCormick,
%    "A nonliner programming model for launch vehicle design and costing",
%    Operations Research, pp. 185-210, 1967.
% 
%    SIF input: P. Driscoll, Virginia Tech., April 1993,
%               corrected and simplified by Ph. L. Toint, May 1993.
% 
%    classification = 'C-COOR2-MY-25-28'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LAUNCH';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','AW1',ix_);
        pb.xnames{iv} = 'AW1';
        [iv,ix_] = s2mpjlib('ii','IW1',ix_);
        pb.xnames{iv} = 'IW1';
        [iv,ix_] = s2mpjlib('ii','MF1',ix_);
        pb.xnames{iv} = 'MF1';
        [iv,ix_] = s2mpjlib('ii','TT1',ix_);
        pb.xnames{iv} = 'TT1';
        [iv,ix_] = s2mpjlib('ii','PW1',ix_);
        pb.xnames{iv} = 'PW1';
        [iv,ix_] = s2mpjlib('ii','ET1',ix_);
        pb.xnames{iv} = 'ET1';
        [iv,ix_] = s2mpjlib('ii','S1L',ix_);
        pb.xnames{iv} = 'S1L';
        [iv,ix_] = s2mpjlib('ii','AW2',ix_);
        pb.xnames{iv} = 'AW2';
        [iv,ix_] = s2mpjlib('ii','IW2',ix_);
        pb.xnames{iv} = 'IW2';
        [iv,ix_] = s2mpjlib('ii','MF2',ix_);
        pb.xnames{iv} = 'MF2';
        [iv,ix_] = s2mpjlib('ii','TT2',ix_);
        pb.xnames{iv} = 'TT2';
        [iv,ix_] = s2mpjlib('ii','PW2',ix_);
        pb.xnames{iv} = 'PW2';
        [iv,ix_] = s2mpjlib('ii','ET2',ix_);
        pb.xnames{iv} = 'ET2';
        [iv,ix_] = s2mpjlib('ii','S2L',ix_);
        pb.xnames{iv} = 'S2L';
        [iv,ix_] = s2mpjlib('ii','AW3',ix_);
        pb.xnames{iv} = 'AW3';
        [iv,ix_] = s2mpjlib('ii','IW3',ix_);
        pb.xnames{iv} = 'IW3';
        [iv,ix_] = s2mpjlib('ii','MF3',ix_);
        pb.xnames{iv} = 'MF3';
        [iv,ix_] = s2mpjlib('ii','TT3',ix_);
        pb.xnames{iv} = 'TT3';
        [iv,ix_] = s2mpjlib('ii','PW3',ix_);
        pb.xnames{iv} = 'PW3';
        [iv,ix_] = s2mpjlib('ii','ET3',ix_);
        pb.xnames{iv} = 'ET3';
        [iv,ix_] = s2mpjlib('ii','S3L',ix_);
        pb.xnames{iv} = 'S3L';
        [iv,ix_] = s2mpjlib('ii','INW',ix_);
        pb.xnames{iv} = 'INW';
        [iv,ix_] = s2mpjlib('ii','BT1',ix_);
        pb.xnames{iv} = 'BT1';
        [iv,ix_] = s2mpjlib('ii','BT2',ix_);
        pb.xnames{iv} = 'BT2';
        [iv,ix_] = s2mpjlib('ii','BT3',ix_);
        pb.xnames{iv} = 'BT3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','STA1',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET1');
        valA(end+1) = 0.0002587;
        pbm.gscale(ig,1) = 1.0e+8;
        [ig,ig_] = s2mpjlib('ii','STA2',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET2');
        valA(end+1) = 0.0002587;
        pbm.gscale(ig,1) = 1.0e+8;
        [ig,ig_] = s2mpjlib('ii','STA3',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET3');
        valA(end+1) = 0.001958;
        pbm.gscale(ig,1) = 1.0e+8;
        [ig,ig_] = s2mpjlib('ii','INST',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = 47.040096;
        pbm.gscale(ig,1) = 1.0e+8;
        [ig,ig_] = s2mpjlib('ii','LAUN',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = 0.003;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = 0.003;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = 0.003;
        pbm.gscale(ig,1) = 39215686.0;
        [ig,ig_] = s2mpjlib('ii','SGTH1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('AW1');
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGTH3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = 0.6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('AW2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGTH5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = 0.7;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('AW3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGTH2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET1');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGTH4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET2');
        valA(end+1) = 5.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGTH6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SGTH6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('ET3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SGSI1A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SGSI1A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -12.0;
        [ig,ig_] = s2mpjlib('ii','SGSI1B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SGSI1B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -17.0;
        [ig,ig_] = s2mpjlib('ii','SGSI2A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SGSI2A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -10.0;
        [ig,ig_] = s2mpjlib('ii','SGSI2B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SGSI2B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -13.0;
        [ig,ig_] = s2mpjlib('ii','SGSI3A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SGSI3A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -7.0;
        [ig,ig_] = s2mpjlib('ii','SGSI3B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SGSI3B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -10.0;
        [ig,ig_] = s2mpjlib('ii','TTIW1A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'TTIW1A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -1.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -1.2;
        [ig,ig_] = s2mpjlib('ii','TTIW1B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'TTIW1B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -1.4;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -1.4;
        [ig,ig_] = s2mpjlib('ii','TTIW2A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'TTIW2A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -0.6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -0.6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -0.6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -0.6;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -0.6;
        [ig,ig_] = s2mpjlib('ii','TTIW2B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'TTIW2B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -0.75;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -0.75;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -0.75;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -0.75;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -0.75;
        [ig,ig_] = s2mpjlib('ii','TTIW3A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'TTIW3A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -0.7;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -0.7;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -0.7;
        [ig,ig_] = s2mpjlib('ii','TTIW3B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'TTIW3B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TT3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -0.9;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -0.9;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -0.9;
        [ig,ig_] = s2mpjlib('ii','SMF1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SMF1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MF1');
        valA(end+1) = 20.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SMF2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SMF2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MF2');
        valA(end+1) = 20.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SMF3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SMF3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MF3');
        valA(end+1) = 20.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('IW3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('INW');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','SI1A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SI1A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = -240.0;
        [ig,ig_] = s2mpjlib('ii','SI1B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SI1B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW1');
        valA(end+1) = -290.0;
        [ig,ig_] = s2mpjlib('ii','SI2A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SI2A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -240.0;
        [ig,ig_] = s2mpjlib('ii','SI2B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SI2B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW2');
        valA(end+1) = -290.0;
        [ig,ig_] = s2mpjlib('ii','SI3A',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'SI3A';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -340.0;
        [ig,ig_] = s2mpjlib('ii','SI3B',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SI3B';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('PW3');
        valA(end+1) = -375.0;
        [ig,ig_] = s2mpjlib('ii','GLGCON',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'GLGCON';
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
        pbm.gconst(ig_('STA1')) = 247.963;
        pbm.gconst(ig_('STA2')) = 247.963;
        pbm.gconst(ig_('STA3')) = 32.591;
        pbm.gconst(ig_('INST')) = 35.5;
        pbm.gconst(ig_('TTIW1A')) = 24.0;
        pbm.gconst(ig_('TTIW1B')) = 28.0;
        pbm.gconst(ig_('TTIW2A')) = 12.0;
        pbm.gconst(ig_('TTIW2B')) = 15.0;
        pbm.gconst(ig_('TTIW3A')) = 14.0;
        pbm.gconst(ig_('TTIW3B')) = 18.0;
        pbm.gconst(ig_('SMF1')) = 20.0;
        pbm.gconst(ig_('SMF2')) = 20.0;
        pbm.gconst(ig_('SMF3')) = 20.0;
        pbm.gconst(ig_('GLGCON')) = 35000.0;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        grange(ig_('GLGCON')) = 15000.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 1.0e-8*ones(pb.n,1);
        pb.xupper = 1.0e+4*ones(pb.n,1);
        pb.xlower(ix_('S1L'),1) = 125.0;
        pb.xupper(ix_('S1L')) = 150.0;
        pb.xlower(ix_('S2L'),1) = 75.0;
        pb.xupper(ix_('S2L')) = 100.0;
        pb.xlower(ix_('S3L'),1) = 50.0;
        pb.xupper(ix_('S3L')) = 70.0;
        pb.xlower(ix_('MF1'),1) = 0.25;
        pb.xupper(ix_('MF1')) = 0.30;
        pb.xlower(ix_('MF2'),1) = 0.24;
        pb.xupper(ix_('MF2')) = 0.29;
        pb.xlower(ix_('MF3'),1) = 0.16;
        pb.xupper(ix_('MF3')) = 0.21;
        pb.xlower(ix_('INW'),1) = 2.5;
        pb.xupper(ix_('INW')) = 4.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('AW1'),1) = 68.0;
        pb.x0(ix_('IW1'),1) = 136.0;
        pb.x0(ix_('MF1'),1) = 0.29988744;
        pb.x0(ix_('TT1'),1) = 3733.0;
        pb.x0(ix_('PW1'),1) = 2177.0;
        pb.x0(ix_('ET1'),1) = 746.6;
        pb.x0(ix_('S1L'),1) = 125.0;
        pb.x0(ix_('AW2'),1) = 28.2;
        pb.x0(ix_('IW2'),1) = 47.0;
        pb.x0(ix_('MF2'),1) = 0.28939109;
        pb.x0(ix_('TT2'),1) = 480.0;
        pb.x0(ix_('PW2'),1) = 566.0;
        pb.x0(ix_('ET2'),1) = 96.0;
        pb.x0(ix_('S2L'),1) = 75.0;
        pb.x0(ix_('AW3'),1) = 11.2;
        pb.x0(ix_('IW3'),1) = 16.0;
        pb.x0(ix_('MF3'),1) = 0.20980926;
        pb.x0(ix_('TT3'),1) = 129.0;
        pb.x0(ix_('PW3'),1) = 145.0;
        pb.x0(ix_('ET3'),1) = 129.0;
        pb.x0(ix_('S3L'),1) = 50.0;
        pb.x0(ix_('INW'),1) = 2.5;
        pb.x0(ix_('BT1'),1) = 155.0;
        pb.x0(ix_('BT2'),1) = 314.0;
        pb.x0(ix_('BT3'),1) = 403.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD1',iet_);
        elftv{it}{1} = 'VA';
        elftv{it}{2} = 'VB';
        elftv{it}{3} = 'VC';
        elftv{it}{4} = 'VD';
        elftv{it}{5} = 'VE';
        [it,iet_] = s2mpjlib( 'ii', 'ePOWER',iet_);
        elftv{it}{1} = 'XX';
        elftp{it}{1} = 'PWR';
        elftp{it}{2} = 'SC';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'VA';
        elftv{it}{2} = 'VB';
        elftv{it}{3} = 'VC';
        elftv{it}{4} = 'VD';
        [it,iet_] = s2mpjlib( 'ii', 'eX7Y',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y1';
        elftv{it}{3} = 'Y2';
        elftv{it}{4} = 'Y3';
        elftv{it}{5} = 'Y4';
        elftv{it}{6} = 'Y5';
        elftv{it}{7} = 'Y6';
        elftv{it}{8} = 'Y7';
        [it,iet_] = s2mpjlib( 'ii', 'eX5Y',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y1';
        elftv{it}{3} = 'Y2';
        elftv{it}{4} = 'Y3';
        elftv{it}{5} = 'Y4';
        elftv{it}{6} = 'Y5';
        [it,iet_] = s2mpjlib( 'ii', 'eX3Y',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y1';
        elftv{it}{3} = 'Y2';
        elftv{it}{4} = 'Y3';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eBIG1',iet_);
        elftv{it}{1} = 'LH';
        elftv{it}{2} = 'TH';
        elftv{it}{3} = 'LL';
        elftv{it}{4} = 'V1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'XPROD1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD1';
        ielftype(ie) = iet_('ePROD1');
        vname = 'AW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TT1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VE',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XPF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.146;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPG';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.648;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPROD2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        vname = 'AW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'S1L';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XPL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.736;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.229;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPROD3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD1';
        ielftype(ie) = iet_('ePROD1');
        vname = 'AW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TT2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VE',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X2PF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.146;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'X2PG';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.648;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPROD4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        vname = 'AW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'S2L';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X2PL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.736;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'X2PM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.229;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPROD5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD1';
        ielftype(ie) = iet_('ePROD1');
        vname = 'AW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TT3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VE',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XQF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.539;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XQG';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.772;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1000.0;
        ename = 'XPROD6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        vname = 'AW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'S3L';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('VD',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XQL';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.33;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 100.0;
        ename = 'XQM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'ET3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PWR',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.498;
        [~,posep] = ismember('SC',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 100.0;
        ename = 'SMFE1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX7Y';
        ielftype(ie) = iet_('eX7Y');
        vname = 'MF1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'INW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'SMFE2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX5Y';
        ielftype(ie) = iet_('eX5Y');
        vname = 'MF2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'INW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'SMFE3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX3Y';
        ielftype(ie) = iet_('eX3Y');
        vname = 'MF3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'INW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TT1BT1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'TT1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TT2BT2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'TT2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TT3BT3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'TT3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XBIG11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBIG1';
        ielftype(ie) = iet_('eBIG1');
        vname = 'TT1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('TH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LL',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XBIG12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBIG1';
        ielftype(ie) = iet_('eBIG1');
        vname = 'TT2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('TH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LL',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'XBIG13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBIG1';
        ielftype(ie) = iet_('eBIG1');
        vname = 'TT3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'BT3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('TH',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PW3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('LL',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'MF3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,1.0e+4,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSUMM',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('STA1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5272.77;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 160.909;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XPG');
        pbm.grelw{ig}(posel) = 282.874;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.64570846;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPL');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 31.136196;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XPM');
        pbm.grelw{ig}(posel) = 12.092112;
        ig = ig_('STA2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5272.77;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X2PF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 160.909;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X2PG');
        pbm.grelw{ig}(posel) = 282.874;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.64570846;
        ig = ig_('STA1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X2PL');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 31.136196;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X2PM');
        pbm.grelw{ig}(posel) = 12.092112;
        ig = ig_('STA3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5272.77;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XQF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 181.806;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XQG');
        pbm.grelw{ig}(posel) = 232.57;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XPROD6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.49783215;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XQL');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.22424514;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XQM');
        pbm.grelw{ig}(posel) = 20.708238;
        ig = ig_('LAUN');
        pbm.grftype{ig} = 'gSUMM';
        ig = ig_('SMF1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('SMFE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SMF2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('SMFE2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SMF3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('SMFE3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI1A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT1BT1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI1B');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT1BT1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI2A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT2BT2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI2B');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT2BT2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI3A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT3BT3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SI3B');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TT3BT3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('GLGCON');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XBIG11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -32.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('XBIG12');
        pbm.grelw{ig}(posel) = -32.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('XBIG13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -32.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-MY-25-28';
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

    case 'eBIG1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LG = log(EV_(4));
        varargout{1} = (EV_(1)*EV_(2)*LG)/EV_(3);
        if(nargout>1)
            g_(4,1) = (EV_(1)*EV_(2))/(EV_(4)*EV_(3));
            g_(1,1) = (EV_(2)*LG)/EV_(3);
            g_(2,1) = (EV_(1)*LG)/EV_(3);
            g_(3,1) = -(EV_(1)*EV_(2)*LG)/(EV_(3)^2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(4,4) = -(EV_(1)*EV_(2))/(EV_(3)*EV_(4)^2);
                H_(4,1) = EV_(2)/(EV_(4)*EV_(3));
                H_(1,4) = H_(4,1);
                H_(4,2) = EV_(1)/(EV_(4)*EV_(3));
                H_(2,4) = H_(4,2);
                H_(4,3) = -(EV_(1)*EV_(2))/(EV_(3)^2*EV_(4));
                H_(3,4) = H_(4,3);
                H_(1,2) = LG/EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = -EV_(2)*LG/EV_(3)^2;
                H_(3,1) = H_(1,3);
                H_(2,3) = -EV_(1)*LG/EV_(3)^2;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*(EV_(1)*EV_(2)*LG)/(EV_(3)^3.0);
                varargout{3} = H_;
            end
        end

    case 'ePROD1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EA = 1.2781;
        VA0 = EV_(1)^EA;
        VA1 = EA*EV_(1)^(EA-1.0);
        VA2 = EA*(EA-1.0)*EV_(1)^(EA-2.0);
        EB = -0.1959;
        VB0 = EV_(2)^EB;
        VB1 = EB*EV_(2)^(EB-1.0);
        VB2 = EB*(EB-1.0)*EV_(2)^(EB-2.0);
        EC = 2.4242;
        VC0 = EV_(3)^EC;
        VC1 = EC*EV_(3)^(EC-1.0);
        VC2 = EC*(EC-1.0)*EV_(3)^(EC-2.0);
        ED = 0.38745;
        VD0 = EV_(4)^ED;
        VD1 = ED*EV_(4)^(ED-1.0);
        VD2 = ED*(ED-1.0)*EV_(4)^(ED-2.0);
        EE = 0.9904;
        VE0 = EV_(5)^EE;
        VE1 = EE*EV_(5)^(EE-1.0);
        VE2 = EE*(EE-1.0)*EV_(5)^(EE-2.0);
        varargout{1} = VA0*VB0*VC0*VD0*VE0;
        if(nargout>1)
            g_(1,1) = VA1*VB0*VC0*VD0*VE0;
            g_(2,1) = VA0*VB1*VC0*VD0*VE0;
            g_(3,1) = VA0*VB0*VC1*VD0*VE0;
            g_(4,1) = VA0*VB0*VC0*VD1*VE0;
            g_(5,1) = VA0*VB0*VC0*VD0*VE1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = VA2*VB0*VC0*VD0*VE0;
                H_(1,2) = VA1*VB1*VC0*VD0*VE0;
                H_(2,1) = H_(1,2);
                H_(1,3) = VA1*VB0*VC1*VD0*VE0;
                H_(3,1) = H_(1,3);
                H_(1,4) = VA1*VB0*VC0*VD1*VE0;
                H_(4,1) = H_(1,4);
                H_(1,5) = VA1*VB0*VC0*VD0*VE1;
                H_(5,1) = H_(1,5);
                H_(2,2) = VA0*VB2*VC0*VD0*VE0;
                H_(2,3) = VA0*VB1*VC1*VD0*VE0;
                H_(3,2) = H_(2,3);
                H_(2,4) = VA0*VB1*VC0*VD1*VE0;
                H_(4,2) = H_(2,4);
                H_(2,5) = VA0*VB1*VC0*VD0*VE1;
                H_(5,2) = H_(2,5);
                H_(3,3) = VA0*VB0*VC2*VD0*VE0;
                H_(3,4) = VA0*VB0*VC1*VD1*VE0;
                H_(4,3) = H_(3,4);
                H_(3,5) = VA0*VB0*VC1*VD0*VE1;
                H_(5,3) = H_(3,5);
                H_(4,4) = VA0*VB0*VC0*VD2*VE0;
                H_(4,5) = VA0*VB0*VC0*VD1*VE1;
                H_(5,4) = H_(4,5);
                H_(5,5) = VA0*VB0*VC0*VD0*VE2;
                varargout{3} = H_;
            end
        end

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EA = 0.3322;
        VA0 = EV_(1)^EA;
        VA1 = EA*EV_(1)^(EA-1.0);
        VA2 = EA*(EA-1.0)*EV_(1)^(EA-2.0);
        EB = -1.5935;
        VB0 = EV_(2)^EB;
        VB1 = EB*EV_(2)^(EB-1.0);
        VB2 = EB*(EB-1.0)*EV_(2)^(EB-2.0);
        EC = 0.2363;
        VC0 = EV_(3)^EC;
        VC1 = EC*EV_(3)^(EC-1.0);
        VC2 = EC*(EC-1.0)*EV_(3)^(EC-2.0);
        ED = 0.1079;
        VD0 = EV_(4)^ED;
        VD1 = ED*EV_(4)^(ED-1.0);
        VD2 = ED*(ED-1.0)*EV_(4)^(ED-2.0);
        varargout{1} = VA0*VB0*VC0*VD0;
        if(nargout>1)
            g_(1,1) = VA1*VB0*VC0*VD0;
            g_(2,1) = VA0*VB1*VC0*VD0;
            g_(3,1) = VA0*VB0*VC1*VD0;
            g_(4,1) = VA0*VB0*VC0*VD1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = VA2*VB0*VC0*VD0;
                H_(1,2) = VA1*VB1*VC0*VD0;
                H_(2,1) = H_(1,2);
                H_(1,3) = VA1*VB0*VC1*VD0;
                H_(3,1) = H_(1,3);
                H_(1,4) = VA1*VB0*VC0*VD1;
                H_(4,1) = H_(1,4);
                H_(2,2) = VA0*VB2*VC0*VD0;
                H_(2,3) = VA0*VB1*VC1*VD0;
                H_(3,2) = H_(2,3);
                H_(2,4) = VA0*VB1*VC0*VD1;
                H_(4,2) = H_(2,4);
                H_(3,3) = VA0*VB0*VC2*VD0;
                H_(3,4) = VA0*VB0*VC1*VD1;
                H_(4,3) = H_(3,4);
                H_(4,4) = VA0*VB0*VC0*VD2;
                varargout{3} = H_;
            end
        end

    case 'ePOWER'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SCPWR = pbm.elpar{iel_}(1)/(pbm.elpar{iel_}(2)^pbm.elpar{iel_}(1));
        varargout{1} = (EV_(1)/pbm.elpar{iel_}(2))^pbm.elpar{iel_}(1);
        if(nargout>1)
            g_(1,1) = SCPWR*EV_(1)^(pbm.elpar{iel_}(1)-1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = SCPWR*(pbm.elpar{iel_}(1)-1.0)*EV_(1)^(pbm.elpar{iel_}(1)-2.0);
                varargout{3} = H_;
            end
        end

    case 'eX7Y'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,8);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)+1;
        U_(2,7) = U_(2,7)+1;
        U_(2,8) = U_(2,8)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eX5Y'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,6);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eX3Y'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

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
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSUMM'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^0.460;
        if(nargout>1)
            g_ = 0.460*GVAR_^(-0.540);
            varargout{2} = g_;
            if(nargout>2)
                H_ = -0.2484*GVAR_^(-1.540);
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

