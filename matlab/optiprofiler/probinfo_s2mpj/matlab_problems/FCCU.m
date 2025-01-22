function varargout = FCCU(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : FCCU
%    *********
% 
%    A simple data reconciliation for a fluid catalytic cracker.
% 
%                      +--------------------------+
%                      | FCCU data reconciliation |
%     +-----------+    +--------------------------+
%  1) | Flowsheet |                      Off_gas                      |------->
%     +-----------+                     |----------->                 | Propane
%                              MF_ohd   |                  |-------->F7
%                              |------>F4<-------|         |          |
%                              |        |        |         |DC3_feed  | Butane
%  Feed      Effluent          ^        |        |         |          |------->
%  ------->F1-------->F2......>|        |---------------->F5                  
%           ^                  v       DC4_feed  |         |DC4_btms  |------->
%           |                  |                 |         |          |  LCN
%           |                  |                 |<-------F6-------->F8
%           |                  | HCN              Lean_oil   C8spl_fd |  MCN
%           |                  |-------->                             |------->
%           |                  | LCO          
%           |                  |-------->              
%           |                  | HCO          
%           |                  |-------->              
%           |                  | MF_btms     
%           |                  v 
%           |<----------------F3-------->
%              Dec_recy         Decant
%     +--------------------+
%  2) | Objective function |
%     +--------------------+
%    Obj = sum 1->i [W_i*(C_flow_i - M_flow_i)**2]
%                           
%     Where: W_i       Weight on term i of objective function
%            C_flow_i  Computed flow i (a variable for this problem)
%            M_flow_i  Measrued flow i (a constant for this problem)
%     +-------------+
%  3) | Constraints |
%     +-------------+
%     These represent the linear mass balances around each
%     node, where a node (Fx) represents a single unit operation
%     in a fluid catalytics cracker.
%     +---------------+
%  4) | Initial point |
%     +---------------+
%     Feed       1.0
%     Effluent   1.0
%     MF_ohd     1.0
%     HCN        1.0
%     LCO        1.0
%     HCO        1.0
%     MF_btms    1.0
%     Decant     1.0
%     Dec_recy   1.0
%     Off_gas    1.0
%     DC4_feed   1.0
%     DC3_feed   1.0
%     DC4_btms   1.0
%     Lean_oil   1.0
%     Propane    1.0
%     Butane     1.0
%     C8spl_fd   1.0
%     LCN        1.0
%     MCN        1.0
%     Obj        7.36259000271320D+03
%     +------------------+
%  5) | Optimal solution |
%     +------------------+
%     Feed       3.11639D+01
%     Effluent   3.53528D+01
%     MF_ohd     1.94669D+01
%     HCN        2.94255D+00
%     LCO        4.94255D+00
%     HCO        3.44255D+00
%     MF_btms    4.55828D+00
%     Decant     3.69371D-01
%     Dec_recy   4.18891D+00
%     Off_gas    2.56075D+00
%     DC4_feed   2.41207D+01
%     DC3_feed   5.15601D+00
%     DC4_btms   1.89647D+01
%     Lean_oil   7.21458D+00
%     Propane    2.42801D+00
%     Butane     2.72801D+00
%     C8spl_fd   1.17501D+01
%     LCN        5.87506D+00
%     MCN        5.87506D+00
%     Obj        1.11491D+01
%     +-----------------------------------------------+
%  6) | SPEC.SPC (remove 1st * of every line to use). |
%     +-----------------------------------------------+
% BEGIN
% * maximizer-sought
% *  check-derivatives
%   ignore-derivative-bugs
% * use-scalings
% * print-scalings
%   finite-difference-gradients
% *  exact-second-derivatives-used
% * bfgs-approximate-second-derivatives-used
% * sr1-approximate-second-derivatives-used
%   bandsolver-preconditioned-cg-solver-used   5
% * diagonal-preconditioned-cg-solver-used
% * gill-murray-ponceleon-saunders-preconditioned-cg-solver-used
% * schnabel-eskow-preconditioned-cg-solver-used
% * munksgaards-preconditioned-cg-solver-used
%   exact-cauchy-point-required
% * inexact-cauchy-point-required
% * solve-bqp-accurately
% * two-norm-trust-region
% * gradient-tolerance    1.0D-5
% * constraint-tolerance  1.0D-5
%   trust-region-radius   1.0D+0
%   maximum-number-of-iterations   1000
%   print-level                    1
%   start-printing-at-iteration    0
%   stop-printing-at-iteration     1000
% END
% 
%    Source:
%    W. J. Korchinski, Profimatics, Inc,
%    325 Rolling Oaks Drive, Thousand Oaks, California, USA 91361-1200
%    Telephone: 1-805 496 6661, Fax: 1-805 373 5108
% 
%    SIF input: W. Korchinski, Spring 1993.
% 
%    classification = 'C-CSLR2-MN-19-8'
% 
% ***************************************************************
%  PROBLEM SPECIFICATION BEGINS HERE.
%  **********************************
%  **********************************           
% *************************************
%  Define objective function weights. *
% *************************************
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'FCCU';

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
        v_('W1') = 0.2;
        v_('W2') = 1.0;
        v_('W3') = 1.0;
        v_('W4') = 0.33333333;
        v_('W5') = 0.33333333;
        v_('W6') = 0.33333333;
        v_('W7') = 1.0;
        v_('W8') = 1.0;
        v_('W9') = 1.0;
        v_('W10') = 1.0;
        v_('W11') = 1.0;
        v_('W12') = 1.0;
        v_('W13') = 1.0;
        v_('W14') = 1.0;
        v_('W15') = 0.33333333;
        v_('W16') = 0.33333333;
        v_('W17') = 1.0;
        v_('W18') = 0.33333333;
        v_('W19') = 0.33333333;
        v_('M1') = 31.0;
        v_('M2') = 36.0;
        v_('M3') = 20.0;
        v_('M4') = 3.0;
        v_('M5') = 5.0;
        v_('M6') = 3.5;
        v_('M7') = 4.2;
        v_('M8') = 0.9;
        v_('M9') = 3.9;
        v_('M10') = 2.2;
        v_('M11') = 22.8;
        v_('M12') = 6.8;
        v_('M13') = 19.0;
        v_('M14') = 8.5;
        v_('M15') = 2.2;
        v_('M16') = 2.5;
        v_('M17') = 10.8;
        v_('M18') = 6.5;
        v_('M19') = 6.5;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','Feed',ix_);
        pb.xnames{iv} = 'Feed';
        [iv,ix_] = s2mpjlib('ii','Effluent',ix_);
        pb.xnames{iv} = 'Effluent';
        [iv,ix_] = s2mpjlib('ii','MFuohd',ix_);
        pb.xnames{iv} = 'MFuohd';
        [iv,ix_] = s2mpjlib('ii','HCN',ix_);
        pb.xnames{iv} = 'HCN';
        [iv,ix_] = s2mpjlib('ii','LCO',ix_);
        pb.xnames{iv} = 'LCO';
        [iv,ix_] = s2mpjlib('ii','HCO',ix_);
        pb.xnames{iv} = 'HCO';
        [iv,ix_] = s2mpjlib('ii','MFubtms',ix_);
        pb.xnames{iv} = 'MFubtms';
        [iv,ix_] = s2mpjlib('ii','Decant',ix_);
        pb.xnames{iv} = 'Decant';
        [iv,ix_] = s2mpjlib('ii','Decurecy',ix_);
        pb.xnames{iv} = 'Decurecy';
        [iv,ix_] = s2mpjlib('ii','Offugas',ix_);
        pb.xnames{iv} = 'Offugas';
        [iv,ix_] = s2mpjlib('ii','DC4ufeed',ix_);
        pb.xnames{iv} = 'DC4ufeed';
        [iv,ix_] = s2mpjlib('ii','DC3ufeed',ix_);
        pb.xnames{iv} = 'DC3ufeed';
        [iv,ix_] = s2mpjlib('ii','DC4ubtms',ix_);
        pb.xnames{iv} = 'DC4ubtms';
        [iv,ix_] = s2mpjlib('ii','Leanuoil',ix_);
        pb.xnames{iv} = 'Leanuoil';
        [iv,ix_] = s2mpjlib('ii','Propane',ix_);
        pb.xnames{iv} = 'Propane';
        [iv,ix_] = s2mpjlib('ii','Butane',ix_);
        pb.xnames{iv} = 'Butane';
        [iv,ix_] = s2mpjlib('ii','C8splufd',ix_);
        pb.xnames{iv} = 'C8splufd';
        [iv,ix_] = s2mpjlib('ii','LCN',ix_);
        pb.xnames{iv} = 'LCN';
        [iv,ix_] = s2mpjlib('ii','MCN',ix_);
        pb.xnames{iv} = 'MCN';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','F1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Feed');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Decurecy');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Effluent');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Effluent');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFuohd');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('HCN');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LCO');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('HCO');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFubtms');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFubtms');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Decant');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Decurecy');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFuohd');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Leanuoil');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Offugas');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ufeed');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ufeed');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC3ufeed');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ubtms');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ubtms');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Leanuoil');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('C8splufd');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC3ufeed');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Propane');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Butane');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','F8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'F8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('C8splufd');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LCN');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MCN');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','Obj1',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Feed');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W1');
        [ig,ig_] = s2mpjlib('ii','Obj2',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Effluent');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W2');
        [ig,ig_] = s2mpjlib('ii','Obj3',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFuohd');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W3');
        [ig,ig_] = s2mpjlib('ii','Obj4',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('HCN');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W4');
        [ig,ig_] = s2mpjlib('ii','Obj5',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LCO');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W5');
        [ig,ig_] = s2mpjlib('ii','Obj6',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('HCO');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W6');
        [ig,ig_] = s2mpjlib('ii','Obj7',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MFubtms');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W7');
        [ig,ig_] = s2mpjlib('ii','Obj8',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Decant');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W8');
        [ig,ig_] = s2mpjlib('ii','Obj9',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Decurecy');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W9');
        [ig,ig_] = s2mpjlib('ii','Obj10',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Offugas');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W10');
        [ig,ig_] = s2mpjlib('ii','Obj11',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ufeed');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W11');
        [ig,ig_] = s2mpjlib('ii','Obj12',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC3ufeed');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W12');
        [ig,ig_] = s2mpjlib('ii','Obj13',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DC4ubtms');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W13');
        [ig,ig_] = s2mpjlib('ii','Obj14',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Leanuoil');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W14');
        [ig,ig_] = s2mpjlib('ii','Obj15',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Propane');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W15');
        [ig,ig_] = s2mpjlib('ii','Obj16',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Butane');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W16');
        [ig,ig_] = s2mpjlib('ii','Obj17',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('C8splufd');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W17');
        [ig,ig_] = s2mpjlib('ii','Obj18',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LCN');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W18');
        [ig,ig_] = s2mpjlib('ii','Obj19',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('MCN');
        valA(end+1) = 1.0;
        pbm.gscale(ig,1) = v_('W19');
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
        pbm.gconst(ig_('Obj1')) = v_('M1');
        pbm.gconst(ig_('Obj2')) = v_('M2');
        pbm.gconst(ig_('Obj3')) = v_('M3');
        pbm.gconst(ig_('Obj4')) = v_('M4');
        pbm.gconst(ig_('Obj5')) = v_('M5');
        pbm.gconst(ig_('Obj6')) = v_('M6');
        pbm.gconst(ig_('Obj7')) = v_('M7');
        pbm.gconst(ig_('Obj8')) = v_('M8');
        pbm.gconst(ig_('Obj9')) = v_('M9');
        pbm.gconst(ig_('Obj10')) = v_('M10');
        pbm.gconst(ig_('Obj11')) = v_('M11');
        pbm.gconst(ig_('Obj12')) = v_('M12');
        pbm.gconst(ig_('Obj13')) = v_('M13');
        pbm.gconst(ig_('Obj14')) = v_('M14');
        pbm.gconst(ig_('Obj15')) = v_('M15');
        pbm.gconst(ig_('Obj16')) = v_('M16');
        pbm.gconst(ig_('Obj17')) = v_('M17');
        pbm.gconst(ig_('Obj18')) = v_('M18');
        pbm.gconst(ig_('Obj19')) = v_('M19');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('Feed'),1) = 1.0;
        pb.x0(ix_('Effluent'),1) = 1.0;
        pb.x0(ix_('MFuohd'),1) = 1.0;
        pb.x0(ix_('HCN'),1) = 1.0;
        pb.x0(ix_('LCO'),1) = 1.0;
        pb.x0(ix_('HCO'),1) = 1.0;
        pb.x0(ix_('MFubtms'),1) = 1.0;
        pb.x0(ix_('Decant'),1) = 1.0;
        pb.x0(ix_('Decurecy'),1) = 1.0;
        pb.x0(ix_('Offugas'),1) = 1.0;
        pb.x0(ix_('DC4ufeed'),1) = 1.0;
        pb.x0(ix_('DC3ufeed'),1) = 1.0;
        pb.x0(ix_('DC4ubtms'),1) = 1.0;
        pb.x0(ix_('Leanuoil'),1) = 1.0;
        pb.x0(ix_('Propane'),1) = 1.0;
        pb.x0(ix_('Butane'),1) = 1.0;
        pb.x0(ix_('C8splufd'),1) = 1.0;
        pb.x0(ix_('LCN'),1) = 1.0;
        pb.x0(ix_('MCN'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSQUARE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('Obj1');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj2');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj3');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj4');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj5');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj6');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj7');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj8');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj9');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj10');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj11');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj12');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj13');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj14');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj15');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj16');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj17');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj18');
        pbm.grftype{ig} = 'gSQUARE';
        ig = ig_('Obj19');
        pbm.grftype{ig} = 'gSQUARE';
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
        pb.pbclass = 'C-CSLR2-MN-19-8';
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

