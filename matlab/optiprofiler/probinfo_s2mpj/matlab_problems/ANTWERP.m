function varargout = ANTWERP(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ANTWERP
%    *********
% 
%    This problem arises in the determination of a synthetic population for
%    Belgian municipalities. The question is to estimate the distribution in 
%    Antwerp that households of the following types:
%       type F (a couple + 1 to 5 children + 0 to 2 additional adults)
%       type W (a woman  + 1 to 3 children + 0 to 2 additional adults)
%       type M (a man    + 1 to 3 children + 0 to 2 additional adults).
%    The data consists in 
%       - the number of individuals in households with 3 to 8 members, 
%       - the number of F, W and N households according to their number of
%         children
%       - and the total number of adults and children.
%    If we define the variables by
%    p1F: probability for a F_household to have 1 child,
%    p2F: probability for a F_household to have 2 children,
%    p3F: probability for a F_household to have 3 children,
%    p4F: probability for a F_household to have 4 children,
%    p5F: probability for a F_household to have 5 children,
%    p1W: probability for a W_household to have 1 child,
%    p2W: probability for a W_household to have 2 children,
%    p3W: probability for a W_household to have 3 children,
%    p1M: probability for a M_household to have 1 child,
%    p2M: probability for a M_household to have 2 children,
%    p3M: probability for a M_household to have 3 children,
%    q0F: probability for a F_household to have 1 additional adult,
%    q1F: probability for a F_household to have 2 additional adults,
%    q2F: probability for a F_household to have 3 additional adults,
%    q0W: probability for a W_household to have 1 additional adult,
%    q1W: probability for a W_household to have 2 additional adults,
%    q2W: probability for a W_household to have 3 additional adults,
%    q0M: probability for a M_household to have 1 additional adult,
%    q1M: probability for a M_household to have 2 additional adults,
%    q2M: probability for a M_household to have 3 additional adults,
%    nF : number of F-households,
%    nW : number of W-households,
%    nM : number of N-households,
%    nC2: number of individuals considered as children in age class 2
%    nC3: number of individuals considered as children in age class 3,
%    nA2: number of individuals considered as adults in age class 2,
%    nA3: number of individuals considered as adults in age class 3,
%    the derived predictions for the observed quantities are then given by
%    1) prediction of the number of individuals in household of size 3:
%       p1F*q0F*nF + (p1W*q1W+p2W*q0W)*nW + (p1M*q1M+p2M*q0M)*nM = M3
%    2) prediction of the number of individuals in household of size 4:
%       (p2F*q0F+p1F*q1F)*nF + (p1W*q2W+p2W*q1W+p3W*q0W)*nW 
%                            + (p1M*q2M+p2M*q1M+p3M*q0M)*nM = M4
%    3) prediction of the number of individuals in household of size 5:
%       (p3F*q0F+p2F*q1F+p1F*q2F)*nF + (p2W*q2W+p3W*q1W)*nW 
%                                    + (p2M*q2M+p3M*q1M)*nM = M5
%    4) prediction of the number of individuals in household of size 6:
%       (p4F*q0F+p3F*q1F+p2F*q2F)*nF + p3W*q2W*nW + p3M*q2M*nM = M6
%    5) prediction of the number of individuals in household of size 7:
%       (p5F*q0F+p4F*q1F+p3F*q2F)*nF = M7
%    6) prediction of the number of individuals in household of size 8:
%       (p5F*q1F+p4F*q2F)*nF = M8
%    7) prediction of the number of F-households with 1 child
%       p1F*nF*(M1F+M2F+M3F) = M1F
%    8) prediction of the number of F-households with 2 children
%       p2F*nF*(M1F+M2F+M3F) = M2F
%    9) prediction of the number of F-households with 3 children or more
%       (p3F+p4F+p5F)*nF*(M1F+M2F+M3F) = M3F
%    10) prediction of the number of W-households with 1 child
%        p1W*nW*(M1W+M2W+M3W) = M1W
%    11) prediction of the number of W-households with 2 children
%        p2W*nW*(M1W+M2W+M3W) = M2W
%    12) prediction of the number of W-households with 3 children
%        p3W*nW*(M1W+M2W+M3W) = M3W
%    13) prediction of the number of M-households with 1 child
%        p1M*nM*(M1M+M2M+M3M) = M1M
%    14) prediction of the number of M-households with 2 children
%        p2M*nM*(M1M+M2M+M3M) = M2M
%  
%    14) prediction of the number of M-households with 3 children
%        p3M*nM*(M1M+M2M+M3M) = M3M
%    15) prediction of the number of children
%        (p1F+2*p2F+3*p3F+4*P4F+5*p5F)*NF + (p1W+2*p2W+3*p3W)*nW
%                  + (p1W+2*p2W+3*p3W)*nW - nC2 -nC3 = N0 + N1
%    16) prediction of the total number of adults
%        (2*q0F+3*q1F+4*q2F)*nF + (q0W+2*q1W+3*q2W)*nW 
%                  + (q0W+2*q1W+3*q2W)*nW - nA2 -nA3 = N4
%    17) composition of age class 2
%        nC2 + nA2 = N2
%    18) composition of age class 3
%        nC3 + nA3 = N3
%    19) prediction of the total number of individuals in F-households
%        (p1F+2*p2F+3*p3F+4*P4F+5*p5F)*NF + (2*q0F+3*q1F+4*q2F)*nF = NINF
%    20) the piF are probabilities and sum up to 1
%        p1F + p2F + p3F + p4F + p5F = 1
%    21) the qiF are probabilities and sum up to 1
%        
%        q0F + q1F + q2F = 1
%    22) the piW are probabilities and sum up to 1
%        p1W + p2W + p3W = 1
%    23) the qiW are probabilities and sum up to 1
%        
%        q0W + q1W + q2W = 1
%    24) the piM are probabilities and sum up to 1
%        p1M + p2M + p3M = 1
%    25) the qiM are probabilities and sum up to 1
%        
%        q0M + q1M + q2M = 1
%    In addition, the following inequalities are imposed
%    26) the fraction of children in age class 2 exceeds that in age class 3
%         nC2/N2 >= nC3/N3
%    27) there are more adults in age class 2 than children
%         nA2 >= nC2
%    28) there are more adults in age class 3 than children
%         nA3 >= nC3
%    and the bounds on the variables are
%        0 <= piF <= 1      ( i = 1, 2, 3, 4, 5 )
%  
%        0 <= qiF <= 1      ( i = 0, 1, 2 ) 
%        0 <= piW <= 1      ( i = 1, 2, 3 )
%  
%        0 <= qiW <= 1      ( i = 0, 1, 2 ) 
%        0 <= piM <= 1      ( i = 1, 2, 3 )
%  
%        0 <= qiM <= 1      ( i = 0, 1, 2 ) 
%        nF >= 0,  nW >= 0,  nM >= 0
%        0 <= nC2 <= N2,   0 <= nA2 <= N2
%        0 <= nC3 <= N3,   0 <= nA3 <= N3
% 
%    The problem is solved as a linearly/bound  constrained nonlinear least-squares
%    problem in 27 variables.  In the least-squares formulation,
%    each equation is scaled in a proportion inverse to its right-hand side. 
% 
%    The problem appears to be very ill-conditioned.
%    SIF input: Ph. Toint, Apr 2006.
% 
%    classification = 'C-CSLR2-RN-27-8-0-3-24-0-2-0-8-0-0-0'
% 
%    Problem initial data
% 
%    Number of households according to their sizes
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ANTWERP';

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
        v_('M3') = 23844.0;
        v_('M4') = 16323.0;
        v_('M5') = 6613.0;
        v_('M6') = 2535.0;
        v_('M7') = 1109.0;
        v_('M8') = 1667.0;
        v_('M1F') = 16405.0;
        v_('M2F') = 13647.0;
        v_('M3F') = 9895.0;
        v_('M1M') = 4041.0;
        v_('M2M') = 1634.0;
        v_('M3M') = 637.0;
        v_('M1W') = 10966.0;
        v_('M2W') = 4566.0;
        v_('M3W') = 1921.0;
        v_('N0') = 15866.0;
        v_('N1') = 59832.0;
        v_('N2') = 61929.0;
        v_('N3') = 32321.0;
        v_('N4') = 73650.0;
        v_('NINF') = 180055.0;
        v_('NINN') = 47677.0;
        v_('TMP') = v_('M1F')+v_('M2F');
        v_('SNF') = v_('TMP')+v_('M3F');
        v_('TMP') = v_('M1M')+v_('M2M');
        v_('SNM') = v_('TMP')+v_('M3M');
        v_('TMP') = v_('M1W')+v_('M2W');
        v_('SNW') = v_('TMP')+v_('M3W');
        v_('N2/2') = 0.5*v_('N2');
        v_('N3/2') = 0.5*v_('N3');
        v_('N23/2') = v_('N2/2')+v_('N3/2');
        v_('N01') = v_('N0')+v_('N1');
        v_('N0123') = v_('N01')+v_('N23/2');
        v_('N234') = v_('N4')+v_('N23/2');
        v_('SP1F') = v_('M1F')/v_('SNF');
        v_('SP2F') = v_('M2F')/v_('SNF');
        v_('SP1W') = v_('M1W')/v_('SNW');
        v_('SP2W') = v_('M2W')/v_('SNW');
        v_('SP3W') = v_('M3W')/v_('SNW');
        v_('SP1M') = v_('M1M')/v_('SNM');
        v_('SP2M') = v_('M2M')/v_('SNM');
        v_('SP3M') = v_('M3M')/v_('SNM');
        v_('1/N3') = 1.0/v_('N3');
        v_('-1/N2') = -1.0/v_('N2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','P1F',ix_);
        pb.xnames{iv} = 'P1F';
        [iv,ix_] = s2mpjlib('ii','P2F',ix_);
        pb.xnames{iv} = 'P2F';
        [iv,ix_] = s2mpjlib('ii','P3F',ix_);
        pb.xnames{iv} = 'P3F';
        [iv,ix_] = s2mpjlib('ii','P4F',ix_);
        pb.xnames{iv} = 'P4F';
        [iv,ix_] = s2mpjlib('ii','P5F',ix_);
        pb.xnames{iv} = 'P5F';
        [iv,ix_] = s2mpjlib('ii','P1W',ix_);
        pb.xnames{iv} = 'P1W';
        [iv,ix_] = s2mpjlib('ii','P2W',ix_);
        pb.xnames{iv} = 'P2W';
        [iv,ix_] = s2mpjlib('ii','P3W',ix_);
        pb.xnames{iv} = 'P3W';
        [iv,ix_] = s2mpjlib('ii','P1M',ix_);
        pb.xnames{iv} = 'P1M';
        [iv,ix_] = s2mpjlib('ii','P2M',ix_);
        pb.xnames{iv} = 'P2M';
        [iv,ix_] = s2mpjlib('ii','P3M',ix_);
        pb.xnames{iv} = 'P3M';
        [iv,ix_] = s2mpjlib('ii','Q0F',ix_);
        pb.xnames{iv} = 'Q0F';
        [iv,ix_] = s2mpjlib('ii','Q1F',ix_);
        pb.xnames{iv} = 'Q1F';
        [iv,ix_] = s2mpjlib('ii','Q2F',ix_);
        pb.xnames{iv} = 'Q2F';
        [iv,ix_] = s2mpjlib('ii','Q0W',ix_);
        pb.xnames{iv} = 'Q0W';
        [iv,ix_] = s2mpjlib('ii','Q1W',ix_);
        pb.xnames{iv} = 'Q1W';
        [iv,ix_] = s2mpjlib('ii','Q2W',ix_);
        pb.xnames{iv} = 'Q2W';
        [iv,ix_] = s2mpjlib('ii','Q0M',ix_);
        pb.xnames{iv} = 'Q0M';
        [iv,ix_] = s2mpjlib('ii','Q1M',ix_);
        pb.xnames{iv} = 'Q1M';
        [iv,ix_] = s2mpjlib('ii','Q2M',ix_);
        pb.xnames{iv} = 'Q2M';
        [iv,ix_] = s2mpjlib('ii','NF',ix_);
        pb.xnames{iv} = 'NF';
        [iv,ix_] = s2mpjlib('ii','NW',ix_);
        pb.xnames{iv} = 'NW';
        [iv,ix_] = s2mpjlib('ii','NM',ix_);
        pb.xnames{iv} = 'NM';
        [iv,ix_] = s2mpjlib('ii','NC2',ix_);
        pb.xnames{iv} = 'NC2';
        [iv,ix_] = s2mpjlib('ii','NA2',ix_);
        pb.xnames{iv} = 'NA2';
        [iv,ix_] = s2mpjlib('ii','NC3',ix_);
        pb.xnames{iv} = 'NC3';
        [iv,ix_] = s2mpjlib('ii','NA3',ix_);
        pb.xnames{iv} = 'NA3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','HSZ3',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M3');
        [ig,ig_] = s2mpjlib('ii','HSZ4',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M4');
        [ig,ig_] = s2mpjlib('ii','HSZ5',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M5');
        [ig,ig_] = s2mpjlib('ii','HSZ6',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M6');
        [ig,ig_] = s2mpjlib('ii','HSZ7',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M7');
        [ig,ig_] = s2mpjlib('ii','HSZ8',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M8');
        [ig,ig_] = s2mpjlib('ii','HST1F',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M1F');
        [ig,ig_] = s2mpjlib('ii','HST2F',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M2F');
        [ig,ig_] = s2mpjlib('ii','HST3F',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M3F');
        [ig,ig_] = s2mpjlib('ii','HST1W',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M1W');
        [ig,ig_] = s2mpjlib('ii','HST2W',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M2W');
        [ig,ig_] = s2mpjlib('ii','HST3W',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M3W');
        [ig,ig_] = s2mpjlib('ii','HST1M',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M1M');
        [ig,ig_] = s2mpjlib('ii','HST2M',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M2M');
        [ig,ig_] = s2mpjlib('ii','HST3M',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('M3M');
        v_('WHCH') = 100.0*v_('N0123');
        [ig,ig_] = s2mpjlib('ii','HCH',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('WHCH');
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC3');
        valA(end+1) = -1.0;
        v_('WHAD') = 100.0*v_('N234');
        [ig,ig_] = s2mpjlib('ii','HAD',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('WHAD');
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','AGE2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'AGE2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA2');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','AGE3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'AGE3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','HINF',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','HINN',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','PSF',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'PSF';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P1F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P2F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P3F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P4F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P5F');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','PSW',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'PSW';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P1W');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P2W');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P3W');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','PSM',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'PSM';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P1M');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P2M');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P3M');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','QSF',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'QSF';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q0F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q1F');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q2F');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','QSM',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'QSM';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q0M');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q1M');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q2M');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','QSW',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'QSW';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q0W');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q1W');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q2W');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','INEQ2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'INEQ2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','INEQ3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'INEQ3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NC3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NA3');
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
        pbm.gconst(ig_('HSZ3')) = v_('M3');
        pbm.gconst(ig_('HSZ4')) = v_('M4');
        pbm.gconst(ig_('HSZ5')) = v_('M5');
        pbm.gconst(ig_('HSZ6')) = v_('M6');
        pbm.gconst(ig_('HSZ7')) = v_('M7');
        pbm.gconst(ig_('HSZ8')) = v_('M8');
        pbm.gconst(ig_('HST1F')) = v_('M1F');
        pbm.gconst(ig_('HST2F')) = v_('M2F');
        pbm.gconst(ig_('HST3F')) = v_('M3F');
        pbm.gconst(ig_('HST1W')) = v_('M1W');
        pbm.gconst(ig_('HST2W')) = v_('M2W');
        pbm.gconst(ig_('HST3W')) = v_('M3W');
        pbm.gconst(ig_('HST1M')) = v_('M1M');
        pbm.gconst(ig_('HST2M')) = v_('M2M');
        pbm.gconst(ig_('HST3M')) = v_('M3M');
        pbm.gconst(ig_('HCH')) = v_('N01');
        pbm.gconst(ig_('HAD')) = v_('N4');
        pbm.gconst(ig_('HINF')) = v_('NINF');
        pbm.gconst(ig_('HINN')) = v_('NINN');
        pbm.gconst(ig_('AGE2')) = v_('N2');
        pbm.gconst(ig_('AGE3')) = v_('N3');
        pbm.gconst(ig_('PSF')) = 1.0;
        pbm.gconst(ig_('PSM')) = 1.0;
        pbm.gconst(ig_('PSW')) = 1.0;
        pbm.gconst(ig_('QSF')) = 1.0;
        pbm.gconst(ig_('QSM')) = 1.0;
        pbm.gconst(ig_('QSW')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        pb.xupper(ix_('NF')) = +Inf;
        pb.xupper(ix_('NM')) = +Inf;
        pb.xupper(ix_('NW')) = +Inf;
        pb.xlower(ix_('NC2'),1) = 0.0;
        pb.xupper(ix_('NC2')) = v_('N2');
        pb.xlower(ix_('NA2'),1) = 0.0;
        pb.xupper(ix_('NA2')) = v_('N2');
        pb.xlower(ix_('NC3'),1) = 0.0;
        pb.xupper(ix_('NC3')) = v_('N3');
        pb.xlower(ix_('NA3'),1) = 0.0;
        pb.xupper(ix_('NA3')) = v_('N3');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'P1F'))
            pb.x0(ix_('P1F'),1) = v_('SP1F');
        else
            pb.y0(find(pbm.congrps==ig_('P1F')),1) = v_('SP1F');
        end
        if(isKey(ix_,'P2F'))
            pb.x0(ix_('P2F'),1) = v_('SP2F');
        else
            pb.y0(find(pbm.congrps==ig_('P2F')),1) = v_('SP2F');
        end
        if(isKey(ix_,'P3F'))
            pb.x0(ix_('P3F'),1) = 0.15;
        else
            pb.y0(find(pbm.congrps==ig_('P3F')),1) = 0.15;
        end
        if(isKey(ix_,'P4F'))
            pb.x0(ix_('P4F'),1) = 0.10;
        else
            pb.y0(find(pbm.congrps==ig_('P4F')),1) = 0.10;
        end
        if(isKey(ix_,'P5F'))
            pb.x0(ix_('P5F'),1) = 0.05;
        else
            pb.y0(find(pbm.congrps==ig_('P5F')),1) = 0.05;
        end
        if(isKey(ix_,'P1W'))
            pb.x0(ix_('P1W'),1) = v_('SP1W');
        else
            pb.y0(find(pbm.congrps==ig_('P1W')),1) = v_('SP1W');
        end
        if(isKey(ix_,'P2W'))
            pb.x0(ix_('P2W'),1) = v_('SP2W');
        else
            pb.y0(find(pbm.congrps==ig_('P2W')),1) = v_('SP2W');
        end
        if(isKey(ix_,'P3W'))
            pb.x0(ix_('P3W'),1) = v_('SP3W');
        else
            pb.y0(find(pbm.congrps==ig_('P3W')),1) = v_('SP3W');
        end
        if(isKey(ix_,'P1M'))
            pb.x0(ix_('P1M'),1) = v_('SP1M');
        else
            pb.y0(find(pbm.congrps==ig_('P1M')),1) = v_('SP1M');
        end
        if(isKey(ix_,'P2M'))
            pb.x0(ix_('P2M'),1) = v_('SP2M');
        else
            pb.y0(find(pbm.congrps==ig_('P2M')),1) = v_('SP2M');
        end
        if(isKey(ix_,'P3M'))
            pb.x0(ix_('P3M'),1) = v_('SP3M');
        else
            pb.y0(find(pbm.congrps==ig_('P3M')),1) = v_('SP3M');
        end
        if(isKey(ix_,'Q0F'))
            pb.x0(ix_('Q0F'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('Q0F')),1) = 0.6;
        end
        if(isKey(ix_,'Q1F'))
            pb.x0(ix_('Q1F'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('Q1F')),1) = 0.3;
        end
        if(isKey(ix_,'Q2F'))
            pb.x0(ix_('Q2F'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('Q2F')),1) = 0.1;
        end
        if(isKey(ix_,'Q0M'))
            pb.x0(ix_('Q0M'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('Q0M')),1) = 0.6;
        end
        if(isKey(ix_,'Q1M'))
            pb.x0(ix_('Q1M'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('Q1M')),1) = 0.3;
        end
        if(isKey(ix_,'Q2M'))
            pb.x0(ix_('Q2M'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('Q2M')),1) = 0.1;
        end
        if(isKey(ix_,'Q0W'))
            pb.x0(ix_('Q0W'),1) = 0.6;
        else
            pb.y0(find(pbm.congrps==ig_('Q0W')),1) = 0.6;
        end
        if(isKey(ix_,'Q1W'))
            pb.x0(ix_('Q1W'),1) = 0.3;
        else
            pb.y0(find(pbm.congrps==ig_('Q1W')),1) = 0.3;
        end
        if(isKey(ix_,'Q2W'))
            pb.x0(ix_('Q2W'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('Q2W')),1) = 0.1;
        end
        if(isKey(ix_,'NF'))
            pb.x0(ix_('NF'),1) = v_('SNF');
        else
            pb.y0(find(pbm.congrps==ig_('NF')),1) = v_('SNF');
        end
        if(isKey(ix_,'NW'))
            pb.x0(ix_('NW'),1) = v_('SNW');
        else
            pb.y0(find(pbm.congrps==ig_('NW')),1) = v_('SNW');
        end
        if(isKey(ix_,'NM'))
            pb.x0(ix_('NM'),1) = v_('SNM');
        else
            pb.y0(find(pbm.congrps==ig_('NM')),1) = v_('SNM');
        end
        if(isKey(ix_,'NC2'))
            pb.x0(ix_('NC2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('NC2')),1) = 0.0;
        end
        if(isKey(ix_,'NC3'))
            pb.x0(ix_('NC3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('NC3')),1) = 0.0;
        end
        if(isKey(ix_,'NA2'))
            pb.x0(ix_('NA2'),1) = v_('N2');
        else
            pb.y0(find(pbm.congrps==ig_('NA2')),1) = v_('N2');
        end
        if(isKey(ix_,'NA3'))
            pb.x0(ix_('NA3'),1) = v_('N3');
        else
            pb.y0(find(pbm.congrps==ig_('NA3')),1) = v_('N3');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'en3PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'P1FQ0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1WQ1WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2WQ0WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1MQ1MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2MQ0MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2FQ0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1FQ1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1WQ2WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2WQ1WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3WQ0WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1MQ2MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2MQ1MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3MQ0MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3FQ0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2FQ1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1FQ2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2WQ2WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3WQ1WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2MQ2MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3MQ1MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P4FQ0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P4F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3FQ1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2FQ2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3WQ2WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3MQ2MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P5FQ0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P5F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P4FQ1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P4F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3FQ2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P3F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P5FQ1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P5F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P4FQ2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'P4F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Q2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P3F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P4FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P4F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P5FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P5F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P3W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P1MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P2MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'P3MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'P3M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q0FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q1FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q1F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q2FNF';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q2F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NF';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q0WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q0W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q1WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q1W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q2WNW';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q2W';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NW';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q0MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q0M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q1MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q1M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'Q2MNM';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'Q2M';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('HSZ3');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1FQ0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P1WQ1WNW');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2WQ0WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P1MQ1MNM');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2MQ0MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HSZ4');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2FQ0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P1FQ1FNF');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1WQ2WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2WQ1WNW');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3WQ0WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P1MQ2MNM');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2MQ1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P3MQ0MNM');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HSZ5');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3FQ0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2FQ1FNF');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1FQ2FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2WQ2WNW');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3WQ1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2MQ2MNM');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3MQ1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HSZ6');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P4FQ0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P3FQ1FNF');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2FQ2FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P3WQ2WNW');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3MQ2MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HSZ7');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P5FQ0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P4FQ1FNF');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3FQ2FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HSZ8');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P5FQ1FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P4FQ2FNF');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST1F');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST2F');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST3F');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P4FNF');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P5FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST1W');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST2W');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST3W');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST1M');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST2M');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HST3M');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('HCH');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2FNF');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P4FNF');
        pbm.grelw{ig}(posel) = 4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P5FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P1MNM');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P2MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P3MNM');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2WNW');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        ig = ig_('HAD');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q0FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q1FNF');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q2FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q0MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q2MNM');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q0WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q2WNW');
        pbm.grelw{ig}(posel) = 3.0;
        ig = ig_('HINF');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2FNF');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P4FNF');
        pbm.grelw{ig}(posel) = 4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P5FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q0FNF');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q1FNF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q2FNF');
        pbm.grelw{ig}(posel) = 4.0;
        ig = ig_('HINN');
        pbm.grftype{ig} = 'gL2';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2WNW');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q0WNW');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q1WNW');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q2WNW');
        pbm.grelw{ig}(posel) = 3.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('P2MNM');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('P3MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q0MNM');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Q1MNM');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Q2MNM');
        pbm.grelw{ig}(posel) = 3.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CSLR2-RN-27-8-0-3-24-0-2-0-8-0-0-0';
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

    case 'en3PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

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

