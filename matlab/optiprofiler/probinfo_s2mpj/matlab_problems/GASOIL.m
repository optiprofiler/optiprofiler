function varargout = GASOIL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : GASOIL
%    *********
% 
%    Determine the reaction coefficients for the catalytic cracking of gas oil
%    and other byproducts. The nonlinear model that describes the process is
% 
%      y_1' = - (theta_1 + theta_3 ) y_1^2
%      y_2' = theta_1 y_1^2 + theta_2 y_2
% 
%    with given initial conditions. The problem is to minimize
% 
%     sum{i=1,20} || y(tau_i,theta) - z_i||^2
% 
%    where the z_i are concentration measurements for y at times tau_i (i=1,20)
% 
%    This is problem 12 in the COPS (Version 2) collection of 
%    E. Dolan and J. More'
%    see "Benchmarking Optimization Software with COPS"
%    Argonne National Labs Technical Report ANL/MCS-246 (2000)
% 
%    SIF input: Nick Gould, November 2000
% 
%    classification = 'C-COOR2-AN-V-V'
% 
%  The number of differential equations
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GASOIL';

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
        v_('NE') = 2;
        if(nargs<1)
            v_('NH') = 10;  %  SIF file default value
        else
            v_('NH') = varargin{1};
        end
%       Alternative values for the SIF file parameters:
% IE NH                  50             $-PARAMETER
% IE NH                  100            $-PARAMETER
% IE NH                  200            $-PARAMETER
% IE NH                  400            $-PARAMETER
        v_('NP') = 3;
        v_('NM') = 21;
        v_('NC') = 4;
        v_('RHO1') = 0.0694318442;
        v_('RHO2') = 0.3300094782;
        v_('RHO3') = 0.6699905218;
        v_('RHO4') = 0.9305681558;
        v_('TAU1') = 0.0;
        v_('TAU2') = 0.025;
        v_('TAU3') = 0.05;
        v_('TAU4') = 0.075;
        v_('TAU5') = 0.10;
        v_('TAU6') = 0.125;
        v_('TAU7') = 0.150;
        v_('TAU8') = 0.175;
        v_('TAU9') = 0.20;
        v_('TAU10') = 0.225;
        v_('TAU11') = 0.250;
        v_('TAU12') = 0.30;
        v_('TAU13') = 0.35;
        v_('TAU14') = 0.40;
        v_('TAU15') = 0.45;
        v_('TAU16') = 0.50;
        v_('TAU17') = 0.55;
        v_('TAU18') = 0.65;
        v_('TAU19') = 0.75;
        v_('TAU20') = 0.85;
        v_('TAU21') = 0.95;
        v_('TF') = v_(['TAU',int2str(round(v_('NM')))]);
        v_('RNH') = v_('NH');
        v_('H') = v_('TF')/v_('RNH');
        v_('Z1,1') = 1.0000;
        v_('Z1,2') = 0.0000;
        v_('Z2,1') = 0.8105;
        v_('Z2,2') = 0.2000;
        v_('Z3,1') = 0.6208;
        v_('Z3,2') = 0.2886;
        v_('Z4,1') = 0.5258;
        v_('Z4,2') = 0.3010;
        v_('Z5,1') = 0.4345;
        v_('Z5,2') = 0.3215;
        v_('Z6,1') = 0.3903;
        v_('Z6,2') = 0.3123;
        v_('Z7,1') = 0.3342;
        v_('Z7,2') = 0.2716;
        v_('Z8,1') = 0.3034;
        v_('Z8,2') = 0.2551;
        v_('Z9,1') = 0.2735;
        v_('Z9,2') = 0.2258;
        v_('Z10,1') = 0.2405;
        v_('Z10,2') = 0.1959;
        v_('Z11,1') = 0.2283;
        v_('Z11,2') = 0.1789;
        v_('Z12,1') = 0.2071;
        v_('Z12,2') = 0.1457;
        v_('Z13,1') = 0.1669;
        v_('Z13,2') = 0.1198;
        v_('Z14,1') = 0.1530;
        v_('Z14,2') = 0.0909;
        v_('Z15,1') = 0.1339;
        v_('Z15,2') = 0.0719;
        v_('Z16,1') = 0.1265;
        v_('Z16,2') = 0.0561;
        v_('Z17,1') = 0.1200;
        v_('Z17,2') = 0.0460;
        v_('Z18,1') = 0.0990;
        v_('Z18,2') = 0.0280;
        v_('Z19,1') = 0.0870;
        v_('Z19,2') = 0.0190;
        v_('Z20,1') = 0.0770;
        v_('Z20,2') = 0.0140;
        v_('Z21,1') = 0.0690;
        v_('Z21,2') = 0.0100;
        v_('BC1') = 1.0;
        v_('BC2') = 0.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('NH-1') = -1+v_('NH');
        v_('FACT0') = 1.0;
        for I=v_('1'):v_('NC')
            v_('RI') = I;
            v_('I-1') = -1+I;
            v_(['FACT',int2str(I)]) = v_(['FACT',int2str(round(v_('I-1')))])*v_('RI');
        end
        for I=v_('1'):v_('NM')
            v_('TAU/H') = v_(['TAU',int2str(I)])/v_('H');
            v_('IT/H') = fix(v_('TAU/H'));
            v_('IT/H+1') = 1+v_('IT/H');
            v_('A') = v_('IT/H+1');
            v_('B') = v_('NH');
            v_('A') = -1*v_('A');
            v_('B') = -1*v_('B');
            v_('A') = v_('A');
            v_('ABSA') = abs(v_('A'));
            v_('ABSA') = fix(v_('ABSA'));
            v_('B') = v_('B');
            v_('ABSB') = abs(v_('B'));
            v_('ABSB') = fix(v_('ABSB'));
            v_('ABSA+B') = v_('ABSA')+v_('ABSB');
            v_('A') = v_('A')+v_('ABSA+B');
            v_('B') = v_('B')+v_('ABSA+B');
            v_('A/B') = fix(v_('A')/v_('B'));
            v_('B/A') = fix(v_('B')/v_('A'));
            v_('SUM') = v_('A/B')+v_('B/A');
            v_('A') = v_('A')*v_('A/B');
            v_('B') = v_('B')*v_('B/A');
            v_('MAXA,B') = v_('A')+v_('B');
            v_('MAXA,B') = fix(v_('MAXA,B')/v_('SUM'));
            v_('MINA,B') = v_('ABSA+B')-v_('MAXA,B');
            v_(['ITAU',int2str(I)]) = v_('MINA,B');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('NP')
            [iv,ix_] = s2mpjlib('ii',['THETA',int2str(I)],ix_);
            pb.xnames{iv} = ['THETA',int2str(I)];
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NE')
                [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
            for K=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [iv,ix_] =...
                          s2mpjlib('ii',['W',int2str(I),',',int2str(K),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['W',int2str(I),',',int2str(K),',',int2str(S)];
                end
            end
            for J=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [iv,ix_] =...
                          s2mpjlib('ii',['U',int2str(I),',',int2str(J),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['U',int2str(I),',',int2str(J),',',int2str(S)];
                    [iv,ix_] =...
                          s2mpjlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['DU',int2str(I),',',int2str(J),',',int2str(S)];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('NM')
            v_('RITAU') = v_(['ITAU',int2str(J)]);
            v_('I') = fix(v_('RITAU'));
            v_('T') = -1.0+v_('RITAU');
            v_('T') = v_('T')*v_('H');
            v_('DIFF') = v_(['TAU',int2str(J)])-v_('T');
            for S=v_('1'):v_('NE')
                v_('RATIO') = v_('DIFF');
                [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(J),',',int2str(S)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I'))),',',int2str(S)]);
                valA(end+1) = 1.0;
                for K=v_('1'):v_('NC')
                    v_('COEF') = v_('RATIO')/v_(['FACT',int2str(K)]);
                    [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(J),',',int2str(S)],ig_);
                    gtype{ig} = '<>';
                    irA(end+1)  = ig;
                    icA(end+1)  =...
                          ix_(['W',int2str(round(v_('I'))),',',int2str(K),',',int2str(S)]);
                    valA(end+1) = v_('COEF');
                    v_('RATIO') = v_('RATIO')*v_('DIFF');
                    v_('RATIO') = v_('RATIO')/v_('H');
                end
            end
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                v_('RH') = v_(['RHO',int2str(J)]);
                [ig,ig_] =...
                      s2mpjlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                irA(end+1)  = ig;
                icA(end+1)  =...
                      ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))]);
                valA(end+1) = -1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('1')))]);
                valA(end+1) = 1.0;
                [ig,ig_] =...
                      s2mpjlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                irA(end+1)  = ig;
                icA(end+1)  =...
                      ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))]);
                valA(end+1) = -1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(round(v_('2')))]);
                valA(end+1) = 1.0;
                v_('PROD') = v_('RH')*v_('H');
                for K=v_('1'):v_('NC')
                    v_('COEF') = v_('PROD')/v_(['FACT',int2str(K)]);
                    [ig,ig_] =...
                          s2mpjlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                    irA(end+1)  = ig;
                    icA(end+1)  =...
                          ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('1')))]);
                    valA(end+1) = v_('COEF');
                    [ig,ig_] =...
                          s2mpjlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                    irA(end+1)  = ig;
                    icA(end+1)  =...
                          ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('2')))]);
                    valA(end+1) = v_('COEF');
                    v_('PROD') = v_('PROD')*v_('RH');
                end
                [ig,ig_] =...
                      s2mpjlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                irA(end+1)  = ig;
                icA(end+1)  =...
                      ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))]);
                valA(end+1) = -1.0;
                [ig,ig_] =...
                      s2mpjlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                irA(end+1)  = ig;
                icA(end+1)  =...
                      ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))]);
                valA(end+1) = -1.0;
                v_('PROD') = 1.0;
                for K=v_('1'):v_('NC')
                    v_('K-1') = -1+K;
                    v_('COEF') = v_('PROD')/v_(['FACT',int2str(round(v_('K-1')))]);
                    [ig,ig_] =...
                          s2mpjlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                    irA(end+1)  = ig;
                    icA(end+1)  =...
                          ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('1')))]);
                    valA(end+1) = v_('COEF');
                    [ig,ig_] =...
                          s2mpjlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                    irA(end+1)  = ig;
                    icA(end+1)  =...
                          ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('2')))]);
                    valA(end+1) = v_('COEF');
                    v_('PROD') = v_('PROD')*v_('RH');
                end
            end
        end
        for I=v_('1'):v_('NH-1')
            v_('I+1') = 1+I;
            for S=v_('1'):v_('NE')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(S)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['C',int2str(I),',',int2str(S)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(I),',',int2str(S)]);
                valA(end+1) = 1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['V',int2str(round(v_('I+1'))),',',int2str(S)]);
                valA(end+1) = -1.0;
                for J=v_('1'):v_('NC')
                    v_('COEF') = v_('H')/v_(['FACT',int2str(J)]);
                    [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(S)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['C',int2str(I),',',int2str(S)];
                    irA(end+1)  = ig;
                    icA(end+1)  = ix_(['W',int2str(I),',',int2str(J),',',int2str(S)]);
                    valA(end+1) = v_('COEF');
                end
            end
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [ig,ig_] =...
                          s2mpjlib('ii',['CO',int2str(I),',',int2str(J),',',int2str(S)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['CO',int2str(I),',',int2str(J),',',int2str(S)];
                    irA(end+1)  = ig;
                    icA(end+1)  = ix_(['DU',int2str(I),',',int2str(J),',',int2str(S)]);
                    valA(end+1) = 1.0;
                end
            end
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
        for J=v_('1'):v_('NM')
            for S=v_('1'):v_('NE')
                pbm.gconst(ig_(['OBJ',int2str(J),',',int2str(S)])) =...
                      v_(['Z',int2str(J),',',int2str(S)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('NP')
            pb.xlower(ix_(['THETA',int2str(I)]),1) = 0.0;
        end
        for S=v_('1'):v_('NE')
            pb.xlower(ix_(['V',int2str(round(v_('1'))),',',int2str(S)]),1) =...
                  v_(['BC',int2str(S)]);
            pb.xupper(ix_(['V',int2str(round(v_('1'))),',',int2str(S)]),1) =...
                  v_(['BC',int2str(S)]);
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('NP')
            if(isKey(ix_,['THETA',int2str(I)]))
                pb.x0(ix_(['THETA',int2str(I)]),1) = 0.0;
            else
                pb.y0(find(pbm.congrps==ig_(['THETA',int2str(I)])),1) = 0.0;
            end
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NE')
                if(isKey(ix_,['V',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['V',int2str(I),',',int2str(J)]),1) = 0.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(I),',',int2str(J)])),1) = 0.0;
                end
            end
        end
        v_('I1') = 1;
        v_('RITAU') = v_(['ITAU',int2str(round(v_('1')))]);
        v_('I2') = fix(v_('RITAU'));
        for I=v_('I1'):v_('I2')
            for S=v_('1'):v_('NE')
                if(isKey(ix_,['V',int2str(I),',',int2str(S)]))
                    pb.x0(ix_(['V',int2str(I),',',int2str(S)]),1) = v_(['BC',int2str(S)]);
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(I),',',int2str(S)])),1) =...
                          v_(['BC',int2str(S)]);
                end
                for J=v_('1'):v_('NC')
                    if(isKey(ix_,['W',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['W',int2str(I),',',int2str(J),',',int2str(S)]),1) = 0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['W',int2str(I),',',int2str(J),',',int2str(S)])),1) = 0.0;
                    end
                    if(isKey(ix_,['U',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['U',int2str(I),',',int2str(J),',',int2str(S)]),1) =...
                              v_(['BC',int2str(S)]);
                    else
                        pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(J),',',int2str(S)])),1) = v_(['BC',int2str(S)]);
                    end
                    if(isKey(ix_,['DU',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['DU',int2str(I),',',int2str(J),',',int2str(S)]),1) = 0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['DU',int2str(I),',',int2str(J),',',int2str(S)])),1) = 0.0;
                    end
                end
            end
        end
        for K=v_('2'):v_('NM')
            v_('I1') = 1+v_('I2');
            v_('RITAU') = v_(['ITAU',int2str(K)]);
            v_('I2') = fix(v_('RITAU'));
            for I=v_('I1'):v_('I2')
                v_('S') = 1;
                if(isKey(ix_,['V',int2str(I),',',int2str(round(v_('S')))]))
                    pb.x0(ix_(['V',int2str(I),',',int2str(round(v_('S')))]),1) =...
                          v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(I),',',int2str(round(v_('S')))])),1) = v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                end
                for J=v_('1'):v_('NC')
                    if(isKey(ix_,['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = 0.0;
                    end
                    if(isKey(ix_,['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                    else
                        pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                    end
                    if(isKey(ix_,['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = 0.0;
                    end
                end
                v_('S') = 2;
                if(isKey(ix_,['V',int2str(I),',',int2str(round(v_('S')))]))
                    pb.x0(ix_(['V',int2str(I),',',int2str(round(v_('S')))]),1) =...
                          v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(I),',',int2str(round(v_('S')))])),1) = v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                end
                for J=v_('1'):v_('NC')
                    if(isKey(ix_,['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['W',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = 0.0;
                    end
                    if(isKey(ix_,['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                    else
                        pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = v_(['Z',int2str(K),',',int2str(round(v_('S')))]);
                    end
                    if(isKey(ix_,['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]))
                        pb.x0(ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))]),1) =...
                              0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('S')))])),1) = 0.0;
                    end
                end
            end
        end
        v_('I1') = 1+v_('I2');
        v_('I2') = v_('NH');
        for I=v_('I1'):v_('I2')
            for S=v_('1'):v_('NE')
                if(isKey(ix_,['V',int2str(I),',',int2str(S)]))
                    pb.x0(ix_(['V',int2str(I),',',int2str(S)]),1) =...
                          v_(['Z',int2str(round(v_('NM'))),',',int2str(S)]);
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(I),',',int2str(S)])),1) =...
                          v_(['Z',int2str(round(v_('NM'))),',',int2str(S)]);
                end
                for J=v_('1'):v_('NC')
                    if(isKey(ix_,['W',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['W',int2str(I),',',int2str(J),',',int2str(S)]),1) = 0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['W',int2str(I),',',int2str(J),',',int2str(S)])),1) = 0.0;
                    end
                    if(isKey(ix_,['U',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['U',int2str(I),',',int2str(J),',',int2str(S)]),1) =...
                              v_(['Z',int2str(round(v_('NM'))),',',int2str(S)]);
                    else
                        pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(J),',',int2str(S)])),1) = v_(['Z',int2str(round(v_('NM'))),',',int2str(S)]);
                    end
                    if(isKey(ix_,['DU',int2str(I),',',int2str(J),',',int2str(S)]))
                        pb.x0(ix_(['DU',int2str(I),',',int2str(J),',',int2str(S)]),1) = 0.0;
                    else
                        pb.y0(find(pbm.congrps==ig_(['DU',int2str(I),',',int2str(J),',',int2str(S)])),1) = 0.0;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD1',iet_);
        elftv{it}{1} = 'THETA1';
        elftv{it}{2} = 'THETA3';
        elftv{it}{3} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'THETA';
        elftv{it}{2} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD3',iet_);
        elftv{it}{1} = 'THETA';
        elftv{it}{2} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                ename = ['P1',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD1';
                ielftype(ie) = iet_('ePROD1');
                vname = ['THETA',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETA1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['THETA',int2str(round(v_('3')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETA3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['P2',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD2';
                ielftype(ie) = iet_('ePROD2');
                vname = ['THETA',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETA',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['P3',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD3';
                ielftype(ie) = iet_('ePROD3');
                vname = ['THETA',int2str(round(v_('2')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETA',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                ig = ig_(['CO',int2str(I),',',int2str(J),',',int2str(round(v_('1')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P1',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['CO',int2str(I),',',int2str(J),',',int2str(round(v_('2')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P2',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['P3',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for J=v_('1'):v_('NM')
            for S=v_('1'):v_('NE')
                ig = ig_(['OBJ',int2str(J),',',int2str(S)]);
                pbm.grftype{ig} = 'gL2';
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLUTION             5.23664D-03   $ (NH=50)
% LO SOLUTION             5.23659D-03   $ (NH=100)
% LO SOLUTION             5.23659D-03   $ (NH=200)
% LO SOLUTION             5.23659D-03   $ (NH=400)
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-V-V';
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

    case 'ePROD1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)^2;
        if(nargout>1)
            g_(1,1) = IV_(2)^2;
            g_(2,1) = 2.0*IV_(1)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 2.0*IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*IV_(1);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(1)*EV_(2)^2;
        if(nargout>1)
            g_(1,1) = -EV_(2)^2;
            g_(2,1) = -2.0*EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -2.0*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -2.0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'ePROD3'

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

