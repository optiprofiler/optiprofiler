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
%    classification = 'OOR2-AN-V-V'
% 
%  The number of differential equations
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GASOIL';

switch(action)

    case 'setup'

    pb.name      = 'GASOIL';
    pb.sifpbname = 'GASOIL';
    pbm.name     = 'GASOIL';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('NE') = 2;
        if(nargin<2)
            v_('NH') = 10;  %  SIF file default value
        else
            v_('NH') = varargin{1};
        end
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
        for I=v_('1'):v_('NP')
            [iv,ix_] = s2xlib('ii',['THETA',int2str(I)],ix_);
            pb.xnames{iv} = ['THETA',int2str(I)];
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NE')
                [iv,ix_] = s2xlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
            for K=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [iv,ix_] = s2xlib('ii',['W',int2str(I),',',int2str(K),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['W',int2str(I),',',int2str(K),',',int2str(S)];
                end
            end
            for J=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [iv,ix_] = s2xlib('ii',['U',int2str(I),',',int2str(J),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['U',int2str(I),',',int2str(J),',',int2str(S)];
                    [iv,ix_] =...
                          s2xlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(S)],ix_);
                    pb.xnames{iv} = ['DU',int2str(I),',',int2str(J),',',int2str(S)];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('NM')
            v_('RITAU') = v_(['ITAU',int2str(J)]);
            v_('I') = fix(v_('RITAU'));
            v_('T') = -1.0+v_('RITAU');
            v_('T') = v_('T')*v_('H');
            v_('DIFF') = v_(['TAU',int2str(J)])-v_('T');
            for S=v_('1'):v_('NE')
                v_('RATIO') = v_('DIFF');
                [ig,ig_] = s2xlib('ii',['OBJ',int2str(J),',',int2str(S)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['V',int2str(round(v_('I'))),',',int2str(S)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                for K=v_('1'):v_('NC')
                    v_('COEF') = v_('RATIO')/v_(['FACT',int2str(K)]);
                    [ig,ig_] = s2xlib('ii',['OBJ',int2str(J),',',int2str(S)],ig_);
                    gtype{ig} = '<>';
                    iv = ix_(['W',int2str(round(v_('I'))),',',int2str(K),',',int2str(S)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                    v_('RATIO') = v_('RATIO')*v_('DIFF');
                    v_('RATIO') = v_('RATIO')/v_('H');
                end
            end
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                v_('RH') = v_(['RHO',int2str(J)]);
                [ig,ig_] =...
                      s2xlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                iv = ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['V',int2str(I),',',int2str(round(v_('1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                [ig,ig_] =...
                      s2xlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                iv = ix_(['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['V',int2str(I),',',int2str(round(v_('2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                v_('PROD') = v_('RH')*v_('H');
                for K=v_('1'):v_('NC')
                    v_('COEF') = v_('PROD')/v_(['FACT',int2str(K)]);
                    [ig,ig_] =...
                          s2xlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                    iv = ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('1')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                    [ig,ig_] =...
                          s2xlib('ii',['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['U',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                    iv = ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('2')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                    v_('PROD') = v_('PROD')*v_('RH');
                end
                [ig,ig_] =...
                      s2xlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                iv = ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] =...
                      s2xlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                iv = ix_(['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                v_('PROD') = 1.0;
                for K=v_('1'):v_('NC')
                    v_('K-1') = -1+K;
                    v_('COEF') = v_('PROD')/v_(['FACT',int2str(round(v_('K-1')))]);
                    [ig,ig_] =...
                          s2xlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('1')))];
                    iv = ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('1')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                    [ig,ig_] =...
                          s2xlib('ii',['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['DU',int2str(I),',',int2str(J),',',int2str(round(v_('2')))];
                    iv = ix_(['W',int2str(I),',',int2str(K),',',int2str(round(v_('2')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                    v_('PROD') = v_('PROD')*v_('RH');
                end
            end
        end
        for I=v_('1'):v_('NH-1')
            v_('I+1') = 1+I;
            for S=v_('1'):v_('NE')
                [ig,ig_] = s2xlib('ii',['C',int2str(I),',',int2str(S)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['C',int2str(I),',',int2str(S)];
                iv = ix_(['V',int2str(I),',',int2str(S)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['V',int2str(round(v_('I+1'))),',',int2str(S)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                for J=v_('1'):v_('NC')
                    v_('COEF') = v_('H')/v_(['FACT',int2str(J)]);
                    [ig,ig_] = s2xlib('ii',['C',int2str(I),',',int2str(S)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['C',int2str(I),',',int2str(S)];
                    iv = ix_(['W',int2str(I),',',int2str(J),',',int2str(S)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('COEF')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('COEF');
                    end
                end
            end
        end
        for I=v_('1'):v_('NH')
            for J=v_('1'):v_('NC')
                for S=v_('1'):v_('NE')
                    [ig,ig_] =...
                          s2xlib('ii',['CO',int2str(I),',',int2str(J),',',int2str(S)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['CO',int2str(I),',',int2str(J),',',int2str(S)];
                    iv = ix_(['DU',int2str(I),',',int2str(J),',',int2str(S)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1.0;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
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
