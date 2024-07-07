function varargout = OPTPRLOC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OPTPRLOC
%    *********
% 
%    Optimal positioning of a new product in a multiattribute space.
%    Consider a market of M existing products, a set of N consumers
%    in a multiattribute (dim K) space.
% 
%    Source: Test problem 4 in M. Duran & I.E. Grossmann,
%    "An outer approximation algorithm for a class of mixed integer nonlinear
%     programs", Mathematical Programming 36, pp. 307-339, 1986.
% 
%    SIF input: S. Leyffer, October 1997
% 
%    classification = 'QQR2-AN-30-30'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OPTPRLOC';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('1') = 1;
        v_('K') = 5;
        v_('M') = 10;
        v_('N') = 25;
        v_('H') = 1000.0;
        v_('Z1,1') = 2.26;
        v_('Z1,2') = 5.15;
        v_('Z1,3') = 4.03;
        v_('Z1,4') = 1.74;
        v_('Z1,5') = 4.74;
        v_('Z2,1') = 5.51;
        v_('Z2,2') = 9.01;
        v_('Z2,3') = 3.84;
        v_('Z2,4') = 1.47;
        v_('Z2,5') = 9.92;
        v_('Z3,1') = 4.06;
        v_('Z3,2') = 1.80;
        v_('Z3,3') = 0.71;
        v_('Z3,4') = 9.09;
        v_('Z3,5') = 8.13;
        v_('Z4,1') = 6.30;
        v_('Z4,2') = 0.11;
        v_('Z4,3') = 4.08;
        v_('Z4,4') = 7.29;
        v_('Z4,5') = 4.24;
        v_('Z5,1') = 2.81;
        v_('Z5,2') = 1.65;
        v_('Z5,3') = 8.08;
        v_('Z5,4') = 3.99;
        v_('Z5,5') = 3.51;
        v_('Z6,1') = 4.29;
        v_('Z6,2') = 9.49;
        v_('Z6,3') = 2.24;
        v_('Z6,4') = 9.78;
        v_('Z6,5') = 1.52;
        v_('Z7,1') = 9.76;
        v_('Z7,2') = 3.64;
        v_('Z7,3') = 6.62;
        v_('Z7,4') = 3.66;
        v_('Z7,5') = 9.08;
        v_('Z8,1') = 1.37;
        v_('Z8,2') = 6.99;
        v_('Z8,3') = 7.19;
        v_('Z8,4') = 3.03;
        v_('Z8,5') = 3.39;
        v_('Z9,1') = 8.89;
        v_('Z9,2') = 8.29;
        v_('Z9,3') = 6.05;
        v_('Z9,4') = 7.48;
        v_('Z9,5') = 4.09;
        v_('Z10,1') = 7.42;
        v_('Z10,2') = 4.60;
        v_('Z10,3') = 0.30;
        v_('Z10,4') = 0.97;
        v_('Z10,5') = 8.77;
        v_('Z11,1') = 1.54;
        v_('Z11,2') = 7.06;
        v_('Z11,3') = 0.01;
        v_('Z11,4') = 1.23;
        v_('Z11,5') = 3.11;
        v_('Z12,1') = 7.74;
        v_('Z12,2') = 4.40;
        v_('Z12,3') = 7.93;
        v_('Z12,4') = 5.95;
        v_('Z12,5') = 4.88;
        v_('Z13,1') = 9.94;
        v_('Z13,2') = 5.21;
        v_('Z13,3') = 8.58;
        v_('Z13,4') = 0.13;
        v_('Z13,5') = 4.57;
        v_('Z14,1') = 9.54;
        v_('Z14,2') = 1.57;
        v_('Z14,3') = 9.66;
        v_('Z14,4') = 5.24;
        v_('Z14,5') = 7.90;
        v_('Z15,1') = 7.46;
        v_('Z15,2') = 8.81;
        v_('Z15,3') = 1.67;
        v_('Z15,4') = 6.47;
        v_('Z15,5') = 1.81;
        v_('Z16,1') = 0.56;
        v_('Z16,2') = 8.10;
        v_('Z16,3') = 0.19;
        v_('Z16,4') = 6.11;
        v_('Z16,5') = 6.40;
        v_('Z17,1') = 3.86;
        v_('Z17,2') = 6.68;
        v_('Z17,3') = 6.42;
        v_('Z17,4') = 7.29;
        v_('Z17,5') = 4.66;
        v_('Z18,1') = 2.98;
        v_('Z18,2') = 2.98;
        v_('Z18,3') = 3.03;
        v_('Z18,4') = 0.02;
        v_('Z18,5') = 0.67;
        v_('Z19,1') = 3.61;
        v_('Z19,2') = 7.62;
        v_('Z19,3') = 1.79;
        v_('Z19,4') = 7.80;
        v_('Z19,5') = 9.81;
        v_('Z20,1') = 5.68;
        v_('Z20,2') = 4.24;
        v_('Z20,3') = 4.17;
        v_('Z20,4') = 6.75;
        v_('Z20,5') = 1.08;
        v_('Z21,1') = 5.48;
        v_('Z21,2') = 3.74;
        v_('Z21,3') = 3.34;
        v_('Z21,4') = 6.22;
        v_('Z21,5') = 7.94;
        v_('Z22,1') = 8.13;
        v_('Z22,2') = 8.72;
        v_('Z22,3') = 3.93;
        v_('Z22,4') = 8.80;
        v_('Z22,5') = 8.56;
        v_('Z23,1') = 1.37;
        v_('Z23,2') = 0.54;
        v_('Z23,3') = 1.55;
        v_('Z23,4') = 5.56;
        v_('Z23,5') = 5.85;
        v_('Z24,1') = 8.79;
        v_('Z24,2') = 5.04;
        v_('Z24,3') = 4.83;
        v_('Z24,4') = 6.94;
        v_('Z24,5') = 0.38;
        v_('Z25,1') = 2.66;
        v_('Z25,2') = 4.19;
        v_('Z25,3') = 6.49;
        v_('Z25,4') = 8.04;
        v_('Z25,5') = 1.66;
        v_('W1,1') = 9.57;
        v_('W1,2') = 2.74;
        v_('W1,3') = 9.75;
        v_('W1,4') = 3.96;
        v_('W1,5') = 8.67;
        v_('W2,1') = 8.38;
        v_('W2,2') = 3.93;
        v_('W2,3') = 5.18;
        v_('W2,4') = 5.20;
        v_('W2,5') = 7.82;
        v_('W3,1') = 9.81;
        v_('W3,2') = 0.04;
        v_('W3,3') = 4.21;
        v_('W3,4') = 7.38;
        v_('W3,5') = 4.11;
        v_('W4,1') = 7.41;
        v_('W4,2') = 6.08;
        v_('W4,3') = 5.46;
        v_('W4,4') = 4.86;
        v_('W4,5') = 1.48;
        v_('W5,1') = 9.96;
        v_('W5,2') = 9.13;
        v_('W5,3') = 2.95;
        v_('W5,4') = 8.25;
        v_('W5,5') = 3.58;
        v_('W6,1') = 9.39;
        v_('W6,2') = 4.27;
        v_('W6,3') = 5.09;
        v_('W6,4') = 1.81;
        v_('W6,5') = 7.58;
        v_('W7,1') = 1.88;
        v_('W7,2') = 7.20;
        v_('W7,3') = 6.65;
        v_('W7,4') = 1.74;
        v_('W7,5') = 2.86;
        v_('W8,1') = 4.01;
        v_('W8,2') = 2.67;
        v_('W8,3') = 4.86;
        v_('W8,4') = 2.55;
        v_('W8,5') = 6.91;
        v_('W9,1') = 4.18;
        v_('W9,2') = 1.92;
        v_('W9,3') = 2.60;
        v_('W9,4') = 7.15;
        v_('W9,5') = 2.86;
        v_('W10,1') = 7.81;
        v_('W10,2') = 2.14;
        v_('W10,3') = 9.63;
        v_('W10,4') = 7.61;
        v_('W10,5') = 9.17;
        v_('W11,1') = 8.96;
        v_('W11,2') = 3.47;
        v_('W11,3') = 5.49;
        v_('W11,4') = 4.73;
        v_('W11,5') = 9.43;
        v_('W12,1') = 9.94;
        v_('W12,2') = 1.63;
        v_('W12,3') = 1.23;
        v_('W12,4') = 4.33;
        v_('W12,5') = 7.08;
        v_('W13,1') = 0.31;
        v_('W13,2') = 5.00;
        v_('W13,3') = 0.16;
        v_('W13,4') = 2.52;
        v_('W13,5') = 3.08;
        v_('W14,1') = 6.02;
        v_('W14,2') = 0.92;
        v_('W14,3') = 7.47;
        v_('W14,4') = 9.74;
        v_('W14,5') = 1.76;
        v_('W15,1') = 5.06;
        v_('W15,2') = 4.52;
        v_('W15,3') = 1.89;
        v_('W15,4') = 1.22;
        v_('W15,5') = 9.05;
        v_('W16,1') = 5.92;
        v_('W16,2') = 2.56;
        v_('W16,3') = 7.74;
        v_('W16,4') = 6.96;
        v_('W16,5') = 5.18;
        v_('W17,1') = 6.45;
        v_('W17,2') = 1.52;
        v_('W17,3') = 0.06;
        v_('W17,4') = 5.34;
        v_('W17,5') = 8.47;
        v_('W18,1') = 1.04;
        v_('W18,2') = 1.36;
        v_('W18,3') = 5.99;
        v_('W18,4') = 8.10;
        v_('W18,5') = 5.22;
        v_('W19,1') = 1.40;
        v_('W19,2') = 1.35;
        v_('W19,3') = 0.59;
        v_('W19,4') = 8.58;
        v_('W19,5') = 1.21;
        v_('W20,1') = 6.68;
        v_('W20,2') = 9.48;
        v_('W20,3') = 1.60;
        v_('W20,4') = 6.74;
        v_('W20,5') = 8.92;
        v_('W21,1') = 1.95;
        v_('W21,2') = 0.46;
        v_('W21,3') = 2.90;
        v_('W21,4') = 1.79;
        v_('W21,5') = 0.99;
        v_('W22,1') = 5.18;
        v_('W22,2') = 5.10;
        v_('W22,3') = 8.81;
        v_('W22,4') = 3.27;
        v_('W22,5') = 9.63;
        v_('W23,1') = 1.47;
        v_('W23,2') = 5.71;
        v_('W23,3') = 6.95;
        v_('W23,4') = 1.42;
        v_('W23,5') = 3.49;
        v_('W24,1') = 5.40;
        v_('W24,2') = 3.12;
        v_('W24,3') = 5.37;
        v_('W24,4') = 6.10;
        v_('W24,5') = 3.71;
        v_('W25,1') = 6.32;
        v_('W25,2') = 0.81;
        v_('W25,3') = 6.12;
        v_('W25,4') = 6.73;
        v_('W25,5') = 7.93;
        v_('DEL1,1') = 0.62;
        v_('DEL1,2') = 5.06;
        v_('DEL1,3') = 7.82;
        v_('DEL1,4') = 0.22;
        v_('DEL1,5') = 4.42;
        v_('DEL2,1') = 5.21;
        v_('DEL2,2') = 2.66;
        v_('DEL2,3') = 9.54;
        v_('DEL2,4') = 5.03;
        v_('DEL2,5') = 8.01;
        v_('DEL3,1') = 5.27;
        v_('DEL3,2') = 7.72;
        v_('DEL3,3') = 7.97;
        v_('DEL3,4') = 3.31;
        v_('DEL3,5') = 6.56;
        v_('DEL4,1') = 1.02;
        v_('DEL4,2') = 8.89;
        v_('DEL4,3') = 8.77;
        v_('DEL4,4') = 3.10;
        v_('DEL4,5') = 6.66;
        v_('DEL5,1') = 1.26;
        v_('DEL5,2') = 6.80;
        v_('DEL5,3') = 2.30;
        v_('DEL5,4') = 1.75;
        v_('DEL5,5') = 6.65;
        v_('DEL6,1') = 3.74;
        v_('DEL6,2') = 9.06;
        v_('DEL6,3') = 9.80;
        v_('DEL6,4') = 3.01;
        v_('DEL6,5') = 9.52;
        v_('DEL7,1') = 4.64;
        v_('DEL7,2') = 7.99;
        v_('DEL7,3') = 6.69;
        v_('DEL7,4') = 5.88;
        v_('DEL7,5') = 8.23;
        v_('DEL8,1') = 8.35;
        v_('DEL8,2') = 3.79;
        v_('DEL8,3') = 1.19;
        v_('DEL8,4') = 1.96;
        v_('DEL8,5') = 5.88;
        v_('DEL9,1') = 6.44;
        v_('DEL9,2') = 0.17;
        v_('DEL9,3') = 9.93;
        v_('DEL9,4') = 6.80;
        v_('DEL9,5') = 9.75;
        v_('DEL10,1') = 6.49;
        v_('DEL10,2') = 1.92;
        v_('DEL10,3') = 0.05;
        v_('DEL10,4') = 4.89;
        v_('DEL10,5') = 6.43;
        v_('R1') = 77.83985;
        v_('R2') = 175.9710;
        v_('R3') = 201.8226;
        v_('R4') = 143.9533;
        v_('R5') = 154.3895;
        v_('R6') = 433.3177;
        v_('R7') = 109.0764;
        v_('R8') = 41.59592;
        v_('R9') = 144.0623;
        v_('R10') = 99.83416;
        v_('R11') = 149.1791;
        v_('R12') = 123.8074;
        v_('R13') = 27.22197;
        v_('R14') = 89.92683;
        v_('R15') = 293.0766;
        v_('R16') = 174.3170;
        v_('R17') = 125.1028;
        v_('R18') = 222.8417;
        v_('R19') = 50.48593;
        v_('R20') = 361.1973;
        v_('R21') = 40.32642;
        v_('R22') = 161.8518;
        v_('R23') = 66.85827;
        v_('R24') = 340.5807;
        v_('R25') = 407.5200;
        for I=v_('1'):v_('N')
            v_(['RPH',int2str(I)]) = v_(['R',int2str(I)])+v_('H');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('K')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','PROFIT',ig_);
        gtype{ig} = '<>';
        iv = ix_('Y1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.2;
        end
        iv = ix_('Y3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.2;
        end
        iv = ix_('Y5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('Y6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('Y7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.1;
        end
        iv = ix_('Y8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.8+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.8;
        end
        iv = ix_('Y9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.4;
        end
        iv = ix_('Y11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.3;
        end
        iv = ix_('Y13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.1;
        end
        iv = ix_('Y14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.3;
        end
        iv = ix_('Y15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.5;
        end
        iv = ix_('Y16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('Y17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.8+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.8;
        end
        iv = ix_('Y18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.1;
        end
        iv = ix_('Y19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('Y20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('Y23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.2;
        end
        iv = ix_('Y24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.7+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.7;
        end
        iv = ix_('Y25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.7+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.7;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.5;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['ELLI',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['ELLI',int2str(I)];
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('H');
            end
        end
        [ig,ig_] = s2mpjlib('ii','LIN1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'LIN1';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','LIN2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'LIN2';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.6+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.6;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.9;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.5;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.1;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','LIN3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'LIN3';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','LIN4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'LIN4';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.157+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.157;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.05;
        end
        [ig,ig_] = s2mpjlib('ii','LIN5',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'LIN5';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.25+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.25;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.05;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.3;
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['ELLI',int2str(I)])) = v_(['RPH',int2str(I)]);
        end
        pbm.gconst(ig_('LIN1')) = 10.0;
        pbm.gconst(ig_('LIN2')) = -0.64;
        pbm.gconst(ig_('LIN3')) = 0.69;
        pbm.gconst(ig_('LIN4')) = 1.5;
        pbm.gconst(ig_('LIN5')) = 4.5;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 2.0;
        pb.xupper(ix_('X1')) = 4.5;
        pb.xupper(ix_('X2')) = 8.0;
        pb.xlower(ix_('X3'),1) = 3.0;
        pb.xupper(ix_('X3')) = 9.0;
        pb.xupper(ix_('X4')) = 5.0;
        pb.xlower(ix_('X5'),1) = 4.0;
        pb.xupper(ix_('X5')) = 10.0;
        for I=v_('1'):v_('N')
            pb.xupper(ix_(['Y',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXMBS',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'B';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OBJSQ1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXMBS';
        ielftype(ie) = iet_('eXMBS');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('B',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'OBJSQ2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXMBS';
        ielftype(ie) = iet_('eXMBS');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('B',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        for I=v_('1'):v_('N')
            for L=v_('1'):v_('K')
                ename = ['SQ',int2str(I),',',int2str(L)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXMBS';
                ielftype(ie) = iet_('eXMBS');
                vname = ['X',int2str(L)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('B',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['Z',int2str(I),',',int2str(L)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('PROFIT');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OBJSQ1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.6;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OBJSQ2');
        pbm.grelw{ig}(posel) = 0.1;
        for I=v_('1'):v_('N')
            for L=v_('1'):v_('K')
                ig = ig_(['ELLI',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['SQ',int2str(I),',',int2str(L)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['W',int2str(I),',',int2str(L)]);
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QQR2-AN-30-30';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eXMBS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-pbm.elpar{iel_}(1))^2;
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-pbm.elpar{iel_}(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

