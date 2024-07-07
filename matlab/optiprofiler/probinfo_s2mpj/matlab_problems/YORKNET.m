function varargout = YORKNET(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A problem arising in the modelling of the Yorkshire water system.
% 
%    Source:
%    an problem submitted for the LANCELOT licence.
% 
%    SIF input: B. Ulanicki, Water Software Systems,De Montfort University,
%               The Gateway, Leicester LE1 9BH, UK.
%               e-mail: bul@uk.ac.dmu * Tel no.0533 577070
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'SOR2-AY-312-256'
% 
% DECLARE CONSTANTS DESCRIBING NETWORK
% 
% STANDARD DECLARATIONS
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'YORKNET';

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
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_('15') = 15;
        v_('16') = 16;
        v_('17') = 17;
        v_('18') = 18;
        v_('19') = 19;
        v_('ONE') = 1;
        v_('NSTEP') = 8;
        v_('NSTEP+1') = 8+v_('NSTEP');
        v_('ST1') = 0.125;
        v_('ST2') = 0.125;
        v_('ST3') = 0.125;
        v_('ST4') = 0.125;
        v_('ST5') = 0.125;
        v_('ST6') = 0.125;
        v_('ST7') = 0.125;
        v_('ST8') = 0.125;
        v_('MST1') = 0.0-v_('ST1');
        v_('MST2') = 0.0-v_('ST2');
        v_('MST3') = 0.0-v_('ST3');
        v_('MST4') = 0.0-v_('ST4');
        v_('MST5') = 0.0-v_('ST5');
        v_('MST6') = 0.0-v_('ST6');
        v_('MST7') = 0.0-v_('ST7');
        v_('MST8') = 0.0-v_('ST8');
        v_('TFL1') = 1.87;
        v_('TFH1') = 4.03;
        v_('TF1') = 0.0+v_('TFL1');
        v_('TF2') = 0.0+v_('TFL1');
        v_('TF3') = 0.0+v_('TFH1');
        v_('TF4') = 0.0+v_('TFH1');
        v_('TF5') = 0.0+v_('TFH1');
        v_('TF6') = 0.0+v_('TFH1');
        v_('TF7') = 0.0+v_('TFH1');
        v_('TF8') = 0.0+v_('TFH1');
        v_('BP1') = 24.0*v_('ST1');
        v_('BP2') = 24.0*v_('ST2');
        v_('BP3') = 24.0*v_('ST3');
        v_('BP4') = 24.0*v_('ST4');
        v_('BP5') = 24.0*v_('ST5');
        v_('BP6') = 24.0*v_('ST6');
        v_('BP7') = 24.0*v_('ST7');
        v_('BP8') = 24.0*v_('ST8');
        v_('CP1') = v_('BP1')*v_('TF1');
        v_('CP2') = v_('BP2')*v_('TF2');
        v_('CP3') = v_('BP3')*v_('TF3');
        v_('CP4') = v_('BP4')*v_('TF4');
        v_('CP5') = v_('BP5')*v_('TF5');
        v_('CP6') = v_('BP6')*v_('TF6');
        v_('CP7') = v_('BP7')*v_('TF7');
        v_('CP8') = v_('BP8')*v_('TF8');
        v_('NEL') = 15;
        v_('SuVAL') = 9;
        v_('EuVAL') = 10;
        v_('SuPIPE') = 1;
        v_('EuPIPE') = 8;
        v_('MGINV1') = -3.365170e-3;
        v_('MGINV2') = -2.314284e-2;
        v_('MGINV3') = -6.631203e-3;
        v_('MGINV4') = -1.702093e-3;
        v_('MGINV5') = -1.205983e-2;
        v_('MGINV6') = -9.617776e-4;
        v_('MGINV7') = -1.392046e-5;
        v_('MGINV8') = -4.411625e-3;
        v_('MGINV9') = -2.019250e-3;
        v_('MGINV10') = -2.288437e-3;
        v_('NPMP') = 5;
        v_('SuPMP') = 11;
        v_('EuPMP') = 15;
        v_('CONA11') = -.035520;
        v_('CONB11') = -.054720;
        v_('CONC11') = 99.80;
        v_('CONA12') = -0.07475;
        v_('CONB12') = -9.05;
        v_('CONC12') = 110;
        v_('CONA13') = -.042420;
        v_('CONB13') = -.005370;
        v_('CONC13') = 175.29;
        v_('CONA14') = -.040733;
        v_('CONB14') = -.032036;
        v_('CONC14') = 139.6;
        v_('CONA15') = -.167495;
        v_('CONB15') = -.0019;
        v_('CONC15') = 139.6;
        v_('NND') = 13;
        v_('D1') = 0.0;
        v_('D4') = -33.0;
        v_('D2') = 0.0;
        v_('D3') = -55.0;
        v_('D5') = 0.0;
        v_('D6') = 0.0;
        v_('D7') = 0.0;
        v_('D8') = -25.0;
        v_('D9') = 0.0;
        v_('D10') = -17.0;
        v_('D11') = 0.0;
        v_('D12') = 0.0;
        v_('D13') = 0.0;
        v_('SuRES') = 1;
        v_('EuRES') = 4;
        v_('EuRES+1') = 1+v_('EuRES');
        v_('HGT1') = 5.77;
        v_('HGT2') = 3.00;
        v_('HGT3') = 131.08;
        v_('HGT4') = 44.0;
        v_('MXHGT1') = 9.60;
        v_('MXHGT2') = 7.89;
        v_('MXHGT3') = 138.76;
        v_('MXHGT4') = 53.34;
        v_('XAR1') = 1.599;
        v_('XAR2') = 4.6421;
        v_('XAR3') = 30.2307;
        v_('XAR4') = 5.3938;
        v_('MXAR1') = 0.0-v_('XAR1');
        v_('MXAR2') = 0.0-v_('XAR2');
        v_('MXAR3') = 0.0-v_('XAR3');
        v_('MXAR4') = 0.0-v_('XAR4');
        v_('RXAR1') = 1.0/v_('XAR1');
        v_('MRXAR1') = 0.0-v_('RXAR1');
        v_('RXAR2') = 1.0/v_('XAR2');
        v_('MRXAR2') = 0.0-v_('RXAR2');
        v_('RXAR3') = 1.0/v_('XAR3');
        v_('MRXAR3') = 0.0-v_('RXAR3');
        v_('RXAR4') = 1.0/v_('XAR4');
        v_('MRXAR4') = 0.0-v_('RXAR4');
        v_('HTXAR1') = v_('XAR1')*v_('HGT1');
        v_('HTXAR2') = v_('XAR2')*v_('HGT2');
        v_('HTXAR3') = v_('XAR3')*v_('HGT3');
        v_('HTXAR4') = v_('XAR4')*v_('HGT4');
        v_('STHGT1') = 8.5;
        v_('STHGT2') = 6.0;
        v_('STHGT3') = 135.6;
        v_('STHGT4') = 48.5;
        v_('STVOL1') = v_('STHGT1')*v_('XAR1');
        v_('STVOL2') = v_('STHGT2')*v_('XAR2');
        v_('STVOL3') = v_('STHGT3')*v_('XAR3');
        v_('STVOL4') = v_('STHGT4')*v_('XAR4');
        v_('MSTVOL1') = 0.0-v_('STVOL1');
        v_('MSTVOL2') = 0.0-v_('STVOL2');
        v_('MSTVOL3') = 0.0-v_('STVOL3');
        v_('MSTVOL4') = 0.0-v_('STVOL4');
        v_('WMN1') = 3.764;
        v_('WMN2') = 11.35;
        v_('WMN3') = 156.648;
        v_('WMN4') = 45.929;
        v_('WMX1') = 5.646;
        v_('WMX2') = 22.133;
        v_('WMX3') = 223.489;
        v_('WMX4') = 61.876;
        v_('H1') = 8.99;
        v_('H2') = 52.84;
        v_('H3') = 138.31;
        v_('H4') = 5.67;
        v_('W1') = 4.728;
        v_('W2') = 15.601;
        v_('W3') = 190.648;
        v_('W4') = 55.00;
        v_('SuTW') = 1;
        v_('EuTW') = 2;
        v_('BCTW1') = 28.34;
        v_('BCTW2') = 18.86;
        v_('BCTW3') = 1.0;
        v_('BCTW4') = 1.0;
        v_('BCTW5') = 28.34;
        v_('BCTW6') = 18.86;
        v_('BCTW7') = 1.0;
        v_('BCTW8') = 1.0;
        v_('CTW1,1') = v_('BCTW1')*v_('ST1');
        v_('CTW1,2') = v_('BCTW1')*v_('ST2');
        v_('CTW1,3') = v_('BCTW1')*v_('ST3');
        v_('CTW1,4') = v_('BCTW1')*v_('ST4');
        v_('CTW1,5') = v_('BCTW1')*v_('ST5');
        v_('CTW1,6') = v_('BCTW1')*v_('ST6');
        v_('CTW1,7') = v_('BCTW1')*v_('ST7');
        v_('CTW1,8') = v_('BCTW1')*v_('ST8');
        v_('CTW2,1') = v_('BCTW2')*v_('ST1');
        v_('CTW2,2') = v_('BCTW2')*v_('ST2');
        v_('CTW2,3') = v_('BCTW2')*v_('ST3');
        v_('CTW2,4') = v_('BCTW2')*v_('ST4');
        v_('CTW2,5') = v_('BCTW2')*v_('ST5');
        v_('CTW2,6') = v_('BCTW2')*v_('ST6');
        v_('CTW2,7') = v_('BCTW2')*v_('ST7');
        v_('CTW2,8') = v_('BCTW2')*v_('ST8');
        v_('WS1') = 94.67;
        v_('WS2') = 30.87;
        v_('QSMX1') = 500.0;
        v_('QSMX2') = 40.0;
        v_('SCuQ') = 1.0;
        v_('SCuH') = 1.0;
        v_('SCuV') = 1.0;
        v_('PIPEuSC') = 1.0;
        v_('VALuSC') = 1.0;
        v_('NDuSC') = 1.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('ONE'):v_('NSTEP')
            for J=v_('ONE'):v_('NEL')
                [iv,ix_] = s2mpjlib('ii',['q',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['q',int2str(J),',',int2str(I)];
            end
            for J=v_('ONE'):v_('NND')
                [iv,ix_] = s2mpjlib('ii',['h',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['h',int2str(J),',',int2str(I)];
            end
            for J=v_('SuRES'):v_('EuRES')
                [iv,ix_] = s2mpjlib('ii',['qr',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['qr',int2str(J),',',int2str(I)];
            end
            for J=v_('SuTW'):v_('EuTW')
                [iv,ix_] = s2mpjlib('ii',['qs',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['qs',int2str(J),',',int2str(I)];
            end
            for J=v_('SuPMP'):v_('EuPMP')
                [iv,ix_] = s2mpjlib('ii',['u',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['u',int2str(J),',',int2str(I)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('ONE'):v_('NSTEP')
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('1'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('1'))),',',int2str(I)];
            iv = ix_(['qr',int2str(round(v_('1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['qs',int2str(round(v_('1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('1'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('1'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('13'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('2'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('2'))),',',int2str(I)];
            iv = ix_(['qr',int2str(round(v_('2'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['qs',int2str(round(v_('2'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('2'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('2'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('14'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('15'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('3'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('3'))),',',int2str(I)];
            iv = ix_(['qr',int2str(round(v_('3'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('5'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('3'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('3'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('6'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('4'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('4'))),',',int2str(I)];
            iv = ix_(['qr',int2str(round(v_('4'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('12'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('4'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('4'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('8'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('9'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('5'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('5'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('13'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('6'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('6'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('2'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('3'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('6'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('6'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('12'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('7'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('7'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('2'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('7'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('7'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('7'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('11'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('8'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('8'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('3'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('7'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('9'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('9'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('4'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('5'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ND',int2str(round(v_('9'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('9'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('11'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['ND',int2str(round(v_('10'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('10'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('4'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['q',int2str(round(v_('6'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['ND',int2str(round(v_('11'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('11'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('10'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('14'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['ND',int2str(round(v_('11'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('11'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('15'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['ND',int2str(round(v_('12'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('12'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('8'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('10'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['ND',int2str(round(v_('13'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ND',int2str(round(v_('13'))),',',int2str(I)];
            iv = ix_(['q',int2str(round(v_('9'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['q',int2str(round(v_('1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for I=v_('1'):v_('NSTEP')
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('1'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('1'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('13'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('5'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('2'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('2'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('7'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('6'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('3'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('3'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('8'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('6'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('4'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('4'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('10'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('9'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('5'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('5'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('3'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('9'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('6'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('6'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('10'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('3'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('7'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('7'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('7'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('8'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('8'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('8'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('4'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('12'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EL',int2str(round(v_('9'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('9'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('4'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('13'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['EL',int2str(round(v_('10'))),',',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EL',int2str(round(v_('10'))),',',int2str(I)];
            iv = ix_(['h',int2str(round(v_('12'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['h',int2str(round(v_('11'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            for J=v_('SuPMP'):v_('EuPMP')
                [ig,ig_] = s2mpjlib('ii',['EL',int2str(J),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['EL',int2str(J),',',int2str(I)];
            end
        end
        for I=v_('ONE'):v_('NSTEP')
            for J=v_('SuTW'):v_('EuTW')
                [ig,ig_] = s2mpjlib('ii',['TWC',int2str(J),',',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['qs',int2str(J),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['CTW',int2str(J),',',int2str(I)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['CTW',int2str(J),',',int2str(I)]);
                end
            end
            for J=v_('SuPMP'):v_('EuPMP')
                [ig,ig_] = s2mpjlib('ii',['PC',int2str(J),',',int2str(I)],ig_);
                gtype{ig} = '<>';
            end
        end
        for J=v_('SuRES'):v_('EuRES')
            [ig,ig_] = s2mpjlib('ii',['RD',int2str(J),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['RD',int2str(J),',',int2str(round(v_('1')))];
            iv = ix_(['h',int2str(J),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['MXAR',int2str(J)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['MXAR',int2str(J)]);
            end
            [ig,ig_] = s2mpjlib('ii',['RD',int2str(J),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['RD',int2str(J),',',int2str(round(v_('1')))];
            iv = ix_(['qr',int2str(J),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['MST',int2str(round(v_('1')))])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['MST',int2str(round(v_('1')))]);
            end
        end
        for I=v_('2'):v_('NSTEP')
            for J=v_('SuRES'):v_('EuRES')
                v_('A') = -1+I;
                [ig,ig_] = s2mpjlib('ii',['RD',int2str(J),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['RD',int2str(J),',',int2str(I)];
                iv = ix_(['h',int2str(J),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['MXAR',int2str(J)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['MXAR',int2str(J)]);
                end
                iv = ix_(['h',int2str(J),',',int2str(round(v_('A')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['XAR',int2str(J)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['XAR',int2str(J)]);
                end
                iv = ix_(['qr',int2str(J),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['MST',int2str(I)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['MST',int2str(I)]);
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
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('ONE'):v_('NSTEP')
            for J=v_('ONE'):v_('NND')
                pbm.gconst(ig_(['ND',int2str(J),',',int2str(I)])) = v_(['D',int2str(J)]);
            end
        end
        for J=v_('SuRES'):v_('EuRES')
            pbm.gconst(ig_(['RD',int2str(J),',',int2str(round(v_('1')))])) =...
                  v_(['MSTVOL',int2str(J)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('ONE'):v_('NSTEP')
            for J=v_('SuRES'):v_('EuRES')
                pb.xlower(ix_(['h',int2str(J),',',int2str(I)]),1) = v_(['HGT',int2str(J)]);
                pb.xupper(ix_(['h',int2str(J),',',int2str(I)])) = v_(['MXHGT',int2str(J)]);
                pb.xlower(ix_(['qr',int2str(J),',',int2str(I)])) = -Inf;
                pb.xupper(ix_(['qr',int2str(J),',',int2str(I)]),1) = +Inf;
            end
            for J=v_('EuRES+1'):v_('NND')
                pb.xlower(ix_(['h',int2str(J),',',int2str(I)]),1) = 0.0;
            end
            for J=v_('SuTW'):v_('EuTW')
                pb.xlower(ix_(['qs',int2str(J),',',int2str(I)]),1) = 0.0;
                pb.xupper(ix_(['qs',int2str(J),',',int2str(I)])) = v_(['QSMX',int2str(J)]);
            end
            for J=v_('ONE'):v_('NEL')
                pb.xlower(ix_(['q',int2str(J),',',int2str(I)])) = -Inf;
                pb.xupper(ix_(['q',int2str(J),',',int2str(I)]),1) = +Inf;
            end
            for J=v_('SuPMP'):v_('EuPMP')
                pb.xlower(ix_(['u',int2str(J),',',int2str(I)]),1) = 0.0;
                pb.xupper(ix_(['u',int2str(J),',',int2str(I)])) = 7.0;
                pb.xlower(ix_(['q',int2str(J),',',int2str(I)]),1) = 0.0;
            end
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        pb.y0 = 0.0*ones(pb.m,1);
        for I=v_('1'):v_('NSTEP')
            for J=v_('EuRES+1'):v_('NND')
                pb.x0(ix_(['h',int2str(J),',',int2str(I)]),1) = 0.0;
            end
            for J=v_('ONE'):v_('NEL')
                pb.x0(ix_(['q',int2str(J),',',int2str(I)]),1) = 20.0;
            end
            for J=v_('SuPMP'):v_('EuPMP')
                pb.x0(ix_(['u',int2str(J),',',int2str(I)]),1) = 3.0;
            end
            pb.x0(ix_(['qs',int2str(round(v_('1'))),',',int2str(I)]),1) = 25.0;
            pb.x0(ix_(['qs',int2str(round(v_('2'))),',',int2str(I)]),1) = 50.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXZ',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eXPOW',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXMYZSQ',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'ePMP1',iet_);
        elftv{it}{1} = 'Q';
        [it,iet_] = s2mpjlib( 'ii', 'ePMP2',iet_);
        elftv{it}{1} = 'Q';
        [it,iet_] = s2mpjlib( 'ii', 'ePMP3',iet_);
        elftv{it}{1} = 'Q';
        [it,iet_] = s2mpjlib( 'ii', 'ePMP4',iet_);
        elftv{it}{1} = 'Q';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NSTEP')
            for J=v_('SuPIPE'):v_('EuPIPE')
                ename = ['EA',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXPOW';
                ielftype(ie) = iet_('eXPOW');
                vname = ['q',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for J=v_('SuVAL'):v_('EuVAL')
                ename = ['EA',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXPOW';
                ielftype(ie) = iet_('eXPOW');
                vname = ['q',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for J=v_('SuPMP'):v_('EuPMP')
                ename = ['EA',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXSQ';
                ielftype(ie) = iet_('eXSQ');
                vname = ['q',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['EB',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXZ';
                ielftype(ie) = iet_('eXZ');
                vname = ['q',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['u',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['EC',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXSQ';
                ielftype(ie) = iet_('eXSQ');
                vname = ['u',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['EH',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXMYZSQ';
                ielftype(ie) = iet_('eXMYZSQ');
                vname = ['u',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            ename = ['EH',int2str(round(v_('11'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('7'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('11'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('9'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('12'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('4'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('12'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('6'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('13'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('1'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('13'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('5'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('14'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('2'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('14'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('11'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('15'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('2'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(round(v_('15'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['h',int2str(round(v_('11'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('NSTEP')
            ename = ['PPW',int2str(round(v_('11'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePMP4';
            ielftype(ie) = iet_('ePMP4');
            ename = ['PPW',int2str(round(v_('12'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePMP3';
            ielftype(ie) = iet_('ePMP3');
            ename = ['PPW',int2str(round(v_('13'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePMP1';
            ielftype(ie) = iet_('ePMP1');
            ename = ['PPW',int2str(round(v_('14'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePMP2';
            ielftype(ie) = iet_('ePMP2');
            ename = ['PPW',int2str(round(v_('15'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePMP2';
            ielftype(ie) = iet_('ePMP2');
            for J=v_('SuPMP'):v_('EuPMP')
                ename = ['PPW',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                vname = ['q',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('Q',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NSTEP')
            for J=v_('SuPIPE'):v_('EuPIPE')
                ig = ig_(['EL',int2str(J),',',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EA',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['MGINV',int2str(J)]);
            end
            for J=v_('SuVAL'):v_('EuVAL')
                ig = ig_(['EL',int2str(J),',',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EA',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['MGINV',int2str(J)]);
            end
            for J=v_('SuPMP'):v_('EuPMP')
                ig = ig_(['EL',int2str(J),',',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EA',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['CONA',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EB',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['CONB',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EC',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['CONC',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['EH',int2str(J),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.0;
            end
        end
        for J=v_('SuPMP'):v_('EuPMP')
            ig = ig_(['PC',int2str(J),',',int2str(round(v_('2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PPW',int2str(J),',',int2str(round(v_('2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 100.0;
            ig = ig_(['PC',int2str(J),',',int2str(round(v_('3')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PPW',int2str(J),',',int2str(round(v_('3')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 100.0;
            ig = ig_(['PC',int2str(J),',',int2str(round(v_('4')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PPW',int2str(J),',',int2str(round(v_('4')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 100.0;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'SOR2-AY-312-256';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eXSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eXZ'

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

    case 'eXMYZSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
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

    case 'eXPOW'

        EV_  = varargin{1};
        iel_ = varargin{2};
        MODX = abs(EV_(1));
        ISNEG = EV_(1)<0.0;
        ISPOS = EV_(1)>=0.0;
        if(ISNEG)
            SGN = -1.0;
        end
        if(ISPOS)
            SGN = +1.0;
        end
        varargout{1} = SGN*MODX^1.852;
        if(nargout>1)
            g_(1,1) = 1.852*MODX^0.852;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = SGN*1.577904*MODX^(-0.148);
                varargout{3} = H_;
            end
        end

    case 'ePMP1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.074*EV_(1)*EV_(1)+3.062*EV_(1)+50.357;
        if(nargout>1)
            g_(1,1) = 0.148*EV_(1)+3.062;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 0.148;
                varargout{3} = H_;
            end
        end

    case 'ePMP2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.747*EV_(1)*EV_(1)-10.287*EV_(1);
        if(nargout>1)
            g_(1,1) = 1.494*EV_(1)-10.287;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.494;
                varargout{3} = H_;
            end
        end

    case 'ePMP3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.034*EV_(1)*EV_(1)+0.220*EV_(1)+6.685;
        if(nargout>1)
            g_(1,1) = 0.068*EV_(1)+0.220;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 0.068;
                varargout{3} = H_;
            end
        end

    case 'ePMP4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.079*EV_(1)*EV_(1)-2.761*EV_(1)+35.014;
        if(nargout>1)
            g_(1,1) = 0.158*EV_(1)-2.761;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 0.158;
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

