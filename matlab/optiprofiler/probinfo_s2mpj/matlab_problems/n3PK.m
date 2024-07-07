function varargout = n3PK(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : n3PK
%    *********
% 
%    A problem arising in the estimation of structured O/D matrix
% 
%    Source:  
%    M. Bierlaire, private communication
%    see also
%    M. Bierlaire and Ph. L. Toint,
%    MEUSE: an origin-destination estimator that exploits structure,
%    Transportation Research B, 29, 1, 47--60, 1995.
% 
%    SIF input: Ph. Toint, Dec 1989, Corrected July 1993.
% 
%    classification = 'SBR2-MN-30-0'
% 
%  Parameters
% 
%  Number of parking columns
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'n3PK';

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
        v_('NPKC') = 3;
        v_('NPKC-1') = -1+v_('NPKC');
        v_('NPKC+1') = 1+v_('NPKC');
        v_('NCENT') = 6;
        v_('NCENT-1') = -1+v_('NCENT');
        v_('RNCENT-1') = v_('NCENT-1');
        v_('GAMMA') = 1.0000e+04;
        v_('FT0') = 0.500000;
        v_('FT1') = 0.500000;
        v_('FT2') = 0.500000;
        v_('WFT0') = 1.000000;
        v_('WFT1') = 1.000000;
        v_('WFT2') = 1.000000;
        v_('COUNT') = 9;
        v_('COUNT-1') = -1+v_('COUNT');
        v_('DEFW') = 999.999953;
        v_('0') = 0;
        v_('1') = 1;
        v_('COU0') = 910.000000;
        v_('COU1') = 175.000000;
        v_('COU2') = 1915.000000;
        v_('COU3') = 450.000000;
        v_('COU4') = 260.000000;
        v_('COU5') = 80.000000;
        v_('COU6') = 670.000000;
        v_('COU7') = 1450.000000;
        v_('COU8') = 990.000000;
        v_('PHI0') = 1.000000;
        v_('PHI1') = 1.000000;
        v_('PHI2') = 1.000000;
        v_('PHI3') = 1.000000;
        v_('PHI4') = 1.000000;
        v_('PHI5') = 1.000000;
        v_('PHI6') = 1.000000;
        v_('PHI7') = 1.000000;
        v_('PHI8') = 1.000000;
        for I=v_('0'):v_('COUNT-1')
            v_(['PHI',int2str(I)]) = v_(['PHI',int2str(I)])/v_('GAMMA');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','A1,0',ix_);
        pb.xnames{iv} = 'A1,0';
        [iv,ix_] = s2mpjlib('ii','A2,0',ix_);
        pb.xnames{iv} = 'A2,0';
        [iv,ix_] = s2mpjlib('ii','A3,0',ix_);
        pb.xnames{iv} = 'A3,0';
        [iv,ix_] = s2mpjlib('ii','A4,0',ix_);
        pb.xnames{iv} = 'A4,0';
        [iv,ix_] = s2mpjlib('ii','A5,0',ix_);
        pb.xnames{iv} = 'A5,0';
        [iv,ix_] = s2mpjlib('ii','A0,1',ix_);
        pb.xnames{iv} = 'A0,1';
        [iv,ix_] = s2mpjlib('ii','A2,1',ix_);
        pb.xnames{iv} = 'A2,1';
        [iv,ix_] = s2mpjlib('ii','A3,1',ix_);
        pb.xnames{iv} = 'A3,1';
        [iv,ix_] = s2mpjlib('ii','A4,1',ix_);
        pb.xnames{iv} = 'A4,1';
        [iv,ix_] = s2mpjlib('ii','A5,1',ix_);
        pb.xnames{iv} = 'A5,1';
        [iv,ix_] = s2mpjlib('ii','A0,2',ix_);
        pb.xnames{iv} = 'A0,2';
        [iv,ix_] = s2mpjlib('ii','A1,2',ix_);
        pb.xnames{iv} = 'A1,2';
        [iv,ix_] = s2mpjlib('ii','A3,2',ix_);
        pb.xnames{iv} = 'A3,2';
        [iv,ix_] = s2mpjlib('ii','A4,2',ix_);
        pb.xnames{iv} = 'A4,2';
        [iv,ix_] = s2mpjlib('ii','A5,2',ix_);
        pb.xnames{iv} = 'A5,2';
        for J=v_('NPKC'):v_('NCENT-1')
            v_('J+1') = 1+J;
            v_('J-1') = -1+J;
            for I=v_('0'):v_('J-1')
                [iv,ix_] = s2mpjlib('ii',['T',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['T',int2str(I),',',int2str(J)];
            end
            for I=v_('J+1'):v_('NCENT-1')
                [iv,ix_] = s2mpjlib('ii',['T',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['T',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','G0,3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.010000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.010000;
        end
        [ig,ig_] = s2mpjlib('ii','G1,3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.007143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.007143;
        end
        [ig,ig_] = s2mpjlib('ii','G2,3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.008333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.008333;
        end
        [ig,ig_] = s2mpjlib('ii','G4,3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T4,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.050000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.050000;
        end
        [ig,ig_] = s2mpjlib('ii','G5,3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T5,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.050000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.050000;
        end
        [ig,ig_] = s2mpjlib('ii','G0,4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.005000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.005000;
        end
        [ig,ig_] = s2mpjlib('ii','G1,4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.005556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.005556;
        end
        [ig,ig_] = s2mpjlib('ii','G2,4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.050000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.050000;
        end
        [ig,ig_] = s2mpjlib('ii','G3,4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.001667+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.001667;
        end
        [ig,ig_] = s2mpjlib('ii','G5,4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T5,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.025000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.025000;
        end
        [ig,ig_] = s2mpjlib('ii','G0,5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.020000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.020000;
        end
        [ig,ig_] = s2mpjlib('ii','G1,5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.033333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.033333;
        end
        [ig,ig_] = s2mpjlib('ii','G2,5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.014286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.014286;
        end
        [ig,ig_] = s2mpjlib('ii','G3,5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.006667+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.006667;
        end
        [ig,ig_] = s2mpjlib('ii','G4,5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T4,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.050000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.050000;
        end
        v_('TMP') = 5.000000*v_('FT0');
        v_('TMP1') = 1.0/v_('TMP');
        [ig,ig_] = s2mpjlib('ii','H0',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('WFT0');
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        v_('TMP') = 5.000000*v_('FT1');
        v_('TMP1') = 1.0/v_('TMP');
        [ig,ig_] = s2mpjlib('ii','H1',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('WFT1');
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        v_('TMP') = 5.000000*v_('FT2');
        v_('TMP1') = 1.0/v_('TMP');
        [ig,ig_] = s2mpjlib('ii','H2',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('WFT2');
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP1');
        end
        for I=v_('0'):v_('COUNT-1')
            [ig,ig_] = s2mpjlib('ii',['K',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_(['PHI',int2str(I)]);
        end
        v_('TMP1') = 200.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 200.000000/v_('COU4');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K4',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 200.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 480.000000/v_('COU8');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 480.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 480.000000/v_('COU6');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 480.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 120.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 360.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 360.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 560.000000/v_('COU8');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 560.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 560.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 240.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 400.000000/v_('COU8');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 400.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 400.000000/v_('COU6');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 400.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 400.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 420.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 420.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 180.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 180.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 180.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 320.000000/v_('COU8');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 320.000000/v_('COU7');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 320.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 320.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 20.000000/v_('COU1');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 20.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 60.000000/v_('COU1');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 40.000000/v_('COU2');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 40.000000/v_('COU1');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 40.000000/v_('COU0');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 120.000000/v_('COU5');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K5',ig_);
        gtype{ig} = '<>';
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 20.000000/v_('COU8');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP1') = 20.000000/v_('COU5');
        v_('TMP') = 1.000000*v_('TMP1');
        [ig,ig_] = s2mpjlib('ii','K5',ig_);
        gtype{ig} = '<>';
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU7');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU3');
        [ig,ig_] = s2mpjlib('ii','K3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU7');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU4');
        [ig,ig_] = s2mpjlib('ii','K4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU8');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU7');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU7');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('T4,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU8');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('T5,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU7');
        [ig,ig_] = s2mpjlib('ii','K7',ig_);
        gtype{ig} = '<>';
        iv = ix_('T5,3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU3');
        [ig,ig_] = s2mpjlib('ii','K3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU4');
        [ig,ig_] = s2mpjlib('ii','K4',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU8');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU3');
        [ig,ig_] = s2mpjlib('ii','K3',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU2');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU8');
        [ig,ig_] = s2mpjlib('ii','K8',ig_);
        gtype{ig} = '<>';
        iv = ix_('T5,4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU0');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('T0,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('T1,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T2,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        iv = ix_('T3,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU2');
        [ig,ig_] = s2mpjlib('ii','K2',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU1');
        [ig,ig_] = s2mpjlib('ii','K1',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU0');
        [ig,ig_] = s2mpjlib('ii','K0',ig_);
        gtype{ig} = '<>';
        iv = ix_('T3,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU6');
        [ig,ig_] = s2mpjlib('ii','K6',ig_);
        gtype{ig} = '<>';
        iv = ix_('T4,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        v_('TMP') = 1.000000/v_('COU5');
        [ig,ig_] = s2mpjlib('ii','K5',ig_);
        gtype{ig} = '<>';
        iv = ix_('T4,5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('TMP')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('TMP');
        end
        [ig,ig_] = s2mpjlib('ii','L1,0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L2,0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L3,0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L4,0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L5,0',ig_);
        gtype{ig} = '<>';
        iv = ix_('A1,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,0');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        [ig,ig_] = s2mpjlib('ii','L0,1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L2,1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L3,1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L4,1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L5,1',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A2,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        [ig,ig_] = s2mpjlib('ii','L0,2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L1,2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L3,2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L4,2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        [ig,ig_] = s2mpjlib('ii','L5,2',ig_);
        gtype{ig} = '<>';
        iv = ix_('A0,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A1,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A3,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A4,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.200000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.200000;
        end
        iv = ix_('A5,2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.800000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.800000;
        end
        pbm.gscale(ig,1) = 0.500000;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for J=v_('NPKC'):v_('NCENT-1')
            v_('J+1') = 1+J;
            v_('J-1') = -1+J;
            for I=v_('0'):v_('J-1')
                pbm.gconst(ig_(['G',int2str(I),',',int2str(J)])) = 1.0;
            end
            for I=v_('J+1'):v_('NCENT-1')
                pbm.gconst(ig_(['G',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        for J=v_('0'):v_('NPKC-1')
            pbm.gconst(ig_(['H',int2str(J)])) = 1.0;
        end
        for J=v_('0'):v_('COUNT-1')
            pbm.gconst(ig_(['K',int2str(J)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('A1,0'),1) = v_('FT0');
        pb.x0(ix_('A2,0'),1) = v_('FT0');
        pb.x0(ix_('A3,0'),1) = v_('FT0');
        pb.x0(ix_('A4,0'),1) = v_('FT0');
        pb.x0(ix_('A5,0'),1) = v_('FT0');
        pb.x0(ix_('A0,1'),1) = v_('FT1');
        pb.x0(ix_('A2,1'),1) = v_('FT1');
        pb.x0(ix_('A3,1'),1) = v_('FT1');
        pb.x0(ix_('A4,1'),1) = v_('FT1');
        pb.x0(ix_('A5,1'),1) = v_('FT1');
        pb.x0(ix_('A0,2'),1) = v_('FT2');
        pb.x0(ix_('A1,2'),1) = v_('FT2');
        pb.x0(ix_('A3,2'),1) = v_('FT2');
        pb.x0(ix_('A4,2'),1) = v_('FT2');
        pb.x0(ix_('A5,2'),1) = v_('FT2');
        pb.x0(ix_('T0,3'),1) = 100.000000;
        pb.x0(ix_('T1,3'),1) = 140.000000;
        pb.x0(ix_('T2,3'),1) = 120.000000;
        pb.x0(ix_('T4,3'),1) = 20.000000;
        pb.x0(ix_('T5,3'),1) = 20.000000;
        pb.x0(ix_('T0,4'),1) = 200.000000;
        pb.x0(ix_('T1,4'),1) = 180.000000;
        pb.x0(ix_('T2,4'),1) = 20.000000;
        pb.x0(ix_('T3,4'),1) = 600.000000;
        pb.x0(ix_('T5,4'),1) = 40.000000;
        pb.x0(ix_('T0,5'),1) = 50.000000;
        pb.x0(ix_('T1,5'),1) = 30.000000;
        pb.x0(ix_('T2,5'),1) = 70.000000;
        pb.x0(ix_('T3,5'),1) = 150.000000;
        pb.x0(ix_('T4,5'),1) = 20.000000;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gSQUARE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-MN-30-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
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

