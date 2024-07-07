function varargout = CmRELOAD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CmRELOAD
%    *********
% 
%    Source: Nuclear Reactor Core Reload Pattern Optimization
%    A.J. Quist et.al., draft paper, September 1997.
%    (2nd data set implemented here)
%    SIF input: S. Leyffer, November 1997
% 
%    classification = 'LOR2-MN-342-284'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CmRELOAD';

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
        v_('N') = 14;
        v_('M') = 3;
        v_('L') = 4;
        v_('T') = 6;
        v_('D11') = 1;
        v_('D12') = 7;
        v_('D21') = 11;
        v_('D22') = 14;
        v_('KFRESH') = 1.25;
        v_('FLIM') = 1.8;
        v_('KEFFuINI') = 0.956145;
        v_('ALPHA') = 6E-6;
        v_('CONSPW') = 364.0;
        v_('CYTIME') = 350.0;
        v_('TT') = v_('T');
        v_('T-1') = -1.0+v_('TT');
        v_('DELTAT') = v_('CYTIME')/v_('T-1');
        v_('ACC') = v_('ALPHA')*v_('CONSPW');
        v_('ACC') = v_('ACC')*v_('DELTAT');
        v_('-ACC') = -1.0*v_('ACC');
        for I=v_('1'):v_('N')
            v_(['V',int2str(I)]) = 1.0;
        end
        v_(['V',int2str(round(v_('D11')))]) = 0.5;
        v_(['V',int2str(round(v_('D12')))]) = 0.5;
        v_(['V',int2str(round(v_('D21')))]) = 0.5;
        v_(['V',int2str(round(v_('D22')))]) = 0.5;
        v_('NROW1') = 2;
        v_(['ROWu',int2str(round(v_('1'))),',',int2str(round(v_('1')))]) = 1;
        v_(['ROWu',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 2;
        v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('1')))]) = 0.705;
        v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 0.25;
        v_('NROW2') = 5;
        v_(['ROWu',int2str(round(v_('2'))),',',int2str(round(v_('1')))]) = 1;
        v_(['ROWu',int2str(round(v_('2'))),',',int2str(round(v_('2')))]) = 2;
        v_(['ROWu',int2str(round(v_('2'))),',',int2str(round(v_('3')))]) = 3;
        v_(['ROWu',int2str(round(v_('2'))),',',int2str(round(v_('4')))]) = 7;
        v_(['ROWu',int2str(round(v_('2'))),',',int2str(round(v_('5')))]) = 8;
        v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('1')))]) = 0.125;
        v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('2')))]) = 0.625;
        v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('3')))]) = 0.125;
        v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('7')))]) = 0.08;
        v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('8')))]) = 0.045;
        v_('NROW3') = 6;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('1')))]) = 2;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('2')))]) = 3;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('3')))]) = 4;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('4')))]) = 7;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('5')))]) = 8;
        v_(['ROWu',int2str(round(v_('3'))),',',int2str(round(v_('6')))]) = 9;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('2')))]) = 0.125;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('3')))]) = 0.58;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('4')))]) = 0.125;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('7')))]) = 0.045;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('8')))]) = 0.08;
        v_(['G',int2str(round(v_('3'))),',',int2str(round(v_('9')))]) = 0.045;
        v_('NROW4') = 6;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('1')))]) = 3;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('2')))]) = 4;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('3')))]) = 5;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('4')))]) = 8;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('5')))]) = 9;
        v_(['ROWu',int2str(round(v_('4'))),',',int2str(round(v_('6')))]) = 10;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('3')))]) = 0.125;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('4')))]) = 0.58;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('5')))]) = 0.125;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('8')))]) = 0.045;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('9')))]) = 0.08;
        v_(['G',int2str(round(v_('4'))),',',int2str(round(v_('10')))]) = 0.045;
        v_('NROW5') = 5;
        v_(['ROWu',int2str(round(v_('5'))),',',int2str(round(v_('1')))]) = 4;
        v_(['ROWu',int2str(round(v_('5'))),',',int2str(round(v_('2')))]) = 5;
        v_(['ROWu',int2str(round(v_('5'))),',',int2str(round(v_('3')))]) = 6;
        v_(['ROWu',int2str(round(v_('5'))),',',int2str(round(v_('4')))]) = 9;
        v_(['ROWu',int2str(round(v_('5'))),',',int2str(round(v_('5')))]) = 10;
        v_(['G',int2str(round(v_('5'))),',',int2str(round(v_('4')))]) = 0.125;
        v_(['G',int2str(round(v_('5'))),',',int2str(round(v_('5')))]) = 0.58;
        v_(['G',int2str(round(v_('5'))),',',int2str(round(v_('6')))]) = 0.125;
        v_(['G',int2str(round(v_('5'))),',',int2str(round(v_('9')))]) = 0.045;
        v_(['G',int2str(round(v_('5'))),',',int2str(round(v_('10')))]) = 0.08;
        v_('NROW6') = 3;
        v_(['ROWu',int2str(round(v_('6'))),',',int2str(round(v_('1')))]) = 5;
        v_(['ROWu',int2str(round(v_('6'))),',',int2str(round(v_('2')))]) = 6;
        v_(['ROWu',int2str(round(v_('6'))),',',int2str(round(v_('3')))]) = 10;
        v_(['G',int2str(round(v_('6'))),',',int2str(round(v_('5')))]) = 0.125;
        v_(['G',int2str(round(v_('6'))),',',int2str(round(v_('6')))]) = 0.58;
        v_(['G',int2str(round(v_('6'))),',',int2str(round(v_('10')))]) = 0.045;
        v_('NROW7') = 5;
        v_(['ROWu',int2str(round(v_('7'))),',',int2str(round(v_('1')))]) = 1;
        v_(['ROWu',int2str(round(v_('7'))),',',int2str(round(v_('2')))]) = 2;
        v_(['ROWu',int2str(round(v_('7'))),',',int2str(round(v_('3')))]) = 7;
        v_(['ROWu',int2str(round(v_('7'))),',',int2str(round(v_('4')))]) = 8;
        v_(['ROWu',int2str(round(v_('7'))),',',int2str(round(v_('5')))]) = 11;
        v_(['G',int2str(round(v_('7'))),',',int2str(round(v_('1')))]) = 0.045;
        v_(['G',int2str(round(v_('7'))),',',int2str(round(v_('2')))]) = 0.16;
        v_(['G',int2str(round(v_('7'))),',',int2str(round(v_('7')))]) = 0.5;
        v_(['G',int2str(round(v_('7'))),',',int2str(round(v_('8')))]) = 0.16;
        v_(['G',int2str(round(v_('7'))),',',int2str(round(v_('11')))]) = 0.045;
        v_('NROW8') = 8;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('1')))]) = 2;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('2')))]) = 3;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('3')))]) = 4;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('4')))]) = 7;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('5')))]) = 8;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('6')))]) = 9;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('7')))]) = 11;
        v_(['ROWu',int2str(round(v_('8'))),',',int2str(round(v_('8')))]) = 12;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('2')))]) = 0.045;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('3')))]) = 0.08;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('4')))]) = 0.045;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('7')))]) = 0.08;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('8')))]) = 0.545;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('9')))]) = 0.08;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('11')))]) = 0.08;
        v_(['G',int2str(round(v_('8'))),',',int2str(round(v_('12')))]) = 0.045;
        v_('NROW9') = 9;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('1')))]) = 3;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('2')))]) = 4;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('3')))]) = 5;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('4')))]) = 8;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('5')))]) = 9;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('6')))]) = 10;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('7')))]) = 11;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('8')))]) = 12;
        v_(['ROWu',int2str(round(v_('9'))),',',int2str(round(v_('9')))]) = 13;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('3')))]) = 0.045;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('4')))]) = 0.08;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('5')))]) = 0.045;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('8')))]) = 0.08;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('9')))]) = 0.5;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('10')))]) = 0.08;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('11')))]) = 0.045;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('12')))]) = 0.08;
        v_(['G',int2str(round(v_('9'))),',',int2str(round(v_('13')))]) = 0.045;
        v_('NROW10') = 7;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('1')))]) = 4;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('2')))]) = 5;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('3')))]) = 6;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('4')))]) = 9;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('5')))]) = 10;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('6')))]) = 12;
        v_(['ROWu',int2str(round(v_('10'))),',',int2str(round(v_('7')))]) = 13;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('4')))]) = 0.045;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('5')))]) = 0.08;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('6')))]) = 0.045;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('9')))]) = 0.08;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('10')))]) = 0.5;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('12')))]) = 0.045;
        v_(['G',int2str(round(v_('10'))),',',int2str(round(v_('13')))]) = 0.08;
        v_('NROW11') = 6;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('1')))]) = 7;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('2')))]) = 8;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('3')))]) = 9;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('4')))]) = 11;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('5')))]) = 12;
        v_(['ROWu',int2str(round(v_('11'))),',',int2str(round(v_('6')))]) = 14;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('7')))]) = 0.045;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('8')))]) = 0.125;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('9')))]) = 0.045;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('11')))]) = 0.5;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('12')))]) = 0.16;
        v_(['G',int2str(round(v_('11'))),',',int2str(round(v_('14')))]) = 0.045;
        v_('NROW12') = 7;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('1')))]) = 8;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('2')))]) = 9;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('3')))]) = 10;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('4')))]) = 11;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('5')))]) = 12;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('6')))]) = 13;
        v_(['ROWu',int2str(round(v_('12'))),',',int2str(round(v_('7')))]) = 14;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('8')))]) = 0.045;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('9')))]) = 0.08;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('10')))]) = 0.045;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('11')))]) = 0.08;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('12')))]) = 0.545;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('13')))]) = 0.08;
        v_(['G',int2str(round(v_('12'))),',',int2str(round(v_('14')))]) = 0.08;
        v_('NROW13') = 5;
        v_(['ROWu',int2str(round(v_('13'))),',',int2str(round(v_('1')))]) = 9;
        v_(['ROWu',int2str(round(v_('13'))),',',int2str(round(v_('2')))]) = 10;
        v_(['ROWu',int2str(round(v_('13'))),',',int2str(round(v_('3')))]) = 12;
        v_(['ROWu',int2str(round(v_('13'))),',',int2str(round(v_('4')))]) = 13;
        v_(['ROWu',int2str(round(v_('13'))),',',int2str(round(v_('5')))]) = 14;
        v_(['G',int2str(round(v_('13'))),',',int2str(round(v_('9')))]) = 0.045;
        v_(['G',int2str(round(v_('13'))),',',int2str(round(v_('10')))]) = 0.08;
        v_(['G',int2str(round(v_('13'))),',',int2str(round(v_('12')))]) = 0.08;
        v_(['G',int2str(round(v_('13'))),',',int2str(round(v_('13')))]) = 0.5;
        v_(['G',int2str(round(v_('13'))),',',int2str(round(v_('14')))]) = 0.045;
        v_('NROW14') = 3;
        v_(['ROWu',int2str(round(v_('14'))),',',int2str(round(v_('1')))]) = 11;
        v_(['ROWu',int2str(round(v_('14'))),',',int2str(round(v_('2')))]) = 12;
        v_(['ROWu',int2str(round(v_('14'))),',',int2str(round(v_('3')))]) = 14;
        v_(['G',int2str(round(v_('14'))),',',int2str(round(v_('11')))]) = 0.045;
        v_(['G',int2str(round(v_('14'))),',',int2str(round(v_('12')))]) = 0.125;
        v_(['G',int2str(round(v_('14'))),',',int2str(round(v_('14')))]) = 0.5;
        v_('T-1') = -1+v_('T');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            for K=v_('1'):v_('L')
                for J=v_('1'):v_('M')
                    [iv,ix_] =...
                          s2mpjlib('ii',['X',int2str(I),',',int2str(K),',',int2str(J)],ix_);
                    pb.xnames{iv} = ['X',int2str(I),',',int2str(K),',',int2str(J)];
                end
            end
        end
        for S=v_('1'):v_('T')
            for I=v_('1'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['KINF',int2str(I),',',int2str(S)],ix_);
                pb.xnames{iv} = ['KINF',int2str(I),',',int2str(S)];
                [iv,ix_] = s2mpjlib('ii',['PHI',int2str(I),',',int2str(S)],ix_);
                pb.xnames{iv} = ['PHI',int2str(I),',',int2str(S)];
            end
            [iv,ix_] = s2mpjlib('ii',['KEFF',int2str(S)],ix_);
            pb.xnames{iv} = ['KEFF',int2str(S)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_(['KEFF',int2str(round(v_('T')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for K=v_('1'):v_('L')
            for J=v_('1'):v_('M')
                for I=v_('1'):v_('N')
                    [ig,ig_] = s2mpjlib('ii',['SUMI',int2str(K),',',int2str(J)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['SUMI',int2str(K),',',int2str(J)];
                    iv = ix_(['X',int2str(I),',',int2str(K),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_(['V',int2str(I)])+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_(['V',int2str(I)]);
                    end
                end
            end
        end
        for I=v_('1'):v_('N')
            for K=v_('1'):v_('L')
                for J=v_('1'):v_('M')
                    [ig,ig_] = s2mpjlib('ii',['SUMLM',int2str(I)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['SUMLM',int2str(I)];
                    iv = ix_(['X',int2str(I),',',int2str(K),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1.0;
                    end
                end
            end
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['PLAC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PLAC',int2str(I)];
            iv = ix_(['KINF',int2str(I),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['PLAC',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['PLAC',int2str(I)];
                iv = ix_(['X',int2str(I),',',int2str(round(v_('1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KFRESH')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KFRESH');
                end
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                [ig,ig_] = s2mpjlib('ii',['KERN',int2str(I),',',int2str(S)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['KERN',int2str(I),',',int2str(S)];
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T-1')
                v_('R') = 1+S;
                [ig,ig_] = s2mpjlib('ii',['KINFF',int2str(I),',',int2str(S)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['KINFF',int2str(I),',',int2str(S)];
                iv = ix_(['KINF',int2str(I),',',int2str(round(v_('R')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['KINF',int2str(I),',',int2str(S)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for S=v_('1'):v_('T')
            [ig,ig_] = s2mpjlib('ii',['CPOW',int2str(S)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CPOW',int2str(S)];
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                [ig,ig_] = s2mpjlib('ii',['PEAK',int2str(I),',',int2str(S)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['PEAK',int2str(I),',',int2str(S)];
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
        for K=v_('1'):v_('L')
            for J=v_('1'):v_('M')
                pbm.gconst(ig_(['SUMI',int2str(K),',',int2str(J)])) = 1.0;
            end
        end
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['SUMLM',int2str(I)])) = 1.0;
        end
        for S=v_('1'):v_('T')
            pbm.gconst(ig_(['CPOW',int2str(S)])) = 1.0;
        end
        v_('TEMP') = 0.0;
        for I=v_('1'):v_('N')
            v_('TEMP') = v_('TEMP')+v_(['V',int2str(I)]);
        end
        v_('TEMP') = v_('FLIM')/v_('TEMP');
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                pbm.gconst(ig_(['PEAK',int2str(I),',',int2str(S)])) = v_('TEMP');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            for K=v_('1'):v_('L')
                for J=v_('1'):v_('M')
                    pb.xupper(ix_(['X',int2str(I),',',int2str(K),',',int2str(J)])) = 1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                pb.xupper(ix_(['KINF',int2str(I),',',int2str(S)])) = v_('KFRESH');
            end
        end
        v_('LOuKEFF') = 0.9;
        v_('UPuKEFF') = 1.5;
        v_('TEMP') = v_('T');
        v_('TEMP') = -0.015*v_('TEMP');
        v_('LOuKEFF') = v_('LOuKEFF')+v_('TEMP');
        v_('UPuKEFF') = v_('UPuKEFF')+v_('TEMP');
        for S=v_('1'):v_('T')
            pb.xlower(ix_(['KEFF',int2str(S)]),1) = v_('LOuKEFF');
            pb.xupper(ix_(['KEFF',int2str(S)])) = v_('UPuKEFF');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('R14') = 14.0;
        v_('TEMP') = v_('KFRESH')/v_('R14');
        for S=v_('1'):v_('T')
            if(isKey(ix_,['KEFF',int2str(S)]))
                pb.x0(ix_(['KEFF',int2str(S)]),1) = v_('KEFFuINI');
            else
                pb.y0(find(pbm.congrps==ig_(['KEFF',int2str(S)])),1) = v_('KEFFuINI');
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                if(isKey(ix_,['KINF',int2str(I),',',int2str(S)]))
                    pb.x0(ix_(['KINF',int2str(I),',',int2str(S)]),1) = v_('KFRESH');
                else
                    pb.y0(find(pbm.congrps==ig_(['KINF',int2str(I),',',int2str(S)])),1) =...
                          v_('KFRESH');
                end
                if(isKey(ix_,['PHI',int2str(I),',',int2str(S)]))
                    pb.x0(ix_(['PHI',int2str(I),',',int2str(S)]),1) = v_('TEMP');
                else
                    pb.y0(find(pbm.congrps==ig_(['PHI',int2str(I),',',int2str(S)])),1) =...
                          v_('TEMP');
                end
            end
        end
        for I=v_('1'):v_('N')
            for K=v_('1'):v_('L')
                for J=v_('1'):v_('M')
                    if(isKey(ix_,['X',int2str(I),',',int2str(K),',',int2str(J)]))
                        pb.x0(ix_(['X',int2str(I),',',int2str(K),',',int2str(J)]),1) = 0.5;
                    else
                        pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(K),',',int2str(J)])),1) = 0.5;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'en3PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        v_('K') = 2;
        v_('K1') = -1+v_('K');
        for J=v_('1'):v_('M')
            for I=v_('1'):v_('N')
                for II=v_('1'):v_('N')
                    ename = ['Au',int2str(I),',',int2str(II),',',int2str(J)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en3PROD';
                    ielftype(ie) = iet_('en3PROD');
                    vname = ['X',int2str(I),',',int2str(round(v_('K'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(II),',',int2str(round(v_('K1'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['KINF',int2str(II),',',int2str(round(v_('T')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        v_('K') = 3;
        v_('K1') = -1+v_('K');
        for J=v_('1'):v_('M')
            for I=v_('1'):v_('N')
                for II=v_('1'):v_('N')
                    ename = ['Bu',int2str(I),',',int2str(II),',',int2str(J)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en3PROD';
                    ielftype(ie) = iet_('en3PROD');
                    vname = ['X',int2str(I),',',int2str(round(v_('K'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(II),',',int2str(round(v_('K1'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['KINF',int2str(II),',',int2str(round(v_('T')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        v_('K') = 4;
        v_('K1') = -1+v_('K');
        for J=v_('1'):v_('M')
            for I=v_('1'):v_('N')
                for II=v_('1'):v_('N')
                    ename = ['Cu',int2str(I),',',int2str(II),',',int2str(J)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en3PROD';
                    ielftype(ie) = iet_('en3PROD');
                    vname = ['X',int2str(I),',',int2str(round(v_('K'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(II),',',int2str(round(v_('K1'))),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['KINF',int2str(II),',',int2str(round(v_('T')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                ename = ['KTP',int2str(I),',',int2str(S)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PROD';
                ielftype(ie) = iet_('en2PROD');
                vname = ['KEFF',int2str(S)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PHI',int2str(I),',',int2str(S)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                ename = ['P',int2str(I),',',int2str(S)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PROD';
                ielftype(ie) = iet_('en2PROD');
                vname = ['KINF',int2str(I),',',int2str(S)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PHI',int2str(I),',',int2str(S)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('M')
                for II=v_('1'):v_('N')
                    ig = ig_(['PLAC',int2str(I)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['Au',int2str(I),',',int2str(II),',',int2str(J)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['V',int2str(II)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['Bu',int2str(I),',',int2str(II),',',int2str(J)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['V',int2str(II)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['Cu',int2str(I),',',int2str(II),',',int2str(J)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['V',int2str(II)]);
                end
            end
        end
        for S=v_('1'):v_('T')
            for I=v_('1'):v_('N')
                ig = ig_(['KERN',int2str(I),',',int2str(S)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['KTP',int2str(I),',',int2str(S)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = -1.0;
                v_('TEMP') = v_(['NROW',int2str(I)]);
                v_('NuROW') = fix(v_('TEMP'));
                for II=v_('1'):v_('NuROW')
                    v_('TEMP') = v_(['ROWu',int2str(I),',',int2str(II)]);
                    v_('III') = fix(v_('TEMP'));
                    ig = ig_(['KERN',int2str(I),',',int2str(S)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('III'))),',',int2str(S)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['G',int2str(I),',',int2str(round(v_('III')))]);
                end
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T-1')
                ig = ig_(['KINFF',int2str(I),',',int2str(S)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(S)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-ACC');
            end
        end
        for S=v_('1'):v_('T')
            for I=v_('1'):v_('N')
                ig = ig_(['CPOW',int2str(S)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(S)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['V',int2str(I)]);
            end
        end
        for I=v_('1'):v_('N')
            for S=v_('1'):v_('T')
                ig = ig_(['PEAK',int2str(I),',',int2str(S)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(S)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MN-342-284';
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

    case 'en2PROD'

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

    case 'en3PROD'

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

