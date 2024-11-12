function varargout = LAKES(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A problem of water resource management in Canada, which may be 
%    formulated as
%    Min  SUM   SUM  (T(i,j)- R(i,j))^2 + (O(i,j)-R(N+i,j))^2)
%        i=1,N j=1,5 
%    subject to
%    T(i+1,1)-T(i,1)+O(i,1)        =  G(i,1)
%    T(i+1,2)-T(i,2)-O(i,1)+O(i,2) =  G(i,2)
%    T(i+1,3)-T(i,3)-O(i,2)+O(i,3) =  G(i,3)
%    T(i+1,4)-T(i,4)-O(i,3)+O(i,4) =  G(i,4)
%    T(i+1,5)-T(i,5)-O(i,4)+O(i,5) =  G(i,5) 
%    i=1,N and T(N+1,j) = T(1,j)  for j=1,5
%    O(i,2)-a*((T(i,2)/480.8+T(i,3)/4.6)/2-543.4)^2 * 
%    (T(i,2)/480.8-T(i,3)/4.6)^.5=0
%    O(i,3)-b*((T(i,3)/4.6-543.4)^2*(T(i,3)/4.6-T(i,4)/105.15)^0.5) = 0
%    O(i,4)-c*(T(i,4)/105.15-550.11)^2.2 = 0
%    where T(i,j) and O(i,j) are variables, R(i,j) are given and
%    a=.0841168  b=.1280849 and c=0.2605.
%    Extra variables 
%    
%    v(i,2) = T(i,2) / 961.6 + T(i,3) / 9.2 - 543.4
%    w(i,2) = T(i,2) / 480.8 - T(i,3) / 4.6
%    v(i,3) = T(i,3) / 4.6 - 543.4
%    w(i,3) = T(i,3) / 4.6 - T(i,4) / 105.15
%    v(i,4) = T(i,4) / 105.15 - 550.11
%    are introduced so that the nonlinear constraints may be rewritten as
%    O(i,2)-a*v(i,2)^2 * w(i,2)^0.5 = 0 ; w(i,2) > 0
%    O(i,3)-b*v(i,3)^2 * w(i,3)^0.5 = 0 ; w(i,3) > 0
%    O(i,4)-c*v(i,4)^2.2 = 0 ; v(i,4) > 0
%    Source:
%    S Jafar Sadjadi
%    Dept. of Systems Design Engineering
%    University of Waterloo
%    Ontario, N2L 3G1 Canada
% 
%    SIF input: Nick Gould and Jafar Sadjadi, November 1995
% 
%    classification = 'C-CQOR2-RN-90-78'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LAKES';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 6;
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        v_('NN') = 2*v_('N');
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
        v_('S1s1') = 202761.072;
        v_('S1s2') = 277791.816;
        v_('S1s3') = 2636.996;
        v_('S1s4') = 59987.0235;
        v_('S1s5') = 19490.4;
        v_('S2s1') = 202703.646;
        v_('S2s2') = 277849.512;
        v_('S2s3') = 2638.1;
        v_('S2s4') = 59998.59;
        v_('S2s5') = 19555.2;
        v_('S3s1') = 202720.536;
        v_('S3s2') = 277955.288;
        v_('S3s3') = 2639.894;
        v_('S3s4') = 60046.959;
        v_('S3s5') = 19597.6;
        v_('S4s1') = 202808.364;
        v_('S4s2') = 278104.336;
        v_('S4s3') = 2640.906;
        v_('S4s4') = 60074.298;
        v_('S4s5') = 19652.8;
        v_('S5s1') = 202916.46;
        v_('S5s2') = 278224.536;
        v_('S5s3') = 2641.458;
        v_('S5s4') = 60091.122;
        v_('S5s5') = 19708.8;
        v_('S6s1') = 202953.618;
        v_('S6s2') = 278277.424;
        v_('S6s3') = 2641.458;
        v_('S6s4') = 60082.71;
        v_('S6s5') = 19706.4;
        v_('O1o1') = 83.728;
        v_('O1o2') = 174.665;
        v_('O1o3') = 180.539;
        v_('O1o4') = 211.558;
        v_('O1o5') = 232.252;
        v_('O2o1') = 83.789;
        v_('O2o2') = 173.255;
        v_('O2o3') = 179.917;
        v_('O2o4') = 210.585;
        v_('O2o5') = 215.254;
        v_('O3o1') = 82.9160;
        v_('O3o2') = 173.721;
        v_('O3o3') = 182.676;
        v_('O3o4') = 207.838;
        v_('O3o5') = 203.855;
        v_('O4o1') = 80.134;
        v_('O4o2') = 178.654;
        v_('O4o3') = 185.917;
        v_('O4o4') = 206.416;
        v_('O4o5') = 186.308;
        v_('O5o1') = 65.345;
        v_('O5o2') = 188.01;
        v_('O5o3') = 192.568;
        v_('O5o4') = 204.3;
        v_('O5o5') = 201.1;
        v_('O6o1') = 72.005;
        v_('O6o2') = 193.833;
        v_('O6o3') = 196.651;
        v_('O6o4') = 204.25;
        v_('O6o5') = 241.079;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('5')
                [iv,ix_] = s2mpjlib('ii',['T',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['T',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['O',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['O',int2str(I),',',int2str(J)];
            end
            [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(round(v_('2')))],ix_);
            pb.xnames{iv} = ['V',int2str(I),',',int2str(round(v_('2')))];
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I),',',int2str(round(v_('2')))],ix_);
            pb.xnames{iv} = ['W',int2str(I),',',int2str(round(v_('2')))];
            [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(round(v_('3')))],ix_);
            pb.xnames{iv} = ['V',int2str(I),',',int2str(round(v_('3')))];
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I),',',int2str(round(v_('3')))],ix_);
            pb.xnames{iv} = ['W',int2str(I),',',int2str(round(v_('3')))];
            [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(round(v_('4')))],ix_);
            pb.xnames{iv} = ['V',int2str(I),',',int2str(round(v_('4')))];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for i=v_('1'):v_('N')
            v_('n+i') = v_('N')+i;
            for j=v_('1'):v_('5')
                [ig,ig_] = s2mpjlib('ii',['R',int2str(i),',',int2str(j)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['T',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                [ig,ig_] =...
                      s2mpjlib('ii',['R',int2str(round(v_('n+i'))),',',int2str(j)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['O',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for i=v_('1'):v_('N-1')
            v_('k+1') = 1+i;
            for j=v_('1'):v_('5')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(i),',',int2str(j)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(i),',',int2str(j)];
                iv = ix_(['T',int2str(round(v_('k+1'))),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['T',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['O',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
            for j=v_('2'):v_('5')
                v_('j-1') = -1+j;
                [ig,ig_] = s2mpjlib('ii',['G',int2str(i),',',int2str(j)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(i),',',int2str(j)];
                iv = ix_(['O',int2str(i),',',int2str(round(v_('j-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for j=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N'))),',',int2str(j)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(round(v_('N'))),',',int2str(j)];
            iv = ix_(['T',int2str(round(v_('1'))),',',int2str(j)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['T',int2str(round(v_('N'))),',',int2str(j)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N'))),',',int2str(j)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(round(v_('N'))),',',int2str(j)];
            iv = ix_(['O',int2str(round(v_('N'))),',',int2str(j)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for j=v_('2'):v_('5')
            v_('j-1') = -1+j;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N'))),',',int2str(j)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(round(v_('N'))),',',int2str(j)];
            iv = ix_(['O',int2str(round(v_('N'))),',',int2str(round(v_('j-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for i=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['A',int2str(i),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(i),',',int2str(round(v_('1')))];
            iv = ix_(['O',int2str(i),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['A',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['O',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['A',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['O',int2str(i),',',int2str(round(v_('4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['V',int2str(i),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('C') = 961.6;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            v_('C') = 9.2;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['W',int2str(i),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('C') = 480.8;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            v_('C') = -4.6;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('2')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['V',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('C') = 4.6;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['W',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('C') = 4.6;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            v_('C') = -105.15;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['W',int2str(i),',',int2str(round(v_('3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['W',int2str(i),',',int2str(round(v_('3')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('4')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('4')))];
            iv = ix_(['V',int2str(i),',',int2str(round(v_('4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('C') = 105.15;
            v_('C') = 1.0/v_('C');
            [ig,ig_] = s2mpjlib('ii',['V',int2str(i),',',int2str(round(v_('4')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(i),',',int2str(round(v_('4')))];
            iv = ix_(['T',int2str(i),',',int2str(round(v_('4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('C')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('C');
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
        pbm.gconst(ig_(['R',int2str(round(v_('1'))),',',int2str(round(v_('1')))])) =...
              v_('S1s1');
        pbm.gconst(ig_(['R',int2str(round(v_('1'))),',',int2str(round(v_('2')))])) =...
              v_('S1s2');
        pbm.gconst(ig_(['R',int2str(round(v_('1'))),',',int2str(round(v_('3')))])) =...
              v_('S1s3');
        pbm.gconst(ig_(['R',int2str(round(v_('1'))),',',int2str(round(v_('4')))])) =...
              v_('S1s4');
        pbm.gconst(ig_(['R',int2str(round(v_('1'))),',',int2str(round(v_('5')))])) =...
              v_('S1s5');
        pbm.gconst(ig_(['R',int2str(round(v_('2'))),',',int2str(round(v_('1')))])) =...
              v_('S2s1');
        pbm.gconst(ig_(['R',int2str(round(v_('2'))),',',int2str(round(v_('2')))])) =...
              v_('S2s2');
        pbm.gconst(ig_(['R',int2str(round(v_('2'))),',',int2str(round(v_('3')))])) =...
              v_('S2s3');
        pbm.gconst(ig_(['R',int2str(round(v_('2'))),',',int2str(round(v_('4')))])) =...
              v_('S2s4');
        pbm.gconst(ig_(['R',int2str(round(v_('2'))),',',int2str(round(v_('5')))])) =...
              v_('S2s5');
        pbm.gconst(ig_(['R',int2str(round(v_('3'))),',',int2str(round(v_('1')))])) =...
              v_('S3s1');
        pbm.gconst(ig_(['R',int2str(round(v_('3'))),',',int2str(round(v_('2')))])) =...
              v_('S3s2');
        pbm.gconst(ig_(['R',int2str(round(v_('3'))),',',int2str(round(v_('3')))])) =...
              v_('S3s3');
        pbm.gconst(ig_(['R',int2str(round(v_('3'))),',',int2str(round(v_('4')))])) =...
              v_('S3s4');
        pbm.gconst(ig_(['R',int2str(round(v_('3'))),',',int2str(round(v_('5')))])) =...
              v_('S3s5');
        pbm.gconst(ig_(['R',int2str(round(v_('4'))),',',int2str(round(v_('1')))])) =...
              v_('S4s1');
        pbm.gconst(ig_(['R',int2str(round(v_('4'))),',',int2str(round(v_('2')))])) =...
              v_('S4s2');
        pbm.gconst(ig_(['R',int2str(round(v_('4'))),',',int2str(round(v_('3')))])) =...
              v_('S4s3');
        pbm.gconst(ig_(['R',int2str(round(v_('4'))),',',int2str(round(v_('4')))])) =...
              v_('S4s4');
        pbm.gconst(ig_(['R',int2str(round(v_('4'))),',',int2str(round(v_('5')))])) =...
              v_('S4s5');
        pbm.gconst(ig_(['R',int2str(round(v_('5'))),',',int2str(round(v_('1')))])) =...
              v_('S5s1');
        pbm.gconst(ig_(['R',int2str(round(v_('5'))),',',int2str(round(v_('2')))])) =...
              v_('S5s2');
        pbm.gconst(ig_(['R',int2str(round(v_('5'))),',',int2str(round(v_('3')))])) =...
              v_('S5s3');
        pbm.gconst(ig_(['R',int2str(round(v_('5'))),',',int2str(round(v_('4')))])) =...
              v_('S5s4');
        pbm.gconst(ig_(['R',int2str(round(v_('5'))),',',int2str(round(v_('5')))])) =...
              v_('S5s5');
        pbm.gconst(ig_(['R',int2str(round(v_('6'))),',',int2str(round(v_('1')))])) =...
              v_('S6s1');
        pbm.gconst(ig_(['R',int2str(round(v_('6'))),',',int2str(round(v_('2')))])) =...
              v_('S6s2');
        pbm.gconst(ig_(['R',int2str(round(v_('6'))),',',int2str(round(v_('3')))])) =...
              v_('S6s3');
        pbm.gconst(ig_(['R',int2str(round(v_('6'))),',',int2str(round(v_('4')))])) =...
              v_('S6s4');
        pbm.gconst(ig_(['R',int2str(round(v_('6'))),',',int2str(round(v_('5')))])) =...
              v_('S6s5');
        pbm.gconst(ig_(['R',int2str(round(v_('7'))),',',int2str(round(v_('1')))])) =...
              v_('O1o5');
        pbm.gconst(ig_(['R',int2str(round(v_('7'))),',',int2str(round(v_('2')))])) =...
              v_('O1o2');
        pbm.gconst(ig_(['R',int2str(round(v_('7'))),',',int2str(round(v_('3')))])) =...
              v_('O1o3');
        pbm.gconst(ig_(['R',int2str(round(v_('7'))),',',int2str(round(v_('4')))])) =...
              v_('O1o4');
        pbm.gconst(ig_(['R',int2str(round(v_('7'))),',',int2str(round(v_('5')))])) =...
              v_('O1o5');
        pbm.gconst(ig_(['R',int2str(round(v_('8'))),',',int2str(round(v_('1')))])) =...
              v_('O2o1');
        pbm.gconst(ig_(['R',int2str(round(v_('8'))),',',int2str(round(v_('2')))])) =...
              v_('O2o2');
        pbm.gconst(ig_(['R',int2str(round(v_('8'))),',',int2str(round(v_('3')))])) =...
              v_('O2o3');
        pbm.gconst(ig_(['R',int2str(round(v_('8'))),',',int2str(round(v_('4')))])) =...
              v_('O2o4');
        pbm.gconst(ig_(['R',int2str(round(v_('8'))),',',int2str(round(v_('5')))])) =...
              v_('O2o5');
        pbm.gconst(ig_(['R',int2str(round(v_('9'))),',',int2str(round(v_('1')))])) =...
              v_('O3o1');
        pbm.gconst(ig_(['R',int2str(round(v_('9'))),',',int2str(round(v_('2')))])) =...
              v_('O3o2');
        pbm.gconst(ig_(['R',int2str(round(v_('9'))),',',int2str(round(v_('3')))])) =...
              v_('O3o3');
        pbm.gconst(ig_(['R',int2str(round(v_('9'))),',',int2str(round(v_('4')))])) =...
              v_('O3o4');
        pbm.gconst(ig_(['R',int2str(round(v_('9'))),',',int2str(round(v_('5')))])) =...
              v_('O3o5');
        pbm.gconst(ig_(['R',int2str(round(v_('10'))),',',int2str(round(v_('1')))])) = v_('O4o1');
        pbm.gconst(ig_(['R',int2str(round(v_('10'))),',',int2str(round(v_('2')))])) = v_('O4o2');
        pbm.gconst(ig_(['R',int2str(round(v_('10'))),',',int2str(round(v_('3')))])) = v_('O4o3');
        pbm.gconst(ig_(['R',int2str(round(v_('10'))),',',int2str(round(v_('4')))])) = v_('O4o4');
        pbm.gconst(ig_(['R',int2str(round(v_('10'))),',',int2str(round(v_('5')))])) = v_('O4o5');
        pbm.gconst(ig_(['R',int2str(round(v_('11'))),',',int2str(round(v_('1')))])) = v_('O5o1');
        pbm.gconst(ig_(['R',int2str(round(v_('11'))),',',int2str(round(v_('2')))])) = v_('O5o2');
        pbm.gconst(ig_(['R',int2str(round(v_('11'))),',',int2str(round(v_('3')))])) = v_('O5o3');
        pbm.gconst(ig_(['R',int2str(round(v_('11'))),',',int2str(round(v_('4')))])) = v_('O5o4');
        pbm.gconst(ig_(['R',int2str(round(v_('11'))),',',int2str(round(v_('5')))])) = v_('O5o5');
        pbm.gconst(ig_(['R',int2str(round(v_('12'))),',',int2str(round(v_('1')))])) = v_('O6o1');
        pbm.gconst(ig_(['R',int2str(round(v_('12'))),',',int2str(round(v_('2')))])) = v_('O6o2');
        pbm.gconst(ig_(['R',int2str(round(v_('12'))),',',int2str(round(v_('3')))])) = v_('O6o3');
        pbm.gconst(ig_(['R',int2str(round(v_('12'))),',',int2str(round(v_('4')))])) = v_('O6o4');
        pbm.gconst(ig_(['R',int2str(round(v_('12'))),',',int2str(round(v_('5')))])) = v_('O6o5');
        pbm.gconst(ig_(['G',int2str(round(v_('1'))),',',int2str(round(v_('1')))])) = -...
             22.0;
        pbm.gconst(ig_(['G',int2str(round(v_('1'))),',',int2str(round(v_('2')))])) = -...
             1.0;
        pbm.gconst(ig_(['G',int2str(round(v_('1'))),',',int2str(round(v_('3')))])) =...
              3.0;
        pbm.gconst(ig_(['G',int2str(round(v_('1'))),',',int2str(round(v_('4')))])) = -...
             27.2;
        pbm.gconst(ig_(['G',int2str(round(v_('1'))),',',int2str(round(v_('5')))])) =...
              51.5;
        pbm.gconst(ig_(['G',int2str(round(v_('2'))),',',int2str(round(v_('1')))])) =...
              44.0;
        pbm.gconst(ig_(['G',int2str(round(v_('2'))),',',int2str(round(v_('2')))])) =...
              162.0;
        pbm.gconst(ig_(['G',int2str(round(v_('2'))),',',int2str(round(v_('3')))])) =...
              8.0;
        pbm.gconst(ig_(['G',int2str(round(v_('2'))),',',int2str(round(v_('4')))])) =...
              12.5;
        pbm.gconst(ig_(['G',int2str(round(v_('2'))),',',int2str(round(v_('5')))])) =...
              53.5;
        pbm.gconst(ig_(['G',int2str(round(v_('3'))),',',int2str(round(v_('1')))])) = -...
             11.0;
        pbm.gconst(ig_(['G',int2str(round(v_('3'))),',',int2str(round(v_('2')))])) =...
              60.0;
        pbm.gconst(ig_(['G',int2str(round(v_('3'))),',',int2str(round(v_('3')))])) =...
              10.0;
        pbm.gconst(ig_(['G',int2str(round(v_('3'))),',',int2str(round(v_('4')))])) =...
              18.0;
        pbm.gconst(ig_(['G',int2str(round(v_('3'))),',',int2str(round(v_('5')))])) =...
              39.0;
        pbm.gconst(ig_(['G',int2str(round(v_('4'))),',',int2str(round(v_('1')))])) =...
              124.0;
        pbm.gconst(ig_(['G',int2str(round(v_('4'))),',',int2str(round(v_('2')))])) =...
              246.0;
        pbm.gconst(ig_(['G',int2str(round(v_('4'))),',',int2str(round(v_('3')))])) =...
              6.0;
        pbm.gconst(ig_(['G',int2str(round(v_('4'))),',',int2str(round(v_('4')))])) =...
              9.7;
        pbm.gconst(ig_(['G',int2str(round(v_('4'))),',',int2str(round(v_('5')))])) =...
              17.2;
        pbm.gconst(ig_(['G',int2str(round(v_('5'))),',',int2str(round(v_('1')))])) =...
              127.0;
        pbm.gconst(ig_(['G',int2str(round(v_('5'))),',',int2str(round(v_('2')))])) =...
              175.0;
        pbm.gconst(ig_(['G',int2str(round(v_('5'))),',',int2str(round(v_('3')))])) =...
              3.0;
        pbm.gconst(ig_(['G',int2str(round(v_('5'))),',',int2str(round(v_('4')))])) =...
              10.0;
        pbm.gconst(ig_(['G',int2str(round(v_('5'))),',',int2str(round(v_('5')))])) =...
              30.2;
        pbm.gconst(ig_(['G',int2str(round(v_('6'))),',',int2str(round(v_('1')))])) =...
              78.0;
        pbm.gconst(ig_(['G',int2str(round(v_('6'))),',',int2str(round(v_('2')))])) =...
              156.0;
        pbm.gconst(ig_(['G',int2str(round(v_('6'))),',',int2str(round(v_('3')))])) =...
              3.0;
        pbm.gconst(ig_(['G',int2str(round(v_('6'))),',',int2str(round(v_('4')))])) =...
              14.0;
        pbm.gconst(ig_(['G',int2str(round(v_('6'))),',',int2str(round(v_('5')))])) =...
              23.2;
        for i=v_('1'):v_('N')
            pbm.gconst(ig_(['V',int2str(i),',',int2str(round(v_('2')))])) = 543.4;
            pbm.gconst(ig_(['W',int2str(i),',',int2str(round(v_('2')))])) = 0.0;
            pbm.gconst(ig_(['V',int2str(i),',',int2str(round(v_('3')))])) = 543.4;
            pbm.gconst(ig_(['W',int2str(i),',',int2str(round(v_('3')))])) = 0.0;
            pbm.gconst(ig_(['V',int2str(i),',',int2str(round(v_('4')))])) = 550.11;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for i=v_('1'):v_('N')
            pb.xlower(ix_(['W',int2str(i),',',int2str(round(v_('2')))]),1) = 0.0001;
            pb.xlower(ix_(['W',int2str(i),',',int2str(round(v_('3')))]),1) = 0.0001;
            pb.xlower(ix_(['V',int2str(i),',',int2str(round(v_('4')))]),1) = 0.0001;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.y0 = 1.0*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en1VAR',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'en2VAR',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for i=v_('1'):v_('N')
            ename = ['B',int2str(i),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2VAR';
            ielftype(ie) = iet_('en2VAR');
            ename = ['B',int2str(i),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(i),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(i),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['W',int2str(i),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(i),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 0.0841168;
            ename = ['B',int2str(i),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2VAR';
            ielftype(ie) = iet_('en2VAR');
            ename = ['B',int2str(i),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(i),',',int2str(round(v_('3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(i),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['W',int2str(i),',',int2str(round(v_('3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(i),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 0.1280849;
            ename = ['B',int2str(i),',',int2str(round(v_('3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en1VAR';
            ielftype(ie) = iet_('en1VAR');
            ename = ['B',int2str(i),',',int2str(round(v_('3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(i),',',int2str(round(v_('4')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['B',int2str(i),',',int2str(round(v_('3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 0.2605;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for i=v_('1'):v_('N')
            ig = ig_(['A',int2str(i),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(i),',',int2str(round(v_('1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['A',int2str(i),',',int2str(round(v_('2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(i),',',int2str(round(v_('2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['A',int2str(i),',',int2str(round(v_('3')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(i),',',int2str(round(v_('3')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for i=v_('1'):v_('N')
            v_('n+i') = v_('N')+i;
            for j=v_('1'):v_('5')
                ig = ig_(['R',int2str(i),',',int2str(j)]);
                pbm.grftype{ig} = 'gL2';
                ig = ig_(['R',int2str(round(v_('n+i'))),',',int2str(j)]);
                pbm.grftype{ig} = 'gL2';
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-RN-90-78';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en1VAR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)^2.2;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*2.2*EV_(1)^1.2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = pbm.elpar{iel_}(1)*2.64*EV_(1)^0.2;
                varargout{3} = H_;
            end
        end

    case 'en2VAR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)^2*EV_(2)^0.5;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*2.0*EV_(1)*EV_(2)^0.5;
            g_(2,1) = pbm.elpar{iel_}(1)*0.5*EV_(1)^2/EV_(2)^0.5;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.elpar{iel_}(1)*2.0*EV_(2)^0.5;
                H_(1,2) = pbm.elpar{iel_}(1)*EV_(1)/EV_(2)^0.5;
                H_(2,1) = H_(1,2);
                H_(2,2) = -pbm.elpar{iel_}(1)*0.25*EV_(1)^2/EV_(2)^1.5;
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

