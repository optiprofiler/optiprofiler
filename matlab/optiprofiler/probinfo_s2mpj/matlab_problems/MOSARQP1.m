function varargout = MOSARQP1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MOSARQP1
%    *********
% 
%    A convex quadratic problem, with variable dimensions.
%    In this problem, half the linear constraints are active at the solution.
% 
%    Source:
%    J.L. Morales-Perez and R.W.H. Sargent,
%    "On the implementation and performance of an interior point method for
%    large sparse convex quadratic programming",
%    Centre for Process Systems Engineering, Imperial College, London,
%    November 1991.
% 
%    SIF input: Ph. Toint, August 1993.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-CQLR2-AN-V-V'
% 
%    Problem variants: these are distinguished by the triplet ( N, M, COND ),
%    where: - N (nb of variables) must be even and have an integer square root
%           - M (nb of constraints) must be at least sqrt(N) 
%             and at most N - sqrt(N)
%           - COND (problem conditioning) is a positive integer
%    Except for the first, the instances suggested are those used by Morales
%    and Sargent.
% 
%       Alternative values for the SIF file parameters:
% IE N                   36             $-PARAMETER
% IE M                   10             $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   100            $-PARAMETER     original value
% IE M                   10             $-PARAMETER     original value
% RE COND                3.0            $-PARAMETER     original value
% 
% IE N                   900            $-PARAMETER
% IE M                   30             $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   30             $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   30             $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   60             $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   60             $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   60             $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   90             $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   90             $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   90             $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   120            $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   120            $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   120            $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   300            $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   300            $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   300            $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   600            $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   600            $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   900            $-PARAMETER
% IE M                   600            $-PARAMETER
% RE COND                3.0            $-PARAMETER
% 
% IE N                   2500           $-PARAMETER
% IE M                   700            $-PARAMETER
% RE COND                1.0            $-PARAMETER
% 
% IE N                   2500           $-PARAMETER
% IE M                   700            $-PARAMETER
% RE COND                2.0            $-PARAMETER
% 
% IE N                   2500           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MOSARQP1';

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
        if(nargs<1)
            v_('N') = 36;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE M                   700            $-PARAMETER
        if(nargs<2)
            v_('M') = 10;  %  SIF file default value
        else
            v_('M') = varargin{2};
        end
% RE COND                3.0            $-PARAMETER
        if(nargs<3)
            v_('COND') = 2.0;  %  SIF file default value
        else
            v_('COND') = varargin{3};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('N-1') = -1+v_('N');
        v_('RN-1') = v_('N-1');
        v_('RN') = v_('N');
        v_('M-1') = -1+v_('M');
        v_('RNP') = 0.1+v_('RN');
        v_('RRTN') = sqrt(v_('RNP'));
        v_('RTN') = fix(v_('RRTN'));
        v_('RTN-1') = -1+v_('RTN');
        v_('RTN-2') = -2+v_('RTN');
        v_('RTN+1') = 1+v_('RTN');
        v_('2RTN') = v_('RTN')+v_('RTN');
        v_('M-RTN+1') = v_('M')-v_('RTN-1');
        for I=v_('1'):v_('2'):v_('N-1')
            v_('I+1') = 1+I;
            v_(['XC',int2str(I)]) = -1.0;
            v_(['XC',int2str(round(v_('I+1')))]) = 1.0;
        end
        v_(['XC',int2str(round(v_('N')))]) = 1.0;
        v_('NNZ') = 10;
        v_('Y1') = -0.3569732;
        v_('Y2') = 0.9871576;
        v_('Y3') = 0.5619363;
        v_('Y4') = -0.1984624;
        v_('Y5') = 0.4653328;
        v_('Y6') = 0.7364367;
        v_('Y7') = -0.4560378;
        v_('Y8') = -0.6457813;
        v_('Y9') = -0.0601357;
        v_('Y10') = 0.1035624;
        v_('NZ1') = 0.68971452;
        v_('NZ2') = 0.13452678;
        v_('NZ3') = 0.51234678;
        v_('NZ4') = 0.76591423;
        v_('NZ5') = 0.20857854;
        v_('NZ6') = 0.85672348;
        v_('NZ7') = 0.04356789;
        v_('NZ8') = 0.44692743;
        v_('NZ9') = 0.30136413;
        v_('NZ10') = 0.91367489;
        v_('YN2') = 0.0;
        for I=v_('1'):v_('NNZ')
            v_('RKI') = v_(['NZ',int2str(I)])*v_('RN');
            v_(['K',int2str(I)]) = 1.1+v_('RKI');
            v_('TMP') = v_(['Y',int2str(I)])*v_(['Y',int2str(I)]);
            v_('YN2') = v_('YN2')+v_('TMP');
        end
        v_('-2/YN2') = -2.0/v_('YN2');
        v_('4/YN4') = v_('-2/YN2')*v_('-2/YN2');
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TMP') = v_('RI-1')/v_('RN-1');
            v_('TMP') = v_('TMP')*v_('COND');
            v_(['D',int2str(I)]) = exp(v_('TMP'));
        end
        v_('YDY') = 0.0;
        v_('YXC') = 0.0;
        v_('YDXC') = 0.0;
        for I=v_('1'):v_('NNZ')
            v_('RKI') = v_(['K',int2str(I)]);
            v_('KI') = fix(v_('RKI'));
            v_(['DY',int2str(I)]) =...
                  v_(['Y',int2str(I)])*v_(['D',int2str(round(v_('KI')))]);
            v_('TMP') = v_(['DY',int2str(I)])*v_(['Y',int2str(I)]);
            v_('YDY') = v_('YDY')+v_('TMP');
            v_('TMP') = v_(['Y',int2str(I)])*v_(['XC',int2str(round(v_('KI')))]);
            v_('YXC') = v_('YXC')+v_('TMP');
            v_('TMP') = v_(['DY',int2str(I)])*v_(['XC',int2str(round(v_('KI')))]);
            v_('YDXC') = v_('YDXC')+v_('TMP');
        end
        v_('AA') = v_('-2/YN2')*v_('YXC');
        v_('DD') = v_('4/YN4')*v_('YDY');
        v_('BB') = v_('DD')*v_('YXC');
        v_('CC') = v_('-2/YN2')*v_('YDXC');
        v_('BB+CC') = v_('BB')+v_('CC');
        v_('DD/2') = 0.5*v_('DD');
        for I=v_('1'):v_('N')
            v_(['C',int2str(I)]) = v_(['D',int2str(I)])*v_(['XC',int2str(I)]);
        end
        for I=v_('1'):v_('NNZ')
            v_('RKI') = v_(['K',int2str(I)]);
            v_('KI') = fix(v_('RKI'));
            v_('TMP') = v_(['DY',int2str(I)])*v_('AA');
            v_(['C',int2str(round(v_('KI')))]) = v_(['C',int2str(round(v_('KI')))])+...
                 v_('TMP');
            v_('TMP') = v_(['Y',int2str(I)])*v_('BB+CC');
            v_(['C',int2str(round(v_('KI')))]) = v_(['C',int2str(round(v_('KI')))])+...
                 v_('TMP');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['C',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['C',int2str(I)]);
            end
        end
        [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['CS',int2str(round(v_('1')))];
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['CS',int2str(round(v_('1')))];
        iv = ix_(['X',int2str(round(v_('RTN+1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_(['X',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('2'):v_('RTN-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('I+RTN') = I+v_('RTN');
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 4.0;
            end
            iv = ix_(['X',int2str(round(v_('I+RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('RTN')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['CS',int2str(round(v_('RTN')))];
        iv = ix_(['X',int2str(round(v_('RTN')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('RTN')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['CS',int2str(round(v_('RTN')))];
        iv = ix_(['X',int2str(round(v_('RTN-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_(['X',int2str(round(v_('2RTN')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        v_('JS') = v_('RTN');
        for J=v_('RTN+1'):v_('RTN'):v_('M-RTN+1')
            v_('J+1') = 1+J;
            v_('JS') = J+v_('RTN-1');
            v_('JS-1') = -1+v_('JS');
            v_('J-RTN') = J-v_('RTN');
            v_('J+RTN') = J+v_('RTN');
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(J)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(J)];
            iv = ix_(['X',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 4.0;
            end
            iv = ix_(['X',int2str(round(v_('J+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('J-RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('J+RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            for I=v_('J+1'):v_('JS-1')
                v_('I+1') = 1+I;
                v_('I-1') = -1+I;
                v_('I+RTN') = I+v_('RTN');
                v_('I-RTN') = I-v_('RTN');
                [ig,ig_] = s2mpjlib('ii',['CS',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['CS',int2str(I)];
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 4.0;
                end
                iv = ix_(['X',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('I+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('I-RTN')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('I+RTN')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
            v_('JS+RTN') = v_('JS')+v_('RTN');
            v_('JS-RTN') = v_('JS')-v_('RTN');
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('JS')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(round(v_('JS')))];
            iv = ix_(['X',int2str(round(v_('JS')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 4.0;
            end
            iv = ix_(['X',int2str(round(v_('JS-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('JS')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(round(v_('JS')))];
            iv = ix_(['X',int2str(round(v_('JS-RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('JS+RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        v_('K') = 1+v_('JS');
        for I=v_('K'):v_('M'):v_('M')
            v_('K+1') = 1+v_('K');
            v_('K+RTN') = v_('K')+v_('RTN');
            v_('K-RTN') = v_('K')-v_('RTN');
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('K')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(round(v_('K')))];
            iv = ix_(['X',int2str(round(v_('K')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 4.0;
            end
            iv = ix_(['X',int2str(round(v_('K+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(round(v_('K')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(round(v_('K')))];
            iv = ix_(['X',int2str(round(v_('K-RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('K+RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        v_('K') = 1+v_('K');
        for I=v_('K'):v_('M')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('I+RTN') = I+v_('RTN');
            v_('I-RTN') = I-v_('RTN');
            [ig,ig_] = s2mpjlib('ii',['CS',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CS',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 4.0;
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('I-RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+RTN')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
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
        pbm.gconst(ig_(['CS',int2str(round(v_('1')))])) = 0.5;
        pbm.gconst(ig_(['CS',int2str(round(v_('RTN')))])) = 0.5;
        v_('K') = v_('RTN+1');
        for J=v_('RTN+1'):v_('RTN'):v_('M-RTN+1')
            v_('K') = 1+v_('K');
            for I=v_('1'):v_('RTN-2')
                v_('K') = 1+v_('K');
                pbm.gconst(ig_(['CS',int2str(round(v_('K')))])) = -0.5;
            end
            v_('K') = 1+v_('K');
        end
        v_('K') = 1+v_('K');
        for J=v_('K'):v_('M')
            pbm.gconst(ig_(['CS',int2str(J)])) = -0.5;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['XSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('NNZ')
            v_('RKI') = v_(['K',int2str(I)]);
            v_('KI') = fix(v_('RKI'));
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('RKJ') = v_(['K',int2str(J)]);
                v_('KJ') = fix(v_('RKJ'));
                ename = ['P',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['X',int2str(round(v_('KI')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('KJ')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            v_('TMP') = 0.5*v_(['D',int2str(I)]);
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('TMP');
        end
        for I=v_('1'):v_('NNZ')
            v_('RKI') = v_(['K',int2str(I)]);
            v_('KI') = fix(v_('RKI'));
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                v_('TMP') = v_(['DY',int2str(I)])*v_(['Y',int2str(J)]);
                v_('WIJ') = v_('TMP')*v_('-2/YN2');
                v_('TMP') = v_(['DY',int2str(J)])*v_(['Y',int2str(I)]);
                v_('TMP') = v_('TMP')*v_('-2/YN2');
                v_('WIJ') = v_('WIJ')+v_('TMP');
                v_('TMP') = v_(['Y',int2str(I)])*v_(['Y',int2str(J)]);
                v_('TMP') = v_('TMP')*v_('DD');
                v_('WIJ') = v_('WIJ')+v_('TMP');
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('WIJ');
            end
            v_('TMP') = v_(['DY',int2str(I)])*v_(['Y',int2str(I)]);
            v_('WII') = v_('TMP')*v_('-2/YN2');
            v_('TMP') = v_(['Y',int2str(I)])*v_(['Y',int2str(I)]);
            v_('TMP') = v_('TMP')*v_('DD/2');
            v_('WII') = v_('WII')+v_('TMP');
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('KI')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('WII');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(  36, 10,2)   -24.13768932
% LO SOLTN( 100, 10,3)   -154.2001028
% LO SOLTN( 900, 30,1)   -380.0891288
% LO SOLTN( 900, 30,2)   -711.7109010
% LO SOLTN( 900, 30,3)   -1424.280649
% LO SOLTN( 900, 60,1)   -374.7138829
% LO SOLTN( 900, 60,2)   -706.1411506
% LO SOLTN( 900, 60,3)   -1418.592879
% LO SOLTN( 900, 90,1)   -369.8384609
% LO SOLTN( 900, 90,2)   -700.8243599
% LO SOLTN( 900, 90,3)   -1412.776689
% LO SOLTN( 900,120,1)   -364.8603691
% LO SOLTN( 900,120,2)   -695.2431416
% LO SOLTN( 900,120,3)   -1406.503648
% LO SOLTN( 900,300,1)   -331.0120280
% LO SOLTN( 900,300,2)   -652.2778434
% LO SOLTN( 900,300,3)   -1351.831332
% LO SOLTN( 900,600,1)   -257.4400842
% LO SOLTN( 900,600,2)   -529.6445809
% LO SOLTN( 900,600,3)   -1145.403000
% LO SOLTN(2500,700,1)   -952.8754378
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-AN-V-V';
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

    case 'eSQ'

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

