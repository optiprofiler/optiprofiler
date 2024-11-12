function varargout = LINVERSE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    Problem : LINVERSE
%    *********
% 
%    The problem is to find the positive definite lower bidiagonal
%    matrix L such that the matrix L(inv)L(inv-transp) best approximates,
%    in the Frobenius norm, a given symmetric target matrix T.
%    More precisely, one is  interested in the positive definite lower
%    bidiagonal L such that
% 
%         || L T L(transp) - I ||     is minimum.
%                                F
% 
%    The positive definite character of L is imposed by requiring
%    that all its diagonal entries to be at least equal to EPSILON,
%    a strictly positive real number.
% 
%    Many variants of the problem can be obtained by varying the target
%    matrix T and the scalar EPSILON.  In the present problem,
%    a) T is chosen to be pentadiagonal with T(i,j) = sin(i)cos(j) (j .leq. i)
%    b) EPSILON = 1.D-8
% 
%    Source:
%    Ph. Toint, private communication, 1991.
% 
%    SIF input: Ph. Toint, March 1991.
% 
%    classification = 'C-CSBR2-AN-V-0'
% 
%    Dimension of the matrix
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER  n = 19    original value
% IE N                   100            $-PARAMETER  n = 199
% IE N                   500            $-PARAMETER  n = 999
% IE N                   1000           $-PARAMETER  n = 1999
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LINVERSE';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('EPSILON') = 1.0e-8;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        for J=v_('1'):v_('N-2')
            v_('J+2') = 2+J;
            v_('RJ') = J;
            v_('COSJ') = cos(v_('RJ'));
            for I=J:v_('J+2')
                v_('RI') = I;
                v_('SINI') = sin(v_('RI'));
                v_(['T',int2str(I),',',int2str(J)]) = v_('SINI')*v_('COSJ');
            end
        end
        v_('RN-1') = v_('N-1');
        v_('SINI') = sin(v_('RN-1'));
        v_('COSJ') = cos(v_('RN-1'));
        v_(['T',int2str(round(v_('N-1'))),',',int2str(round(v_('N-1')))]) = v_('SINI')*...
             v_('COSJ');
        v_('RN') = v_('N');
        v_('SINI') = sin(v_('RN'));
        v_(['T',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))]) = v_('SINI')*...
             v_('COSJ');
        v_('COSJ') = cos(v_('RN'));
        v_(['T',int2str(round(v_('N'))),',',int2str(round(v_('N')))]) = v_('SINI')*...
             v_('COSJ');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['A',int2str(I)],ix_);
            pb.xnames{iv} = ['A',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii',['A',int2str(round(v_('N')))],ix_);
        pb.xnames{iv} = ['A',int2str(round(v_('N')))];
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('N-2')
            v_('J+1') = 1+J;
            v_('J+2') = 2+J;
            [ig,ig_] = s2mpjlib('ii',['O',int2str(J),',',int2str(J)],ig_);
            gtype{ig} = '<>';
            [ig,ig_] =...
                  s2mpjlib('ii',['O',int2str(round(v_('J+1'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.5;
            [ig,ig_] =...
                  s2mpjlib('ii',['O',int2str(round(v_('J+2'))),',',int2str(J)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.5;
        end
        [ig,ig_] =...
              s2mpjlib('ii',['O',int2str(round(v_('N-1'))),',',int2str(round(v_('N-1')))],ig_);
        gtype{ig} = '<>';
        [ig,ig_] =...
              s2mpjlib('ii',['O',int2str(round(v_('N'))),',',int2str(round(v_('N-1')))],ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 0.5;
        [ig,ig_] =...
              s2mpjlib('ii',['O',int2str(round(v_('N'))),',',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['O',int2str(I),',',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['A',int2str(I)]),1) = v_('EPSILON');
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = -1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = ['S',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('2'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('1')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('2'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('3'))),',',int2str(round(v_('2')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['S',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['U',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['A',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['V',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W',int2str(round(v_('3'))),',',int2str(round(v_('3')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['B',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('4'):v_('N')
            v_('I-1') = -1+I;
            v_('I-2') = -2+I;
            ename = ['S',int2str(I),',',int2str(round(v_('I-2')))];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
            end
            vname = ['A',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['S',int2str(I),',',int2str(round(v_('I-2')))];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
            end
            vname = ['A',int2str(round(v_('I-2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['V',int2str(I),',',int2str(round(v_('I-2')))];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
            end
            vname = ['B',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['V',int2str(I),',',int2str(round(v_('I-2')))];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
            end
            vname = ['A',int2str(round(v_('I-2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            for J=v_('I-1'):I
                v_('J-1') = -1+J;
                ename = ['S',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                end
                vname = ['A',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['A',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['U',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                end
                vname = ['A',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['B',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['V',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                end
                vname = ['B',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['A',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['W',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                end
                vname = ['B',int2str(round(v_('I-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['B',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_(['O',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        ig = ig_(['O',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['V',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        ig = ig_(['O',int2str(round(v_('3'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('3'))),',',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['V',int2str(round(v_('3'))),',',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        ig = ig_(['O',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['U',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['V',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['W',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        ig = ig_(['O',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['U',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['V',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['W',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
        ig = ig_(['O',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['S',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['U',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['V',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('3'))),',',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['W',int2str(round(v_('3'))),',',int2str(round(v_('3')))]);
        pbm.grelw{ig}(posel) =...
              v_(['T',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        for I=v_('4'):v_('N')
            v_('I-1') = -1+I;
            v_('I-2') = -2+I;
            ig = ig_(['O',int2str(I),',',int2str(round(v_('I-2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(round(v_('I-2')))]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(round(v_('I-2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['V',int2str(I),',',int2str(round(v_('I-2')))]);
            pbm.grelw{ig}(posel) =...
                  v_(['T',int2str(round(v_('I-1'))),',',int2str(round(v_('I-2')))]);
            ig = ig_(['O',int2str(I),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(round(v_('I-2')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['V',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) =...
                  v_(['T',int2str(round(v_('I-1'))),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['W',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) =...
                  v_(['T',int2str(round(v_('I-1'))),',',int2str(round(v_('I-2')))]);
            ig = ig_(['O',int2str(I),',',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(I)]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U',int2str(I),',',int2str(I)]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['V',int2str(I),',',int2str(I)]);
            pbm.grelw{ig}(posel) = v_(['T',int2str(I),',',int2str(round(v_('I-1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['W',int2str(I),',',int2str(I)]);
            pbm.grelw{ig}(posel) =...
                  v_(['T',int2str(round(v_('I-1'))),',',int2str(round(v_('I-1')))]);
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           6.00000000
% LO SOLTN(100)          68.0000000
% LO SOLTN(500)          340.000000
% LO SOLTN(1000)         ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSBR2-AN-V-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

