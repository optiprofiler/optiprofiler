function varargout = MESH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    The goodness of a finite element grid is characterized by the
%    smallest angle of all triangles. Given a triangulation of a domain 
%    in R**2. Find a topological equivalent triangulation so, that the 
%    smallest angel becomes as large as possible. Topological equivalent 
%    means shifting the edges of the grid only in such a way that 
%    neighbouring triangles remain neighbours. 
% 
%    Source: Prof. Dr. Michael Kraetzschmar, Institut fuer Angewandte 
%            Mathematik der Fachhochschule Flensburg, Kanzleistrasse 91-93, 
%            D-24943 FLENSBURG, GERMANY
% 
%    SIF input: Prof. Dr. Michael Kraetzschmar
% 
%    classification = 'C-COOR2-AY-41-48'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MESH';

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
        v_('omega1') = 1.0e3;
        v_('omega2') = -1.0e3;
        v_('omega3') = -1.0e5;
        v_('s') = 0.700000;
        v_('pi') = acos(-1.0);
        v_('sqrt3/2') = sqrt(0.75);
        v_('h') = v_('sqrt3/2')*v_('s');
        v_('drei') = 3.0;
        v_('pi/3') = v_('pi')/v_('drei');
        v_('1') = 1;
        v_('np') = 5;
        v_('nk') = 8;
        v_('nd') = 4;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for i=v_('1'):v_('np')
            [iv,ix_] = s2mpjlib('ii',['x',int2str(i)],ix_);
            pb.xnames{iv} = ['x',int2str(i)];
            [iv,ix_] = s2mpjlib('ii',['y',int2str(i)],ix_);
            pb.xnames{iv} = ['y',int2str(i)];
        end
        for i=v_('1'):v_('nk')
            [iv,ix_] = s2mpjlib('ii',['l',int2str(i)],ix_);
            pb.xnames{iv} = ['l',int2str(i)];
        end
        for i=v_('1'):v_('nd')
            [iv,ix_] = s2mpjlib('ii',['alpha',int2str(i)],ix_);
            pb.xnames{iv} = ['alpha',int2str(i)];
            [iv,ix_] = s2mpjlib('ii',['beta',int2str(i)],ix_);
            pb.xnames{iv} = ['beta',int2str(i)];
            [iv,ix_] = s2mpjlib('ii',['gamma',int2str(i)],ix_);
            pb.xnames{iv} = ['gamma',int2str(i)];
            [iv,ix_] = s2mpjlib('ii',['delta',int2str(i)],ix_);
            pb.xnames{iv} = ['delta',int2str(i)];
            [iv,ix_] = s2mpjlib('ii',['f',int2str(i)],ix_);
            pb.xnames{iv} = ['f',int2str(i)];
        end
        [iv,ix_] = s2mpjlib('ii','deltamin',ix_);
        pb.xnames{iv} = 'deltamin';
        [iv,ix_] = s2mpjlib('ii','fmin',ix_);
        pb.xnames{iv} = 'fmin';
        [iv,ix_] = s2mpjlib('ii','fmax',ix_);
        pb.xnames{iv} = 'fmax';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','obj1',ig_);
        gtype{ig} = '<>';
        iv = ix_('deltamin');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = v_('omega1');
        [ig,ig_] = s2mpjlib('ii','obj2',ig_);
        gtype{ig} = '<>';
        iv = ix_('fmax');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('fmin');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('omega2');
        [ig,ig_] = s2mpjlib('ii','obj3',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('omega3');
        for i=v_('1'):v_('nk')
            [ig,ig_] = s2mpjlib('ii',['seit',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['seit',int2str(i)];
        end
        for i=v_('1'):v_('nd')
            [ig,ig_] = s2mpjlib('ii',['skal',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['skal',int2str(i)];
            [ig,ig_] = s2mpjlib('ii',['skbe',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['skbe',int2str(i)];
        end
        for i=v_('1'):v_('nd')
            [ig,ig_] = s2mpjlib('ii',['doppf',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['doppf',int2str(i)];
            iv = ix_(['f',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for i=v_('1'):v_('nd')
            [ig,ig_] = s2mpjlib('ii',['wisum',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['wisum',int2str(i)];
            iv = ix_(['alpha',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['beta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['gamma',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for i=v_('1'):v_('nd')
            [ig,ig_] = s2mpjlib('ii',['alphd',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['alphd',int2str(i)];
            iv = ix_(['alpha',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['delta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['betad',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['betad',int2str(i)];
            iv = ix_(['beta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['delta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['gammd',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['gammd',int2str(i)];
            iv = ix_(['gamma',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['delta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['deltd',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['deltd',int2str(i)];
            iv = ix_(['delta',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('deltamin');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for i=v_('1'):v_('nd')
            [ig,ig_] = s2mpjlib('ii',['fmind',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['fmind',int2str(i)];
            iv = ix_(['f',int2str(i)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('fmin');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['fmaxd',int2str(i)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['fmaxd',int2str(i)];
            iv = ix_('fmax');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['f',int2str(i)]);
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
        for i=v_('1'):v_('nd')
            pbm.gconst(ig_(['wisum',int2str(i)])) = v_('pi');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('deltamin')) = v_('pi');
        pb.xlower(ix_('x1'),1) = 0.000000;
        pb.xupper(ix_('x1'),1) = 0.000000;
        pb.xlower(ix_('y1'),1) = 0.000000;
        pb.xupper(ix_('y1'),1) = 0.000000;
        pb.xlower(ix_('x2'),1) = 0.000000;
        pb.xupper(ix_('x2'),1) = 0.000000;
        pb.xlower(ix_('y2'),1) = 1.000000;
        pb.xupper(ix_('y2'),1) = 1.000000;
        pb.xlower(ix_('x3'),1) = 1.000000;
        pb.xupper(ix_('x3'),1) = 1.000000;
        pb.xlower(ix_('y3'),1) = 1.000000;
        pb.xupper(ix_('y3'),1) = 1.000000;
        pb.xlower(ix_('x4'),1) = 1.000000;
        pb.xupper(ix_('x4'),1) = 1.000000;
        pb.xlower(ix_('y4'),1) = 0.000000;
        pb.xupper(ix_('y4'),1) = 0.000000;
        pb.xlower(ix_('x5')) = -Inf;
        pb.xupper(ix_('x5'),1) = +Inf;
        pb.xlower(ix_('y5')) = -Inf;
        pb.xupper(ix_('y5'),1) = +Inf;
        for i=v_('1'):v_('nd')
            pb.xupper(ix_(['alpha',int2str(i)])) = v_('pi');
            pb.xupper(ix_(['beta',int2str(i)])) = v_('pi');
            pb.xupper(ix_(['gamma',int2str(i)])) = v_('pi');
            pb.xupper(ix_(['delta',int2str(i)])) = v_('pi');
        end
        pb.xupper(ix_('deltamin')) = v_('pi');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('x1'),1) = 0.000000;
        pb.x0(ix_('y1'),1) = 0.000000;
        pb.x0(ix_('x2'),1) = 0.000000;
        pb.x0(ix_('y2'),1) = 1.000000;
        pb.x0(ix_('x3'),1) = 1.000000;
        pb.x0(ix_('y3'),1) = 1.000000;
        pb.x0(ix_('x4'),1) = 1.000000;
        pb.x0(ix_('y4'),1) = 0.000000;
        pb.x0(ix_('x5'),1) = 0.350000;
        pb.x0(ix_('y5'),1) = 0.606218;
        pb.x0(ix_('l1'),1) = 0.700000;
        pb.x0(ix_('l2'),1) = 0.526844;
        pb.x0(ix_('l3'),1) = 1.000000;
        pb.x0(ix_('l4'),1) = 0.759977;
        pb.x0(ix_('l5'),1) = 1.000000;
        pb.x0(ix_('l6'),1) = 0.888819;
        pb.x0(ix_('l7'),1) = 1.000000;
        pb.x0(ix_('l8'),1) = 1.000000;
        pb.x0(ix_('alpha1'),1) = 0.523599;
        pb.x0(ix_('beta1'),1) = 1.891392;
        pb.x0(ix_('gamma1'),1) = 0.726602;
        pb.x0(ix_('delta1'),1) = 0.523599;
        pb.x0(ix_('f1'),1) = 0.350000;
        pb.x0(ix_('alpha2'),1) = 0.844195;
        pb.x0(ix_('beta2'),1) = 1.752711;
        pb.x0(ix_('gamma2'),1) = 0.544687;
        pb.x0(ix_('delta2'),1) = 0.544687;
        pb.x0(ix_('f2'),1) = 0.393782;
        pb.x0(ix_('alpha3'),1) = 1.026109;
        pb.x0(ix_('beta3'),1) = 1.295247;
        pb.x0(ix_('gamma3'),1) = 0.820236;
        pb.x0(ix_('delta3'),1) = 0.820236;
        pb.x0(ix_('f3'),1) = 0.650000;
        pb.x0(ix_('alpha4'),1) = 0.750560;
        pb.x0(ix_('beta4'),1) = 1.343835;
        pb.x0(ix_('gamma4'),1) = 1.047198;
        pb.x0(ix_('delta4'),1) = 0.750560;
        pb.x0(ix_('f4'),1) = 0.606218;
        pb.x0(ix_('deltamin'),1) = 0.523599;
        pb.x0(ix_('fmin'),1) = 0.350000;
        pb.x0(ix_('fmax'),1) = 0.650000;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ediffsq',iet_);
        elftv{it}{1} = 'winkel';
        elftv{it}{2} = 'minwi';
        [it,iet_] = s2mpjlib( 'ii', 'elaenge',iet_);
        elftv{it}{1} = 'lp1x';
        elftv{it}{2} = 'lp1y';
        elftv{it}{3} = 'lp2x';
        elftv{it}{4} = 'lp2y';
        elftv{it}{5} = 'llaenge';
        [it,iet_] = s2mpjlib( 'ii', 'evekprod',iet_);
        elftv{it}{1} = 'vp1x';
        elftv{it}{2} = 'vp1y';
        elftv{it}{3} = 'vp2x';
        elftv{it}{4} = 'vp2y';
        elftv{it}{5} = 'vp3x';
        elftv{it}{6} = 'vp3y';
        [it,iet_] = s2mpjlib( 'ii', 'esklprod',iet_);
        elftv{it}{1} = 'sp1x';
        elftv{it}{2} = 'sp1y';
        elftv{it}{3} = 'sp2x';
        elftv{it}{4} = 'sp2y';
        elftv{it}{5} = 'sp3x';
        elftv{it}{6} = 'sp3y';
        [it,iet_] = s2mpjlib( 'ii', 'ecosprod',iet_);
        elftv{it}{1} = 'cl1';
        elftv{it}{2} = 'cl2';
        elftv{it}{3} = 'cw';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for i=v_('1'):v_('nd')
            ename = ['aldsq',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ediffsq';
            ielftype(ie) = iet_('ediffsq');
            vname = ['alpha',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('winkel',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['delta',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('minwi',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['bedsq',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ediffsq';
            ielftype(ie) = iet_('ediffsq');
            vname = ['beta',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('winkel',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['delta',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('minwi',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['gadsq',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ediffsq';
            ielftype(ie) = iet_('ediffsq');
            vname = ['gamma',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('winkel',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['delta',int2str(i)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('minwi',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for i=v_('1'):v_('nk')
            ename = ['laeng',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'elaenge';
            ielftype(ie) = iet_('elaenge');
        end
        ename = 'laeng1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'laeng8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('lp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('llaenge',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for i=v_('1'):v_('nd')
            ename = ['sal',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'esklprod';
            ielftype(ie) = iet_('esklprod');
        end
        ename = 'sal1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sal2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sal3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sal4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for i=v_('1'):v_('nd')
            ename = ['cal',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ecosprod';
            ielftype(ie) = iet_('ecosprod');
        end
        ename = 'cal1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'alpha1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cal2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'alpha2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cal3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'alpha3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cal4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'alpha4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for i=v_('1'):v_('nd')
            ename = ['sbe',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'esklprod';
            ielftype(ie) = iet_('esklprod');
        end
        ename = 'sbe1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sbe2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sbe3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'sbe4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('sp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for i=v_('1'):v_('nd')
            ename = ['cbe',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ecosprod';
            ielftype(ie) = iet_('ecosprod');
        end
        ename = 'cbe1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'beta1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cbe2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'beta2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cbe3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'beta3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'cbe4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'l6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'l1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cl2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'beta4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('cw',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for i=v_('1'):v_('nd')
            ename = ['flae',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'evekprod';
            ielftype(ie) = iet_('evekprod');
        end
        ename = 'flae1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'flae2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'flae3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'flae4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'x4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp1y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp2y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'x1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3x',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'y1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('vp3y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gsquare',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('obj1');
        pbm.grftype{ig} = 'gsquare';
        ig = ig_('obj2');
        pbm.grftype{ig} = 'gsquare';
        for i=v_('1'):v_('nd')
            ig = ig_('obj3');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['aldsq',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['bedsq',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['gadsq',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
        end
        for i=v_('1'):v_('nk')
            ig = ig_(['seit',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['laeng',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
        end
        for i=v_('1'):v_('nd')
            ig = ig_(['skal',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['sal',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['cal',int2str(i)]);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['skbe',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['sbe',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['cbe',int2str(i)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        for i=v_('1'):v_('nd')
            ig = ig_(['doppf',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['flae',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN              5.9213448D-4
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AY-41-48';
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

    case 'ediffsq'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'elaenge'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,5);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        U_(3,5) = U_(3,5)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(3)*IV_(3)-IV_(1)*IV_(1)-IV_(2)*IV_(2);
        if(nargout>1)
            g_(3,1) = IV_(3)+IV_(3);
            g_(1,1) = -(IV_(1)+IV_(1));
            g_(2,1) = -(IV_(2)+IV_(2));
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(3,3) = 2.0;
                H_(3,1) = 0.0;
                H_(1,3) = H_(3,1);
                H_(3,2) = 0.0;
                H_(2,3) = H_(3,2);
                H_(1,2) = 0.0;
                H_(2,1) = H_(1,2);
                H_(1,1) = -2.0;
                H_(2,2) = -2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'evekprod'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(4,6);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        U_(3,5) = U_(3,5)+1;
        U_(3,3) = U_(3,3)-1;
        U_(4,6) = U_(4,6)+1;
        U_(4,4) = U_(4,4)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        IV_(4) = U_(4,:)*EV_;
        varargout{1} = IV_(3)*IV_(2)-IV_(1)*IV_(4);
        if(nargout>1)
            g_(1,1) = -IV_(4);
            g_(2,1) = IV_(3);
            g_(3,1) = IV_(2);
            g_(4,1) = -IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(3,2) = 1.0;
                H_(2,3) = H_(3,2);
                H_(1,4) = -1.0;
                H_(4,1) = H_(1,4);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'esklprod'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(4,6);
        U_(1,1) = U_(1,1)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,2) = U_(2,2)+1;
        U_(2,4) = U_(2,4)-1;
        U_(3,5) = U_(3,5)+1;
        U_(3,3) = U_(3,3)-1;
        U_(4,6) = U_(4,6)+1;
        U_(4,4) = U_(4,4)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        IV_(4) = U_(4,:)*EV_;
        varargout{1} = IV_(1)*IV_(3)+IV_(2)*IV_(4);
        if(nargout>1)
            g_(1,1) = IV_(3);
            g_(2,1) = IV_(4);
            g_(3,1) = IV_(1);
            g_(4,1) = IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,3) = 1.0;
                H_(3,1) = H_(1,3);
                H_(2,4) = 1.0;
                H_(4,2) = H_(2,4);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'ecosprod'

        EV_  = varargin{1};
        iel_ = varargin{2};
        cosa = cos(EV_(3));
        sina = sin(EV_(3));
        prod2 = EV_(1)*EV_(2);
        prod3 = prod2*cosa;
        varargout{1} = prod3;
        if(nargout>1)
            g_(1,1) = EV_(2)*cosa;
            g_(2,1) = EV_(1)*cosa;
            g_(3,1) = -sina*prod2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = cosa;
                H_(2,1) = H_(1,2);
                H_(1,3) = -sina*EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = -sina*EV_(1);
                H_(3,2) = H_(2,3);
                H_(3,3) = -prod3;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gsquare'

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

