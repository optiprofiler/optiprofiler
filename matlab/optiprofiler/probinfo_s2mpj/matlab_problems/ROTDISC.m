function varargout = ROTDISC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Optimal design of a rotating disk of minimal weight, with constraints on
%    stress and profile.
%    The problem arise in the mechanical design of turbine where several disks 
%    are assembled, as in jet engines and steam turbines in power generation
%    systems. The data correspond to the real design problem for the engine of 
%    a small civil jet.
%    The problem has a linear objective function, linear constraints and 
%    quadratic equality constraints.  The solution lies at a vertex of the
%    feasible set.
% 
%    Source:
%    B. Apraxine and E. Loute,
%    "The optimal design of a rotating disk: a test case for nonlinear
%    programming codes",
%    Facultes Universitaires Saint Louis, Brussels, 1993.
%    See also:
%    J. P. Nigoghossian,
%    "Problem: the optimization of jet engine discs",
%    in "Optimisation and Design", M. Avriel, M. J. Rijckaert and D. J. Wilde,
%    eds., Prentice Hall, Englewood Cliffs, 1973.
% 
%    SIF input : E. Loute and Ph. L. Toint, April 1993
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-LQR2-RN-905-1081'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 6 X 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ROTDISC';

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
        v_('K') = 180;
        v_('rotspeed') = 22000.0;
        v_('ri') = 41.0;
        v_('ro') = 131.0;
        v_('rim') = 127.0;
        v_('rho') = 8200.0;
        v_('E') = 18525.0;
        v_('epsc') = 14.0e-6;
        v_('nu') = 0.3;
        v_('wro') = 13.0;
        v_('wmin') = 4.0;
        v_('wmax') = 34.0;
        v_('sigmaro') = 18.5;
        v_('sigmari') = 0.0;
        v_('sigmati') = 72.0;
        v_('sigmato') = 47.5;
        v_('sigmatA') = 66.0;
        v_('sigmaru') = 55.0;
        v_('tempo') = 450.0;
        v_('tempn') = 150.0;
        v_('n') = 4;
        v_('ech1') = 100.0;
        v_('ech2') = 1.0E-13;
        v_('ech3') = 1.0E-9;
        v_('1') = 1;
        v_('0') = 0;
        v_('K-1') = -1+v_('K');
        v_('RK') = v_('K');
        v_('Dr') = v_('ro')-v_('ri');
        v_('dr') = v_('Dr')/v_('RK');
        v_('pi') = 3.1415926535;
        v_('2pi') = 2.0*v_('pi');
        v_('aux1') = v_('2pi')*v_('rotspeed');
        v_('60.0') = 60.0;
        v_('omega') = v_('aux1')/v_('60.0');
        v_('omega2') = v_('omega')*v_('omega');
        v_('romg2') = v_('rho')*v_('omega2');
        v_('romg2/2') = 0.5*v_('romg2');
        v_('aux2') = v_('romg2/2')*v_('ech2');
        v_('1+nu') = 1.0+v_('nu');
        v_('3nu') = 3.0*v_('nu');
        v_('1+3nu') = 1.0+v_('3nu');
        v_('3+nu') = 3.0+v_('nu');
        v_('(1+nu)/2') = 0.5*v_('1+nu');
        v_('(1+3nu)/2') = 0.5*v_('1+3nu');
        v_('(3+nu)/2') = 0.5*v_('3+nu');
        v_('pirho') = v_('pi')*v_('rho');
        v_('dr/2') = 0.5*v_('dr');
        v_('aux3') = v_('aux2')*v_('dr');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for k=v_('0'):v_('K')
            [iv,ix_] = s2mpjlib('ii',['w',int2str(k)],ix_);
            pb.xnames{iv} = ['w',int2str(k)];
            [iv,ix_] = s2mpjlib('ii',['sigt',int2str(k)],ix_);
            pb.xnames{iv} = ['sigt',int2str(k)];
            [iv,ix_] = s2mpjlib('ii',['sigr',int2str(k)],ix_);
            pb.xnames{iv} = ['sigr',int2str(k)];
            [iv,ix_] = s2mpjlib('ii',['x',int2str(k)],ix_);
            pb.xnames{iv} = ['x',int2str(k)];
            [iv,ix_] = s2mpjlib('ii',['y',int2str(k)],ix_);
            pb.xnames{iv} = ['y',int2str(k)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        v_('rk') = v_('ri');
        v_('-dr/2') = -1.0*v_('dr/2');
        v_('rk2') = v_('rk')*v_('rk');
        for k=v_('0'):v_('K-1')
            v_('k+1') = 1+k;
            v_('rk+1') = v_('rk')+v_('dr');
            v_('coef1') = v_('aux3')*v_('rk2');
            v_('rk+1sq') = v_('rk+1')*v_('rk+1');
            v_('coef2') = v_('aux3')*v_('rk+1sq');
            [ig,ig_] = s2mpjlib('ii',['SR',int2str(k)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['SR',int2str(k)];
            pbm.gscale(ig,1) = v_('ech1');
            iv = ix_(['w',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef1');
            end
            iv = ix_(['w',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef2');
            end
            v_('tmp1') = v_('(1+3nu)/2')*v_('rk');
            v_('tmp2') = v_('(1+nu)/2')*v_('rk+1');
            v_('tmp3') = v_('tmp1')-v_('tmp2');
            v_('coef3') = v_('tmp3')/v_('rk');
            [ig,ig_] = s2mpjlib('ii',['ST',int2str(k)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ST',int2str(k)];
            iv = ix_(['sigr',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef3')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef3');
            end
            v_('tmp4') = v_('(3+nu)/2')*v_('rk');
            v_('tmp5') = v_('(1+nu)/2')*v_('rk+1');
            v_('tmp6') = v_('tmp5')-v_('tmp4');
            v_('coef4') = v_('tmp6')/v_('rk');
            iv = ix_(['sigt',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef4')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef4');
            end
            v_('tmp7') = v_('(1+3nu)/2')*v_('rk+1');
            v_('tmp8') = v_('(1+nu)/2')*v_('rk');
            v_('tmp9') = v_('tmp8')-v_('tmp7');
            v_('coef5') = v_('tmp9')/v_('rk+1');
            iv = ix_(['sigr',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef5')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef5');
            end
            v_('tmp10') = v_('(3+nu)/2')*v_('rk+1');
            v_('tmp11') = v_('(1+nu)/2')*v_('rk');
            v_('tmp12') = v_('tmp10')-v_('tmp11');
            v_('coef6') = v_('tmp12')/v_('rk+1');
            iv = ix_(['sigt',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('coef6')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('coef6');
            end
            [ig,ig_] = s2mpjlib('ii',['STAy',int2str(k)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['STAy',int2str(k)];
            iv = ix_(['y',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['y',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['STAx',int2str(k)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['STAx',int2str(k)];
            iv = ix_(['x',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['x',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['w',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-dr/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-dr/2');
            end
            iv = ix_(['w',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-dr/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-dr/2');
            end
            v_('rk') = v_('rk+1');
            v_('rk2') = v_('rk+1sq');
        end
        v_('rk-1') = v_('ri');
        v_('rk') = v_('rk-1')+v_('dr');
        v_('rk-1sq') = v_('rk-1')*v_('rk-1');
        v_('rk2') = v_('rk')*v_('rk');
        v_('aux3') = v_('rk-1sq')-v_('rk2');
        v_('aux4') = v_('aux3')*v_('pirho');
        v_('coef1') = v_('aux4')*v_('ech3');
        v_('-coef1') = -1.0*v_('coef1');
        [ig,ig_] = s2mpjlib('ii','WEIGHT',ig_);
        gtype{ig} = '<>';
        iv = ix_(['w',int2str(round(v_('0')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-coef1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-coef1');
        end
        for k=v_('1'):v_('K-1')
            v_('k-1') = -1+k;
            v_('rk+1') = v_('rk')+v_('dr');
            v_('rk+1sq') = v_('rk+1')*v_('rk+1');
            v_('aux3') = v_('rk-1sq')-v_('rk+1sq');
            v_('aux4') = v_('aux3')*v_('pirho');
            v_('coef1') = v_('aux4')*v_('ech3');
            v_('-coef1') = -1.0*v_('coef1');
            [ig,ig_] = s2mpjlib('ii','WEIGHT',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'WEIGHT';
            iv = ix_(['w',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-coef1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-coef1');
            end
            v_('rk-1sq') = v_('rk2');
            v_('rk') = v_('rk+1');
            v_('rk2') = v_('rk+1sq');
        end
        v_('aux3') = v_('rk-1sq')-v_('rk2');
        v_('aux4') = v_('pirho')*v_('aux3');
        v_('coef1') = v_('aux4')*v_('ech3');
        v_('-coef1') = -1.0*v_('coef1');
        [ig,ig_] = s2mpjlib('ii','WEIGHT',ig_);
        gtype{ig} = '<>';
        iv = ix_(['w',int2str(round(v_('K')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-coef1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-coef1');
        end
        for k=v_('0'):v_('K-1')
            v_('k+1') = 1+k;
            [ig,ig_] = s2mpjlib('ii',['SLOP',int2str(k)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['SLOP',int2str(k)];
            iv = ix_(['w',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['w',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['SLOM',int2str(k)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['SLOM',int2str(k)];
            iv = ix_(['w',int2str(round(v_('k+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['w',int2str(k)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        v_('-sigmatA') = -1.0*v_('sigmatA');
        [ig,ig_] = s2mpjlib('ii','AVsigt',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'AVsigt';
        iv = ix_(['y',int2str(round(v_('K')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_(['x',int2str(round(v_('K')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-sigmatA')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-sigmatA');
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
        v_('Eepsc') = v_('E')*v_('epsc');
        v_('rk') = v_('ri');
        v_('aux1') = v_('rk')/v_('ro');
        v_('aux2') = v_('aux1')*v_('aux1');
        v_('aux2') = v_('aux2')*v_('aux2');
        v_('tmp1') = v_('aux2')*v_('tempn');
        v_('tk') = v_('tmp1')+v_('tempo');
        for k=v_('0'):v_('K-1')
            v_('rk+1') = v_('rk')+v_('dr');
            v_('aux1') = v_('rk+1')/v_('ro');
            v_('aux2') = v_('aux1')*v_('aux1');
            v_('aux2') = v_('aux2')*v_('aux2');
            v_('tmp1') = v_('aux2')*v_('tempn');
            v_('tk+1') = v_('tmp1')+v_('tempo');
            v_('tmp2') = v_('tk')-v_('tk+1');
            v_('coef1') = v_('tmp2')*v_('Eepsc');
            pbm.gconst(ig_(['ST',int2str(k)])) = v_('coef1');
            v_('rk') = v_('rk+1');
            v_('tk') = v_('tk+1');
        end
        v_('4dr') = 4.0*v_('dr');
        for k=v_('0'):v_('K-1')
            pbm.gconst(ig_(['SLOP',int2str(k)])) = v_('4dr');
            pbm.gconst(ig_(['SLOM',int2str(k)])) = v_('4dr');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_(['sigr',int2str(round(v_('0')))]),1) = v_('sigmari');
        pb.xupper(ix_(['sigr',int2str(round(v_('0')))]),1) = v_('sigmari');
        pb.xupper(ix_(['sigt',int2str(round(v_('0')))])) = v_('sigmati');
        pb.xlower(ix_(['sigt',int2str(round(v_('0')))]),1) = -100.0;
        pb.xlower(ix_(['x',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['x',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['y',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['y',int2str(round(v_('0')))]),1) = 0.0;
        for k=v_('1'):v_('K')
            pb.xlower(ix_(['x',int2str(k)]),1) = 0.0;
        end
        for k=v_('1'):v_('K-1')
            pb.xlower(ix_(['sigr',int2str(k)]),1) = 0.0;
            pb.xupper(ix_(['sigr',int2str(k)])) = v_('sigmaru');
            pb.xlower(ix_(['sigt',int2str(k)]),1) = -100.0;
            pb.xupper(ix_(['sigt',int2str(k)])) = 100.0;
        end
        v_('K-9') = -9+v_('K');
        for k=v_('0'):v_('K-9')
            pb.xlower(ix_(['w',int2str(k)]),1) = v_('wmin');
            pb.xupper(ix_(['w',int2str(k)])) = v_('wmax');
        end
        v_('K-8') = -8+v_('K');
        for k=v_('K-8'):v_('K')
            pb.xlower(ix_(['w',int2str(k)]),1) = v_('wro');
            pb.xupper(ix_(['w',int2str(k)]),1) = v_('wro');
        end
        pb.xupper(ix_(['sigt',int2str(round(v_('K')))])) = v_('sigmato');
        pb.xlower(ix_(['sigr',int2str(round(v_('K')))]),1) = v_('sigmaro');
        pb.xupper(ix_(['sigr',int2str(round(v_('K')))]),1) = v_('sigmaro');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'w0'))
            pb.x0(ix_('w0'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w0')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt0'))
            pb.x0(ix_('sigt0'),1) = 70.224043123;
        else
            pb.y0(find(pbm.congrps==ig_('sigt0')),1) = 70.224043123;
        end
        if(isKey(ix_,'sigr0'))
            pb.x0(ix_('sigr0'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('sigr0')),1) = 0.0;
        end
        if(isKey(ix_,'x0'))
            pb.x0(ix_('x0'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('x0')),1) = 0.0;
        end
        if(isKey(ix_,'y0'))
            pb.x0(ix_('y0'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('y0')),1) = 0.0;
        end
        if(isKey(ix_,'w1'))
            pb.x0(ix_('w1'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w1')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt1'))
            pb.x0(ix_('sigt1'),1) = 69.337178701;
        else
            pb.y0(find(pbm.congrps==ig_('sigt1')),1) = 69.337178701;
        end
        if(isKey(ix_,'sigr1'))
            pb.x0(ix_('sigr1'),1) = .75150203489;
        else
            pb.y0(find(pbm.congrps==ig_('sigr1')),1) = .75150203489;
        end
        if(isKey(ix_,'x1'))
            pb.x0(ix_('x1'),1) = 15.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x1')),1) = 15.000000000;
        end
        if(isKey(ix_,'y1'))
            pb.x0(ix_('y1'),1) = 1046.7091637;
        else
            pb.y0(find(pbm.congrps==ig_('y1')),1) = 1046.7091637;
        end
        if(isKey(ix_,'w2'))
            pb.x0(ix_('w2'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w2')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt2'))
            pb.x0(ix_('sigt2'),1) = 68.478656583;
        else
            pb.y0(find(pbm.congrps==ig_('sigt2')),1) = 68.478656583;
        end
        if(isKey(ix_,'sigr2'))
            pb.x0(ix_('sigr2'),1) = 1.4725717269;
        else
            pb.y0(find(pbm.congrps==ig_('sigr2')),1) = 1.4725717269;
        end
        if(isKey(ix_,'x2'))
            pb.x0(ix_('x2'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x2')),1) = 30.000000000;
        end
        if(isKey(ix_,'y2'))
            pb.x0(ix_('y2'),1) = 2080.3279283;
        else
            pb.y0(find(pbm.congrps==ig_('y2')),1) = 2080.3279283;
        end
        if(isKey(ix_,'w3'))
            pb.x0(ix_('w3'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w3')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt3'))
            pb.x0(ix_('sigt3'),1) = 67.647086003;
        else
            pb.y0(find(pbm.congrps==ig_('sigt3')),1) = 67.647086003;
        end
        if(isKey(ix_,'sigr3'))
            pb.x0(ix_('sigr3'),1) = 2.1645828149;
        else
            pb.y0(find(pbm.congrps==ig_('sigr3')),1) = 2.1645828149;
        end
        if(isKey(ix_,'x3'))
            pb.x0(ix_('x3'),1) = 45.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x3')),1) = 45.000000000;
        end
        if(isKey(ix_,'y3'))
            pb.x0(ix_('y3'),1) = 3101.2709977;
        else
            pb.y0(find(pbm.congrps==ig_('y3')),1) = 3101.2709977;
        end
        if(isKey(ix_,'w4'))
            pb.x0(ix_('w4'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w4')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt4'))
            pb.x0(ix_('sigt4'),1) = 66.841155637;
        else
            pb.y0(find(pbm.congrps==ig_('sigt4')),1) = 66.841155637;
        end
        if(isKey(ix_,'sigr4'))
            pb.x0(ix_('sigr4'),1) = 2.8288294334;
        else
            pb.y0(find(pbm.congrps==ig_('sigr4')),1) = 2.8288294334;
        end
        if(isKey(ix_,'x4'))
            pb.x0(ix_('x4'),1) = 60.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x4')),1) = 60.000000000;
        end
        if(isKey(ix_,'y4'))
            pb.x0(ix_('y4'),1) = 4109.9328100;
        else
            pb.y0(find(pbm.congrps==ig_('y4')),1) = 4109.9328100;
        end
        if(isKey(ix_,'w5'))
            pb.x0(ix_('w5'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w5')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt5'))
            pb.x0(ix_('sigt5'),1) = 66.059628152;
        else
            pb.y0(find(pbm.congrps==ig_('sigt5')),1) = 66.059628152;
        end
        if(isKey(ix_,'sigr5'))
            pb.x0(ix_('sigr5'),1) = 3.4665315686;
        else
            pb.y0(find(pbm.congrps==ig_('sigr5')),1) = 3.4665315686;
        end
        if(isKey(ix_,'x5'))
            pb.x0(ix_('x5'),1) = 75.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x5')),1) = 75.000000000;
        end
        if(isKey(ix_,'y5'))
            pb.x0(ix_('y5'),1) = 5106.6886884;
        else
            pb.y0(find(pbm.congrps==ig_('y5')),1) = 5106.6886884;
        end
        if(isKey(ix_,'w6'))
            pb.x0(ix_('w6'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w6')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt6'))
            pb.x0(ix_('sigt6'),1) = 65.301335165;
        else
            pb.y0(find(pbm.congrps==ig_('sigt6')),1) = 65.301335165;
        end
        if(isKey(ix_,'sigr6'))
            pb.x0(ix_('sigr6'),1) = 4.0788400842;
        else
            pb.y0(find(pbm.congrps==ig_('sigr6')),1) = 4.0788400842;
        end
        if(isKey(ix_,'x6'))
            pb.x0(ix_('x6'),1) = 90.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('x6')),1) = 90.000000000;
        end
        if(isKey(ix_,'y6'))
            pb.x0(ix_('y6'),1) = 6091.8959133;
        else
            pb.y0(find(pbm.congrps==ig_('y6')),1) = 6091.8959133;
        end
        if(isKey(ix_,'w7'))
            pb.x0(ix_('w7'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w7')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt7'))
            pb.x0(ix_('sigt7'),1) = 64.565172618;
        else
            pb.y0(find(pbm.congrps==ig_('sigt7')),1) = 64.565172618;
        end
        if(isKey(ix_,'sigr7'))
            pb.x0(ix_('sigr7'),1) = 4.6668413532;
        else
            pb.y0(find(pbm.congrps==ig_('sigr7')),1) = 4.6668413532;
        end
        if(isKey(ix_,'x7'))
            pb.x0(ix_('x7'),1) = 105.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x7')),1) = 105.00000000;
        end
        if(isKey(ix_,'y7'))
            pb.x0(ix_('y7'),1) = 7065.8947217;
        else
            pb.y0(find(pbm.congrps==ig_('y7')),1) = 7065.8947217;
        end
        if(isKey(ix_,'w8'))
            pb.x0(ix_('w8'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w8')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt8'))
            pb.x0(ix_('sigt8'),1) = 63.850096500;
        else
            pb.y0(find(pbm.congrps==ig_('sigt8')),1) = 63.850096500;
        end
        if(isKey(ix_,'sigr8'))
            pb.x0(ix_('sigr8'),1) = 5.2315615316;
        else
            pb.y0(find(pbm.congrps==ig_('sigr8')),1) = 5.2315615316;
        end
        if(isKey(ix_,'x8'))
            pb.x0(ix_('x8'),1) = 120.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x8')),1) = 120.00000000;
        end
        if(isKey(ix_,'y8'))
            pb.x0(ix_('y8'),1) = 8029.0092401;
        else
            pb.y0(find(pbm.congrps==ig_('y8')),1) = 8029.0092401;
        end
        if(isKey(ix_,'w9'))
            pb.x0(ix_('w9'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w9')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt9'))
            pb.x0(ix_('sigt9'),1) = 63.155118892;
        else
            pb.y0(find(pbm.congrps==ig_('sigt9')),1) = 63.155118892;
        end
        if(isKey(ix_,'sigr9'))
            pb.x0(ix_('sigr9'),1) = 5.7739705049;
        else
            pb.y0(find(pbm.congrps==ig_('sigr9')),1) = 5.7739705049;
        end
        if(isKey(ix_,'x9'))
            pb.x0(ix_('x9'),1) = 135.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x9')),1) = 135.00000000;
        end
        if(isKey(ix_,'y9'))
            pb.x0(ix_('y9'),1) = 8981.5483555;
        else
            pb.y0(find(pbm.congrps==ig_('y9')),1) = 8981.5483555;
        end
        if(isKey(ix_,'w10'))
            pb.x0(ix_('w10'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w10')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt10'))
            pb.x0(ix_('sigt10'),1) = 62.479304325;
        else
            pb.y0(find(pbm.congrps==ig_('sigt10')),1) = 62.479304325;
        end
        if(isKey(ix_,'sigr10'))
            pb.x0(ix_('sigr10'),1) = 6.2949855370;
        else
            pb.y0(find(pbm.congrps==ig_('sigr10')),1) = 6.2949855370;
        end
        if(isKey(ix_,'x10'))
            pb.x0(ix_('x10'),1) = 150.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x10')),1) = 150.00000000;
        end
        if(isKey(ix_,'y10'))
            pb.x0(ix_('y10'),1) = 9923.8065296;
        else
            pb.y0(find(pbm.congrps==ig_('y10')),1) = 9923.8065296;
        end
        if(isKey(ix_,'w11'))
            pb.x0(ix_('w11'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w11')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt11'))
            pb.x0(ix_('sigt11'),1) = 61.821766398;
        else
            pb.y0(find(pbm.congrps==ig_('sigt11')),1) = 61.821766398;
        end
        if(isKey(ix_,'sigr11'))
            pb.x0(ix_('sigr11'),1) = 6.7954746441;
        else
            pb.y0(find(pbm.congrps==ig_('sigr11')),1) = 6.7954746441;
        end
        if(isKey(ix_,'x11'))
            pb.x0(ix_('x11'),1) = 165.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x11')),1) = 165.00000000;
        end
        if(isKey(ix_,'y11'))
            pb.x0(ix_('y11'),1) = 10856.064560;
        else
            pb.y0(find(pbm.congrps==ig_('y11')),1) = 10856.064560;
        end
        if(isKey(ix_,'w12'))
            pb.x0(ix_('w12'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w12')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt12'))
            pb.x0(ix_('sigt12'),1) = 61.181664651;
        else
            pb.y0(find(pbm.congrps==ig_('sigt12')),1) = 61.181664651;
        end
        if(isKey(ix_,'sigr12'))
            pb.x0(ix_('sigr12'),1) = 7.2762597204;
        else
            pb.y0(find(pbm.congrps==ig_('sigr12')),1) = 7.2762597204;
        end
        if(isKey(ix_,'x12'))
            pb.x0(ix_('x12'),1) = 180.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x12')),1) = 180.00000000;
        end
        if(isKey(ix_,'y12'))
            pb.x0(ix_('y12'),1) = 11778.590293;
        else
            pb.y0(find(pbm.congrps==ig_('y12')),1) = 11778.590293;
        end
        if(isKey(ix_,'w13'))
            pb.x0(ix_('w13'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w13')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt13'))
            pb.x0(ix_('sigt13'),1) = 60.558201672;
        else
            pb.y0(find(pbm.congrps==ig_('sigt13')),1) = 60.558201672;
        end
        if(isKey(ix_,'sigr13'))
            pb.x0(ix_('sigr13'),1) = 7.7381194336;
        else
            pb.y0(find(pbm.congrps==ig_('sigr13')),1) = 7.7381194336;
        end
        if(isKey(ix_,'x13'))
            pb.x0(ix_('x13'),1) = 195.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x13')),1) = 195.00000000;
        end
        if(isKey(ix_,'y13'))
            pb.x0(ix_('y13'),1) = 12691.639290;
        else
            pb.y0(find(pbm.congrps==ig_('y13')),1) = 12691.639290;
        end
        if(isKey(ix_,'w14'))
            pb.x0(ix_('w14'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w14')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt14'))
            pb.x0(ix_('sigt14'),1) = 59.950620407;
        else
            pb.y0(find(pbm.congrps==ig_('sigt14')),1) = 59.950620407;
        end
        if(isKey(ix_,'sigr14'))
            pb.x0(ix_('sigr14'),1) = 8.1817919107;
        else
            pb.y0(find(pbm.congrps==ig_('sigr14')),1) = 8.1817919107;
        end
        if(isKey(ix_,'x14'))
            pb.x0(ix_('x14'),1) = 210.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x14')),1) = 210.00000000;
        end
        if(isKey(ix_,'y14'))
            pb.x0(ix_('y14'),1) = 13595.455456;
        else
            pb.y0(find(pbm.congrps==ig_('y14')),1) = 13595.455456;
        end
        if(isKey(ix_,'w15'))
            pb.x0(ix_('w15'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w15')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt15'))
            pb.x0(ix_('sigt15'),1) = 59.358201664;
        else
            pb.y0(find(pbm.congrps==ig_('sigt15')),1) = 59.358201664;
        end
        if(isKey(ix_,'sigr15'))
            pb.x0(ix_('sigr15'),1) = 8.6079772309;
        else
            pb.y0(find(pbm.congrps==ig_('sigr15')),1) = 8.6079772309;
        end
        if(isKey(ix_,'x15'))
            pb.x0(ix_('x15'),1) = 225.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x15')),1) = 225.00000000;
        end
        if(isKey(ix_,'y15'))
            pb.x0(ix_('y15'),1) = 14490.271621;
        else
            pb.y0(find(pbm.congrps==ig_('y15')),1) = 14490.271621;
        end
        if(isKey(ix_,'w16'))
            pb.x0(ix_('w16'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w16')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt16'))
            pb.x0(ix_('sigt16'),1) = 58.780261800;
        else
            pb.y0(find(pbm.congrps==ig_('sigt16')),1) = 58.780261800;
        end
        if(isKey(ix_,'sigr16'))
            pb.x0(ix_('sigr16'),1) = 9.0173397415;
        else
            pb.y0(find(pbm.congrps==ig_('sigr16')),1) = 9.0173397415;
        end
        if(isKey(ix_,'x16'))
            pb.x0(ix_('x16'),1) = 240.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x16')),1) = 240.00000000;
        end
        if(isKey(ix_,'y16'))
            pb.x0(ix_('y16'),1) = 15376.310097;
        else
            pb.y0(find(pbm.congrps==ig_('y16')),1) = 15376.310097;
        end
        if(isKey(ix_,'w17'))
            pb.x0(ix_('w17'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w17')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt17'))
            pb.x0(ix_('sigt17'),1) = 58.216150564;
        else
            pb.y0(find(pbm.congrps==ig_('sigt17')),1) = 58.216150564;
        end
        if(isKey(ix_,'sigr17'))
            pb.x0(ix_('sigr17'),1) = 9.4105102106;
        else
            pb.y0(find(pbm.congrps==ig_('sigr17')),1) = 9.4105102106;
        end
        if(isKey(ix_,'x17'))
            pb.x0(ix_('x17'),1) = 255.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x17')),1) = 255.00000000;
        end
        if(isKey(ix_,'y17'))
            pb.x0(ix_('y17'),1) = 16253.783190;
        else
            pb.y0(find(pbm.congrps==ig_('y17')),1) = 16253.783190;
        end
        if(isKey(ix_,'w18'))
            pb.x0(ix_('w18'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w18')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt18'))
            pb.x0(ix_('sigt18'),1) = 57.665249095;
        else
            pb.y0(find(pbm.congrps==ig_('sigt18')),1) = 57.665249095;
        end
        if(isKey(ix_,'sigr18'))
            pb.x0(ix_('sigr18'),1) = 9.7880878301;
        else
            pb.y0(find(pbm.congrps==ig_('sigr18')),1) = 9.7880878301;
        end
        if(isKey(ix_,'x18'))
            pb.x0(ix_('x18'),1) = 270.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x18')),1) = 270.00000000;
        end
        if(isKey(ix_,'y18'))
            pb.x0(ix_('y18'),1) = 17122.893688;
        else
            pb.y0(find(pbm.congrps==ig_('y18')),1) = 17122.893688;
        end
        if(isKey(ix_,'w19'))
            pb.x0(ix_('w19'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w19')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt19'))
            pb.x0(ix_('sigt19'),1) = 57.126968056;
        else
            pb.y0(find(pbm.congrps==ig_('sigt19')),1) = 57.126968056;
        end
        if(isKey(ix_,'sigr19'))
            pb.x0(ix_('sigr19'),1) = 10.150642080;
        else
            pb.y0(find(pbm.congrps==ig_('sigr19')),1) = 10.150642080;
        end
        if(isKey(ix_,'x19'))
            pb.x0(ix_('x19'),1) = 285.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x19')),1) = 285.00000000;
        end
        if(isKey(ix_,'y19'))
            pb.x0(ix_('y19'),1) = 17983.835316;
        else
            pb.y0(find(pbm.congrps==ig_('y19')),1) = 17983.835316;
        end
        if(isKey(ix_,'w20'))
            pb.x0(ix_('w20'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w20')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt20'))
            pb.x0(ix_('sigt20'),1) = 56.600745894;
        else
            pb.y0(find(pbm.congrps==ig_('sigt20')),1) = 56.600745894;
        end
        if(isKey(ix_,'sigr20'))
            pb.x0(ix_('sigr20'),1) = 10.498714467;
        else
            pb.y0(find(pbm.congrps==ig_('sigr20')),1) = 10.498714467;
        end
        if(isKey(ix_,'x20'))
            pb.x0(ix_('x20'),1) = 300.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x20')),1) = 300.00000000;
        end
        if(isKey(ix_,'y20'))
            pb.x0(ix_('y20'),1) = 18836.793171;
        else
            pb.y0(find(pbm.congrps==ig_('y20')),1) = 18836.793171;
        end
        if(isKey(ix_,'w21'))
            pb.x0(ix_('w21'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w21')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt21'))
            pb.x0(ix_('sigt21'),1) = 56.086047224;
        else
            pb.y0(find(pbm.congrps==ig_('sigt21')),1) = 56.086047224;
        end
        if(isKey(ix_,'sigr21'))
            pb.x0(ix_('sigr21'),1) = 10.832820143;
        else
            pb.y0(find(pbm.congrps==ig_('sigr21')),1) = 10.832820143;
        end
        if(isKey(ix_,'x21'))
            pb.x0(ix_('x21'),1) = 315.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x21')),1) = 315.00000000;
        end
        if(isKey(ix_,'y21'))
            pb.x0(ix_('y21'),1) = 19681.944119;
        else
            pb.y0(find(pbm.congrps==ig_('y21')),1) = 19681.944119;
        end
        if(isKey(ix_,'w22'))
            pb.x0(ix_('w22'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w22')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt22'))
            pb.x0(ix_('sigt22'),1) = 55.582361311;
        else
            pb.y0(find(pbm.congrps==ig_('sigt22')),1) = 55.582361311;
        end
        if(isKey(ix_,'sigr22'))
            pb.x0(ix_('sigr22'),1) = 11.153449416;
        else
            pb.y0(find(pbm.congrps==ig_('sigr22')),1) = 11.153449416;
        end
        if(isKey(ix_,'x22'))
            pb.x0(ix_('x22'),1) = 330.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x22')),1) = 330.00000000;
        end
        if(isKey(ix_,'y22'))
            pb.x0(ix_('y22'),1) = 20519.457183;
        else
            pb.y0(find(pbm.congrps==ig_('y22')),1) = 20519.457183;
        end
        if(isKey(ix_,'w23'))
            pb.x0(ix_('w23'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w23')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt23'))
            pb.x0(ix_('sigt23'),1) = 55.089200663;
        else
            pb.y0(find(pbm.congrps==ig_('sigt23')),1) = 55.089200663;
        end
        if(isKey(ix_,'sigr23'))
            pb.x0(ix_('sigr23'),1) = 11.461069163;
        else
            pb.y0(find(pbm.congrps==ig_('sigr23')),1) = 11.461069163;
        end
        if(isKey(ix_,'x23'))
            pb.x0(ix_('x23'),1) = 345.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x23')),1) = 345.00000000;
        end
        if(isKey(ix_,'y23'))
            pb.x0(ix_('y23'),1) = 21349.493898;
        else
            pb.y0(find(pbm.congrps==ig_('y23')),1) = 21349.493898;
        end
        if(isKey(ix_,'w24'))
            pb.x0(ix_('w24'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w24')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt24'))
            pb.x0(ix_('sigt24'),1) = 54.606099709;
        else
            pb.y0(find(pbm.congrps==ig_('sigt24')),1) = 54.606099709;
        end
        if(isKey(ix_,'sigr24'))
            pb.x0(ix_('sigr24'),1) = 11.756124148;
        else
            pb.y0(find(pbm.congrps==ig_('sigr24')),1) = 11.756124148;
        end
        if(isKey(ix_,'x24'))
            pb.x0(ix_('x24'),1) = 360.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x24')),1) = 360.00000000;
        end
        if(isKey(ix_,'y24'))
            pb.x0(ix_('y24'),1) = 22172.208651;
        else
            pb.y0(find(pbm.congrps==ig_('y24')),1) = 22172.208651;
        end
        if(isKey(ix_,'w25'))
            pb.x0(ix_('w25'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w25')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt25'))
            pb.x0(ix_('sigt25'),1) = 54.132613569;
        else
            pb.y0(find(pbm.congrps==ig_('sigt25')),1) = 54.132613569;
        end
        if(isKey(ix_,'sigr25'))
            pb.x0(ix_('sigr25'),1) = 12.039038253;
        else
            pb.y0(find(pbm.congrps==ig_('sigr25')),1) = 12.039038253;
        end
        if(isKey(ix_,'x25'))
            pb.x0(ix_('x25'),1) = 375.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x25')),1) = 375.00000000;
        end
        if(isKey(ix_,'y25'))
            pb.x0(ix_('y25'),1) = 22987.749000;
        else
            pb.y0(find(pbm.congrps==ig_('y25')),1) = 22987.749000;
        end
        if(isKey(ix_,'w26'))
            pb.x0(ix_('w26'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w26')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt26'))
            pb.x0(ix_('sigt26'),1) = 53.668316898;
        else
            pb.y0(find(pbm.congrps==ig_('sigt26')),1) = 53.668316898;
        end
        if(isKey(ix_,'sigr26'))
            pb.x0(ix_('sigr26'),1) = 12.310215632;
        else
            pb.y0(find(pbm.congrps==ig_('sigr26')),1) = 12.310215632;
        end
        if(isKey(ix_,'x26'))
            pb.x0(ix_('x26'),1) = 390.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x26')),1) = 390.00000000;
        end
        if(isKey(ix_,'y26'))
            pb.x0(ix_('y26'),1) = 23796.255979;
        else
            pb.y0(find(pbm.congrps==ig_('y26')),1) = 23796.255979;
        end
        if(isKey(ix_,'w27'))
            pb.x0(ix_('w27'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w27')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt27'))
            pb.x0(ix_('sigt27'),1) = 53.212802809;
        else
            pb.y0(find(pbm.congrps==ig_('sigt27')),1) = 53.212802809;
        end
        if(isKey(ix_,'sigr27'))
            pb.x0(ix_('sigr27'),1) = 12.570041791;
        else
            pb.y0(find(pbm.congrps==ig_('sigr27')),1) = 12.570041791;
        end
        if(isKey(ix_,'x27'))
            pb.x0(ix_('x27'),1) = 405.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x27')),1) = 405.00000000;
        end
        if(isKey(ix_,'y27'))
            pb.x0(ix_('y27'),1) = 24597.864377;
        else
            pb.y0(find(pbm.congrps==ig_('y27')),1) = 24597.864377;
        end
        if(isKey(ix_,'w28'))
            pb.x0(ix_('w28'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w28')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt28'))
            pb.x0(ix_('sigt28'),1) = 52.765681856;
        else
            pb.y0(find(pbm.congrps==ig_('sigt28')),1) = 52.765681856;
        end
        if(isKey(ix_,'sigr28'))
            pb.x0(ix_('sigr28'),1) = 12.818884596;
        else
            pb.y0(find(pbm.congrps==ig_('sigr28')),1) = 12.818884596;
        end
        if(isKey(ix_,'x28'))
            pb.x0(ix_('x28'),1) = 420.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x28')),1) = 420.00000000;
        end
        if(isKey(ix_,'y28'))
            pb.x0(ix_('y28'),1) = 25392.703012;
        else
            pb.y0(find(pbm.congrps==ig_('y28')),1) = 25392.703012;
        end
        if(isKey(ix_,'w29'))
            pb.x0(ix_('w29'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w29')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt29'))
            pb.x0(ix_('sigt29'),1) = 52.326581096;
        else
            pb.y0(find(pbm.congrps==ig_('sigt29')),1) = 52.326581096;
        end
        if(isKey(ix_,'sigr29'))
            pb.x0(ix_('sigr29'),1) = 13.057095224;
        else
            pb.y0(find(pbm.congrps==ig_('sigr29')),1) = 13.057095224;
        end
        if(isKey(ix_,'x29'))
            pb.x0(ix_('x29'),1) = 435.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x29')),1) = 435.00000000;
        end
        if(isKey(ix_,'y29'))
            pb.x0(ix_('y29'),1) = 26180.894984;
        else
            pb.y0(find(pbm.congrps==ig_('y29')),1) = 26180.894984;
        end
        if(isKey(ix_,'w30'))
            pb.x0(ix_('w30'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w30')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt30'))
            pb.x0(ix_('sigt30'),1) = 51.895143190;
        else
            pb.y0(find(pbm.congrps==ig_('sigt30')),1) = 51.895143190;
        end
        if(isKey(ix_,'sigr30'))
            pb.x0(ix_('sigr30'),1) = 13.285009049;
        else
            pb.y0(find(pbm.congrps==ig_('sigr30')),1) = 13.285009049;
        end
        if(isKey(ix_,'x30'))
            pb.x0(ix_('x30'),1) = 450.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x30')),1) = 450.00000000;
        end
        if(isKey(ix_,'y30'))
            pb.x0(ix_('y30'),1) = 26962.557916;
        else
            pb.y0(find(pbm.congrps==ig_('y30')),1) = 26962.557916;
        end
        if(isKey(ix_,'w31'))
            pb.x0(ix_('w31'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w31')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt31'))
            pb.x0(ix_('sigt31'),1) = 51.471025575;
        else
            pb.y0(find(pbm.congrps==ig_('sigt31')),1) = 51.471025575;
        end
        if(isKey(ix_,'sigr31'))
            pb.x0(ix_('sigr31'),1) = 13.502946477;
        else
            pb.y0(find(pbm.congrps==ig_('sigr31')),1) = 13.502946477;
        end
        if(isKey(ix_,'x31'))
            pb.x0(ix_('x31'),1) = 465.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x31')),1) = 465.00000000;
        end
        if(isKey(ix_,'y31'))
            pb.x0(ix_('y31'),1) = 27737.804182;
        else
            pb.y0(find(pbm.congrps==ig_('y31')),1) = 27737.804182;
        end
        if(isKey(ix_,'w32'))
            pb.x0(ix_('w32'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w32')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt32'))
            pb.x0(ix_('sigt32'),1) = 51.053899680;
        else
            pb.y0(find(pbm.congrps==ig_('sigt32')),1) = 51.053899680;
        end
        if(isKey(ix_,'sigr32'))
            pb.x0(ix_('sigr32'),1) = 13.711213727;
        else
            pb.y0(find(pbm.congrps==ig_('sigr32')),1) = 13.711213727;
        end
        if(isKey(ix_,'x32'))
            pb.x0(ix_('x32'),1) = 480.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x32')),1) = 480.00000000;
        end
        if(isKey(ix_,'y32'))
            pb.x0(ix_('y32'),1) = 28506.741121;
        else
            pb.y0(find(pbm.congrps==ig_('y32')),1) = 28506.741121;
        end
        if(isKey(ix_,'w33'))
            pb.x0(ix_('w33'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w33')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt33'))
            pb.x0(ix_('sigt33'),1) = 50.643450190;
        else
            pb.y0(find(pbm.congrps==ig_('sigt33')),1) = 50.643450190;
        end
        if(isKey(ix_,'sigr33'))
            pb.x0(ix_('sigr33'),1) = 13.910103570;
        else
            pb.y0(find(pbm.congrps==ig_('sigr33')),1) = 13.910103570;
        end
        if(isKey(ix_,'x33'))
            pb.x0(ix_('x33'),1) = 495.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x33')),1) = 495.00000000;
        end
        if(isKey(ix_,'y33'))
            pb.x0(ix_('y33'),1) = 29269.471245;
        else
            pb.y0(find(pbm.congrps==ig_('y33')),1) = 29269.471245;
        end
        if(isKey(ix_,'w34'))
            pb.x0(ix_('w34'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w34')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt34'))
            pb.x0(ix_('sigt34'),1) = 50.239374354;
        else
            pb.y0(find(pbm.congrps==ig_('sigt34')),1) = 50.239374354;
        end
        if(isKey(ix_,'sigr34'))
            pb.x0(ix_('sigr34'),1) = 14.099896013;
        else
            pb.y0(find(pbm.congrps==ig_('sigr34')),1) = 14.099896013;
        end
        if(isKey(ix_,'x34'))
            pb.x0(ix_('x34'),1) = 510.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x34')),1) = 510.00000000;
        end
        if(isKey(ix_,'y34'))
            pb.x0(ix_('y34'),1) = 30026.092429;
        else
            pb.y0(find(pbm.congrps==ig_('y34')),1) = 30026.092429;
        end
        if(isKey(ix_,'w35'))
            pb.x0(ix_('w35'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w35')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt35'))
            pb.x0(ix_('sigt35'),1) = 49.841381338;
        else
            pb.y0(find(pbm.congrps==ig_('sigt35')),1) = 49.841381338;
        end
        if(isKey(ix_,'sigr35'))
            pb.x0(ix_('sigr35'),1) = 14.280858959;
        else
            pb.y0(find(pbm.congrps==ig_('sigr35')),1) = 14.280858959;
        end
        if(isKey(ix_,'x35'))
            pb.x0(ix_('x35'),1) = 525.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x35')),1) = 525.00000000;
        end
        if(isKey(ix_,'y35'))
            pb.x0(ix_('y35'),1) = 30776.698097;
        else
            pb.y0(find(pbm.congrps==ig_('y35')),1) = 30776.698097;
        end
        if(isKey(ix_,'w36'))
            pb.x0(ix_('w36'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w36')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt36'))
            pb.x0(ix_('sigt36'),1) = 49.449191607;
        else
            pb.y0(find(pbm.congrps==ig_('sigt36')),1) = 49.449191607;
        end
        if(isKey(ix_,'sigr36'))
            pb.x0(ix_('sigr36'),1) = 14.453248807;
        else
            pb.y0(find(pbm.congrps==ig_('sigr36')),1) = 14.453248807;
        end
        if(isKey(ix_,'x36'))
            pb.x0(ix_('x36'),1) = 540.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x36')),1) = 540.00000000;
        end
        if(isKey(ix_,'y36'))
            pb.x0(ix_('y36'),1) = 31521.377394;
        else
            pb.y0(find(pbm.congrps==ig_('y36')),1) = 31521.377394;
        end
        if(isKey(ix_,'w37'))
            pb.x0(ix_('w37'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w37')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt37'))
            pb.x0(ix_('sigt37'),1) = 49.062536356;
        else
            pb.y0(find(pbm.congrps==ig_('sigt37')),1) = 49.062536356;
        end
        if(isKey(ix_,'sigr37'))
            pb.x0(ix_('sigr37'),1) = 14.617311038;
        else
            pb.y0(find(pbm.congrps==ig_('sigr37')),1) = 14.617311038;
        end
        if(isKey(ix_,'x37'))
            pb.x0(ix_('x37'),1) = 555.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x37')),1) = 555.00000000;
        end
        if(isKey(ix_,'y37'))
            pb.x0(ix_('y37'),1) = 32260.215354;
        else
            pb.y0(find(pbm.congrps==ig_('y37')),1) = 32260.215354;
        end
        if(isKey(ix_,'w38'))
            pb.x0(ix_('w38'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w38')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt38'))
            pb.x0(ix_('sigt38'),1) = 48.681156965;
        else
            pb.y0(find(pbm.congrps==ig_('sigt38')),1) = 48.681156965;
        end
        if(isKey(ix_,'sigr38'))
            pb.x0(ix_('sigr38'),1) = 14.773280750;
        else
            pb.y0(find(pbm.congrps==ig_('sigr38')),1) = 14.773280750;
        end
        if(isKey(ix_,'x38'))
            pb.x0(ix_('x38'),1) = 570.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x38')),1) = 570.00000000;
        end
        if(isKey(ix_,'y38'))
            pb.x0(ix_('y38'),1) = 32993.293054;
        else
            pb.y0(find(pbm.congrps==ig_('y38')),1) = 32993.293054;
        end
        if(isKey(ix_,'w39'))
            pb.x0(ix_('w39'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w39')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt39'))
            pb.x0(ix_('sigt39'),1) = 48.304804485;
        else
            pb.y0(find(pbm.congrps==ig_('sigt39')),1) = 48.304804485;
        end
        if(isKey(ix_,'sigr39'))
            pb.x0(ix_('sigr39'),1) = 14.921383174;
        else
            pb.y0(find(pbm.congrps==ig_('sigr39')),1) = 14.921383174;
        end
        if(isKey(ix_,'x39'))
            pb.x0(ix_('x39'),1) = 585.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x39')),1) = 585.00000000;
        end
        if(isKey(ix_,'y39'))
            pb.x0(ix_('y39'),1) = 33720.687765;
        else
            pb.y0(find(pbm.congrps==ig_('y39')),1) = 33720.687765;
        end
        if(isKey(ix_,'w40'))
            pb.x0(ix_('w40'),1) = 30.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w40')),1) = 30.000000000;
        end
        if(isKey(ix_,'sigt40'))
            pb.x0(ix_('sigt40'),1) = 47.933239160;
        else
            pb.y0(find(pbm.congrps==ig_('sigt40')),1) = 47.933239160;
        end
        if(isKey(ix_,'sigr40'))
            pb.x0(ix_('sigr40'),1) = 15.061834152;
        else
            pb.y0(find(pbm.congrps==ig_('sigr40')),1) = 15.061834152;
        end
        if(isKey(ix_,'x40'))
            pb.x0(ix_('x40'),1) = 600.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x40')),1) = 600.00000000;
        end
        if(isKey(ix_,'y40'))
            pb.x0(ix_('y40'),1) = 34442.473092;
        else
            pb.y0(find(pbm.congrps==ig_('y40')),1) = 34442.473092;
        end
        if(isKey(ix_,'w41'))
            pb.x0(ix_('w41'),1) = 29.017500000;
        else
            pb.y0(find(pbm.congrps==ig_('w41')),1) = 29.017500000;
        end
        if(isKey(ix_,'sigt41'))
            pb.x0(ix_('sigt41'),1) = 47.721358592;
        else
            pb.y0(find(pbm.congrps==ig_('sigt41')),1) = 47.721358592;
        end
        if(isKey(ix_,'sigr41'))
            pb.x0(ix_('sigr41'),1) = 15.705670256;
        else
            pb.y0(find(pbm.congrps==ig_('sigr41')),1) = 15.705670256;
        end
        if(isKey(ix_,'x41'))
            pb.x0(ix_('x41'),1) = 614.75437500;
        else
            pb.y0(find(pbm.congrps==ig_('x41')),1) = 614.75437500;
        end
        if(isKey(ix_,'y41'))
            pb.x0(ix_('y41'),1) = 35148.161016;
        else
            pb.y0(find(pbm.congrps==ig_('y41')),1) = 35148.161016;
        end
        if(isKey(ix_,'w42'))
            pb.x0(ix_('w42'),1) = 28.035000000;
        else
            pb.y0(find(pbm.congrps==ig_('w42')),1) = 28.035000000;
        end
        if(isKey(ix_,'sigt42'))
            pb.x0(ix_('sigt42'),1) = 47.528871380;
        else
            pb.y0(find(pbm.congrps==ig_('sigt42')),1) = 47.528871380;
        end
        if(isKey(ix_,'sigr42'))
            pb.x0(ix_('sigr42'),1) = 16.379639595;
        else
            pb.y0(find(pbm.congrps==ig_('sigr42')),1) = 16.379639595;
        end
        if(isKey(ix_,'x42'))
            pb.x0(ix_('x42'),1) = 629.01750000;
        else
            pb.y0(find(pbm.congrps==ig_('x42')),1) = 629.01750000;
        end
        if(isKey(ix_,'y42'))
            pb.x0(ix_('y42'),1) = 35827.467624;
        else
            pb.y0(find(pbm.congrps==ig_('y42')),1) = 35827.467624;
        end
        if(isKey(ix_,'w43'))
            pb.x0(ix_('w43'),1) = 27.052500000;
        else
            pb.y0(find(pbm.congrps==ig_('w43')),1) = 27.052500000;
        end
        if(isKey(ix_,'sigt43'))
            pb.x0(ix_('sigt43'),1) = 47.356912245;
        else
            pb.y0(find(pbm.congrps==ig_('sigt43')),1) = 47.356912245;
        end
        if(isKey(ix_,'sigr43'))
            pb.x0(ix_('sigr43'),1) = 17.087816055;
        else
            pb.y0(find(pbm.congrps==ig_('sigr43')),1) = 17.087816055;
        end
        if(isKey(ix_,'x43'))
            pb.x0(ix_('x43'),1) = 642.78937500;
        else
            pb.y0(find(pbm.congrps==ig_('x43')),1) = 642.78937500;
        end
        if(isKey(ix_,'y43'))
            pb.x0(ix_('y43'),1) = 36480.866319;
        else
            pb.y0(find(pbm.congrps==ig_('y43')),1) = 36480.866319;
        end
        if(isKey(ix_,'w44'))
            pb.x0(ix_('w44'),1) = 26.070000000;
        else
            pb.y0(find(pbm.congrps==ig_('w44')),1) = 26.070000000;
        end
        if(isKey(ix_,'sigt44'))
            pb.x0(ix_('sigt44'),1) = 47.206824226;
        else
            pb.y0(find(pbm.congrps==ig_('sigt44')),1) = 47.206824226;
        end
        if(isKey(ix_,'sigr44'))
            pb.x0(ix_('sigr44'),1) = 17.834855648;
        else
            pb.y0(find(pbm.congrps==ig_('sigr44')),1) = 17.834855648;
        end
        if(isKey(ix_,'x44'))
            pb.x0(ix_('x44'),1) = 656.07000000;
        else
            pb.y0(find(pbm.congrps==ig_('x44')),1) = 656.07000000;
        end
        if(isKey(ix_,'y44'))
            pb.x0(ix_('y44'),1) = 37108.817513;
        else
            pb.y0(find(pbm.congrps==ig_('y44')),1) = 37108.817513;
        end
        if(isKey(ix_,'w45'))
            pb.x0(ix_('w45'),1) = 25.087500000;
        else
            pb.y0(find(pbm.congrps==ig_('w45')),1) = 25.087500000;
        end
        if(isKey(ix_,'sigt45'))
            pb.x0(ix_('sigt45'),1) = 47.080196532;
        else
            pb.y0(find(pbm.congrps==ig_('sigt45')),1) = 47.080196532;
        end
        if(isKey(ix_,'sigr45'))
            pb.x0(ix_('sigr45'),1) = 18.626112102;
        else
            pb.y0(find(pbm.congrps==ig_('sigr45')),1) = 18.626112102;
        end
        if(isKey(ix_,'x45'))
            pb.x0(ix_('x45'),1) = 668.85937500;
        else
            pb.y0(find(pbm.congrps==ig_('x45')),1) = 668.85937500;
        end
        if(isKey(ix_,'y45'))
            pb.x0(ix_('y45'),1) = 37711.769097;
        else
            pb.y0(find(pbm.congrps==ig_('y45')),1) = 37711.769097;
        end
        if(isKey(ix_,'w46'))
            pb.x0(ix_('w46'),1) = 24.105000000;
        else
            pb.y0(find(pbm.congrps==ig_('w46')),1) = 24.105000000;
        end
        if(isKey(ix_,'sigt46'))
            pb.x0(ix_('sigt46'),1) = 46.978911596;
        else
            pb.y0(find(pbm.congrps==ig_('sigt46')),1) = 46.978911596;
        end
        if(isKey(ix_,'sigr46'))
            pb.x0(ix_('sigr46'),1) = 19.467780620;
        else
            pb.y0(find(pbm.congrps==ig_('sigr46')),1) = 19.467780620;
        end
        if(isKey(ix_,'x46'))
            pb.x0(ix_('x46'),1) = 681.15750000;
        else
            pb.y0(find(pbm.congrps==ig_('x46')),1) = 681.15750000;
        end
        if(isKey(ix_,'y46'))
            pb.x0(ix_('y46'),1) = 38290.156871;
        else
            pb.y0(find(pbm.congrps==ig_('y46')),1) = 38290.156871;
        end
        if(isKey(ix_,'w47'))
            pb.x0(ix_('w47'),1) = 23.122500000;
        else
            pb.y0(find(pbm.congrps==ig_('w47')),1) = 23.122500000;
        end
        if(isKey(ix_,'sigt47'))
            pb.x0(ix_('sigt47'),1) = 46.905204055;
        else
            pb.y0(find(pbm.congrps==ig_('sigt47')),1) = 46.905204055;
        end
        if(isKey(ix_,'sigr47'))
            pb.x0(ix_('sigr47'),1) = 20.367078186;
        else
            pb.y0(find(pbm.congrps==ig_('sigr47')),1) = 20.367078186;
        end
        if(isKey(ix_,'x47'))
            pb.x0(ix_('x47'),1) = 692.96437500;
        else
            pb.y0(find(pbm.congrps==ig_('x47')),1) = 692.96437500;
        end
        if(isKey(ix_,'y47'))
            pb.x0(ix_('y47'),1) = 38844.404932;
        else
            pb.y0(find(pbm.congrps==ig_('y47')),1) = 38844.404932;
        end
        if(isKey(ix_,'w48'))
            pb.x0(ix_('w48'),1) = 22.140000000;
        else
            pb.y0(find(pbm.congrps==ig_('w48')),1) = 22.140000000;
        end
        if(isKey(ix_,'sigt48'))
            pb.x0(ix_('sigt48'),1) = 46.861735282;
        else
            pb.y0(find(pbm.congrps==ig_('sigt48')),1) = 46.861735282;
        end
        if(isKey(ix_,'sigr48'))
            pb.x0(ix_('sigr48'),1) = 21.332471779;
        else
            pb.y0(find(pbm.congrps==ig_('sigr48')),1) = 21.332471779;
        end
        if(isKey(ix_,'x48'))
            pb.x0(ix_('x48'),1) = 704.28000000;
        else
            pb.y0(find(pbm.congrps==ig_('x48')),1) = 704.28000000;
        end
        if(isKey(ix_,'y48'))
            pb.x0(ix_('y48'),1) = 39374.926032;
        else
            pb.y0(find(pbm.congrps==ig_('y48')),1) = 39374.926032;
        end
        if(isKey(ix_,'w49'))
            pb.x0(ix_('w49'),1) = 21.372500000;
        else
            pb.y0(find(pbm.congrps==ig_('w49')),1) = 21.372500000;
        end
        if(isKey(ix_,'sigt49'))
            pb.x0(ix_('sigt49'),1) = 46.783637776;
        else
            pb.y0(find(pbm.congrps==ig_('sigt49')),1) = 46.783637776;
        end
        if(isKey(ix_,'sigr49'))
            pb.x0(ix_('sigr49'),1) = 22.149717856;
        else
            pb.y0(find(pbm.congrps==ig_('sigr49')),1) = 22.149717856;
        end
        if(isKey(ix_,'x49'))
            pb.x0(ix_('x49'),1) = 715.15812500;
        else
            pb.y0(find(pbm.congrps==ig_('x49')),1) = 715.15812500;
        end
        if(isKey(ix_,'y49'))
            pb.x0(ix_('y49'),1) = 39884.276561;
        else
            pb.y0(find(pbm.congrps==ig_('y49')),1) = 39884.276561;
        end
        if(isKey(ix_,'w50'))
            pb.x0(ix_('w50'),1) = 20.605000000;
        else
            pb.y0(find(pbm.congrps==ig_('w50')),1) = 20.605000000;
        end
        if(isKey(ix_,'sigt50'))
            pb.x0(ix_('sigt50'),1) = 46.729531640;
        else
            pb.y0(find(pbm.congrps==ig_('sigt50')),1) = 46.729531640;
        end
        if(isKey(ix_,'sigr50'))
            pb.x0(ix_('sigr50'),1) = 23.016346353;
        else
            pb.y0(find(pbm.congrps==ig_('sigr50')),1) = 23.016346353;
        end
        if(isKey(ix_,'x50'))
            pb.x0(ix_('x50'),1) = 725.65250000;
        else
            pb.y0(find(pbm.congrps==ig_('x50')),1) = 725.65250000;
        end
        if(isKey(ix_,'y50'))
            pb.x0(ix_('y50'),1) = 40374.962886;
        else
            pb.y0(find(pbm.congrps==ig_('y50')),1) = 40374.962886;
        end
        if(isKey(ix_,'w51'))
            pb.x0(ix_('w51'),1) = 19.837500000;
        else
            pb.y0(find(pbm.congrps==ig_('w51')),1) = 19.837500000;
        end
        if(isKey(ix_,'sigt51'))
            pb.x0(ix_('sigt51'),1) = 46.701411241;
        else
            pb.y0(find(pbm.congrps==ig_('sigt51')),1) = 46.701411241;
        end
        if(isKey(ix_,'sigr51'))
            pb.x0(ix_('sigr51'),1) = 23.938737423;
        else
            pb.y0(find(pbm.congrps==ig_('sigr51')),1) = 23.938737423;
        end
        if(isKey(ix_,'x51'))
            pb.x0(ix_('x51'),1) = 735.76312500;
        else
            pb.y0(find(pbm.congrps==ig_('x51')),1) = 735.76312500;
        end
        if(isKey(ix_,'y51'))
            pb.x0(ix_('y51'),1) = 40847.288197;
        else
            pb.y0(find(pbm.congrps==ig_('y51')),1) = 40847.288197;
        end
        if(isKey(ix_,'w52'))
            pb.x0(ix_('w52'),1) = 19.070000000;
        else
            pb.y0(find(pbm.congrps==ig_('w52')),1) = 19.070000000;
        end
        if(isKey(ix_,'sigt52'))
            pb.x0(ix_('sigt52'),1) = 46.701615132;
        else
            pb.y0(find(pbm.congrps==ig_('sigt52')),1) = 46.701615132;
        end
        if(isKey(ix_,'sigr52'))
            pb.x0(ix_('sigr52'),1) = 24.924274110;
        else
            pb.y0(find(pbm.congrps==ig_('sigr52')),1) = 24.924274110;
        end
        if(isKey(ix_,'x52'))
            pb.x0(ix_('x52'),1) = 745.49000000;
        else
            pb.y0(find(pbm.congrps==ig_('x52')),1) = 745.49000000;
        end
        if(isKey(ix_,'y52'))
            pb.x0(ix_('y52'),1) = 41301.547959;
        else
            pb.y0(find(pbm.congrps==ig_('y52')),1) = 41301.547959;
        end
        if(isKey(ix_,'w53'))
            pb.x0(ix_('w53'),1) = 18.302500000;
        else
            pb.y0(find(pbm.congrps==ig_('w53')),1) = 18.302500000;
        end
        if(isKey(ix_,'sigt53'))
            pb.x0(ix_('sigt53'),1) = 46.732895222;
        else
            pb.y0(find(pbm.congrps==ig_('sigt53')),1) = 46.732895222;
        end
        if(isKey(ix_,'sigr53'))
            pb.x0(ix_('sigr53'),1) = 25.981553727;
        else
            pb.y0(find(pbm.congrps==ig_('sigr53')),1) = 25.981553727;
        end
        if(isKey(ix_,'x53'))
            pb.x0(ix_('x53'),1) = 754.83312500;
        else
            pb.y0(find(pbm.congrps==ig_('x53')),1) = 754.83312500;
        end
        if(isKey(ix_,'y53'))
            pb.x0(ix_('y53'),1) = 41738.030113;
        else
            pb.y0(find(pbm.congrps==ig_('y53')),1) = 41738.030113;
        end
        if(isKey(ix_,'w54'))
            pb.x0(ix_('w54'),1) = 17.535000000;
        else
            pb.y0(find(pbm.congrps==ig_('w54')),1) = 17.535000000;
        end
        if(isKey(ix_,'sigt54'))
            pb.x0(ix_('sigt54'),1) = 46.798503908;
        else
            pb.y0(find(pbm.congrps==ig_('sigt54')),1) = 46.798503908;
        end
        if(isKey(ix_,'sigr54'))
            pb.x0(ix_('sigr54'),1) = 27.120654675;
        else
            pb.y0(find(pbm.congrps==ig_('sigr54')),1) = 27.120654675;
        end
        if(isKey(ix_,'x54'))
            pb.x0(ix_('x54'),1) = 763.79250000;
        else
            pb.y0(find(pbm.congrps==ig_('x54')),1) = 763.79250000;
        end
        if(isKey(ix_,'y54'))
            pb.x0(ix_('y54'),1) = 42157.015258;
        else
            pb.y0(find(pbm.congrps==ig_('y54')),1) = 42157.015258;
        end
        if(isKey(ix_,'w55'))
            pb.x0(ix_('w55'),1) = 16.767500000;
        else
            pb.y0(find(pbm.congrps==ig_('w55')),1) = 16.767500000;
        end
        if(isKey(ix_,'sigt55'))
            pb.x0(ix_('sigt55'),1) = 46.902304859;
        else
            pb.y0(find(pbm.congrps==ig_('sigt55')),1) = 46.902304859;
        end
        if(isKey(ix_,'sigr55'))
            pb.x0(ix_('sigr55'),1) = 28.353476460;
        else
            pb.y0(find(pbm.congrps==ig_('sigr55')),1) = 28.353476460;
        end
        if(isKey(ix_,'x55'))
            pb.x0(ix_('x55'),1) = 772.36812500;
        else
            pb.y0(find(pbm.congrps==ig_('x55')),1) = 772.36812500;
        end
        if(isKey(ix_,'y55'))
            pb.x0(ix_('y55'),1) = 42558.776798;
        else
            pb.y0(find(pbm.congrps==ig_('y55')),1) = 42558.776798;
        end
        if(isKey(ix_,'w56'))
            pb.x0(ix_('w56'),1) = 16.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w56')),1) = 16.000000000;
        end
        if(isKey(ix_,'sigt56'))
            pb.x0(ix_('sigt56'),1) = 47.048915288;
        else
            pb.y0(find(pbm.congrps==ig_('sigt56')),1) = 47.048915288;
        end
        if(isKey(ix_,'sigr56'))
            pb.x0(ix_('sigr56'),1) = 29.694177498;
        else
            pb.y0(find(pbm.congrps==ig_('sigr56')),1) = 29.694177498;
        end
        if(isKey(ix_,'x56'))
            pb.x0(ix_('x56'),1) = 780.56000000;
        else
            pb.y0(find(pbm.congrps==ig_('x56')),1) = 780.56000000;
        end
        if(isKey(ix_,'y56'))
            pb.x0(ix_('y56'),1) = 42943.581059;
        else
            pb.y0(find(pbm.congrps==ig_('y56')),1) = 42943.581059;
        end
        if(isKey(ix_,'w57'))
            pb.x0(ix_('w57'),1) = 15.440000000;
        else
            pb.y0(find(pbm.congrps==ig_('w57')),1) = 15.440000000;
        end
        if(isKey(ix_,'sigt57'))
            pb.x0(ix_('sigt57'),1) = 47.117144627;
        else
            pb.y0(find(pbm.congrps==ig_('sigt57')),1) = 47.117144627;
        end
        if(isKey(ix_,'sigr57'))
            pb.x0(ix_('sigr57'),1) = 30.741797333;
        else
            pb.y0(find(pbm.congrps==ig_('sigr57')),1) = 30.741797333;
        end
        if(isKey(ix_,'x57'))
            pb.x0(ix_('x57'),1) = 788.42000000;
        else
            pb.y0(find(pbm.congrps==ig_('x57')),1) = 788.42000000;
        end
        if(isKey(ix_,'y57'))
            pb.x0(ix_('y57'),1) = 43313.648898;
        else
            pb.y0(find(pbm.congrps==ig_('y57')),1) = 43313.648898;
        end
        if(isKey(ix_,'w58'))
            pb.x0(ix_('w58'),1) = 14.880000000;
        else
            pb.y0(find(pbm.congrps==ig_('w58')),1) = 14.880000000;
        end
        if(isKey(ix_,'sigt58'))
            pb.x0(ix_('sigt58'),1) = 47.215105046;
        else
            pb.y0(find(pbm.congrps==ig_('sigt58')),1) = 47.215105046;
        end
        if(isKey(ix_,'sigr58'))
            pb.x0(ix_('sigr58'),1) = 31.860061245;
        else
            pb.y0(find(pbm.congrps==ig_('sigr58')),1) = 31.860061245;
        end
        if(isKey(ix_,'x58'))
            pb.x0(ix_('x58'),1) = 796.00000000;
        else
            pb.y0(find(pbm.congrps==ig_('x58')),1) = 796.00000000;
        end
        if(isKey(ix_,'y58'))
            pb.x0(ix_('y58'),1) = 43671.161267;
        else
            pb.y0(find(pbm.congrps==ig_('y58')),1) = 43671.161267;
        end
        if(isKey(ix_,'w59'))
            pb.x0(ix_('w59'),1) = 14.320000000;
        else
            pb.y0(find(pbm.congrps==ig_('w59')),1) = 14.320000000;
        end
        if(isKey(ix_,'sigt59'))
            pb.x0(ix_('sigt59'),1) = 47.345666162;
        else
            pb.y0(find(pbm.congrps==ig_('sigt59')),1) = 47.345666162;
        end
        if(isKey(ix_,'sigr59'))
            pb.x0(ix_('sigr59'),1) = 33.057769694;
        else
            pb.y0(find(pbm.congrps==ig_('sigr59')),1) = 33.057769694;
        end
        if(isKey(ix_,'x59'))
            pb.x0(ix_('x59'),1) = 803.30000000;
        else
            pb.y0(find(pbm.congrps==ig_('x59')),1) = 803.30000000;
        end
        if(isKey(ix_,'y59'))
            pb.x0(ix_('y59'),1) = 44016.298943;
        else
            pb.y0(find(pbm.congrps==ig_('y59')),1) = 44016.298943;
        end
        if(isKey(ix_,'w60'))
            pb.x0(ix_('w60'),1) = 13.760000000;
        else
            pb.y0(find(pbm.congrps==ig_('w60')),1) = 13.760000000;
        end
        if(isKey(ix_,'sigt60'))
            pb.x0(ix_('sigt60'),1) = 47.512175218;
        else
            pb.y0(find(pbm.congrps==ig_('sigt60')),1) = 47.512175218;
        end
        if(isKey(ix_,'sigr60'))
            pb.x0(ix_('sigr60'),1) = 34.345138122;
        else
            pb.y0(find(pbm.congrps==ig_('sigr60')),1) = 34.345138122;
        end
        if(isKey(ix_,'x60'))
            pb.x0(ix_('x60'),1) = 810.32000000;
        else
            pb.y0(find(pbm.congrps==ig_('x60')),1) = 810.32000000;
        end
        if(isKey(ix_,'y60'))
            pb.x0(ix_('y60'),1) = 44349.238310;
        else
            pb.y0(find(pbm.congrps==ig_('y60')),1) = 44349.238310;
        end
        if(isKey(ix_,'w61'))
            pb.x0(ix_('w61'),1) = 13.200000000;
        else
            pb.y0(find(pbm.congrps==ig_('w61')),1) = 13.200000000;
        end
        if(isKey(ix_,'sigt61'))
            pb.x0(ix_('sigt61'),1) = 47.718555318;
        else
            pb.y0(find(pbm.congrps==ig_('sigt61')),1) = 47.718555318;
        end
        if(isKey(ix_,'sigr61'))
            pb.x0(ix_('sigr61'),1) = 35.734097806;
        else
            pb.y0(find(pbm.congrps==ig_('sigr61')),1) = 35.734097806;
        end
        if(isKey(ix_,'x61'))
            pb.x0(ix_('x61'),1) = 817.06000000;
        else
            pb.y0(find(pbm.congrps==ig_('x61')),1) = 817.06000000;
        end
        if(isKey(ix_,'y61'))
            pb.x0(ix_('y61'),1) = 44670.151426;
        else
            pb.y0(find(pbm.congrps==ig_('y61')),1) = 44670.151426;
        end
        if(isKey(ix_,'w62'))
            pb.x0(ix_('w62'),1) = 12.640000000;
        else
            pb.y0(find(pbm.congrps==ig_('w62')),1) = 12.640000000;
        end
        if(isKey(ix_,'sigt62'))
            pb.x0(ix_('sigt62'),1) = 47.969429433;
        else
            pb.y0(find(pbm.congrps==ig_('sigt62')),1) = 47.969429433;
        end
        if(isKey(ix_,'sigr62'))
            pb.x0(ix_('sigr62'),1) = 37.238676641;
        else
            pb.y0(find(pbm.congrps==ig_('sigr62')),1) = 37.238676641;
        end
        if(isKey(ix_,'x62'))
            pb.x0(ix_('x62'),1) = 823.52000000;
        else
            pb.y0(find(pbm.congrps==ig_('x62')),1) = 823.52000000;
        end
        if(isKey(ix_,'y62'))
            pb.x0(ix_('y62'),1) = 44979.206055;
        else
            pb.y0(find(pbm.congrps==ig_('y62')),1) = 44979.206055;
        end
        if(isKey(ix_,'w63'))
            pb.x0(ix_('w63'),1) = 12.080000000;
        else
            pb.y0(find(pbm.congrps==ig_('w63')),1) = 12.080000000;
        end
        if(isKey(ix_,'sigt63'))
            pb.x0(ix_('sigt63'),1) = 48.270278469;
        else
            pb.y0(find(pbm.congrps==ig_('sigt63')),1) = 48.270278469;
        end
        if(isKey(ix_,'sigr63'))
            pb.x0(ix_('sigr63'),1) = 38.875485767;
        else
            pb.y0(find(pbm.congrps==ig_('sigr63')),1) = 38.875485767;
        end
        if(isKey(ix_,'x63'))
            pb.x0(ix_('x63'),1) = 829.70000000;
        else
            pb.y0(find(pbm.congrps==ig_('x63')),1) = 829.70000000;
        end
        if(isKey(ix_,'y63'))
            pb.x0(ix_('y63'),1) = 45276.565693;
        else
            pb.y0(find(pbm.congrps==ig_('y63')),1) = 45276.565693;
        end
        if(isKey(ix_,'w64'))
            pb.x0(ix_('w64'),1) = 11.520000000;
        else
            pb.y0(find(pbm.congrps==ig_('w64')),1) = 11.520000000;
        end
        if(isKey(ix_,'sigt64'))
            pb.x0(ix_('sigt64'),1) = 48.627644833;
        else
            pb.y0(find(pbm.congrps==ig_('sigt64')),1) = 48.627644833;
        end
        if(isKey(ix_,'sigr64'))
            pb.x0(ix_('sigr64'),1) = 40.664348076;
        else
            pb.y0(find(pbm.congrps==ig_('sigr64')),1) = 40.664348076;
        end
        if(isKey(ix_,'x64'))
            pb.x0(ix_('x64'),1) = 835.60000000;
        else
            pb.y0(find(pbm.congrps==ig_('x64')),1) = 835.60000000;
        end
        if(isKey(ix_,'y64'))
            pb.x0(ix_('y64'),1) = 45562.389551;
        else
            pb.y0(find(pbm.congrps==ig_('y64')),1) = 45562.389551;
        end
        if(isKey(ix_,'w65'))
            pb.x0(ix_('w65'),1) = 11.175000000;
        else
            pb.y0(find(pbm.congrps==ig_('w65')),1) = 11.175000000;
        end
        if(isKey(ix_,'sigt65'))
            pb.x0(ix_('sigt65'),1) = 48.801073466;
        else
            pb.y0(find(pbm.congrps==ig_('sigt65')),1) = 48.801073466;
        end
        if(isKey(ix_,'sigr65'))
            pb.x0(ix_('sigr65'),1) = 41.809788300;
        else
            pb.y0(find(pbm.congrps==ig_('sigr65')),1) = 41.809788300;
        end
        if(isKey(ix_,'x65'))
            pb.x0(ix_('x65'),1) = 841.27375000;
        else
            pb.y0(find(pbm.congrps==ig_('x65')),1) = 841.27375000;
        end
        if(isKey(ix_,'y65'))
            pb.x0(ix_('y65'),1) = 45838.775167;
        else
            pb.y0(find(pbm.congrps==ig_('y65')),1) = 45838.775167;
        end
        if(isKey(ix_,'w66'))
            pb.x0(ix_('w66'),1) = 10.830000000;
        else
            pb.y0(find(pbm.congrps==ig_('w66')),1) = 10.830000000;
        end
        if(isKey(ix_,'sigt66'))
            pb.x0(ix_('sigt66'),1) = 49.001998851;
        else
            pb.y0(find(pbm.congrps==ig_('sigt66')),1) = 49.001998851;
        end
        if(isKey(ix_,'sigr66'))
            pb.x0(ix_('sigr66'),1) = 43.023368828;
        else
            pb.y0(find(pbm.congrps==ig_('sigr66')),1) = 43.023368828;
        end
        if(isKey(ix_,'x66'))
            pb.x0(ix_('x66'),1) = 846.77500000;
        else
            pb.y0(find(pbm.congrps==ig_('x66')),1) = 846.77500000;
        end
        if(isKey(ix_,'y66'))
            pb.x0(ix_('y66'),1) = 46107.786078;
        else
            pb.y0(find(pbm.congrps==ig_('y66')),1) = 46107.786078;
        end
        if(isKey(ix_,'w67'))
            pb.x0(ix_('w67'),1) = 10.485000000;
        else
            pb.y0(find(pbm.congrps==ig_('w67')),1) = 10.485000000;
        end
        if(isKey(ix_,'sigt67'))
            pb.x0(ix_('sigt67'),1) = 49.232764164;
        else
            pb.y0(find(pbm.congrps==ig_('sigt67')),1) = 49.232764164;
        end
        if(isKey(ix_,'sigr67'))
            pb.x0(ix_('sigr67'),1) = 44.312155780;
        else
            pb.y0(find(pbm.congrps==ig_('sigr67')),1) = 44.312155780;
        end
        if(isKey(ix_,'x67'))
            pb.x0(ix_('x67'),1) = 852.10375000;
        else
            pb.y0(find(pbm.congrps==ig_('x67')),1) = 852.10375000;
        end
        if(isKey(ix_,'y67'))
            pb.x0(ix_('y67'),1) = 46369.510373;
        else
            pb.y0(find(pbm.congrps==ig_('y67')),1) = 46369.510373;
        end
        if(isKey(ix_,'w68'))
            pb.x0(ix_('w68'),1) = 10.140000000;
        else
            pb.y0(find(pbm.congrps==ig_('w68')),1) = 10.140000000;
        end
        if(isKey(ix_,'sigt68'))
            pb.x0(ix_('sigt68'),1) = 49.496036475;
        else
            pb.y0(find(pbm.congrps==ig_('sigt68')),1) = 49.496036475;
        end
        if(isKey(ix_,'sigr68'))
            pb.x0(ix_('sigr68'),1) = 45.684166505;
        else
            pb.y0(find(pbm.congrps==ig_('sigr68')),1) = 45.684166505;
        end
        if(isKey(ix_,'x68'))
            pb.x0(ix_('x68'),1) = 857.26000000;
        else
            pb.y0(find(pbm.congrps==ig_('x68')),1) = 857.26000000;
        end
        if(isKey(ix_,'y68'))
            pb.x0(ix_('y68'),1) = 46624.034209;
        else
            pb.y0(find(pbm.congrps==ig_('y68')),1) = 46624.034209;
        end
        if(isKey(ix_,'w69'))
            pb.x0(ix_('w69'),1) = 9.7950000000;
        else
            pb.y0(find(pbm.congrps==ig_('w69')),1) = 9.7950000000;
        end
        if(isKey(ix_,'sigt69'))
            pb.x0(ix_('sigt69'),1) = 49.794861993;
        else
            pb.y0(find(pbm.congrps==ig_('sigt69')),1) = 49.794861993;
        end
        if(isKey(ix_,'sigr69'))
            pb.x0(ix_('sigr69'),1) = 47.148537492;
        else
            pb.y0(find(pbm.congrps==ig_('sigr69')),1) = 47.148537492;
        end
        if(isKey(ix_,'x69'))
            pb.x0(ix_('x69'),1) = 862.24375000;
        else
            pb.y0(find(pbm.congrps==ig_('x69')),1) = 862.24375000;
        end
        if(isKey(ix_,'y69'))
            pb.x0(ix_('y69'),1) = 46871.441830;
        else
            pb.y0(find(pbm.congrps==ig_('y69')),1) = 46871.441830;
        end
        if(isKey(ix_,'w70'))
            pb.x0(ix_('w70'),1) = 9.4500000000;
        else
            pb.y0(find(pbm.congrps==ig_('w70')),1) = 9.4500000000;
        end
        if(isKey(ix_,'sigt70'))
            pb.x0(ix_('sigt70'),1) = 50.132733249;
        else
            pb.y0(find(pbm.congrps==ig_('sigt70')),1) = 50.132733249;
        end
        if(isKey(ix_,'sigr70'))
            pb.x0(ix_('sigr70'),1) = 48.715729024;
        else
            pb.y0(find(pbm.congrps==ig_('sigr70')),1) = 48.715729024;
        end
        if(isKey(ix_,'x70'))
            pb.x0(ix_('x70'),1) = 867.05500000;
        else
            pb.y0(find(pbm.congrps==ig_('x70')),1) = 867.05500000;
        end
        if(isKey(ix_,'y70'))
            pb.x0(ix_('y70'),1) = 47111.815580;
        else
            pb.y0(find(pbm.congrps==ig_('y70')),1) = 47111.815580;
        end
        if(isKey(ix_,'w71'))
            pb.x0(ix_('w71'),1) = 9.1050000000;
        else
            pb.y0(find(pbm.congrps==ig_('w71')),1) = 9.1050000000;
        end
        if(isKey(ix_,'sigt71'))
            pb.x0(ix_('sigt71'),1) = 50.513671349;
        else
            pb.y0(find(pbm.congrps==ig_('sigt71')),1) = 50.513671349;
        end
        if(isKey(ix_,'sigr71'))
            pb.x0(ix_('sigr71'),1) = 50.397776340;
        else
            pb.y0(find(pbm.congrps==ig_('sigr71')),1) = 50.397776340;
        end
        if(isKey(ix_,'x71'))
            pb.x0(ix_('x71'),1) = 871.69375000;
        else
            pb.y0(find(pbm.congrps==ig_('x71')),1) = 871.69375000;
        end
        if(isKey(ix_,'y71'))
            pb.x0(ix_('y71'),1) = 47345.235907;
        else
            pb.y0(find(pbm.congrps==ig_('y71')),1) = 47345.235907;
        end
        if(isKey(ix_,'w72'))
            pb.x0(ix_('w72'),1) = 8.7600000000;
        else
            pb.y0(find(pbm.congrps==ig_('w72')),1) = 8.7600000000;
        end
        if(isKey(ix_,'sigt72'))
            pb.x0(ix_('sigt72'),1) = 50.942327392;
        else
            pb.y0(find(pbm.congrps==ig_('sigt72')),1) = 50.942327392;
        end
        if(isKey(ix_,'sigr72'))
            pb.x0(ix_('sigr72'),1) = 52.208600106;
        else
            pb.y0(find(pbm.congrps==ig_('sigr72')),1) = 52.208600106;
        end
        if(isKey(ix_,'x72'))
            pb.x0(ix_('x72'),1) = 876.16000000;
        else
            pb.y0(find(pbm.congrps==ig_('x72')),1) = 876.16000000;
        end
        if(isKey(ix_,'y72'))
            pb.x0(ix_('y72'),1) = 47571.781348;
        else
            pb.y0(find(pbm.congrps==ig_('y72')),1) = 47571.781348;
        end
        if(isKey(ix_,'w73'))
            pb.x0(ix_('w73'),1) = 8.6275000000;
        else
            pb.y0(find(pbm.congrps==ig_('w73')),1) = 8.6275000000;
        end
        if(isKey(ix_,'sigt73'))
            pb.x0(ix_('sigt73'),1) = 51.020192955;
        else
            pb.y0(find(pbm.congrps==ig_('sigt73')),1) = 51.020192955;
        end
        if(isKey(ix_,'sigr73'))
            pb.x0(ix_('sigr73'),1) = 52.831000764;
        else
            pb.y0(find(pbm.congrps==ig_('sigr73')),1) = 52.831000764;
        end
        if(isKey(ix_,'x73'))
            pb.x0(ix_('x73'),1) = 880.50687500;
        else
            pb.y0(find(pbm.congrps==ig_('x73')),1) = 880.50687500;
        end
        if(isKey(ix_,'y73'))
            pb.x0(ix_('y73'),1) = 47793.389224;
        else
            pb.y0(find(pbm.congrps==ig_('y73')),1) = 47793.389224;
        end
        if(isKey(ix_,'w74'))
            pb.x0(ix_('w74'),1) = 8.4950000000;
        else
            pb.y0(find(pbm.congrps==ig_('w74')),1) = 8.4950000000;
        end
        if(isKey(ix_,'sigt74'))
            pb.x0(ix_('sigt74'),1) = 51.105469720;
        else
            pb.y0(find(pbm.congrps==ig_('sigt74')),1) = 51.105469720;
        end
        if(isKey(ix_,'sigr74'))
            pb.x0(ix_('sigr74'),1) = 53.471000483;
        else
            pb.y0(find(pbm.congrps==ig_('sigr74')),1) = 53.471000483;
        end
        if(isKey(ix_,'x74'))
            pb.x0(ix_('x74'),1) = 884.78750000;
        else
            pb.y0(find(pbm.congrps==ig_('x74')),1) = 884.78750000;
        end
        if(isKey(ix_,'y74'))
            pb.x0(ix_('y74'),1) = 48011.968644;
        else
            pb.y0(find(pbm.congrps==ig_('y74')),1) = 48011.968644;
        end
        if(isKey(ix_,'w75'))
            pb.x0(ix_('w75'),1) = 8.3625000000;
        else
            pb.y0(find(pbm.congrps==ig_('w75')),1) = 8.3625000000;
        end
        if(isKey(ix_,'sigt75'))
            pb.x0(ix_('sigt75'),1) = 51.198444013;
        else
            pb.y0(find(pbm.congrps==ig_('sigt75')),1) = 51.198444013;
        end
        if(isKey(ix_,'sigr75'))
            pb.x0(ix_('sigr75'),1) = 54.129557211;
        else
            pb.y0(find(pbm.congrps==ig_('sigr75')),1) = 54.129557211;
        end
        if(isKey(ix_,'x75'))
            pb.x0(ix_('x75'),1) = 889.00187500;
        else
            pb.y0(find(pbm.congrps==ig_('x75')),1) = 889.00187500;
        end
        if(isKey(ix_,'y75'))
            pb.x0(ix_('y75'),1) = 48227.540632;
        else
            pb.y0(find(pbm.congrps==ig_('y75')),1) = 48227.540632;
        end
        if(isKey(ix_,'w76'))
            pb.x0(ix_('w76'),1) = 8.2300000000;
        else
            pb.y0(find(pbm.congrps==ig_('w76')),1) = 8.2300000000;
        end
        if(isKey(ix_,'sigt76'))
            pb.x0(ix_('sigt76'),1) = 51.299424731;
        else
            pb.y0(find(pbm.congrps==ig_('sigt76')),1) = 51.299424731;
        end
        if(isKey(ix_,'sigr76'))
            pb.x0(ix_('sigr76'),1) = 54.807687715;
        else
            pb.y0(find(pbm.congrps==ig_('sigr76')),1) = 54.807687715;
        end
        if(isKey(ix_,'x76'))
            pb.x0(ix_('x76'),1) = 893.15000000;
        else
            pb.y0(find(pbm.congrps==ig_('x76')),1) = 893.15000000;
        end
        if(isKey(ix_,'y76'))
            pb.x0(ix_('y76'),1) = 48440.125946;
        else
            pb.y0(find(pbm.congrps==ig_('y76')),1) = 48440.125946;
        end
        if(isKey(ix_,'w77'))
            pb.x0(ix_('w77'),1) = 8.0975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w77')),1) = 8.0975000000;
        end
        if(isKey(ix_,'sigt77'))
            pb.x0(ix_('sigt77'),1) = 51.408745006;
        else
            pb.y0(find(pbm.congrps==ig_('sigt77')),1) = 51.408745006;
        end
        if(isKey(ix_,'sigr77'))
            pb.x0(ix_('sigr77'),1) = 55.506472516;
        else
            pb.y0(find(pbm.congrps==ig_('sigr77')),1) = 55.506472516;
        end
        if(isKey(ix_,'x77'))
            pb.x0(ix_('x77'),1) = 897.23187500;
        else
            pb.y0(find(pbm.congrps==ig_('x77')),1) = 897.23187500;
        end
        if(isKey(ix_,'y77'))
            pb.x0(ix_('y77'),1) = 48649.745090;
        else
            pb.y0(find(pbm.congrps==ig_('y77')),1) = 48649.745090;
        end
        if(isKey(ix_,'w78'))
            pb.x0(ix_('w78'),1) = 7.9650000000;
        else
            pb.y0(find(pbm.congrps==ig_('w78')),1) = 7.9650000000;
        end
        if(isKey(ix_,'sigt78'))
            pb.x0(ix_('sigt78'),1) = 51.526764037;
        else
            pb.y0(find(pbm.congrps==ig_('sigt78')),1) = 51.526764037;
        end
        if(isKey(ix_,'sigr78'))
            pb.x0(ix_('sigr78'),1) = 56.227061303;
        else
            pb.y0(find(pbm.congrps==ig_('sigr78')),1) = 56.227061303;
        end
        if(isKey(ix_,'x78'))
            pb.x0(ix_('x78'),1) = 901.24750000;
        else
            pb.y0(find(pbm.congrps==ig_('x78')),1) = 901.24750000;
        end
        if(isKey(ix_,'y78'))
            pb.x0(ix_('y78'),1) = 48856.418337;
        else
            pb.y0(find(pbm.congrps==ig_('y78')),1) = 48856.418337;
        end
        if(isKey(ix_,'w79'))
            pb.x0(ix_('w79'),1) = 7.8325000000;
        else
            pb.y0(find(pbm.congrps==ig_('w79')),1) = 7.8325000000;
        end
        if(isKey(ix_,'sigt79'))
            pb.x0(ix_('sigt79'),1) = 51.653869098;
        else
            pb.y0(find(pbm.congrps==ig_('sigt79')),1) = 51.653869098;
        end
        if(isKey(ix_,'sigr79'))
            pb.x0(ix_('sigr79'),1) = 56.970678892;
        else
            pb.y0(find(pbm.congrps==ig_('sigr79')),1) = 56.970678892;
        end
        if(isKey(ix_,'x79'))
            pb.x0(ix_('x79'),1) = 905.19687500;
        else
            pb.y0(find(pbm.congrps==ig_('x79')),1) = 905.19687500;
        end
        if(isKey(ix_,'y79'))
            pb.x0(ix_('y79'),1) = 49060.165739;
        else
            pb.y0(find(pbm.congrps==ig_('y79')),1) = 49060.165739;
        end
        if(isKey(ix_,'w80'))
            pb.x0(ix_('w80'),1) = 7.7000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w80')),1) = 7.7000000000;
        end
        if(isKey(ix_,'sigt80'))
            pb.x0(ix_('sigt80'),1) = 51.790477770;
        else
            pb.y0(find(pbm.congrps==ig_('sigt80')),1) = 51.790477770;
        end
        if(isKey(ix_,'sigr80'))
            pb.x0(ix_('sigr80'),1) = 57.738631802;
        else
            pb.y0(find(pbm.congrps==ig_('sigr80')),1) = 57.738631802;
        end
        if(isKey(ix_,'x80'))
            pb.x0(ix_('x80'),1) = 909.08000000;
        else
            pb.y0(find(pbm.congrps==ig_('x80')),1) = 909.08000000;
        end
        if(isKey(ix_,'y80'))
            pb.x0(ix_('y80'),1) = 49261.007141;
        else
            pb.y0(find(pbm.congrps==ig_('y80')),1) = 49261.007141;
        end
        if(isKey(ix_,'w81'))
            pb.x0(ix_('w81'),1) = 7.6725000000;
        else
            pb.y0(find(pbm.congrps==ig_('w81')),1) = 7.6725000000;
        end
        if(isKey(ix_,'sigt81'))
            pb.x0(ix_('sigt81'),1) = 51.694572152;
        else
            pb.y0(find(pbm.congrps==ig_('sigt81')),1) = 51.694572152;
        end
        if(isKey(ix_,'sigr81'))
            pb.x0(ix_('sigr81'),1) = 57.731509685;
        else
            pb.y0(find(pbm.congrps==ig_('sigr81')),1) = 57.731509685;
        end
        if(isKey(ix_,'x81'))
            pb.x0(ix_('x81'),1) = 912.92312500;
        else
            pb.y0(find(pbm.congrps==ig_('x81')),1) = 912.92312500;
        end
        if(isKey(ix_,'y81'))
            pb.x0(ix_('y81'),1) = 49459.860462;
        else
            pb.y0(find(pbm.congrps==ig_('y81')),1) = 49459.860462;
        end
        if(isKey(ix_,'w82'))
            pb.x0(ix_('w82'),1) = 7.6450000000;
        else
            pb.y0(find(pbm.congrps==ig_('w82')),1) = 7.6450000000;
        end
        if(isKey(ix_,'sigt82'))
            pb.x0(ix_('sigt82'),1) = 51.596245048;
        else
            pb.y0(find(pbm.congrps==ig_('sigt82')),1) = 51.596245048;
        end
        if(isKey(ix_,'sigr82'))
            pb.x0(ix_('sigr82'),1) = 57.723691663;
        else
            pb.y0(find(pbm.congrps==ig_('sigr82')),1) = 57.723691663;
        end
        if(isKey(ix_,'x82'))
            pb.x0(ix_('x82'),1) = 916.75250000;
        else
            pb.y0(find(pbm.congrps==ig_('x82')),1) = 916.75250000;
        end
        if(isKey(ix_,'y82'))
            pb.x0(ix_('y82'),1) = 49657.630436;
        else
            pb.y0(find(pbm.congrps==ig_('y82')),1) = 49657.630436;
        end
        if(isKey(ix_,'w83'))
            pb.x0(ix_('w83'),1) = 7.6175000000;
        else
            pb.y0(find(pbm.congrps==ig_('w83')),1) = 7.6175000000;
        end
        if(isKey(ix_,'sigt83'))
            pb.x0(ix_('sigt83'),1) = 51.495471481;
        else
            pb.y0(find(pbm.congrps==ig_('sigt83')),1) = 51.495471481;
        end
        if(isKey(ix_,'sigr83'))
            pb.x0(ix_('sigr83'),1) = 57.715173744;
        else
            pb.y0(find(pbm.congrps==ig_('sigr83')),1) = 57.715173744;
        end
        if(isKey(ix_,'x83'))
            pb.x0(ix_('x83'),1) = 920.56812500;
        else
            pb.y0(find(pbm.congrps==ig_('x83')),1) = 920.56812500;
        end
        if(isKey(ix_,'y83'))
            pb.x0(ix_('y83'),1) = 49854.310448;
        else
            pb.y0(find(pbm.congrps==ig_('y83')),1) = 49854.310448;
        end
        if(isKey(ix_,'w84'))
            pb.x0(ix_('w84'),1) = 7.5900000000;
        else
            pb.y0(find(pbm.congrps==ig_('w84')),1) = 7.5900000000;
        end
        if(isKey(ix_,'sigt84'))
            pb.x0(ix_('sigt84'),1) = 51.392226289;
        else
            pb.y0(find(pbm.congrps==ig_('sigt84')),1) = 51.392226289;
        end
        if(isKey(ix_,'sigr84'))
            pb.x0(ix_('sigr84'),1) = 57.705951941;
        else
            pb.y0(find(pbm.congrps==ig_('sigr84')),1) = 57.705951941;
        end
        if(isKey(ix_,'x84'))
            pb.x0(ix_('x84'),1) = 924.37000000;
        else
            pb.y0(find(pbm.congrps==ig_('x84')),1) = 924.37000000;
        end
        if(isKey(ix_,'y84'))
            pb.x0(ix_('y84'),1) = 50049.893886;
        else
            pb.y0(find(pbm.congrps==ig_('y84')),1) = 50049.893886;
        end
        if(isKey(ix_,'w85'))
            pb.x0(ix_('w85'),1) = 7.5625000000;
        else
            pb.y0(find(pbm.congrps==ig_('w85')),1) = 7.5625000000;
        end
        if(isKey(ix_,'sigt85'))
            pb.x0(ix_('sigt85'),1) = 51.286484124;
        else
            pb.y0(find(pbm.congrps==ig_('sigt85')),1) = 51.286484124;
        end
        if(isKey(ix_,'sigr85'))
            pb.x0(ix_('sigr85'),1) = 57.696022277;
        else
            pb.y0(find(pbm.congrps==ig_('sigr85')),1) = 57.696022277;
        end
        if(isKey(ix_,'x85'))
            pb.x0(ix_('x85'),1) = 928.15812500;
        else
            pb.y0(find(pbm.congrps==ig_('x85')),1) = 928.15812500;
        end
        if(isKey(ix_,'y85'))
            pb.x0(ix_('y85'),1) = 50244.374144;
        else
            pb.y0(find(pbm.congrps==ig_('y85')),1) = 50244.374144;
        end
        if(isKey(ix_,'w86'))
            pb.x0(ix_('w86'),1) = 7.5350000000;
        else
            pb.y0(find(pbm.congrps==ig_('w86')),1) = 7.5350000000;
        end
        if(isKey(ix_,'sigt86'))
            pb.x0(ix_('sigt86'),1) = 51.178219453;
        else
            pb.y0(find(pbm.congrps==ig_('sigt86')),1) = 51.178219453;
        end
        if(isKey(ix_,'sigr86'))
            pb.x0(ix_('sigr86'),1) = 57.685380778;
        else
            pb.y0(find(pbm.congrps==ig_('sigr86')),1) = 57.685380778;
        end
        if(isKey(ix_,'x86'))
            pb.x0(ix_('x86'),1) = 931.93250000;
        else
            pb.y0(find(pbm.congrps==ig_('x86')),1) = 931.93250000;
        end
        if(isKey(ix_,'y86'))
            pb.x0(ix_('y86'),1) = 50437.744624;
        else
            pb.y0(find(pbm.congrps==ig_('y86')),1) = 50437.744624;
        end
        if(isKey(ix_,'w87'))
            pb.x0(ix_('w87'),1) = 7.5075000000;
        else
            pb.y0(find(pbm.congrps==ig_('w87')),1) = 7.5075000000;
        end
        if(isKey(ix_,'sigt87'))
            pb.x0(ix_('sigt87'),1) = 51.067406560;
        else
            pb.y0(find(pbm.congrps==ig_('sigt87')),1) = 51.067406560;
        end
        if(isKey(ix_,'sigr87'))
            pb.x0(ix_('sigr87'),1) = 57.674023474;
        else
            pb.y0(find(pbm.congrps==ig_('sigr87')),1) = 57.674023474;
        end
        if(isKey(ix_,'x87'))
            pb.x0(ix_('x87'),1) = 935.69312500;
        else
            pb.y0(find(pbm.congrps==ig_('x87')),1) = 935.69312500;
        end
        if(isKey(ix_,'y87'))
            pb.x0(ix_('y87'),1) = 50629.998734;
        else
            pb.y0(find(pbm.congrps==ig_('y87')),1) = 50629.998734;
        end
        if(isKey(ix_,'w88'))
            pb.x0(ix_('w88'),1) = 7.4800000000;
        else
            pb.y0(find(pbm.congrps==ig_('w88')),1) = 7.4800000000;
        end
        if(isKey(ix_,'sigt88'))
            pb.x0(ix_('sigt88'),1) = 50.954019546;
        else
            pb.y0(find(pbm.congrps==ig_('sigt88')),1) = 50.954019546;
        end
        if(isKey(ix_,'sigr88'))
            pb.x0(ix_('sigr88'),1) = 57.661946402;
        else
            pb.y0(find(pbm.congrps==ig_('sigr88')),1) = 57.661946402;
        end
        if(isKey(ix_,'x88'))
            pb.x0(ix_('x88'),1) = 939.44000000;
        else
            pb.y0(find(pbm.congrps==ig_('x88')),1) = 939.44000000;
        end
        if(isKey(ix_,'y88'))
            pb.x0(ix_('y88'),1) = 50821.129889;
        else
            pb.y0(find(pbm.congrps==ig_('y88')),1) = 50821.129889;
        end
        if(isKey(ix_,'w89'))
            pb.x0(ix_('w89'),1) = 7.4450000000;
        else
            pb.y0(find(pbm.congrps==ig_('w89')),1) = 7.4450000000;
        end
        if(isKey(ix_,'sigt89'))
            pb.x0(ix_('sigt89'),1) = 50.855607376;
        else
            pb.y0(find(pbm.congrps==ig_('sigt89')),1) = 50.855607376;
        end
        if(isKey(ix_,'sigr89'))
            pb.x0(ix_('sigr89'),1) = 57.707215988;
        else
            pb.y0(find(pbm.congrps==ig_('sigr89')),1) = 57.707215988;
        end
        if(isKey(ix_,'x89'))
            pb.x0(ix_('x89'),1) = 943.17125000;
        else
            pb.y0(find(pbm.congrps==ig_('x89')),1) = 943.17125000;
        end
        if(isKey(ix_,'y89'))
            pb.x0(ix_('y89'),1) = 51011.068905;
        else
            pb.y0(find(pbm.congrps==ig_('y89')),1) = 51011.068905;
        end
        if(isKey(ix_,'w90'))
            pb.x0(ix_('w90'),1) = 7.4100000000;
        else
            pb.y0(find(pbm.congrps==ig_('w90')),1) = 7.4100000000;
        end
        if(isKey(ix_,'sigt90'))
            pb.x0(ix_('sigt90'),1) = 50.755029735;
        else
            pb.y0(find(pbm.congrps==ig_('sigt90')),1) = 50.755029735;
        end
        if(isKey(ix_,'sigr90'))
            pb.x0(ix_('sigr90'),1) = 57.752273609;
        else
            pb.y0(find(pbm.congrps==ig_('sigr90')),1) = 57.752273609;
        end
        if(isKey(ix_,'x90'))
            pb.x0(ix_('x90'),1) = 946.88500000;
        else
            pb.y0(find(pbm.congrps==ig_('x90')),1) = 946.88500000;
        end
        if(isKey(ix_,'y90'))
            pb.x0(ix_('y90'),1) = 51199.747597;
        else
            pb.y0(find(pbm.congrps==ig_('y90')),1) = 51199.747597;
        end
        if(isKey(ix_,'w91'))
            pb.x0(ix_('w91'),1) = 7.3750000000;
        else
            pb.y0(find(pbm.congrps==ig_('w91')),1) = 7.3750000000;
        end
        if(isKey(ix_,'sigt91'))
            pb.x0(ix_('sigt91'),1) = 50.652260905;
        else
            pb.y0(find(pbm.congrps==ig_('sigt91')),1) = 50.652260905;
        end
        if(isKey(ix_,'sigr91'))
            pb.x0(ix_('sigr91'),1) = 57.797128035;
        else
            pb.y0(find(pbm.congrps==ig_('sigr91')),1) = 57.797128035;
        end
        if(isKey(ix_,'x91'))
            pb.x0(ix_('x91'),1) = 950.58125000;
        else
            pb.y0(find(pbm.congrps==ig_('x91')),1) = 950.58125000;
        end
        if(isKey(ix_,'y91'))
            pb.x0(ix_('y91'),1) = 51387.161395;
        else
            pb.y0(find(pbm.congrps==ig_('y91')),1) = 51387.161395;
        end
        if(isKey(ix_,'w92'))
            pb.x0(ix_('w92'),1) = 7.3400000000;
        else
            pb.y0(find(pbm.congrps==ig_('w92')),1) = 7.3400000000;
        end
        if(isKey(ix_,'sigt92'))
            pb.x0(ix_('sigt92'),1) = 50.547275164;
        else
            pb.y0(find(pbm.congrps==ig_('sigt92')),1) = 50.547275164;
        end
        if(isKey(ix_,'sigr92'))
            pb.x0(ix_('sigr92'),1) = 57.841788143;
        else
            pb.y0(find(pbm.congrps==ig_('sigr92')),1) = 57.841788143;
        end
        if(isKey(ix_,'x92'))
            pb.x0(ix_('x92'),1) = 954.26000000;
        else
            pb.y0(find(pbm.congrps==ig_('x92')),1) = 954.26000000;
        end
        if(isKey(ix_,'y92'))
            pb.x0(ix_('y92'),1) = 51573.305751;
        else
            pb.y0(find(pbm.congrps==ig_('y92')),1) = 51573.305751;
        end
        if(isKey(ix_,'w93'))
            pb.x0(ix_('w93'),1) = 7.3050000000;
        else
            pb.y0(find(pbm.congrps==ig_('w93')),1) = 7.3050000000;
        end
        if(isKey(ix_,'sigt93'))
            pb.x0(ix_('sigt93'),1) = 50.440046782;
        else
            pb.y0(find(pbm.congrps==ig_('sigt93')),1) = 50.440046782;
        end
        if(isKey(ix_,'sigr93'))
            pb.x0(ix_('sigr93'),1) = 57.886262921;
        else
            pb.y0(find(pbm.congrps==ig_('sigr93')),1) = 57.886262921;
        end
        if(isKey(ix_,'x93'))
            pb.x0(ix_('x93'),1) = 957.92125000;
        else
            pb.y0(find(pbm.congrps==ig_('x93')),1) = 957.92125000;
        end
        if(isKey(ix_,'y93'))
            pb.x0(ix_('y93'),1) = 51758.176137;
        else
            pb.y0(find(pbm.congrps==ig_('y93')),1) = 51758.176137;
        end
        if(isKey(ix_,'w94'))
            pb.x0(ix_('w94'),1) = 7.2700000000;
        else
            pb.y0(find(pbm.congrps==ig_('w94')),1) = 7.2700000000;
        end
        if(isKey(ix_,'sigt94'))
            pb.x0(ix_('sigt94'),1) = 50.330550025;
        else
            pb.y0(find(pbm.congrps==ig_('sigt94')),1) = 50.330550025;
        end
        if(isKey(ix_,'sigr94'))
            pb.x0(ix_('sigr94'),1) = 57.930561476;
        else
            pb.y0(find(pbm.congrps==ig_('sigr94')),1) = 57.930561476;
        end
        if(isKey(ix_,'x94'))
            pb.x0(ix_('x94'),1) = 961.56500000;
        else
            pb.y0(find(pbm.congrps==ig_('x94')),1) = 961.56500000;
        end
        if(isKey(ix_,'y94'))
            pb.x0(ix_('y94'),1) = 51941.768047;
        else
            pb.y0(find(pbm.congrps==ig_('y94')),1) = 51941.768047;
        end
        if(isKey(ix_,'w95'))
            pb.x0(ix_('w95'),1) = 7.2350000000;
        else
            pb.y0(find(pbm.congrps==ig_('w95')),1) = 7.2350000000;
        end
        if(isKey(ix_,'sigt95'))
            pb.x0(ix_('sigt95'),1) = 50.218759151;
        else
            pb.y0(find(pbm.congrps==ig_('sigt95')),1) = 50.218759151;
        end
        if(isKey(ix_,'sigr95'))
            pb.x0(ix_('sigr95'),1) = 57.974693041;
        else
            pb.y0(find(pbm.congrps==ig_('sigr95')),1) = 57.974693041;
        end
        if(isKey(ix_,'x95'))
            pb.x0(ix_('x95'),1) = 965.19125000;
        else
            pb.y0(find(pbm.congrps==ig_('x95')),1) = 965.19125000;
        end
        if(isKey(ix_,'y95'))
            pb.x0(ix_('y95'),1) = 52124.077002;
        else
            pb.y0(find(pbm.congrps==ig_('y95')),1) = 52124.077002;
        end
        if(isKey(ix_,'w96'))
            pb.x0(ix_('w96'),1) = 7.2000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w96')),1) = 7.2000000000;
        end
        if(isKey(ix_,'sigt96'))
            pb.x0(ix_('sigt96'),1) = 50.104648409;
        else
            pb.y0(find(pbm.congrps==ig_('sigt96')),1) = 50.104648409;
        end
        if(isKey(ix_,'sigr96'))
            pb.x0(ix_('sigr96'),1) = 58.018666983;
        else
            pb.y0(find(pbm.congrps==ig_('sigr96')),1) = 58.018666983;
        end
        if(isKey(ix_,'x96'))
            pb.x0(ix_('x96'),1) = 968.80000000;
        else
            pb.y0(find(pbm.congrps==ig_('x96')),1) = 968.80000000;
        end
        if(isKey(ix_,'y96'))
            pb.x0(ix_('y96'),1) = 52305.098550;
        else
            pb.y0(find(pbm.congrps==ig_('y96')),1) = 52305.098550;
        end
        if(isKey(ix_,'w97'))
            pb.x0(ix_('w97'),1) = 7.1712500000;
        else
            pb.y0(find(pbm.congrps==ig_('w97')),1) = 7.1712500000;
        end
        if(isKey(ix_,'sigt97'))
            pb.x0(ix_('sigt97'),1) = 49.972881021;
        else
            pb.y0(find(pbm.congrps==ig_('sigt97')),1) = 49.972881021;
        end
        if(isKey(ix_,'sigr97'))
            pb.x0(ix_('sigr97'),1) = 58.011883331;
        else
            pb.y0(find(pbm.congrps==ig_('sigr97')),1) = 58.011883331;
        end
        if(isKey(ix_,'x97'))
            pb.x0(ix_('x97'),1) = 972.39281250;
        else
            pb.y0(find(pbm.congrps==ig_('x97')),1) = 972.39281250;
        end
        if(isKey(ix_,'y97'))
            pb.x0(ix_('y97'),1) = 52484.878923;
        else
            pb.y0(find(pbm.congrps==ig_('y97')),1) = 52484.878923;
        end
        if(isKey(ix_,'w98'))
            pb.x0(ix_('w98'),1) = 7.1425000000;
        else
            pb.y0(find(pbm.congrps==ig_('w98')),1) = 7.1425000000;
        end
        if(isKey(ix_,'sigt98'))
            pb.x0(ix_('sigt98'),1) = 49.838337287;
        else
            pb.y0(find(pbm.congrps==ig_('sigt98')),1) = 49.838337287;
        end
        if(isKey(ix_,'sigr98'))
            pb.x0(ix_('sigr98'),1) = 58.004462269;
        else
            pb.y0(find(pbm.congrps==ig_('sigr98')),1) = 58.004462269;
        end
        if(isKey(ix_,'x98'))
            pb.x0(ix_('x98'),1) = 975.97125000;
        else
            pb.y0(find(pbm.congrps==ig_('x98')),1) = 975.97125000;
        end
        if(isKey(ix_,'y98'))
            pb.x0(ix_('y98'),1) = 52663.463510;
        else
            pb.y0(find(pbm.congrps==ig_('y98')),1) = 52663.463510;
        end
        if(isKey(ix_,'w99'))
            pb.x0(ix_('w99'),1) = 7.1137500000;
        else
            pb.y0(find(pbm.congrps==ig_('w99')),1) = 7.1137500000;
        end
        if(isKey(ix_,'sigt99'))
            pb.x0(ix_('sigt99'),1) = 49.700989951;
        else
            pb.y0(find(pbm.congrps==ig_('sigt99')),1) = 49.700989951;
        end
        if(isKey(ix_,'sigr99'))
            pb.x0(ix_('sigr99'),1) = 57.996401618;
        else
            pb.y0(find(pbm.congrps==ig_('sigr99')),1) = 57.996401618;
        end
        if(isKey(ix_,'x99'))
            pb.x0(ix_('x99'),1) = 979.53531250;
        else
            pb.y0(find(pbm.congrps==ig_('x99')),1) = 979.53531250;
        end
        if(isKey(ix_,'y99'))
            pb.x0(ix_('y99'),1) = 52840.846195;
        else
            pb.y0(find(pbm.congrps==ig_('y99')),1) = 52840.846195;
        end
        if(isKey(ix_,'w100'))
            pb.x0(ix_('w100'),1) = 7.0850000000;
        else
            pb.y0(find(pbm.congrps==ig_('w100')),1) = 7.0850000000;
        end
        if(isKey(ix_,'sigt100'))
            pb.x0(ix_('sigt100'),1) = 49.560811600;
        else
            pb.y0(find(pbm.congrps==ig_('sigt100')),1) = 49.560811600;
        end
        if(isKey(ix_,'sigr100'))
            pb.x0(ix_('sigr100'),1) = 57.987699220;
        else
            pb.y0(find(pbm.congrps==ig_('sigr100')),1) = 57.987699220;
        end
        if(isKey(ix_,'x100'))
            pb.x0(ix_('x100'),1) = 983.08500000;
        else
            pb.y0(find(pbm.congrps==ig_('x100')),1) = 983.08500000;
        end
        if(isKey(ix_,'y100'))
            pb.x0(ix_('y100'),1) = 53017.020887;
        else
            pb.y0(find(pbm.congrps==ig_('y100')),1) = 53017.020887;
        end
        if(isKey(ix_,'w101'))
            pb.x0(ix_('w101'),1) = 7.0562500000;
        else
            pb.y0(find(pbm.congrps==ig_('w101')),1) = 7.0562500000;
        end
        if(isKey(ix_,'sigt101'))
            pb.x0(ix_('sigt101'),1) = 49.417774661;
        else
            pb.y0(find(pbm.congrps==ig_('sigt101')),1) = 49.417774661;
        end
        if(isKey(ix_,'sigr101'))
            pb.x0(ix_('sigr101'),1) = 57.978352944;
        else
            pb.y0(find(pbm.congrps==ig_('sigr101')),1) = 57.978352944;
        end
        if(isKey(ix_,'x101'))
            pb.x0(ix_('x101'),1) = 986.62031250;
        else
            pb.y0(find(pbm.congrps==ig_('x101')),1) = 986.62031250;
        end
        if(isKey(ix_,'y101'))
            pb.x0(ix_('y101'),1) = 53191.981517;
        else
            pb.y0(find(pbm.congrps==ig_('y101')),1) = 53191.981517;
        end
        if(isKey(ix_,'w102'))
            pb.x0(ix_('w102'),1) = 7.0275000000;
        else
            pb.y0(find(pbm.congrps==ig_('w102')),1) = 7.0275000000;
        end
        if(isKey(ix_,'sigt102'))
            pb.x0(ix_('sigt102'),1) = 49.271851404;
        else
            pb.y0(find(pbm.congrps==ig_('sigt102')),1) = 49.271851404;
        end
        if(isKey(ix_,'sigr102'))
            pb.x0(ix_('sigr102'),1) = 57.968360681;
        else
            pb.y0(find(pbm.congrps==ig_('sigr102')),1) = 57.968360681;
        end
        if(isKey(ix_,'x102'))
            pb.x0(ix_('x102'),1) = 990.14125000;
        else
            pb.y0(find(pbm.congrps==ig_('x102')),1) = 990.14125000;
        end
        if(isKey(ix_,'y102'))
            pb.x0(ix_('y102'),1) = 53365.722044;
        else
            pb.y0(find(pbm.congrps==ig_('y102')),1) = 53365.722044;
        end
        if(isKey(ix_,'w103'))
            pb.x0(ix_('w103'),1) = 6.9987500000;
        else
            pb.y0(find(pbm.congrps==ig_('w103')),1) = 6.9987500000;
        end
        if(isKey(ix_,'sigt103'))
            pb.x0(ix_('sigt103'),1) = 49.123013943;
        else
            pb.y0(find(pbm.congrps==ig_('sigt103')),1) = 49.123013943;
        end
        if(isKey(ix_,'sigr103'))
            pb.x0(ix_('sigr103'),1) = 57.957720351;
        else
            pb.y0(find(pbm.congrps==ig_('sigr103')),1) = 57.957720351;
        end
        if(isKey(ix_,'x103'))
            pb.x0(ix_('x103'),1) = 993.64781250;
        else
            pb.y0(find(pbm.congrps==ig_('x103')),1) = 993.64781250;
        end
        if(isKey(ix_,'y103'))
            pb.x0(ix_('y103'),1) = 53538.236452;
        else
            pb.y0(find(pbm.congrps==ig_('y103')),1) = 53538.236452;
        end
        if(isKey(ix_,'w104'))
            pb.x0(ix_('w104'),1) = 6.9700000000;
        else
            pb.y0(find(pbm.congrps==ig_('w104')),1) = 6.9700000000;
        end
        if(isKey(ix_,'sigt104'))
            pb.x0(ix_('sigt104'),1) = 48.971234234;
        else
            pb.y0(find(pbm.congrps==ig_('sigt104')),1) = 48.971234234;
        end
        if(isKey(ix_,'sigr104'))
            pb.x0(ix_('sigr104'),1) = 57.946429896;
        else
            pb.y0(find(pbm.congrps==ig_('sigr104')),1) = 57.946429896;
        end
        if(isKey(ix_,'x104'))
            pb.x0(ix_('x104'),1) = 997.14000000;
        else
            pb.y0(find(pbm.congrps==ig_('x104')),1) = 997.14000000;
        end
        if(isKey(ix_,'y104'))
            pb.x0(ix_('y104'),1) = 53709.518751;
        else
            pb.y0(find(pbm.congrps==ig_('y104')),1) = 53709.518751;
        end
        if(isKey(ix_,'w105'))
            pb.x0(ix_('w105'),1) = 6.9412500000;
        else
            pb.y0(find(pbm.congrps==ig_('w105')),1) = 6.9412500000;
        end
        if(isKey(ix_,'sigt105'))
            pb.x0(ix_('sigt105'),1) = 48.816484081;
        else
            pb.y0(find(pbm.congrps==ig_('sigt105')),1) = 48.816484081;
        end
        if(isKey(ix_,'sigr105'))
            pb.x0(ix_('sigr105'),1) = 57.934487286;
        else
            pb.y0(find(pbm.congrps==ig_('sigr105')),1) = 57.934487286;
        end
        if(isKey(ix_,'x105'))
            pb.x0(ix_('x105'),1) = 1000.6178125;
        else
            pb.y0(find(pbm.congrps==ig_('x105')),1) = 1000.6178125;
        end
        if(isKey(ix_,'y105'))
            pb.x0(ix_('y105'),1) = 53879.562982;
        else
            pb.y0(find(pbm.congrps==ig_('y105')),1) = 53879.562982;
        end
        if(isKey(ix_,'w106'))
            pb.x0(ix_('w106'),1) = 6.9125000000;
        else
            pb.y0(find(pbm.congrps==ig_('w106')),1) = 6.9125000000;
        end
        if(isKey(ix_,'sigt106'))
            pb.x0(ix_('sigt106'),1) = 48.658735129;
        else
            pb.y0(find(pbm.congrps==ig_('sigt106')),1) = 48.658735129;
        end
        if(isKey(ix_,'sigr106'))
            pb.x0(ix_('sigr106'),1) = 57.921890519;
        else
            pb.y0(find(pbm.congrps==ig_('sigr106')),1) = 57.921890519;
        end
        if(isKey(ix_,'x106'))
            pb.x0(ix_('x106'),1) = 1004.0812500;
        else
            pb.y0(find(pbm.congrps==ig_('x106')),1) = 1004.0812500;
        end
        if(isKey(ix_,'y106'))
            pb.x0(ix_('y106'),1) = 54048.363213;
        else
            pb.y0(find(pbm.congrps==ig_('y106')),1) = 54048.363213;
        end
        if(isKey(ix_,'w107'))
            pb.x0(ix_('w107'),1) = 6.8837500000;
        else
            pb.y0(find(pbm.congrps==ig_('w107')),1) = 6.8837500000;
        end
        if(isKey(ix_,'sigt107'))
            pb.x0(ix_('sigt107'),1) = 48.497958873;
        else
            pb.y0(find(pbm.congrps==ig_('sigt107')),1) = 48.497958873;
        end
        if(isKey(ix_,'sigr107'))
            pb.x0(ix_('sigr107'),1) = 57.908637618;
        else
            pb.y0(find(pbm.congrps==ig_('sigr107')),1) = 57.908637618;
        end
        if(isKey(ix_,'x107'))
            pb.x0(ix_('x107'),1) = 1007.5303125;
        else
            pb.y0(find(pbm.congrps==ig_('x107')),1) = 1007.5303125;
        end
        if(isKey(ix_,'y107'))
            pb.x0(ix_('y107'),1) = 54215.913546;
        else
            pb.y0(find(pbm.congrps==ig_('y107')),1) = 54215.913546;
        end
        if(isKey(ix_,'w108'))
            pb.x0(ix_('w108'),1) = 6.8550000000;
        else
            pb.y0(find(pbm.congrps==ig_('w108')),1) = 6.8550000000;
        end
        if(isKey(ix_,'sigt108'))
            pb.x0(ix_('sigt108'),1) = 48.334126653;
        else
            pb.y0(find(pbm.congrps==ig_('sigt108')),1) = 48.334126653;
        end
        if(isKey(ix_,'sigr108'))
            pb.x0(ix_('sigr108'),1) = 57.894726636;
        else
            pb.y0(find(pbm.congrps==ig_('sigr108')),1) = 57.894726636;
        end
        if(isKey(ix_,'x108'))
            pb.x0(ix_('x108'),1) = 1010.9650000;
        else
            pb.y0(find(pbm.congrps==ig_('x108')),1) = 1010.9650000;
        end
        if(isKey(ix_,'y108'))
            pb.x0(ix_('y108'),1) = 54382.208112;
        else
            pb.y0(find(pbm.congrps==ig_('y108')),1) = 54382.208112;
        end
        if(isKey(ix_,'w109'))
            pb.x0(ix_('w109'),1) = 6.8262500000;
        else
            pb.y0(find(pbm.congrps==ig_('w109')),1) = 6.8262500000;
        end
        if(isKey(ix_,'sigt109'))
            pb.x0(ix_('sigt109'),1) = 48.167209655;
        else
            pb.y0(find(pbm.congrps==ig_('sigt109')),1) = 48.167209655;
        end
        if(isKey(ix_,'sigr109'))
            pb.x0(ix_('sigr109'),1) = 57.880155655;
        else
            pb.y0(find(pbm.congrps==ig_('sigr109')),1) = 57.880155655;
        end
        if(isKey(ix_,'x109'))
            pb.x0(ix_('x109'),1) = 1014.3853125;
        else
            pb.y0(find(pbm.congrps==ig_('x109')),1) = 1014.3853125;
        end
        if(isKey(ix_,'y109'))
            pb.x0(ix_('y109'),1) = 54547.241075;
        else
            pb.y0(find(pbm.congrps==ig_('y109')),1) = 54547.241075;
        end
        if(isKey(ix_,'w110'))
            pb.x0(ix_('w110'),1) = 6.7975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w110')),1) = 6.7975000000;
        end
        if(isKey(ix_,'sigt110'))
            pb.x0(ix_('sigt110'),1) = 47.997178915;
        else
            pb.y0(find(pbm.congrps==ig_('sigt110')),1) = 47.997178915;
        end
        if(isKey(ix_,'sigr110'))
            pb.x0(ix_('sigr110'),1) = 57.864922786;
        else
            pb.y0(find(pbm.congrps==ig_('sigr110')),1) = 57.864922786;
        end
        if(isKey(ix_,'x110'))
            pb.x0(ix_('x110'),1) = 1017.7912500;
        else
            pb.y0(find(pbm.congrps==ig_('x110')),1) = 1017.7912500;
        end
        if(isKey(ix_,'y110'))
            pb.x0(ix_('y110'),1) = 54711.006635;
        else
            pb.y0(find(pbm.congrps==ig_('y110')),1) = 54711.006635;
        end
        if(isKey(ix_,'w111'))
            pb.x0(ix_('w111'),1) = 6.7687500000;
        else
            pb.y0(find(pbm.congrps==ig_('w111')),1) = 6.7687500000;
        end
        if(isKey(ix_,'sigt111'))
            pb.x0(ix_('sigt111'),1) = 47.824005317;
        else
            pb.y0(find(pbm.congrps==ig_('sigt111')),1) = 47.824005317;
        end
        if(isKey(ix_,'sigr111'))
            pb.x0(ix_('sigr111'),1) = 57.849026171;
        else
            pb.y0(find(pbm.congrps==ig_('sigr111')),1) = 57.849026171;
        end
        if(isKey(ix_,'x111'))
            pb.x0(ix_('x111'),1) = 1021.1828125;
        else
            pb.y0(find(pbm.congrps==ig_('x111')),1) = 1021.1828125;
        end
        if(isKey(ix_,'y111'))
            pb.x0(ix_('y111'),1) = 54873.499025;
        else
            pb.y0(find(pbm.congrps==ig_('y111')),1) = 54873.499025;
        end
        if(isKey(ix_,'w112'))
            pb.x0(ix_('w112'),1) = 6.7400000000;
        else
            pb.y0(find(pbm.congrps==ig_('w112')),1) = 6.7400000000;
        end
        if(isKey(ix_,'sigt112'))
            pb.x0(ix_('sigt112'),1) = 47.647659596;
        else
            pb.y0(find(pbm.congrps==ig_('sigt112')),1) = 47.647659596;
        end
        if(isKey(ix_,'sigr112'))
            pb.x0(ix_('sigr112'),1) = 57.832463983;
        else
            pb.y0(find(pbm.congrps==ig_('sigr112')),1) = 57.832463983;
        end
        if(isKey(ix_,'x112'))
            pb.x0(ix_('x112'),1) = 1024.5600000;
        else
            pb.y0(find(pbm.congrps==ig_('x112')),1) = 1024.5600000;
        end
        if(isKey(ix_,'y112'))
            pb.x0(ix_('y112'),1) = 55034.712515;
        else
            pb.y0(find(pbm.congrps==ig_('y112')),1) = 55034.712515;
        end
        if(isKey(ix_,'w113'))
            pb.x0(ix_('w113'),1) = 6.7187500000;
        else
            pb.y0(find(pbm.congrps==ig_('w113')),1) = 6.7187500000;
        end
        if(isKey(ix_,'sigt113'))
            pb.x0(ix_('sigt113'),1) = 47.448591006;
        else
            pb.y0(find(pbm.congrps==ig_('sigt113')),1) = 47.448591006;
        end
        if(isKey(ix_,'sigr113'))
            pb.x0(ix_('sigr113'),1) = 57.750663878;
        else
            pb.y0(find(pbm.congrps==ig_('sigr113')),1) = 57.750663878;
        end
        if(isKey(ix_,'x113'))
            pb.x0(ix_('x113'),1) = 1027.9246875;
        else
            pb.y0(find(pbm.congrps==ig_('x113')),1) = 1027.9246875;
        end
        if(isKey(ix_,'y113'))
            pb.x0(ix_('y113'),1) = 55194.697627;
        else
            pb.y0(find(pbm.congrps==ig_('y113')),1) = 55194.697627;
        end
        if(isKey(ix_,'w114'))
            pb.x0(ix_('w114'),1) = 6.6975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w114')),1) = 6.6975000000;
        end
        if(isKey(ix_,'sigt114'))
            pb.x0(ix_('sigt114'),1) = 47.245860884;
        else
            pb.y0(find(pbm.congrps==ig_('sigt114')),1) = 47.245860884;
        end
        if(isKey(ix_,'sigr114'))
            pb.x0(ix_('sigr114'),1) = 57.667755969;
        else
            pb.y0(find(pbm.congrps==ig_('sigr114')),1) = 57.667755969;
        end
        if(isKey(ix_,'x114'))
            pb.x0(ix_('x114'),1) = 1031.2787500;
        else
            pb.y0(find(pbm.congrps==ig_('x114')),1) = 1031.2787500;
        end
        if(isKey(ix_,'y114'))
            pb.x0(ix_('y114'),1) = 55353.503720;
        else
            pb.y0(find(pbm.congrps==ig_('y114')),1) = 55353.503720;
        end
        if(isKey(ix_,'w115'))
            pb.x0(ix_('w115'),1) = 6.6762500000;
        else
            pb.y0(find(pbm.congrps==ig_('w115')),1) = 6.6762500000;
        end
        if(isKey(ix_,'sigt115'))
            pb.x0(ix_('sigt115'),1) = 47.039439651;
        else
            pb.y0(find(pbm.congrps==ig_('sigt115')),1) = 47.039439651;
        end
        if(isKey(ix_,'sigr115'))
            pb.x0(ix_('sigr115'),1) = 57.583729120;
        else
            pb.y0(find(pbm.congrps==ig_('sigr115')),1) = 57.583729120;
        end
        if(isKey(ix_,'x115'))
            pb.x0(ix_('x115'),1) = 1034.6221875;
        else
            pb.y0(find(pbm.congrps==ig_('x115')),1) = 1034.6221875;
        end
        if(isKey(ix_,'y115'))
            pb.x0(ix_('y115'),1) = 55511.122773;
        else
            pb.y0(find(pbm.congrps==ig_('y115')),1) = 55511.122773;
        end
        if(isKey(ix_,'w116'))
            pb.x0(ix_('w116'),1) = 6.6550000000;
        else
            pb.y0(find(pbm.congrps==ig_('w116')),1) = 6.6550000000;
        end
        if(isKey(ix_,'sigt116'))
            pb.x0(ix_('sigt116'),1) = 46.829297460;
        else
            pb.y0(find(pbm.congrps==ig_('sigt116')),1) = 46.829297460;
        end
        if(isKey(ix_,'sigr116'))
            pb.x0(ix_('sigr116'),1) = 57.498572198;
        else
            pb.y0(find(pbm.congrps==ig_('sigr116')),1) = 57.498572198;
        end
        if(isKey(ix_,'x116'))
            pb.x0(ix_('x116'),1) = 1037.9550000;
        else
            pb.y0(find(pbm.congrps==ig_('x116')),1) = 1037.9550000;
        end
        if(isKey(ix_,'y116'))
            pb.x0(ix_('y116'),1) = 55667.546782;
        else
            pb.y0(find(pbm.congrps==ig_('y116')),1) = 55667.546782;
        end
        if(isKey(ix_,'w117'))
            pb.x0(ix_('w117'),1) = 6.6337500000;
        else
            pb.y0(find(pbm.congrps==ig_('w117')),1) = 6.6337500000;
        end
        if(isKey(ix_,'sigt117'))
            pb.x0(ix_('sigt117'),1) = 46.615404203;
        else
            pb.y0(find(pbm.congrps==ig_('sigt117')),1) = 46.615404203;
        end
        if(isKey(ix_,'sigr117'))
            pb.x0(ix_('sigr117'),1) = 57.412274068;
        else
            pb.y0(find(pbm.congrps==ig_('sigr117')),1) = 57.412274068;
        end
        if(isKey(ix_,'x117'))
            pb.x0(ix_('x117'),1) = 1041.2771875;
        else
            pb.y0(find(pbm.congrps==ig_('x117')),1) = 1041.2771875;
        end
        if(isKey(ix_,'y117'))
            pb.x0(ix_('y117'),1) = 55822.767760;
        else
            pb.y0(find(pbm.congrps==ig_('y117')),1) = 55822.767760;
        end
        if(isKey(ix_,'w118'))
            pb.x0(ix_('w118'),1) = 6.6125000000;
        else
            pb.y0(find(pbm.congrps==ig_('w118')),1) = 6.6125000000;
        end
        if(isKey(ix_,'sigt118'))
            pb.x0(ix_('sigt118'),1) = 46.397729510;
        else
            pb.y0(find(pbm.congrps==ig_('sigt118')),1) = 46.397729510;
        end
        if(isKey(ix_,'sigr118'))
            pb.x0(ix_('sigr118'),1) = 57.324823593;
        else
            pb.y0(find(pbm.congrps==ig_('sigr118')),1) = 57.324823593;
        end
        if(isKey(ix_,'x118'))
            pb.x0(ix_('x118'),1) = 1044.5887500;
        else
            pb.y0(find(pbm.congrps==ig_('x118')),1) = 1044.5887500;
        end
        if(isKey(ix_,'y118'))
            pb.x0(ix_('y118'),1) = 55976.777741;
        else
            pb.y0(find(pbm.congrps==ig_('y118')),1) = 55976.777741;
        end
        if(isKey(ix_,'w119'))
            pb.x0(ix_('w119'),1) = 6.5912500000;
        else
            pb.y0(find(pbm.congrps==ig_('w119')),1) = 6.5912500000;
        end
        if(isKey(ix_,'sigt119'))
            pb.x0(ix_('sigt119'),1) = 46.176242754;
        else
            pb.y0(find(pbm.congrps==ig_('sigt119')),1) = 46.176242754;
        end
        if(isKey(ix_,'sigr119'))
            pb.x0(ix_('sigr119'),1) = 57.236209626;
        else
            pb.y0(find(pbm.congrps==ig_('sigr119')),1) = 57.236209626;
        end
        if(isKey(ix_,'x119'))
            pb.x0(ix_('x119'),1) = 1047.8896875;
        else
            pb.y0(find(pbm.congrps==ig_('x119')),1) = 1047.8896875;
        end
        if(isKey(ix_,'y119'))
            pb.x0(ix_('y119'),1) = 56129.568777;
        else
            pb.y0(find(pbm.congrps==ig_('y119')),1) = 56129.568777;
        end
        if(isKey(ix_,'w120'))
            pb.x0(ix_('w120'),1) = 6.5700000000;
        else
            pb.y0(find(pbm.congrps==ig_('w120')),1) = 6.5700000000;
        end
        if(isKey(ix_,'sigt120'))
            pb.x0(ix_('sigt120'),1) = 45.950913048;
        else
            pb.y0(find(pbm.congrps==ig_('sigt120')),1) = 45.950913048;
        end
        if(isKey(ix_,'sigr120'))
            pb.x0(ix_('sigr120'),1) = 57.146421013;
        else
            pb.y0(find(pbm.congrps==ig_('sigr120')),1) = 57.146421013;
        end
        if(isKey(ix_,'x120'))
            pb.x0(ix_('x120'),1) = 1051.1800000;
        else
            pb.y0(find(pbm.congrps==ig_('x120')),1) = 1051.1800000;
        end
        if(isKey(ix_,'y120'))
            pb.x0(ix_('y120'),1) = 56281.132942;
        else
            pb.y0(find(pbm.congrps==ig_('y120')),1) = 56281.132942;
        end
        if(isKey(ix_,'w121'))
            pb.x0(ix_('w121'),1) = 6.5412500000;
        else
            pb.y0(find(pbm.congrps==ig_('w121')),1) = 6.5412500000;
        end
        if(isKey(ix_,'sigt121'))
            pb.x0(ix_('sigt121'),1) = 45.741494802;
        else
            pb.y0(find(pbm.congrps==ig_('sigt121')),1) = 45.741494802;
        end
        if(isKey(ix_,'sigr121'))
            pb.x0(ix_('sigr121'),1) = 57.120910879;
        else
            pb.y0(find(pbm.congrps==ig_('sigr121')),1) = 57.120910879;
        end
        if(isKey(ix_,'x121'))
            pb.x0(ix_('x121'),1) = 1054.4578125;
        else
            pb.y0(find(pbm.congrps==ig_('x121')),1) = 1054.4578125;
        end
        if(isKey(ix_,'y121'))
            pb.x0(ix_('y121'),1) = 56431.408955;
        else
            pb.y0(find(pbm.congrps==ig_('y121')),1) = 56431.408955;
        end
        if(isKey(ix_,'w122'))
            pb.x0(ix_('w122'),1) = 6.5125000000;
        else
            pb.y0(find(pbm.congrps==ig_('w122')),1) = 6.5125000000;
        end
        if(isKey(ix_,'sigt122'))
            pb.x0(ix_('sigt122'),1) = 45.528600901;
        else
            pb.y0(find(pbm.congrps==ig_('sigt122')),1) = 45.528600901;
        end
        if(isKey(ix_,'sigr122'))
            pb.x0(ix_('sigr122'),1) = 57.094665862;
        else
            pb.y0(find(pbm.congrps==ig_('sigr122')),1) = 57.094665862;
        end
        if(isKey(ix_,'x122'))
            pb.x0(ix_('x122'),1) = 1057.7212500;
        else
            pb.y0(find(pbm.congrps==ig_('x122')),1) = 1057.7212500;
        end
        if(isKey(ix_,'y122'))
            pb.x0(ix_('y122'),1) = 56580.336846;
        else
            pb.y0(find(pbm.congrps==ig_('y122')),1) = 56580.336846;
        end
        if(isKey(ix_,'w123'))
            pb.x0(ix_('w123'),1) = 6.4837500000;
        else
            pb.y0(find(pbm.congrps==ig_('w123')),1) = 6.4837500000;
        end
        if(isKey(ix_,'sigt123'))
            pb.x0(ix_('sigt123'),1) = 45.312199922;
        else
            pb.y0(find(pbm.congrps==ig_('sigt123')),1) = 45.312199922;
        end
        if(isKey(ix_,'sigr123'))
            pb.x0(ix_('sigr123'),1) = 57.067684218;
        else
            pb.y0(find(pbm.congrps==ig_('sigr123')),1) = 57.067684218;
        end
        if(isKey(ix_,'x123'))
            pb.x0(ix_('x123'),1) = 1060.9703125;
        else
            pb.y0(find(pbm.congrps==ig_('x123')),1) = 1060.9703125;
        end
        if(isKey(ix_,'y123'))
            pb.x0(ix_('y123'),1) = 56727.911344;
        else
            pb.y0(find(pbm.congrps==ig_('y123')),1) = 56727.911344;
        end
        if(isKey(ix_,'w124'))
            pb.x0(ix_('w124'),1) = 6.4550000000;
        else
            pb.y0(find(pbm.congrps==ig_('w124')),1) = 6.4550000000;
        end
        if(isKey(ix_,'sigt124'))
            pb.x0(ix_('sigt124'),1) = 45.092260303;
        else
            pb.y0(find(pbm.congrps==ig_('sigt124')),1) = 45.092260303;
        end
        if(isKey(ix_,'sigr124'))
            pb.x0(ix_('sigr124'),1) = 57.039964224;
        else
            pb.y0(find(pbm.congrps==ig_('sigr124')),1) = 57.039964224;
        end
        if(isKey(ix_,'x124'))
            pb.x0(ix_('x124'),1) = 1064.2050000;
        else
            pb.y0(find(pbm.congrps==ig_('x124')),1) = 1064.2050000;
        end
        if(isKey(ix_,'y124'))
            pb.x0(ix_('y124'),1) = 56874.127223;
        else
            pb.y0(find(pbm.congrps==ig_('y124')),1) = 56874.127223;
        end
        if(isKey(ix_,'w125'))
            pb.x0(ix_('w125'),1) = 6.4262500000;
        else
            pb.y0(find(pbm.congrps==ig_('w125')),1) = 6.4262500000;
        end
        if(isKey(ix_,'sigt125'))
            pb.x0(ix_('sigt125'),1) = 44.868750344;
        else
            pb.y0(find(pbm.congrps==ig_('sigt125')),1) = 44.868750344;
        end
        if(isKey(ix_,'sigr125'))
            pb.x0(ix_('sigr125'),1) = 57.011504190;
        else
            pb.y0(find(pbm.congrps==ig_('sigr125')),1) = 57.011504190;
        end
        if(isKey(ix_,'x125'))
            pb.x0(ix_('x125'),1) = 1067.4253125;
        else
            pb.y0(find(pbm.congrps==ig_('x125')),1) = 1067.4253125;
        end
        if(isKey(ix_,'y125'))
            pb.x0(ix_('y125'),1) = 57018.979310;
        else
            pb.y0(find(pbm.congrps==ig_('y125')),1) = 57018.979310;
        end
        if(isKey(ix_,'w126'))
            pb.x0(ix_('w126'),1) = 6.3975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w126')),1) = 6.3975000000;
        end
        if(isKey(ix_,'sigt126'))
            pb.x0(ix_('sigt126'),1) = 44.641638205;
        else
            pb.y0(find(pbm.congrps==ig_('sigt126')),1) = 44.641638205;
        end
        if(isKey(ix_,'sigr126'))
            pb.x0(ix_('sigr126'),1) = 56.982302451;
        else
            pb.y0(find(pbm.congrps==ig_('sigr126')),1) = 56.982302451;
        end
        if(isKey(ix_,'x126'))
            pb.x0(ix_('x126'),1) = 1070.6312500;
        else
            pb.y0(find(pbm.congrps==ig_('x126')),1) = 1070.6312500;
        end
        if(isKey(ix_,'y126'))
            pb.x0(ix_('y126'),1) = 57162.462482;
        else
            pb.y0(find(pbm.congrps==ig_('y126')),1) = 57162.462482;
        end
        if(isKey(ix_,'w127'))
            pb.x0(ix_('w127'),1) = 6.3687500000;
        else
            pb.y0(find(pbm.congrps==ig_('w127')),1) = 6.3687500000;
        end
        if(isKey(ix_,'sigt127'))
            pb.x0(ix_('sigt127'),1) = 44.410891909;
        else
            pb.y0(find(pbm.congrps==ig_('sigt127')),1) = 44.410891909;
        end
        if(isKey(ix_,'sigr127'))
            pb.x0(ix_('sigr127'),1) = 56.952357376;
        else
            pb.y0(find(pbm.congrps==ig_('sigr127')),1) = 56.952357376;
        end
        if(isKey(ix_,'x127'))
            pb.x0(ix_('x127'),1) = 1073.8228125;
        else
            pb.y0(find(pbm.congrps==ig_('x127')),1) = 1073.8228125;
        end
        if(isKey(ix_,'y127'))
            pb.x0(ix_('y127'),1) = 57304.571669;
        else
            pb.y0(find(pbm.congrps==ig_('y127')),1) = 57304.571669;
        end
        if(isKey(ix_,'w128'))
            pb.x0(ix_('w128'),1) = 6.3400000000;
        else
            pb.y0(find(pbm.congrps==ig_('w128')),1) = 6.3400000000;
        end
        if(isKey(ix_,'sigt128'))
            pb.x0(ix_('sigt128'),1) = 44.176479340;
        else
            pb.y0(find(pbm.congrps==ig_('sigt128')),1) = 44.176479340;
        end
        if(isKey(ix_,'sigr128'))
            pb.x0(ix_('sigr128'),1) = 56.921667364;
        else
            pb.y0(find(pbm.congrps==ig_('sigr128')),1) = 56.921667364;
        end
        if(isKey(ix_,'x128'))
            pb.x0(ix_('x128'),1) = 1077.0000000;
        else
            pb.y0(find(pbm.congrps==ig_('x128')),1) = 1077.0000000;
        end
        if(isKey(ix_,'y128'))
            pb.x0(ix_('y128'),1) = 57445.301855;
        else
            pb.y0(find(pbm.congrps==ig_('y128')),1) = 57445.301855;
        end
        if(isKey(ix_,'w129'))
            pb.x0(ix_('w129'),1) = 6.3162500000;
        else
            pb.y0(find(pbm.congrps==ig_('w129')),1) = 6.3162500000;
        end
        if(isKey(ix_,'sigt129'))
            pb.x0(ix_('sigt129'),1) = 43.924748683;
        else
            pb.y0(find(pbm.congrps==ig_('sigt129')),1) = 43.924748683;
        end
        if(isKey(ix_,'sigr129'))
            pb.x0(ix_('sigr129'),1) = 56.845155310;
        else
            pb.y0(find(pbm.congrps==ig_('sigr129')),1) = 56.845155310;
        end
        if(isKey(ix_,'x129'))
            pb.x0(ix_('x129'),1) = 1080.1640625;
        else
            pb.y0(find(pbm.congrps==ig_('x129')),1) = 1080.1640625;
        end
        if(isKey(ix_,'y129'))
            pb.x0(ix_('y129'),1) = 57584.681499;
        else
            pb.y0(find(pbm.congrps==ig_('y129')),1) = 57584.681499;
        end
        if(isKey(ix_,'w130'))
            pb.x0(ix_('w130'),1) = 6.2925000000;
        else
            pb.y0(find(pbm.congrps==ig_('w130')),1) = 6.2925000000;
        end
        if(isKey(ix_,'sigt130'))
            pb.x0(ix_('sigt130'),1) = 43.668982053;
        else
            pb.y0(find(pbm.congrps==ig_('sigt130')),1) = 43.668982053;
        end
        if(isKey(ix_,'sigr130'))
            pb.x0(ix_('sigr130'),1) = 56.767522016;
        else
            pb.y0(find(pbm.congrps==ig_('sigr130')),1) = 56.767522016;
        end
        if(isKey(ix_,'x130'))
            pb.x0(ix_('x130'),1) = 1083.3162500;
        else
            pb.y0(find(pbm.congrps==ig_('x130')),1) = 1083.3162500;
        end
        if(isKey(ix_,'y130'))
            pb.x0(ix_('y130'),1) = 57722.738189;
        else
            pb.y0(find(pbm.congrps==ig_('y130')),1) = 57722.738189;
        end
        if(isKey(ix_,'w131'))
            pb.x0(ix_('w131'),1) = 6.2687500000;
        else
            pb.y0(find(pbm.congrps==ig_('w131')),1) = 6.2687500000;
        end
        if(isKey(ix_,'sigt131'))
            pb.x0(ix_('sigt131'),1) = 43.409146055;
        else
            pb.y0(find(pbm.congrps==ig_('sigt131')),1) = 43.409146055;
        end
        if(isKey(ix_,'sigr131'))
            pb.x0(ix_('sigr131'),1) = 56.688758490;
        else
            pb.y0(find(pbm.congrps==ig_('sigr131')),1) = 56.688758490;
        end
        if(isKey(ix_,'x131'))
            pb.x0(ix_('x131'),1) = 1086.4565625;
        else
            pb.y0(find(pbm.congrps==ig_('x131')),1) = 1086.4565625;
        end
        if(isKey(ix_,'y131'))
            pb.x0(ix_('y131'),1) = 57859.465228;
        else
            pb.y0(find(pbm.congrps==ig_('y131')),1) = 57859.465228;
        end
        if(isKey(ix_,'w132'))
            pb.x0(ix_('w132'),1) = 6.2450000000;
        else
            pb.y0(find(pbm.congrps==ig_('w132')),1) = 6.2450000000;
        end
        if(isKey(ix_,'sigt132'))
            pb.x0(ix_('sigt132'),1) = 43.145207079;
        else
            pb.y0(find(pbm.congrps==ig_('sigt132')),1) = 43.145207079;
        end
        if(isKey(ix_,'sigr132'))
            pb.x0(ix_('sigr132'),1) = 56.608855703;
        else
            pb.y0(find(pbm.congrps==ig_('sigr132')),1) = 56.608855703;
        end
        if(isKey(ix_,'x132'))
            pb.x0(ix_('x132'),1) = 1089.5850000;
        else
            pb.y0(find(pbm.congrps==ig_('x132')),1) = 1089.5850000;
        end
        if(isKey(ix_,'y132'))
            pb.x0(ix_('y132'),1) = 57994.855954;
        else
            pb.y0(find(pbm.congrps==ig_('y132')),1) = 57994.855954;
        end
        if(isKey(ix_,'w133'))
            pb.x0(ix_('w133'),1) = 6.2212500000;
        else
            pb.y0(find(pbm.congrps==ig_('w133')),1) = 6.2212500000;
        end
        if(isKey(ix_,'sigt133'))
            pb.x0(ix_('sigt133'),1) = 42.877131300;
        else
            pb.y0(find(pbm.congrps==ig_('sigt133')),1) = 42.877131300;
        end
        if(isKey(ix_,'sigr133'))
            pb.x0(ix_('sigr133'),1) = 56.527804595;
        else
            pb.y0(find(pbm.congrps==ig_('sigr133')),1) = 56.527804595;
        end
        if(isKey(ix_,'x133'))
            pb.x0(ix_('x133'),1) = 1092.7015625;
        else
            pb.y0(find(pbm.congrps==ig_('x133')),1) = 1092.7015625;
        end
        if(isKey(ix_,'y133'))
            pb.x0(ix_('y133'),1) = 58128.903746;
        else
            pb.y0(find(pbm.congrps==ig_('y133')),1) = 58128.903746;
        end
        if(isKey(ix_,'w134'))
            pb.x0(ix_('w134'),1) = 6.1975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w134')),1) = 6.1975000000;
        end
        if(isKey(ix_,'sigt134'))
            pb.x0(ix_('sigt134'),1) = 42.604884678;
        else
            pb.y0(find(pbm.congrps==ig_('sigt134')),1) = 42.604884678;
        end
        if(isKey(ix_,'sigr134'))
            pb.x0(ix_('sigr134'),1) = 56.445596072;
        else
            pb.y0(find(pbm.congrps==ig_('sigr134')),1) = 56.445596072;
        end
        if(isKey(ix_,'x134'))
            pb.x0(ix_('x134'),1) = 1095.8062500;
        else
            pb.y0(find(pbm.congrps==ig_('x134')),1) = 1095.8062500;
        end
        if(isKey(ix_,'y134'))
            pb.x0(ix_('y134'),1) = 58261.602028;
        else
            pb.y0(find(pbm.congrps==ig_('y134')),1) = 58261.602028;
        end
        if(isKey(ix_,'w135'))
            pb.x0(ix_('w135'),1) = 6.1737500000;
        else
            pb.y0(find(pbm.congrps==ig_('w135')),1) = 6.1737500000;
        end
        if(isKey(ix_,'sigt135'))
            pb.x0(ix_('sigt135'),1) = 42.328432957;
        else
            pb.y0(find(pbm.congrps==ig_('sigt135')),1) = 42.328432957;
        end
        if(isKey(ix_,'sigr135'))
            pb.x0(ix_('sigr135'),1) = 56.362220999;
        else
            pb.y0(find(pbm.congrps==ig_('sigr135')),1) = 56.362220999;
        end
        if(isKey(ix_,'x135'))
            pb.x0(ix_('x135'),1) = 1098.8990625;
        else
            pb.y0(find(pbm.congrps==ig_('x135')),1) = 1098.8990625;
        end
        if(isKey(ix_,'y135'))
            pb.x0(ix_('y135'),1) = 58392.944262;
        else
            pb.y0(find(pbm.congrps==ig_('y135')),1) = 58392.944262;
        end
        if(isKey(ix_,'w136'))
            pb.x0(ix_('w136'),1) = 6.1500000000;
        else
            pb.y0(find(pbm.congrps==ig_('w136')),1) = 6.1500000000;
        end
        if(isKey(ix_,'sigt136'))
            pb.x0(ix_('sigt136'),1) = 42.047741672;
        else
            pb.y0(find(pbm.congrps==ig_('sigt136')),1) = 42.047741672;
        end
        if(isKey(ix_,'sigr136'))
            pb.x0(ix_('sigr136'),1) = 56.277670207;
        else
            pb.y0(find(pbm.congrps==ig_('sigr136')),1) = 56.277670207;
        end
        if(isKey(ix_,'x136'))
            pb.x0(ix_('x136'),1) = 1101.9800000;
        else
            pb.y0(find(pbm.congrps==ig_('x136')),1) = 1101.9800000;
        end
        if(isKey(ix_,'y136'))
            pb.x0(ix_('y136'),1) = 58522.923955;
        else
            pb.y0(find(pbm.congrps==ig_('y136')),1) = 58522.923955;
        end
        if(isKey(ix_,'w137'))
            pb.x0(ix_('w137'),1) = 6.1150000000;
        else
            pb.y0(find(pbm.congrps==ig_('w137')),1) = 6.1150000000;
        end
        if(isKey(ix_,'sigt137'))
            pb.x0(ix_('sigt137'),1) = 41.794038608;
        else
            pb.y0(find(pbm.congrps==ig_('sigt137')),1) = 41.794038608;
        end
        if(isKey(ix_,'sigr137'))
            pb.x0(ix_('sigr137'),1) = 56.295428092;
        else
            pb.y0(find(pbm.congrps==ig_('sigr137')),1) = 56.295428092;
        end
        if(isKey(ix_,'x137'))
            pb.x0(ix_('x137'),1) = 1105.0462500;
        else
            pb.y0(find(pbm.congrps==ig_('x137')),1) = 1105.0462500;
        end
        if(isKey(ix_,'y137'))
            pb.x0(ix_('y137'),1) = 58651.464995;
        else
            pb.y0(find(pbm.congrps==ig_('y137')),1) = 58651.464995;
        end
        if(isKey(ix_,'w138'))
            pb.x0(ix_('w138'),1) = 6.0800000000;
        else
            pb.y0(find(pbm.congrps==ig_('w138')),1) = 6.0800000000;
        end
        if(isKey(ix_,'sigt138'))
            pb.x0(ix_('sigt138'),1) = 41.536786267;
        else
            pb.y0(find(pbm.congrps==ig_('sigt138')),1) = 41.536786267;
        end
        if(isKey(ix_,'sigr138'))
            pb.x0(ix_('sigr138'),1) = 56.313098515;
        else
            pb.y0(find(pbm.congrps==ig_('sigr138')),1) = 56.313098515;
        end
        if(isKey(ix_,'x138'))
            pb.x0(ix_('x138'),1) = 1108.0950000;
        else
            pb.y0(find(pbm.congrps==ig_('x138')),1) = 1108.0950000;
        end
        if(isKey(ix_,'y138'))
            pb.x0(ix_('y138'),1) = 58778.493546;
        else
            pb.y0(find(pbm.congrps==ig_('y138')),1) = 58778.493546;
        end
        if(isKey(ix_,'w139'))
            pb.x0(ix_('w139'),1) = 6.0450000000;
        else
            pb.y0(find(pbm.congrps==ig_('w139')),1) = 6.0450000000;
        end
        if(isKey(ix_,'sigt139'))
            pb.x0(ix_('sigt139'),1) = 41.275954987;
        else
            pb.y0(find(pbm.congrps==ig_('sigt139')),1) = 41.275954987;
        end
        if(isKey(ix_,'sigr139'))
            pb.x0(ix_('sigr139'),1) = 56.330696252;
        else
            pb.y0(find(pbm.congrps==ig_('sigr139')),1) = 56.330696252;
        end
        if(isKey(ix_,'x139'))
            pb.x0(ix_('x139'),1) = 1111.1262500;
        else
            pb.y0(find(pbm.congrps==ig_('x139')),1) = 1111.1262500;
        end
        if(isKey(ix_,'y139'))
            pb.x0(ix_('y139'),1) = 58904.007748;
        else
            pb.y0(find(pbm.congrps==ig_('y139')),1) = 58904.007748;
        end
        if(isKey(ix_,'w140'))
            pb.x0(ix_('w140'),1) = 6.0100000000;
        else
            pb.y0(find(pbm.congrps==ig_('w140')),1) = 6.0100000000;
        end
        if(isKey(ix_,'sigt140'))
            pb.x0(ix_('sigt140'),1) = 41.011515151;
        else
            pb.y0(find(pbm.congrps==ig_('sigt140')),1) = 41.011515151;
        end
        if(isKey(ix_,'sigr140'))
            pb.x0(ix_('sigr140'),1) = 56.348236444;
        else
            pb.y0(find(pbm.congrps==ig_('sigr140')),1) = 56.348236444;
        end
        if(isKey(ix_,'x140'))
            pb.x0(ix_('x140'),1) = 1114.1400000;
        else
            pb.y0(find(pbm.congrps==ig_('x140')),1) = 1114.1400000;
        end
        if(isKey(ix_,'y140'))
            pb.x0(ix_('y140'),1) = 59028.005837;
        else
            pb.y0(find(pbm.congrps==ig_('y140')),1) = 59028.005837;
        end
        if(isKey(ix_,'w141'))
            pb.x0(ix_('w141'),1) = 5.9750000000;
        else
            pb.y0(find(pbm.congrps==ig_('w141')),1) = 5.9750000000;
        end
        if(isKey(ix_,'sigt141'))
            pb.x0(ix_('sigt141'),1) = 40.743437190;
        else
            pb.y0(find(pbm.congrps==ig_('sigt141')),1) = 40.743437190;
        end
        if(isKey(ix_,'sigr141'))
            pb.x0(ix_('sigr141'),1) = 56.365734613;
        else
            pb.y0(find(pbm.congrps==ig_('sigr141')),1) = 56.365734613;
        end
        if(isKey(ix_,'x141'))
            pb.x0(ix_('x141'),1) = 1117.1362500;
        else
            pb.y0(find(pbm.congrps==ig_('x141')),1) = 1117.1362500;
        end
        if(isKey(ix_,'y141'))
            pb.x0(ix_('y141'),1) = 59150.486148;
        else
            pb.y0(find(pbm.congrps==ig_('y141')),1) = 59150.486148;
        end
        if(isKey(ix_,'w142'))
            pb.x0(ix_('w142'),1) = 5.9400000000;
        else
            pb.y0(find(pbm.congrps==ig_('w142')),1) = 5.9400000000;
        end
        if(isKey(ix_,'sigt142'))
            pb.x0(ix_('sigt142'),1) = 40.471691591;
        else
            pb.y0(find(pbm.congrps==ig_('sigt142')),1) = 40.471691591;
        end
        if(isKey(ix_,'sigr142'))
            pb.x0(ix_('sigr142'),1) = 56.383206675;
        else
            pb.y0(find(pbm.congrps==ig_('sigr142')),1) = 56.383206675;
        end
        if(isKey(ix_,'x142'))
            pb.x0(ix_('x142'),1) = 1120.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x142')),1) = 1120.1150000;
        end
        if(isKey(ix_,'y142'))
            pb.x0(ix_('y142'),1) = 59271.447119;
        else
            pb.y0(find(pbm.congrps==ig_('y142')),1) = 59271.447119;
        end
        if(isKey(ix_,'w143'))
            pb.x0(ix_('w143'),1) = 5.9050000000;
        else
            pb.y0(find(pbm.congrps==ig_('w143')),1) = 5.9050000000;
        end
        if(isKey(ix_,'sigt143'))
            pb.x0(ix_('sigt143'),1) = 40.196248899;
        else
            pb.y0(find(pbm.congrps==ig_('sigt143')),1) = 40.196248899;
        end
        if(isKey(ix_,'sigr143'))
            pb.x0(ix_('sigr143'),1) = 56.400668955;
        else
            pb.y0(find(pbm.congrps==ig_('sigr143')),1) = 56.400668955;
        end
        if(isKey(ix_,'x143'))
            pb.x0(ix_('x143'),1) = 1123.0762500;
        else
            pb.y0(find(pbm.congrps==ig_('x143')),1) = 1123.0762500;
        end
        if(isKey(ix_,'y143'))
            pb.x0(ix_('y143'),1) = 59390.887293;
        else
            pb.y0(find(pbm.congrps==ig_('y143')),1) = 59390.887293;
        end
        if(isKey(ix_,'w144'))
            pb.x0(ix_('w144'),1) = 5.8700000000;
        else
            pb.y0(find(pbm.congrps==ig_('w144')),1) = 5.8700000000;
        end
        if(isKey(ix_,'sigt144'))
            pb.x0(ix_('sigt144'),1) = 39.917079719;
        else
            pb.y0(find(pbm.congrps==ig_('sigt144')),1) = 39.917079719;
        end
        if(isKey(ix_,'sigr144'))
            pb.x0(ix_('sigr144'),1) = 56.418138202;
        else
            pb.y0(find(pbm.congrps==ig_('sigr144')),1) = 56.418138202;
        end
        if(isKey(ix_,'x144'))
            pb.x0(ix_('x144'),1) = 1126.0200000;
        else
            pb.y0(find(pbm.congrps==ig_('x144')),1) = 1126.0200000;
        end
        if(isKey(ix_,'y144'))
            pb.x0(ix_('y144'),1) = 59508.805320;
        else
            pb.y0(find(pbm.congrps==ig_('y144')),1) = 59508.805320;
        end
        if(isKey(ix_,'w145'))
            pb.x0(ix_('w145'),1) = 5.8412500000;
        else
            pb.y0(find(pbm.congrps==ig_('w145')),1) = 5.8412500000;
        end
        if(isKey(ix_,'sigt145'))
            pb.x0(ix_('sigt145'),1) = 39.615894752;
        else
            pb.y0(find(pbm.congrps==ig_('sigt145')),1) = 39.615894752;
        end
        if(isKey(ix_,'sigr145'))
            pb.x0(ix_('sigr145'),1) = 56.375167857;
        else
            pb.y0(find(pbm.congrps==ig_('sigr145')),1) = 56.375167857;
        end
        if(isKey(ix_,'x145'))
            pb.x0(ix_('x145'),1) = 1128.9478125;
        else
            pb.y0(find(pbm.congrps==ig_('x145')),1) = 1128.9478125;
        end
        if(isKey(ix_,'y145'))
            pb.x0(ix_('y145'),1) = 59625.235221;
        else
            pb.y0(find(pbm.congrps==ig_('y145')),1) = 59625.235221;
        end
        if(isKey(ix_,'w146'))
            pb.x0(ix_('w146'),1) = 5.8125000000;
        else
            pb.y0(find(pbm.congrps==ig_('w146')),1) = 5.8125000000;
        end
        if(isKey(ix_,'sigt146'))
            pb.x0(ix_('sigt146'),1) = 39.310443929;
        else
            pb.y0(find(pbm.congrps==ig_('sigt146')),1) = 39.310443929;
        end
        if(isKey(ix_,'sigr146'))
            pb.x0(ix_('sigr146'),1) = 56.331441829;
        else
            pb.y0(find(pbm.congrps==ig_('sigr146')),1) = 56.331441829;
        end
        if(isKey(ix_,'x146'))
            pb.x0(ix_('x146'),1) = 1131.8612500;
        else
            pb.y0(find(pbm.congrps==ig_('x146')),1) = 1131.8612500;
        end
        if(isKey(ix_,'y146'))
            pb.x0(ix_('y146'),1) = 59740.209796;
        else
            pb.y0(find(pbm.congrps==ig_('y146')),1) = 59740.209796;
        end
        if(isKey(ix_,'w147'))
            pb.x0(ix_('w147'),1) = 5.7837500000;
        else
            pb.y0(find(pbm.congrps==ig_('w147')),1) = 5.7837500000;
        end
        if(isKey(ix_,'sigt147'))
            pb.x0(ix_('sigt147'),1) = 39.000692735;
        else
            pb.y0(find(pbm.congrps==ig_('sigt147')),1) = 39.000692735;
        end
        if(isKey(ix_,'sigr147'))
            pb.x0(ix_('sigr147'),1) = 56.286959540;
        else
            pb.y0(find(pbm.congrps==ig_('sigr147')),1) = 56.286959540;
        end
        if(isKey(ix_,'x147'))
            pb.x0(ix_('x147'),1) = 1134.7603125;
        else
            pb.y0(find(pbm.congrps==ig_('x147')),1) = 1134.7603125;
        end
        if(isKey(ix_,'y147'))
            pb.x0(ix_('y147'),1) = 59853.725349;
        else
            pb.y0(find(pbm.congrps==ig_('y147')),1) = 59853.725349;
        end
        if(isKey(ix_,'w148'))
            pb.x0(ix_('w148'),1) = 5.7550000000;
        else
            pb.y0(find(pbm.congrps==ig_('w148')),1) = 5.7550000000;
        end
        if(isKey(ix_,'sigt148'))
            pb.x0(ix_('sigt148'),1) = 38.686606531;
        else
            pb.y0(find(pbm.congrps==ig_('sigt148')),1) = 38.686606531;
        end
        if(isKey(ix_,'sigr148'))
            pb.x0(ix_('sigr148'),1) = 56.241720481;
        else
            pb.y0(find(pbm.congrps==ig_('sigr148')),1) = 56.241720481;
        end
        if(isKey(ix_,'x148'))
            pb.x0(ix_('x148'),1) = 1137.6450000;
        else
            pb.y0(find(pbm.congrps==ig_('x148')),1) = 1137.6450000;
        end
        if(isKey(ix_,'y148'))
            pb.x0(ix_('y148'),1) = 59965.778269;
        else
            pb.y0(find(pbm.congrps==ig_('y148')),1) = 59965.778269;
        end
        if(isKey(ix_,'w149'))
            pb.x0(ix_('w149'),1) = 5.7262500000;
        else
            pb.y0(find(pbm.congrps==ig_('w149')),1) = 5.7262500000;
        end
        if(isKey(ix_,'sigt149'))
            pb.x0(ix_('sigt149'),1) = 38.368150555;
        else
            pb.y0(find(pbm.congrps==ig_('sigt149')),1) = 38.368150555;
        end
        if(isKey(ix_,'sigr149'))
            pb.x0(ix_('sigr149'),1) = 56.195724220;
        else
            pb.y0(find(pbm.congrps==ig_('sigr149')),1) = 56.195724220;
        end
        if(isKey(ix_,'x149'))
            pb.x0(ix_('x149'),1) = 1140.5153125;
        else
            pb.y0(find(pbm.congrps==ig_('x149')),1) = 1140.5153125;
        end
        if(isKey(ix_,'y149'))
            pb.x0(ix_('y149'),1) = 60076.365029;
        else
            pb.y0(find(pbm.congrps==ig_('y149')),1) = 60076.365029;
        end
        if(isKey(ix_,'w150'))
            pb.x0(ix_('w150'),1) = 5.6975000000;
        else
            pb.y0(find(pbm.congrps==ig_('w150')),1) = 5.6975000000;
        end
        if(isKey(ix_,'sigt150'))
            pb.x0(ix_('sigt150'),1) = 38.045289920;
        else
            pb.y0(find(pbm.congrps==ig_('sigt150')),1) = 38.045289920;
        end
        if(isKey(ix_,'sigr150'))
            pb.x0(ix_('sigr150'),1) = 56.148970399;
        else
            pb.y0(find(pbm.congrps==ig_('sigr150')),1) = 56.148970399;
        end
        if(isKey(ix_,'x150'))
            pb.x0(ix_('x150'),1) = 1143.3712500;
        else
            pb.y0(find(pbm.congrps==ig_('x150')),1) = 1143.3712500;
        end
        if(isKey(ix_,'y150'))
            pb.x0(ix_('y150'),1) = 60185.482195;
        else
            pb.y0(find(pbm.congrps==ig_('y150')),1) = 60185.482195;
        end
        if(isKey(ix_,'w151'))
            pb.x0(ix_('w151'),1) = 5.6687500000;
        else
            pb.y0(find(pbm.congrps==ig_('w151')),1) = 5.6687500000;
        end
        if(isKey(ix_,'sigt151'))
            pb.x0(ix_('sigt151'),1) = 37.717989622;
        else
            pb.y0(find(pbm.congrps==ig_('sigt151')),1) = 37.717989622;
        end
        if(isKey(ix_,'sigr151'))
            pb.x0(ix_('sigr151'),1) = 56.101458742;
        else
            pb.y0(find(pbm.congrps==ig_('sigr151')),1) = 56.101458742;
        end
        if(isKey(ix_,'x151'))
            pb.x0(ix_('x151'),1) = 1146.2128125;
        else
            pb.y0(find(pbm.congrps==ig_('x151')),1) = 1146.2128125;
        end
        if(isKey(ix_,'y151'))
            pb.x0(ix_('y151'),1) = 60293.126418;
        else
            pb.y0(find(pbm.congrps==ig_('y151')),1) = 60293.126418;
        end
        if(isKey(ix_,'w152'))
            pb.x0(ix_('w152'),1) = 5.6400000000;
        else
            pb.y0(find(pbm.congrps==ig_('w152')),1) = 5.6400000000;
        end
        if(isKey(ix_,'sigt152'))
            pb.x0(ix_('sigt152'),1) = 37.386214534;
        else
            pb.y0(find(pbm.congrps==ig_('sigt152')),1) = 37.386214534;
        end
        if(isKey(ix_,'sigr152'))
            pb.x0(ix_('sigr152'),1) = 56.053189053;
        else
            pb.y0(find(pbm.congrps==ig_('sigr152')),1) = 56.053189053;
        end
        if(isKey(ix_,'x152'))
            pb.x0(ix_('x152'),1) = 1149.0400000;
        else
            pb.y0(find(pbm.congrps==ig_('x152')),1) = 1149.0400000;
        end
        if(isKey(ix_,'y152'))
            pb.x0(ix_('y152'),1) = 60399.294444;
        else
            pb.y0(find(pbm.congrps==ig_('y152')),1) = 60399.294444;
        end
        if(isKey(ix_,'w153'))
            pb.x0(ix_('w153'),1) = 5.6112500000;
        else
            pb.y0(find(pbm.congrps==ig_('w153')),1) = 5.6112500000;
        end
        if(isKey(ix_,'sigt153'))
            pb.x0(ix_('sigt153'),1) = 37.049929410;
        else
            pb.y0(find(pbm.congrps==ig_('sigt153')),1) = 37.049929410;
        end
        if(isKey(ix_,'sigr153'))
            pb.x0(ix_('sigr153'),1) = 56.004161223;
        else
            pb.y0(find(pbm.congrps==ig_('sigr153')),1) = 56.004161223;
        end
        if(isKey(ix_,'x153'))
            pb.x0(ix_('x153'),1) = 1151.8528125;
        else
            pb.y0(find(pbm.congrps==ig_('x153')),1) = 1151.8528125;
        end
        if(isKey(ix_,'y153'))
            pb.x0(ix_('y153'),1) = 60503.983110;
        else
            pb.y0(find(pbm.congrps==ig_('y153')),1) = 60503.983110;
        end
        if(isKey(ix_,'w154'))
            pb.x0(ix_('w154'),1) = 5.5825000000;
        else
            pb.y0(find(pbm.congrps==ig_('w154')),1) = 5.5825000000;
        end
        if(isKey(ix_,'sigt154'))
            pb.x0(ix_('sigt154'),1) = 36.709098885;
        else
            pb.y0(find(pbm.congrps==ig_('sigt154')),1) = 36.709098885;
        end
        if(isKey(ix_,'sigr154'))
            pb.x0(ix_('sigr154'),1) = 55.954375227;
        else
            pb.y0(find(pbm.congrps==ig_('sigr154')),1) = 55.954375227;
        end
        if(isKey(ix_,'x154'))
            pb.x0(ix_('x154'),1) = 1154.6512500;
        else
            pb.y0(find(pbm.congrps==ig_('x154')),1) = 1154.6512500;
        end
        if(isKey(ix_,'y154'))
            pb.x0(ix_('y154'),1) = 60607.189351;
        else
            pb.y0(find(pbm.congrps==ig_('y154')),1) = 60607.189351;
        end
        if(isKey(ix_,'w155'))
            pb.x0(ix_('w155'),1) = 5.5537500000;
        else
            pb.y0(find(pbm.congrps==ig_('w155')),1) = 5.5537500000;
        end
        if(isKey(ix_,'sigt155'))
            pb.x0(ix_('sigt155'),1) = 36.363687480;
        else
            pb.y0(find(pbm.congrps==ig_('sigt155')),1) = 36.363687480;
        end
        if(isKey(ix_,'sigr155'))
            pb.x0(ix_('sigr155'),1) = 55.903831135;
        else
            pb.y0(find(pbm.congrps==ig_('sigr155')),1) = 55.903831135;
        end
        if(isKey(ix_,'x155'))
            pb.x0(ix_('x155'),1) = 1157.4353125;
        else
            pb.y0(find(pbm.congrps==ig_('x155')),1) = 1157.4353125;
        end
        if(isKey(ix_,'y155'))
            pb.x0(ix_('y155'),1) = 60708.910194;
        else
            pb.y0(find(pbm.congrps==ig_('y155')),1) = 60708.910194;
        end
        if(isKey(ix_,'w156'))
            pb.x0(ix_('w156'),1) = 5.5250000000;
        else
            pb.y0(find(pbm.congrps==ig_('w156')),1) = 5.5250000000;
        end
        if(isKey(ix_,'sigt156'))
            pb.x0(ix_('sigt156'),1) = 36.013659597;
        else
            pb.y0(find(pbm.congrps==ig_('sigt156')),1) = 36.013659597;
        end
        if(isKey(ix_,'sigr156'))
            pb.x0(ix_('sigr156'),1) = 55.852529108;
        else
            pb.y0(find(pbm.congrps==ig_('sigr156')),1) = 55.852529108;
        end
        if(isKey(ix_,'x156'))
            pb.x0(ix_('x156'),1) = 1160.2050000;
        else
            pb.y0(find(pbm.congrps==ig_('x156')),1) = 1160.2050000;
        end
        if(isKey(ix_,'y156'))
            pb.x0(ix_('y156'),1) = 60809.142769;
        else
            pb.y0(find(pbm.congrps==ig_('y156')),1) = 60809.142769;
        end
        if(isKey(ix_,'w157'))
            pb.x0(ix_('w157'),1) = 5.4962500000;
        else
            pb.y0(find(pbm.congrps==ig_('w157')),1) = 5.4962500000;
        end
        if(isKey(ix_,'sigt157'))
            pb.x0(ix_('sigt157'),1) = 35.658979526;
        else
            pb.y0(find(pbm.congrps==ig_('sigt157')),1) = 35.658979526;
        end
        if(isKey(ix_,'sigr157'))
            pb.x0(ix_('sigr157'),1) = 55.800469405;
        else
            pb.y0(find(pbm.congrps==ig_('sigr157')),1) = 55.800469405;
        end
        if(isKey(ix_,'x157'))
            pb.x0(ix_('x157'),1) = 1162.9603125;
        else
            pb.y0(find(pbm.congrps==ig_('x157')),1) = 1162.9603125;
        end
        if(isKey(ix_,'y157'))
            pb.x0(ix_('y157'),1) = 60907.884303;
        else
            pb.y0(find(pbm.congrps==ig_('y157')),1) = 60907.884303;
        end
        if(isKey(ix_,'w158'))
            pb.x0(ix_('w158'),1) = 5.4675000000;
        else
            pb.y0(find(pbm.congrps==ig_('w158')),1) = 5.4675000000;
        end
        if(isKey(ix_,'sigt158'))
            pb.x0(ix_('sigt158'),1) = 35.299611440;
        else
            pb.y0(find(pbm.congrps==ig_('sigt158')),1) = 35.299611440;
        end
        if(isKey(ix_,'sigr158'))
            pb.x0(ix_('sigr158'),1) = 55.747652383;
        else
            pb.y0(find(pbm.congrps==ig_('sigr158')),1) = 55.747652383;
        end
        if(isKey(ix_,'x158'))
            pb.x0(ix_('x158'),1) = 1165.7012500;
        else
            pb.y0(find(pbm.congrps==ig_('x158')),1) = 1165.7012500;
        end
        if(isKey(ix_,'y158'))
            pb.x0(ix_('y158'),1) = 61005.132126;
        else
            pb.y0(find(pbm.congrps==ig_('y158')),1) = 61005.132126;
        end
        if(isKey(ix_,'w159'))
            pb.x0(ix_('w159'),1) = 5.4387500000;
        else
            pb.y0(find(pbm.congrps==ig_('w159')),1) = 5.4387500000;
        end
        if(isKey(ix_,'sigt159'))
            pb.x0(ix_('sigt159'),1) = 34.935519403;
        else
            pb.y0(find(pbm.congrps==ig_('sigt159')),1) = 34.935519403;
        end
        if(isKey(ix_,'sigr159'))
            pb.x0(ix_('sigr159'),1) = 55.694078505;
        else
            pb.y0(find(pbm.congrps==ig_('sigr159')),1) = 55.694078505;
        end
        if(isKey(ix_,'x159'))
            pb.x0(ix_('x159'),1) = 1168.4278125;
        else
            pb.y0(find(pbm.congrps==ig_('x159')),1) = 1168.4278125;
        end
        if(isKey(ix_,'y159'))
            pb.x0(ix_('y159'),1) = 61100.883671;
        else
            pb.y0(find(pbm.congrps==ig_('y159')),1) = 61100.883671;
        end
        if(isKey(ix_,'w160'))
            pb.x0(ix_('w160'),1) = 5.4100000000;
        else
            pb.y0(find(pbm.congrps==ig_('w160')),1) = 5.4100000000;
        end
        if(isKey(ix_,'sigt160'))
            pb.x0(ix_('sigt160'),1) = 34.566667368;
        else
            pb.y0(find(pbm.congrps==ig_('sigt160')),1) = 34.566667368;
        end
        if(isKey(ix_,'sigr160'))
            pb.x0(ix_('sigr160'),1) = 55.639748340;
        else
            pb.y0(find(pbm.congrps==ig_('sigr160')),1) = 55.639748340;
        end
        if(isKey(ix_,'x160'))
            pb.x0(ix_('x160'),1) = 1171.1400000;
        else
            pb.y0(find(pbm.congrps==ig_('x160')),1) = 1171.1400000;
        end
        if(isKey(ix_,'y160'))
            pb.x0(ix_('y160'),1) = 61195.136478;
        else
            pb.y0(find(pbm.congrps==ig_('y160')),1) = 61195.136478;
        end
        if(isKey(ix_,'w161'))
            pb.x0(ix_('w161'),1) = 5.6025000000;
        else
            pb.y0(find(pbm.congrps==ig_('w161')),1) = 5.6025000000;
        end
        if(isKey(ix_,'sigt161'))
            pb.x0(ix_('sigt161'),1) = 33.529237255;
        else
            pb.y0(find(pbm.congrps==ig_('sigt161')),1) = 33.529237255;
        end
        if(isKey(ix_,'sigr161'))
            pb.x0(ix_('sigr161'),1) = 53.385743943;
        else
            pb.y0(find(pbm.congrps==ig_('sigr161')),1) = 53.385743943;
        end
        if(isKey(ix_,'x161'))
            pb.x0(ix_('x161'),1) = 1173.8931250;
        else
            pb.y0(find(pbm.congrps==ig_('x161')),1) = 1173.8931250;
        end
        if(isKey(ix_,'y161'))
            pb.x0(ix_('y161'),1) = 61288.849783;
        else
            pb.y0(find(pbm.congrps==ig_('y161')),1) = 61288.849783;
        end
        if(isKey(ix_,'w162'))
            pb.x0(ix_('w162'),1) = 5.7950000000;
        else
            pb.y0(find(pbm.congrps==ig_('w162')),1) = 5.7950000000;
        end
        if(isKey(ix_,'sigt162'))
            pb.x0(ix_('sigt162'),1) = 32.521949841;
        else
            pb.y0(find(pbm.congrps==ig_('sigt162')),1) = 32.521949841;
        end
        if(isKey(ix_,'sigr162'))
            pb.x0(ix_('sigr162'),1) = 51.273873576;
        else
            pb.y0(find(pbm.congrps==ig_('sigr162')),1) = 51.273873576;
        end
        if(isKey(ix_,'x162'))
            pb.x0(ix_('x162'),1) = 1176.7425000;
        else
            pb.y0(find(pbm.congrps==ig_('x162')),1) = 1176.7425000;
        end
        if(isKey(ix_,'y162'))
            pb.x0(ix_('y162'),1) = 61382.927846;
        else
            pb.y0(find(pbm.congrps==ig_('y162')),1) = 61382.927846;
        end
        if(isKey(ix_,'w163'))
            pb.x0(ix_('w163'),1) = 5.9875000000;
        else
            pb.y0(find(pbm.congrps==ig_('w163')),1) = 5.9875000000;
        end
        if(isKey(ix_,'sigt163'))
            pb.x0(ix_('sigt163'),1) = 31.541203081;
        else
            pb.y0(find(pbm.congrps==ig_('sigt163')),1) = 31.541203081;
        end
        if(isKey(ix_,'sigr163'))
            pb.x0(ix_('sigr163'),1) = 49.290216376;
        else
            pb.y0(find(pbm.congrps==ig_('sigr163')),1) = 49.290216376;
        end
        if(isKey(ix_,'x163'))
            pb.x0(ix_('x163'),1) = 1179.6881250;
        else
            pb.y0(find(pbm.congrps==ig_('x163')),1) = 1179.6881250;
        end
        if(isKey(ix_,'y163'))
            pb.x0(ix_('y163'),1) = 61477.257259;
        else
            pb.y0(find(pbm.congrps==ig_('y163')),1) = 61477.257259;
        end
        if(isKey(ix_,'w164'))
            pb.x0(ix_('w164'),1) = 6.1800000000;
        else
            pb.y0(find(pbm.congrps==ig_('w164')),1) = 6.1800000000;
        end
        if(isKey(ix_,'sigt164'))
            pb.x0(ix_('sigt164'),1) = 30.583856128;
        else
            pb.y0(find(pbm.congrps==ig_('sigt164')),1) = 30.583856128;
        end
        if(isKey(ix_,'sigr164'))
            pb.x0(ix_('sigr164'),1) = 47.422585757;
        else
            pb.y0(find(pbm.congrps==ig_('sigr164')),1) = 47.422585757;
        end
        if(isKey(ix_,'x164'))
            pb.x0(ix_('x164'),1) = 1182.7300000;
        else
            pb.y0(find(pbm.congrps==ig_('x164')),1) = 1182.7300000;
        end
        if(isKey(ix_,'y164'))
            pb.x0(ix_('y164'),1) = 61571.722555;
        else
            pb.y0(find(pbm.congrps==ig_('y164')),1) = 61571.722555;
        end
        if(isKey(ix_,'w165'))
            pb.x0(ix_('w165'),1) = 6.8125000000;
        else
            pb.y0(find(pbm.congrps==ig_('w165')),1) = 6.8125000000;
        end
        if(isKey(ix_,'sigt165'))
            pb.x0(ix_('sigt165'),1) = 28.755024885;
        else
            pb.y0(find(pbm.congrps==ig_('sigt165')),1) = 28.755024885;
        end
        if(isKey(ix_,'sigr165'))
            pb.x0(ix_('sigr165'),1) = 42.704591925;
        else
            pb.y0(find(pbm.congrps==ig_('sigr165')),1) = 42.704591925;
        end
        if(isKey(ix_,'x165'))
            pb.x0(ix_('x165'),1) = 1185.9781250;
        else
            pb.y0(find(pbm.congrps==ig_('x165')),1) = 1185.9781250;
        end
        if(isKey(ix_,'y165'))
            pb.x0(ix_('y165'),1) = 61667.948015;
        else
            pb.y0(find(pbm.congrps==ig_('y165')),1) = 61667.948015;
        end
        if(isKey(ix_,'w166'))
            pb.x0(ix_('w166'),1) = 7.4450000000;
        else
            pb.y0(find(pbm.congrps==ig_('w166')),1) = 7.4450000000;
        end
        if(isKey(ix_,'sigt166'))
            pb.x0(ix_('sigt166'),1) = 27.140969642;
        else
            pb.y0(find(pbm.congrps==ig_('sigt166')),1) = 27.140969642;
        end
        if(isKey(ix_,'sigr166'))
            pb.x0(ix_('sigr166'),1) = 38.769382165;
        else
            pb.y0(find(pbm.congrps==ig_('sigr166')),1) = 38.769382165;
        end
        if(isKey(ix_,'x166'))
            pb.x0(ix_('x166'),1) = 1189.5425000;
        else
            pb.y0(find(pbm.congrps==ig_('x166')),1) = 1189.5425000;
        end
        if(isKey(ix_,'y166'))
            pb.x0(ix_('y166'),1) = 61767.437546;
        else
            pb.y0(find(pbm.congrps==ig_('y166')),1) = 61767.437546;
        end
        if(isKey(ix_,'w167'))
            pb.x0(ix_('w167'),1) = 8.0775000000;
        else
            pb.y0(find(pbm.congrps==ig_('w167')),1) = 8.0775000000;
        end
        if(isKey(ix_,'sigt167'))
            pb.x0(ix_('sigt167'),1) = 25.689067179;
        else
            pb.y0(find(pbm.congrps==ig_('sigt167')),1) = 25.689067179;
        end
        if(isKey(ix_,'sigr167'))
            pb.x0(ix_('sigr167'),1) = 35.432574692;
        else
            pb.y0(find(pbm.congrps==ig_('sigr167')),1) = 35.432574692;
        end
        if(isKey(ix_,'x167'))
            pb.x0(ix_('x167'),1) = 1193.4231250;
        else
            pb.y0(find(pbm.congrps==ig_('x167')),1) = 1193.4231250;
        end
        if(isKey(ix_,'y167'))
            pb.x0(ix_('y167'),1) = 61869.829536;
        else
            pb.y0(find(pbm.congrps==ig_('y167')),1) = 61869.829536;
        end
        if(isKey(ix_,'w168'))
            pb.x0(ix_('w168'),1) = 8.7100000000;
        else
            pb.y0(find(pbm.congrps==ig_('w168')),1) = 8.7100000000;
        end
        if(isKey(ix_,'sigt168'))
            pb.x0(ix_('sigt168'),1) = 24.362140259;
        else
            pb.y0(find(pbm.congrps==ig_('sigt168')),1) = 24.362140259;
        end
        if(isKey(ix_,'sigr168'))
            pb.x0(ix_('sigr168'),1) = 32.563342979;
        else
            pb.y0(find(pbm.congrps==ig_('sigr168')),1) = 32.563342979;
        end
        if(isKey(ix_,'x168'))
            pb.x0(ix_('x168'),1) = 1197.6200000;
        else
            pb.y0(find(pbm.congrps==ig_('x168')),1) = 1197.6200000;
        end
        if(isKey(ix_,'y168'))
            pb.x0(ix_('y168'),1) = 61974.753956;
        else
            pb.y0(find(pbm.congrps==ig_('y168')),1) = 61974.753956;
        end
        if(isKey(ix_,'w169'))
            pb.x0(ix_('w169'),1) = 9.6750000000;
        else
            pb.y0(find(pbm.congrps==ig_('w169')),1) = 9.6750000000;
        end
        if(isKey(ix_,'sigt169'))
            pb.x0(ix_('sigt169'),1) = 22.820212971;
        else
            pb.y0(find(pbm.congrps==ig_('sigt169')),1) = 22.820212971;
        end
        if(isKey(ix_,'sigr169'))
            pb.x0(ix_('sigr169'),1) = 29.029268985;
        else
            pb.y0(find(pbm.congrps==ig_('sigr169')),1) = 29.029268985;
        end
        if(isKey(ix_,'x169'))
            pb.x0(ix_('x169'),1) = 1202.2162500;
        else
            pb.y0(find(pbm.congrps==ig_('x169')),1) = 1202.2162500;
        end
        if(isKey(ix_,'y169'))
            pb.x0(ix_('y169'),1) = 62082.998907;
        else
            pb.y0(find(pbm.congrps==ig_('y169')),1) = 62082.998907;
        end
        if(isKey(ix_,'w170'))
            pb.x0(ix_('w170'),1) = 10.640000000;
        else
            pb.y0(find(pbm.congrps==ig_('w170')),1) = 10.640000000;
        end
        if(isKey(ix_,'sigt170'))
            pb.x0(ix_('sigt170'),1) = 21.448594019;
        else
            pb.y0(find(pbm.congrps==ig_('sigt170')),1) = 21.448594019;
        end
        if(isKey(ix_,'sigr170'))
            pb.x0(ix_('sigr170'),1) = 26.114653928;
        else
            pb.y0(find(pbm.congrps==ig_('sigr170')),1) = 26.114653928;
        end
        if(isKey(ix_,'x170'))
            pb.x0(ix_('x170'),1) = 1207.2950000;
        else
            pb.y0(find(pbm.congrps==ig_('x170')),1) = 1207.2950000;
        end
        if(isKey(ix_,'y170'))
            pb.x0(ix_('y170'),1) = 62195.248557;
        else
            pb.y0(find(pbm.congrps==ig_('y170')),1) = 62195.248557;
        end
        if(isKey(ix_,'w171'))
            pb.x0(ix_('w171'),1) = 11.820000000;
        else
            pb.y0(find(pbm.congrps==ig_('w171')),1) = 11.820000000;
        end
        if(isKey(ix_,'sigt171'))
            pb.x0(ix_('sigt171'),1) = 20.072296567;
        else
            pb.y0(find(pbm.congrps==ig_('sigt171')),1) = 20.072296567;
        end
        if(isKey(ix_,'sigr171'))
            pb.x0(ix_('sigr171'),1) = 23.231954054;
        else
            pb.y0(find(pbm.congrps==ig_('sigr171')),1) = 23.231954054;
        end
        if(isKey(ix_,'x171'))
            pb.x0(ix_('x171'),1) = 1212.9100000;
        else
            pb.y0(find(pbm.congrps==ig_('x171')),1) = 1212.9100000;
        end
        if(isKey(ix_,'y171'))
            pb.x0(ix_('y171'),1) = 62311.615454;
        else
            pb.y0(find(pbm.congrps==ig_('y171')),1) = 62311.615454;
        end
        if(isKey(ix_,'w172'))
            pb.x0(ix_('w172'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w172')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt172'))
            pb.x0(ix_('sigt172'),1) = 18.833069241;
        else
            pb.y0(find(pbm.congrps==ig_('sigt172')),1) = 18.833069241;
        end
        if(isKey(ix_,'sigr172'))
            pb.x0(ix_('sigr172'),1) = 20.850204821;
        else
            pb.y0(find(pbm.congrps==ig_('sigr172')),1) = 20.850204821;
        end
        if(isKey(ix_,'x172'))
            pb.x0(ix_('x172'),1) = 1219.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x172')),1) = 1219.1150000;
        end
        if(isKey(ix_,'y172'))
            pb.x0(ix_('y172'),1) = 62432.136565;
        else
            pb.y0(find(pbm.congrps==ig_('y172')),1) = 62432.136565;
        end
        if(isKey(ix_,'w173'))
            pb.x0(ix_('w173'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w173')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt173'))
            pb.x0(ix_('sigt173'),1) = 18.214200054;
        else
            pb.y0(find(pbm.congrps==ig_('sigt173')),1) = 18.214200054;
        end
        if(isKey(ix_,'sigr173'))
            pb.x0(ix_('sigr173'),1) = 20.564709344;
        else
            pb.y0(find(pbm.congrps==ig_('sigr173')),1) = 20.564709344;
        end
        if(isKey(ix_,'x173'))
            pb.x0(ix_('x173'),1) = 1225.6150000;
        else
            pb.y0(find(pbm.congrps==ig_('x173')),1) = 1225.6150000;
        end
        if(isKey(ix_,'y173'))
            pb.x0(ix_('y173'),1) = 62552.540190;
        else
            pb.y0(find(pbm.congrps==ig_('y173')),1) = 62552.540190;
        end
        if(isKey(ix_,'w174'))
            pb.x0(ix_('w174'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w174')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt174'))
            pb.x0(ix_('sigt174'),1) = 17.589839561;
        else
            pb.y0(find(pbm.congrps==ig_('sigt174')),1) = 17.589839561;
        end
        if(isKey(ix_,'sigr174'))
            pb.x0(ix_('sigr174'),1) = 20.276848479;
        else
            pb.y0(find(pbm.congrps==ig_('sigr174')),1) = 20.276848479;
        end
        if(isKey(ix_,'x174'))
            pb.x0(ix_('x174'),1) = 1232.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x174')),1) = 1232.1150000;
        end
        if(isKey(ix_,'y174'))
            pb.x0(ix_('y174'),1) = 62668.903319;
        else
            pb.y0(find(pbm.congrps==ig_('y174')),1) = 62668.903319;
        end
        if(isKey(ix_,'w175'))
            pb.x0(ix_('w175'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w175')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt175'))
            pb.x0(ix_('sigt175'),1) = 16.959939451;
        else
            pb.y0(find(pbm.congrps==ig_('sigt175')),1) = 16.959939451;
        end
        if(isKey(ix_,'sigr175'))
            pb.x0(ix_('sigr175'),1) = 19.986619910;
        else
            pb.y0(find(pbm.congrps==ig_('sigr175')),1) = 19.986619910;
        end
        if(isKey(ix_,'x175'))
            pb.x0(ix_('x175'),1) = 1238.6150000;
        else
            pb.y0(find(pbm.congrps==ig_('x175')),1) = 1238.6150000;
        end
        if(isKey(ix_,'y175'))
            pb.x0(ix_('y175'),1) = 62781.190101;
        else
            pb.y0(find(pbm.congrps==ig_('y175')),1) = 62781.190101;
        end
        if(isKey(ix_,'w176'))
            pb.x0(ix_('w176'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w176')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt176'))
            pb.x0(ix_('sigt176'),1) = 16.324451370;
        else
            pb.y0(find(pbm.congrps==ig_('sigt176')),1) = 16.324451370;
        end
        if(isKey(ix_,'sigr176'))
            pb.x0(ix_('sigr176'),1) = 19.694021169;
        else
            pb.y0(find(pbm.congrps==ig_('sigr176')),1) = 19.694021169;
        end
        if(isKey(ix_,'x176'))
            pb.x0(ix_('x176'),1) = 1245.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x176')),1) = 1245.1150000;
        end
        if(isKey(ix_,'y176'))
            pb.x0(ix_('y176'),1) = 62889.364371;
        else
            pb.y0(find(pbm.congrps==ig_('y176')),1) = 62889.364371;
        end
        if(isKey(ix_,'w177'))
            pb.x0(ix_('w177'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w177')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt177'))
            pb.x0(ix_('sigt177'),1) = 15.683326910;
        else
            pb.y0(find(pbm.congrps==ig_('sigt177')),1) = 15.683326910;
        end
        if(isKey(ix_,'sigr177'))
            pb.x0(ix_('sigr177'),1) = 19.399049638;
        else
            pb.y0(find(pbm.congrps==ig_('sigr177')),1) = 19.399049638;
        end
        if(isKey(ix_,'x177'))
            pb.x0(ix_('x177'),1) = 1251.6150000;
        else
            pb.y0(find(pbm.congrps==ig_('x177')),1) = 1251.6150000;
        end
        if(isKey(ix_,'y177'))
            pb.x0(ix_('y177'),1) = 62993.389650;
        else
            pb.y0(find(pbm.congrps==ig_('y177')),1) = 62993.389650;
        end
        if(isKey(ix_,'w178'))
            pb.x0(ix_('w178'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w178')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt178'))
            pb.x0(ix_('sigt178'),1) = 15.036517615;
        else
            pb.y0(find(pbm.congrps==ig_('sigt178')),1) = 15.036517615;
        end
        if(isKey(ix_,'sigr178'))
            pb.x0(ix_('sigr178'),1) = 19.101702555;
        else
            pb.y0(find(pbm.congrps==ig_('sigr178')),1) = 19.101702555;
        end
        if(isKey(ix_,'x178'))
            pb.x0(ix_('x178'),1) = 1258.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x178')),1) = 1258.1150000;
        end
        if(isKey(ix_,'y178'))
            pb.x0(ix_('y178'),1) = 63093.229145;
        else
            pb.y0(find(pbm.congrps==ig_('y178')),1) = 63093.229145;
        end
        if(isKey(ix_,'w179'))
            pb.x0(ix_('w179'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w179')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigt179'))
            pb.x0(ix_('sigt179'),1) = 14.383974972;
        else
            pb.y0(find(pbm.congrps==ig_('sigt179')),1) = 14.383974972;
        end
        if(isKey(ix_,'sigr179'))
            pb.x0(ix_('sigr179'),1) = 18.801977014;
        else
            pb.y0(find(pbm.congrps==ig_('sigr179')),1) = 18.801977014;
        end
        if(isKey(ix_,'x179'))
            pb.x0(ix_('x179'),1) = 1264.6150000;
        else
            pb.y0(find(pbm.congrps==ig_('x179')),1) = 1264.6150000;
        end
        if(isKey(ix_,'y179'))
            pb.x0(ix_('y179'),1) = 63188.845746;
        else
            pb.y0(find(pbm.congrps==ig_('y179')),1) = 63188.845746;
        end
        if(isKey(ix_,'w180'))
            pb.x0(ix_('w180'),1) = 13.000000000;
        else
            pb.y0(find(pbm.congrps==ig_('w180')),1) = 13.000000000;
        end
        if(isKey(ix_,'sigr180'))
            pb.x0(ix_('sigr180'),1) = 18.500000000;
        else
            pb.y0(find(pbm.congrps==ig_('sigr180')),1) = 18.500000000;
        end
        if(isKey(ix_,'sigt180'))
            pb.x0(ix_('sigt180'),1) = 13.725650411;
        else
            pb.y0(find(pbm.congrps==ig_('sigt180')),1) = 13.725650411;
        end
        if(isKey(ix_,'x180'))
            pb.x0(ix_('x180'),1) = 1271.1150000;
        else
            pb.y0(find(pbm.congrps==ig_('x180')),1) = 1271.1150000;
        end
        if(isKey(ix_,'y180'))
            pb.x0(ix_('y180'),1) = 63280.202028;
        else
            pb.y0(find(pbm.congrps==ig_('y180')),1) = 63280.202028;
        end
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
        for k=v_('0'):v_('K')
            ename = ['WSR',int2str(k)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['w',int2str(k)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['sigr',int2str(k)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for k=v_('0'):v_('K')
            ename = ['WST',int2str(k)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['w',int2str(k)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['sigt',int2str(k)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        v_('rk') = v_('ri');
        for k=v_('0'):v_('K-1')
            v_('rk+1') = v_('rk')+v_('dr');
            v_('k+1') = 1+k;
            v_('-rk') = -1.0*v_('rk');
            ig = ig_(['SR',int2str(k)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WSR',int2str(k)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-rk');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WSR',int2str(round(v_('k+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('rk+1');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WST',int2str(k)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-dr/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WST',int2str(round(v_('k+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-dr/2');
            ig = ig_(['STAy',int2str(k)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WST',int2str(k)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-dr/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WST',int2str(round(v_('k+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-dr/2');
            v_('rk') = v_('rk+1');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 5.0;
% LO SOLUTION            7.872067544
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-LQR2-RN-905-1081';
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

