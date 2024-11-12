function varargout = MRIBASIS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ******** 
% 
%    An optmization problem arising in the design of medical apparatus.
% 
%    Source:
%    Contribution from a LANCELOT user.
% 
%    SIF input: Arie Quist, TU Delft (NL), 1994.
%    Adaptation for CUTE: Ph. Toint, November 1994.
% 
%    classification = 'C-CLOR2-MY-36-55'
% 
%    useful constants
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MRIBASIS';

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
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('R2') = v_('2');
        v_('k1') = 22.8443;
        v_('k2') = 12.4402;
        v_('k3') = 5.23792;
        v_('k4') = 5.12238;
        v_('k5') = 6.44999;
        v_('k6') = 5.32383;
        v_('k7') = 0.58392;
        v_('k8') = 3.94584;
        v_('k9') = -2.75584;
        v_('k10') = 32.0669;
        v_('k11') = 18.2179;
        v_('k12') = 31.7496;
        v_('S1,1') = -0.377126;
        v_('S1,2') = 0.919679;
        v_('S1,3') = 0.109389;
        v_('S2,1') = 0.634857;
        v_('S2,2') = 0.170704;
        v_('S2,3') = 0.753536;
        v_('S3,1') = 0.674338;
        v_('S3,2') = 0.353624;
        v_('S3,3') = -0.648242;
        v_('xlo1') = 1*v_('k7');
        v_('k10/2') = v_('k10')/v_('R2');
        v_('k8/2') = v_('k8')/v_('R2');
        v_('xup1') = v_('k10/2')-v_('k8/2');
        v_('xlo2') = v_('k10/2')+v_('k8/2');
        v_('k13') = v_('k7')*v_('k5');
        v_('k14') = v_('k8/2')*v_('k6');
        v_('xup2') = 1*v_('k12');
        v_('-k1') = -1*v_('k1');
        v_('-k2') = -1*v_('k2');
        v_('k4-') = v_('k4')-v_('k14');
        v_('-k3') = -1*v_('k3');
        v_('-S1,1') = -1*v_('S1,1');
        v_('-S1,2') = -1*v_('S1,2');
        v_('-S1,3') = -1*v_('S1,3');
        v_('-S2,1') = -1*v_('S2,1');
        v_('-S2,2') = -1*v_('S2,2');
        v_('-S2,3') = -1*v_('S2,3');
        v_('-S3,1') = -1*v_('S3,1');
        v_('-S3,2') = -1*v_('S3,2');
        v_('-S3,3') = -1*v_('S3,3');
        v_('2S1,1') = 2*v_('S1,1');
        v_('2S1,2') = 2*v_('S1,2');
        v_('2S1,3') = 2*v_('S1,3');
        v_('2S2,1') = 2*v_('S2,1');
        v_('2S2,2') = 2*v_('S2,2');
        v_('2S2,3') = 2*v_('S2,3');
        v_('2S3,1') = 2*v_('S3,1');
        v_('2S3,2') = 2*v_('S3,2');
        v_('2S3,3') = 2*v_('S3,3');
        v_('-2S1,1') = -2*v_('S1,1');
        v_('-2S1,2') = -2*v_('S1,2');
        v_('-2S1,3') = -2*v_('S1,3');
        v_('-2S2,1') = -2*v_('S2,1');
        v_('-2S2,2') = -2*v_('S2,2');
        v_('-2S2,3') = -2*v_('S2,3');
        v_('-2S3,1') = -2*v_('S3,1');
        v_('-2S3,2') = -2*v_('S3,2');
        v_('-2S3,3') = -2*v_('S3,3');
        v_('Llo1,1') = v_('S1,1')*v_('k5');
        v_('Llo1,2') = v_('S1,2')*v_('k5');
        v_('Llo1,3') = v_('S1,3')*v_('k5');
        v_('Lup1,1') = v_('S1,1')*v_('k6');
        v_('Lup1,2') = v_('S1,2')*v_('k6');
        v_('Lup1,3') = v_('S1,3')*v_('k6');
        v_('Llo2,1') = v_('S1,1')*v_('k6');
        v_('Llo2,2') = v_('S1,2')*v_('k6');
        v_('Llo2,3') = v_('S1,3')*v_('k6');
        v_('4') = 4;
        v_('xm') = 6;
        v_('Lm') = 4;
        v_('xm-') = -1+v_('xm');
        v_('xm-2') = -2+v_('xm');
        v_('Lm-') = -1+v_('Lm');
        v_('Lm-2') = -2+v_('Lm');
        v_('R12') = 12;
        v_('1/12') = 1/v_('R12');
        v_('-1/12') = -1*v_('1/12');
        v_('R0') = v_('0');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for j=v_('1'):v_('2')
            for k=v_('1'):v_('xm')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(j),',',int2str(k)],ix_);
                pb.xnames{iv} = ['X',int2str(j),',',int2str(k)];
            end
        end
        for i=v_('1'):v_('3')
            for j=v_('1'):v_('2')
                for k=v_('1'):v_('Lm')
                    [iv,ix_] =...
                          s2mpjlib('ii',['L',int2str(i),',',int2str(j),',',int2str(k)],ix_);
                    pb.xnames{iv} = ['L',int2str(i),',',int2str(j),',',int2str(k)];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','Object',ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('xm')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        for j=v_('1'):v_('2')
            for k=v_('2'):v_('xm-2')
                v_('k+') = 1+k;
                [ig,ig_] = s2mpjlib('ii',['PS',int2str(j),',',int2str(k)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['PS',int2str(j),',',int2str(k)];
                iv = ix_(['X',int2str(j),',',int2str(round(v_('k+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1;
                end
                iv = ix_(['X',int2str(j),',',int2str(k)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1;
                end
            end
        end
        [ig,ig_] = s2mpjlib('ii','PL',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'PL';
        iv = ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('xm')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('xm-')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        for i=v_('1'):v_('3')
            for j=v_('1'):v_('2')
                for k=v_('1'):v_('Lm-')
                    v_('k+') = 1+k;
                    v_('2k') = 2*k;
                    v_('2k-') = -1+v_('2k');
                    [ig,ig_] =...
                          s2mpjlib('ii',['SU',int2str(i),',',int2str(j),',',int2str(k)],ig_);
                    gtype{ig}  = '<=';
                    cnames{ig} = ['SU',int2str(i),',',int2str(j),',',int2str(k)];
                    iv = ix_(['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1;
                    end
                    iv = ix_(['L',int2str(i),',',int2str(j),',',int2str(k)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = -1+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = -1;
                    end
                    iv = ix_(['X',int2str(j),',',int2str(round(v_('2k')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-k1')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-k1');
                    end
                    iv = ix_(['X',int2str(j),',',int2str(round(v_('2k-')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('k1')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('k1');
                    end
                    [ig,ig_] =...
                          s2mpjlib('ii',['SL',int2str(i),',',int2str(j),',',int2str(k)],ig_);
                    gtype{ig}  = '>=';
                    cnames{ig} = ['SL',int2str(i),',',int2str(j),',',int2str(k)];
                    iv = ix_(['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = 1+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = 1;
                    end
                    iv = ix_(['L',int2str(i),',',int2str(j),',',int2str(k)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = -1+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = -1;
                    end
                    iv = ix_(['X',int2str(j),',',int2str(round(v_('2k')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('k1')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('k1');
                    end
                    iv = ix_(['X',int2str(j),',',int2str(round(v_('2k-')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-k1')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-k1');
                    end
                end
            end
        end
        for i=v_('1'):v_('3')
            for k=v_('2'):v_('Lm-')
                [ig,ig_] = s2mpjlib('ii','cc1',ig_);
                gtype{ig}  = '==';
                cnames{ig} = 'cc1';
                iv = ix_(['L',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) =...
                          v_(['S',int2str(round(v_('3'))),',',int2str(i)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['S',int2str(round(v_('3'))),',',int2str(i)]);
                end
            end
        end
        for i=v_('1'):v_('3')
            [ig,ig_] = s2mpjlib('ii',['c2const',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['c2const',int2str(i)];
            iv =...
                  ix_(['L',int2str(i),',',int2str(round(v_('2'))),',',int2str(round(v_('Lm')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('k10')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('k10');
            end
        end
        [ig,ig_] = s2mpjlib('ii','c3con1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'c3con1';
        [ig,ig_] = s2mpjlib('ii','c3con2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'c3con2';
        [ig,ig_] = s2mpjlib('ii','c4const',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'c4const';
        for i=v_('1'):v_('3')
            [ig,ig_] = s2mpjlib('ii',['c5con',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['c5con',int2str(i)];
        end
        for i=v_('1'):v_('2')
            [ig,ig_] = s2mpjlib('ii',['c6cn',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['c6cn',int2str(i)];
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
        v_('Opmr1') = v_('k3')*v_('S2,1');
        v_('Opmr2') = v_('k3')*v_('S2,2');
        v_('Opmr3') = v_('k3')*v_('S2,3');
        for i=v_('1'):v_('3')
            pbm.gconst(ig_(['c2const',int2str(i)])) = v_(['Opmr',int2str(i)]);
        end
        pbm.gconst(ig_('c3con1')) = v_('k4-');
        pbm.gconst(ig_('c3con2')) = v_('k13');
        pbm.gconst(ig_('c5con1')) = v_('k13');
        pbm.gconst(ig_('c5con2')) = v_('-k3');
        pbm.gconst(ig_('c5con3')) = v_('k9');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for j=v_('1'):v_('2')
            for k=v_('2'):v_('xm-')
                pb.xlower(ix_(['X',int2str(j),',',int2str(k)]),1) = v_(['xlo',int2str(j)]);
                pb.xupper(ix_(['X',int2str(j),',',int2str(k)])) = v_(['xup',int2str(j)]);
            end
            pb.xlower(ix_(['X',int2str(j),',',int2str(round(v_('1')))]),1) =...
                  v_(['xlo',int2str(j)]);
            pb.xupper(ix_(['X',int2str(j),',',int2str(round(v_('1')))]),1) =...
                  v_(['xlo',int2str(j)]);
        end
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('xm')))]),1) = v_('xup1');
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('xm')))]),1) = v_('xup1');
        pb.xlower(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('xm')))]),1) = v_('k11');
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('xm')))])) =...
              v_('k12');
        for i=v_('1'):v_('3')
            for j=v_('1'):v_('2')
                for k=v_('2'):v_('Lm-')
                    pb.xlower(ix_(['L',int2str(i),',',int2str(j),',',int2str(k)]),1) =...
                          v_('-k2');
                    pb.xupper(ix_(['L',int2str(i),',',int2str(j),',',int2str(k)])) = v_('k2');
                end
                pb.xlower(ix_(['L',int2str(i),',',int2str(j),',',int2str(round(v_('1')))]),1) = v_(['Llo',int2str(j),',',int2str(i)]);
                pb.xupper(ix_(['L',int2str(i),',',int2str(j),',',int2str(round(v_('1')))]),1) = v_(['Llo',int2str(j),',',int2str(i)]);
            end
            pb.xlower(ix_(['L',int2str(i),',',int2str(round(v_('1'))),',',int2str(round(v_('Lm')))]),1) = v_(['Lup',int2str(round(v_('1'))),',',int2str(i)]);
            pb.xupper(ix_(['L',int2str(i),',',int2str(round(v_('1'))),',',int2str(round(v_('Lm')))]),1) = v_(['Lup',int2str(round(v_('1'))),',',int2str(i)]);
            pb.xlower(ix_(['L',int2str(i),',',int2str(round(v_('2'))),',',int2str(round(v_('Lm')))]),1) = v_('-k2');
            pb.xupper(ix_(['L',int2str(i),',',int2str(round(v_('2'))),',',int2str(round(v_('Lm')))])) = v_('k2');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('intlen1') = v_('xup1')-v_('xlo1');
        v_('Rxm-') = v_('xm-');
        v_('dx1') = v_('intlen1')/v_('Rxm-');
        v_('intlen2') = v_('k11')-v_('xlo2');
        v_('dx2') = v_('intlen2')/v_('Rxm-');
        for k=v_('1'):v_('xm-')
            v_('Rk') = k;
            v_('dist1') = v_('dx1')*v_('Rk');
            v_('strtv1') = v_('xlo1')+v_('dist1');
            v_('dist2') = v_('dx2')*v_('Rk');
            v_('strtv2') = v_('xlo2')+v_('dist2');
            v_('k+') = 1+k;
            pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('k+')))]),1) =...
                  v_('strtv1');
            pb.x0(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('k+')))]),1) =...
                  v_('strtv2');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'euv1',iet_);
        elftv{it}{1} = 'v1';
        elftv{it}{2} = 'v2';
        elftv{it}{3} = 'v3';
        elftv{it}{4} = 'v4';
        [it,iet_] = s2mpjlib( 'ii', 'euv2',iet_);
        elftv{it}{1} = 'v1';
        elftv{it}{2} = 'v2';
        elftv{it}{3} = 'v3';
        [it,iet_] = s2mpjlib( 'ii', 'euvw1',iet_);
        elftv{it}{1} = 'v1';
        elftv{it}{2} = 'v2';
        elftv{it}{3} = 'v3';
        elftp{it}{1} = 'p1';
        [it,iet_] = s2mpjlib( 'ii', 'emo',iet_);
        elftv{it}{1} = 's1';
        elftv{it}{2} = 's2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for i=v_('1'):v_('3')
            for j=v_('1'):v_('2')
                for k=v_('1'):v_('Lm-')
                    v_('2k') = 2*k;
                    v_('k+') = 1+k;
                    v_('2k-') = -1+v_('2k');
                    ename = ['e1',int2str(i),',',int2str(j),',',int2str(k)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'euv1';
                    ielftype(ie) = iet_('euv1');
                    vname = ['X',int2str(j),',',int2str(round(v_('2k')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k-')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(k)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v4',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['e3',int2str(i),',',int2str(j),',',int2str(k)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'euvw1';
                    ielftype(ie) = iet_('euvw1');
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(k)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k-')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    [~,posep] = ismember('p1',elftp{ielftype(ie)});
                    pbm.elpar{ie}(posep) = v_('1/12');
                    ename = ['e5',int2str(i),',',int2str(j),',',int2str(k)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'euvw1';
                    ielftype(ie) = iet_('euvw1');
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k-')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    [~,posep] = ismember('p1',elftp{ielftype(ie)});
                    pbm.elpar{ie}(posep) = v_('-1/12');
                end
                for k=v_('1'):v_('Lm-2')
                    v_('2k') = 2*k;
                    v_('k+') = 1+k;
                    v_('2k+') = 1+v_('2k');
                    ename = ['e2',int2str(i),',',int2str(j),',',int2str(k)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'euv2';
                    ielftype(ie) = iet_('euv2');
                    vname = ['X',int2str(j),',',int2str(round(v_('2k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['e4',int2str(i),',',int2str(j),',',int2str(k)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'euvw1';
                    ielftype(ie) = iet_('euvw1');
                    vname = ['L',int2str(i),',',int2str(j),',',int2str(round(v_('k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k+')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(j),',',int2str(round(v_('2k')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('v3',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    [~,posep] = ismember('p1',elftp{ielftype(ie)});
                    pbm.elpar{ie}(posep) = v_('R0');
                end
            end
        end
        for i=v_('1'):v_('3')
            ename = ['factr',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'emo';
            ielftype(ie) = iet_('emo');
            vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('xm')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('s1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname =...
                  ['L',int2str(i),',',int2str(round(v_('2'))),',',int2str(round(v_('Lm')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('s2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for i=v_('1'):v_('3')
            ig = ig_(['c2const',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['factr',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for j=v_('1'):v_('2')
            for i=v_('1'):v_('3')
                for k=v_('1'):v_('Lm-')
                    ig = ig_(['c3con',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e1',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(round(v_('1'))),',',int2str(i)]);
                end
                for k=v_('1'):v_('Lm-2')
                    ig = ig_(['c3con',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e2',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(round(v_('1'))),',',int2str(i)]);
                end
            end
        end
        for i=v_('1'):v_('3')
            for k=v_('1'):v_('Lm-')
                ig = ig_('c4const');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['e1',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['S',int2str(round(v_('2'))),',',int2str(i)]);
            end
            for k=v_('1'):v_('Lm-2')
                ig = ig_('c4const');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['e2',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['S',int2str(round(v_('2'))),',',int2str(i)]);
            end
        end
        for j=v_('1'):v_('3')
            for i=v_('1'):v_('3')
                for k=v_('1'):v_('Lm-')
                    ig = ig_(['c5con',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e1',int2str(i),',',int2str(round(v_('1'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['-S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e1',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(j),',',int2str(i)]);
                end
                for k=v_('1'):v_('Lm-2')
                    ig = ig_(['c5con',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e2',int2str(i),',',int2str(round(v_('1'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['-S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e2',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(j),',',int2str(i)]);
                end
            end
        end
        for j=v_('1'):v_('2')
            for i=v_('1'):v_('3')
                for k=v_('1'):v_('Lm-')
                    ig = ig_(['c6cn',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e3',int2str(i),',',int2str(round(v_('1'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e5',int2str(i),',',int2str(round(v_('1'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e3',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['-S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e5',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['-S',int2str(j),',',int2str(i)]);
                end
                for k=v_('1'):v_('Lm-2')
                    ig = ig_(['c6cn',int2str(j)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e4',int2str(i),',',int2str(round(v_('1'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['2S',int2str(j),',',int2str(i)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['e4',int2str(i),',',int2str(round(v_('2'))),',',int2str(k)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['-2S',int2str(j),',',int2str(i)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               18.2179000000
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-MY-36-55';
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

    case 'euv1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = 0.5e0*IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = 0.5e0*IV_(2);
            g_(2,1) = 0.5e0*IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 0.5e0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'euv2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'euvw1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        d14 = 0.25e0;
        varargout{1} =...
              EV_(1)*(EV_(2)-EV_(3))*((d14+pbm.elpar{iel_}(1))*EV_(2)+(d14-pbm.elpar{iel_}(1))*EV_(3));
        if(nargout>1)
            g_(1,1) =...
                  (EV_(2)-EV_(3))*((d14+pbm.elpar{iel_}(1))*EV_(2)+(d14-pbm.elpar{iel_}(1))*EV_(3));
            g_(2,1) =...
                  EV_(1)*2.0e0*(EV_(2)*(d14+pbm.elpar{iel_}(1))-pbm.elpar{iel_}(1)*EV_(3));
            g_(3,1) =...
                  EV_(1)*2.0e0*(-EV_(3)*(d14-pbm.elpar{iel_}(1))-pbm.elpar{iel_}(1)*EV_(2));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = 2.0e0*(EV_(2)*(d14+pbm.elpar{iel_}(1))-pbm.elpar{iel_}(1)*EV_(3));
                H_(2,1) = H_(1,2);
                H_(1,3) =...
                      2.0e0*(-EV_(3)*(d14-pbm.elpar{iel_}(1))-pbm.elpar{iel_}(1)*EV_(2));
                H_(3,1) = H_(1,3);
                H_(2,2) = EV_(1)*2.0e0*(d14+pbm.elpar{iel_}(1));
                H_(2,3) = -EV_(1)*2.0e0*pbm.elpar{iel_}(1);
                H_(3,2) = H_(2,3);
                H_(3,3) = -EV_(1)*2.0e0*(d14-pbm.elpar{iel_}(1));
                varargout{3} = H_;
            end
        end

    case 'emo'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(2)*EV_(1);
        if(nargout>1)
            g_(1,1) = -EV_(2);
            g_(2,1) = -EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0e0;
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

