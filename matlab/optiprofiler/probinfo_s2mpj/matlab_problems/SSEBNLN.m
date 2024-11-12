function varargout = SSEBNLN(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SSEBNLN
%    *********
%    The Power Generation problem for the SSGB.
% 
%    Source:
%    N. Gould, private communication.
% 
%    SIF input: Nick Gould, 23 October 1989
% 
%    classification = 'C-CLQR2-RN-194-96'
% 
%    period is the number of time periods
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SSEBNLN';

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
        v_('HOURS') = 24;
        v_('DAYS') = 1;
        v_('PERIOD') = v_('HOURS')*v_('DAYS');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('Z1') = 517.0;
        v_('D1,1') = 578.0;
        v_('D1,2') = 517.0;
        v_('D1,3') = 461.0;
        v_('D1,4') = 369.0;
        v_('D1,5') = 299.0;
        v_('D1,6') = 269.0;
        v_('D1,7') = 370.0;
        v_('D1,8') = 559.0;
        v_('D1,9') = 689.0;
        v_('D1,10') = 728.0;
        v_('D1,11') = 683.0;
        v_('D1,12') = 626.0;
        v_('D1,13') = 619.0;
        v_('D1,14') = 586.0;
        v_('D1,15') = 582.0;
        v_('D1,16') = 625.0;
        v_('D1,17') = 821.0;
        v_('D1,18') = 883.0;
        v_('D1,19') = 768.0;
        v_('D1,20') = 711.0;
        v_('D1,21') = 677.0;
        v_('D1,22') = 630.0;
        v_('D1,23') = 545.0;
        v_('D1,24') = 565.0;
        v_('Z2') = 400.0;
        v_('D2,1') = 631.0;
        v_('D2,2') = 574.0;
        v_('D2,3') = 521.0;
        v_('D2,4') = 446.0;
        v_('D2,5') = 359.0;
        v_('D2,6') = 336.0;
        v_('D2,7') = 420.0;
        v_('D2,8') = 588.0;
        v_('D2,9') = 697.0;
        v_('D2,10') = 732.0;
        v_('D2,11') = 713.0;
        v_('D2,12') = 682.0;
        v_('D2,13') = 695.0;
        v_('D2,14') = 651.0;
        v_('D2,15') = 645.0;
        v_('D2,16') = 664.0;
        v_('D2,17') = 816.0;
        v_('D2,18') = 858.0;
        v_('D2,19') = 760.0;
        v_('D2,20') = 700.0;
        v_('D2,21') = 659.0;
        v_('D2,22') = 623.0;
        v_('D2,23') = 517.0;
        v_('D2,24') = 542.0;
        v_('Z3') = 1017.0;
        v_('D3,1') = 582.0;
        v_('D3,2') = 501.0;
        v_('D3,3') = 443.0;
        v_('D3,4') = 367.0;
        v_('D3,5') = 288.0;
        v_('D3,6') = 265.0;
        v_('D3,7') = 349.0;
        v_('D3,8') = 503.0;
        v_('D3,9') = 663.0;
        v_('D3,10') = 651.0;
        v_('D3,11') = 625.0;
        v_('D3,12') = 596.0;
        v_('D3,13') = 608.0;
        v_('D3,14') = 566.0;
        v_('D3,15') = 555.0;
        v_('D3,16') = 584.0;
        v_('D3,17') = 763.0;
        v_('D3,18') = 803.0;
        v_('D3,19') = 710.0;
        v_('D3,20') = 648.0;
        v_('D3,21') = 626.0;
        v_('D3,22') = 590.0;
        v_('D3,23') = 486.0;
        v_('D3,24') = 540.0;
        v_('Z4') = 667.0;
        v_('D4,1') = 602.0;
        v_('D4,2') = 533.0;
        v_('D4,3') = 450.0;
        v_('D4,4') = 378.0;
        v_('D4,5') = 298.0;
        v_('D4,6') = 272.0;
        v_('D4,7') = 369.0;
        v_('D4,8') = 539.0;
        v_('D4,9') = 647.0;
        v_('D4,10') = 652.0;
        v_('D4,11') = 607.0;
        v_('D4,12') = 585.0;
        v_('D4,13') = 587.0;
        v_('D4,14') = 549.0;
        v_('D4,15') = 535.0;
        v_('D4,16') = 564.0;
        v_('D4,17') = 748.0;
        v_('D4,18') = 808.0;
        v_('D4,19') = 710.0;
        v_('D4,20') = 646.0;
        v_('D4,21') = 620.0;
        v_('D4,22') = 581.0;
        v_('D4,23') = 483.0;
        v_('D4,24') = 514.0;
        v_('Z5') = 600.0;
        v_('D5,1') = 579.0;
        v_('D5,2') = 518.0;
        v_('D5,3') = 447.0;
        v_('D5,4') = 355.0;
        v_('D5,5') = 284.0;
        v_('D5,6') = 261.0;
        v_('D5,7') = 348.0;
        v_('D5,8') = 530.0;
        v_('D5,9') = 644.0;
        v_('D5,10') = 648.0;
        v_('D5,11') = 607.0;
        v_('D5,12') = 570.0;
        v_('D5,13') = 577.0;
        v_('D5,14') = 536.0;
        v_('D5,15') = 544.0;
        v_('D5,16') = 554.0;
        v_('D5,17') = 716.0;
        v_('D5,18') = 765.0;
        v_('D5,19') = 676.0;
        v_('D5,20') = 631.0;
        v_('D5,21') = 576.0;
        v_('D5,22') = 528.0;
        v_('D5,23') = 445.0;
        v_('D5,24') = 520.0;
        v_('Z6') = 421.0;
        v_('D6,1') = 618.0;
        v_('D6,2') = 547.0;
        v_('D6,3') = 430.0;
        v_('D6,4') = 327.0;
        v_('D6,5') = 249.0;
        v_('D6,6') = 211.0;
        v_('D6,7') = 227.0;
        v_('D6,8') = 258.0;
        v_('D6,9') = 347.0;
        v_('D6,10') = 491.0;
        v_('D6,11') = 524.0;
        v_('D6,12') = 492.0;
        v_('D6,13') = 467.0;
        v_('D6,14') = 418.0;
        v_('D6,15') = 358.0;
        v_('D6,16') = 378.0;
        v_('D6,17') = 544.0;
        v_('D6,18') = 666.0;
        v_('D6,19') = 589.0;
        v_('D6,20') = 533.0;
        v_('D6,21') = 494.0;
        v_('D6,22') = 460.0;
        v_('D6,23') = 404.0;
        v_('D6,24') = 512.0;
        v_('Z7') = 425.0;
        v_('D7,1') = 615.0;
        v_('D7,2') = 587.0;
        v_('D7,3') = 450.0;
        v_('D7,4') = 320.0;
        v_('D7,5') = 235.0;
        v_('D7,6') = 198.0;
        v_('D7,7') = 195.0;
        v_('D7,8') = 173.0;
        v_('D7,9') = 197.0;
        v_('D7,10') = 349.0;
        v_('D7,11') = 441.0;
        v_('D7,12') = 459.0;
        v_('D7,13') = 485.0;
        v_('D7,14') = 445.0;
        v_('D7,15') = 410.0;
        v_('D7,16') = 421.0;
        v_('D7,17') = 568.0;
        v_('D7,18') = 643.0;
        v_('D7,19') = 596.0;
        v_('D7,20') = 566.0;
        v_('D7,21') = 541.0;
        v_('D7,22') = 532.0;
        v_('D7,23') = 454.0;
        v_('D7,24') = 511.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] =...
              s2mpjlib('ii',['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))],ix_);
        pb.xnames{iv} =...
              ['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))];
        [iv,ix_] =...
              s2mpjlib('ii',['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))],ix_);
        pb.xnames{iv} =...
              ['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))];
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                [iv,ix_] = s2mpjlib('ii',['P1',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['P1',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['P2',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['P2',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['QH',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['QH',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['S',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['S',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['QG',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['QG',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['QP',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['QP',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['V',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['V',int2str(ID),',',int2str(IH)];
                [iv,ix_] = s2mpjlib('ii',['R',int2str(ID),',',int2str(IH)],ix_);
                pb.xnames{iv} = ['R',int2str(ID),',',int2str(IH)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
                gtype{ig} = '<>';
                iv = ix_(['P1',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1000.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1000.0;
                end
                iv = ix_(['P2',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1500.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1500.0;
                end
                iv = ix_(['QH',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1200.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1200.0;
                end
                iv = ix_(['S',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1200.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1200.0;
                end
                iv = ix_(['QG',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1200.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1200.0;
                end
                iv = ix_(['QP',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1200.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1200.0;
                end
            end
        end
        for ID=v_('1'):v_('DAYS')
            v_('P') = -1+ID;
            [ig,ig_] = s2mpjlib('ii',['H',int2str(ID),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['H',int2str(ID),',',int2str(round(v_('1')))];
            iv = ix_(['V',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['V',int2str(round(v_('P'))),',',int2str(round(v_('HOURS')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['H',int2str(ID),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['H',int2str(ID),',',int2str(round(v_('1')))];
            iv = ix_(['S',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['QH',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('2'):v_('HOURS')
                v_('IH-1') = -1+IH;
                [ig,ig_] = s2mpjlib('ii',['H',int2str(ID),',',int2str(IH)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['H',int2str(ID),',',int2str(IH)];
                iv = ix_(['V',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['V',int2str(ID),',',int2str(round(v_('IH-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['S',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['QH',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for ID=v_('1'):v_('DAYS')
            v_('P') = -1+ID;
            [ig,ig_] = s2mpjlib('ii',['R',int2str(ID),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(ID),',',int2str(round(v_('1')))];
            iv = ix_(['R',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['R',int2str(round(v_('P'))),',',int2str(round(v_('HOURS')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['R',int2str(ID),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(ID),',',int2str(round(v_('1')))];
            iv = ix_(['QG',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['QP',int2str(ID),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('2'):v_('HOURS')
                v_('IH-1') = -1+IH;
                [ig,ig_] = s2mpjlib('ii',['R',int2str(ID),',',int2str(IH)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['R',int2str(ID),',',int2str(IH)];
                iv = ix_(['R',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['R',int2str(ID),',',int2str(round(v_('IH-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['QG',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['QP',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                [ig,ig_] = s2mpjlib('ii',['D',int2str(ID),',',int2str(IH)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['D',int2str(ID),',',int2str(IH)];
                iv = ix_(['P1',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['P2',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['QH',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['QG',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['QP',int2str(ID),',',int2str(IH)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.33+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.33;
                end
            end
        end
        for D=v_('1'):v_('DAYS')
            for H=v_('1'):v_('HOURS')
                [ig,ig_] = s2mpjlib('ii',['QG*QP',int2str(D),',',int2str(H)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['QG*QP',int2str(D),',',int2str(H)];
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
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                pbm.gconst(ig_(['H',int2str(ID),',',int2str(IH)])) = v_(['Z',int2str(ID)]);
            end
        end
        for ID=v_('1'):v_('DAYS')
            v_('0.01Z') = 0.01*v_(['Z',int2str(ID)]);
            for IH=v_('1'):v_('HOURS')
                pbm.gconst(ig_(['R',int2str(ID),',',int2str(IH)])) = v_('0.01Z');
            end
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                pbm.gconst(ig_(['D',int2str(ID),',',int2str(IH)])) =...
                      v_(['D',int2str(ID),',',int2str(IH)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 240000.0;
        pb.xupper(ix_(['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 240000.0;
        pb.xlower(ix_(['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 3500.0;
        pb.xupper(ix_(['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 3500.0;
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                pb.xlower(ix_(['P1',int2str(ID),',',int2str(IH)]),1) = 70.0;
                pb.xlower(ix_(['P2',int2str(ID),',',int2str(IH)]),1) = 90.0;
                pb.xlower(ix_(['QH',int2str(ID),',',int2str(IH)]),1) = 25.0;
                pb.xlower(ix_(['V',int2str(ID),',',int2str(IH)]),1) = 180000.0;
                pb.xupper(ix_(['P1',int2str(ID),',',int2str(IH)])) = 325.0;
                pb.xupper(ix_(['P2',int2str(ID),',',int2str(IH)])) = 290.0;
                pb.xupper(ix_(['QH',int2str(ID),',',int2str(IH)])) = 500.0;
                pb.xupper(ix_(['QP',int2str(ID),',',int2str(IH)])) = 225.0;
                pb.xupper(ix_(['QG',int2str(ID),',',int2str(IH)])) = 300.0;
                pb.xupper(ix_(['V',int2str(ID),',',int2str(IH)])) = 280000.0;
                pb.xupper(ix_(['R',int2str(ID),',',int2str(IH)])) = 6000.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]))
            pb.x0(ix_(['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 240000.0;
        else
            pb.y0(find(pbm.congrps==ig_(['V',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))])),1) = 240000.0;
        end
        if(isKey(ix_,['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]))
            pb.x0(ix_(['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))]),1) = 3500.0;
        else
            pb.y0(find(pbm.congrps==ig_(['R',int2str(round(v_('0'))),',',int2str(round(v_('HOURS')))])),1) = 3500.0;
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                if(isKey(ix_,['P1',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['P1',int2str(ID),',',int2str(IH)]),1) = 70.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['P1',int2str(ID),',',int2str(IH)])),1) = 70.0;
                end
                if(isKey(ix_,['P2',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['P2',int2str(ID),',',int2str(IH)]),1) = 90.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['P2',int2str(ID),',',int2str(IH)])),1) = 90.0;
                end
                if(isKey(ix_,['QH',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['QH',int2str(ID),',',int2str(IH)]),1) = 25.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['QH',int2str(ID),',',int2str(IH)])),1) = 25.0;
                end
                if(isKey(ix_,['QP',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['QP',int2str(ID),',',int2str(IH)]),1) = 225.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['QP',int2str(ID),',',int2str(IH)])),1) =...
                          225.0;
                end
                if(isKey(ix_,['V',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['V',int2str(ID),',',int2str(IH)]),1) = 240000.0;
                else
                    pb.y0(find(pbm.congrps==ig_(['V',int2str(ID),',',int2str(IH)])),1) =...
                          240000.0;
                end
                if(isKey(ix_,['R',int2str(ID),',',int2str(IH)]))
                    pb.x0(ix_(['R',int2str(ID),',',int2str(IH)]),1) = 3500;
                else
                    pb.y0(find(pbm.congrps==ig_(['R',int2str(ID),',',int2str(IH)])),1) = 3500;
                end
            end
        end
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'QP';
        elftv{it}{2} = 'QG';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for ID=v_('1'):v_('DAYS')
            for IH=v_('1'):v_('HOURS')
                ename = ['P',int2str(ID),',',int2str(IH)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['QP',int2str(ID),',',int2str(IH)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('QP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['QG',int2str(ID),',',int2str(IH)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('QG',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for D=v_('1'):v_('DAYS')
            for H=v_('1'):v_('HOURS')
                ig = ig_(['QG*QP',int2str(D),',',int2str(H)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(D),',',int2str(H)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               1.617060D+07
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLQR2-RN-194-96';
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

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(2)*EV_(1);
        if(nargout>1)
            g_(2,1) = EV_(1);
            g_(1,1) = EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = 1.0e+0;
                H_(1,2) = H_(2,1);
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

