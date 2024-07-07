function varargout = WATER(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A small nonlinear network problem.
%    The problem is to compute the flows in a water distribution network
%    with 7 nodes and 8 links, subject to known supply/demand at the nodes 
%    and a unique reservoir at node 1.
% 
%    The problem is convex.
% 
%    Source:
%    an exercize for L. Watson course on LANCELOT in the Spring 1993.
% 
%    SIF input: E. P. Smith, Virginia Tech., Spring 1993.
% 
%    classification = 'ONR2-MN-31-10'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'WATER';

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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','obj0102',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 105665.6;
        [ig,ig_] = s2mpjlib('ii','obj0203',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 3613.412;
        [ig,ig_] = s2mpjlib('ii','obj0204',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 105665.6;
        [ig,ig_] = s2mpjlib('ii','obj0305',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 890.1553;
        [ig,ig_] = s2mpjlib('ii','obj0405',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 76.66088;
        [ig,ig_] = s2mpjlib('ii','obj0406',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 55145.82;
        [ig,ig_] = s2mpjlib('ii','obj0607',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 26030.46;
        [ig,ig_] = s2mpjlib('ii','obj0705',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = 890.1553;
        [ig,ig_] = s2mpjlib('ii','obj',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','c1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c1';
        [ig,ig_] = s2mpjlib('ii','c2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c2';
        [ig,ig_] = s2mpjlib('ii','c3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c3';
        [ig,ig_] = s2mpjlib('ii','c4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c4';
        [ig,ig_] = s2mpjlib('ii','c5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c5';
        [ig,ig_] = s2mpjlib('ii','c6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c6';
        [ig,ig_] = s2mpjlib('ii','c7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c7';
        [ig,ig_] = s2mpjlib('ii','c8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c8';
        [ig,ig_] = s2mpjlib('ii','c9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c9';
        [ig,ig_] = s2mpjlib('ii','c10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'c10';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = numEntries(ig_);
        [iv,ix_] = s2mpjlib('ii','Q0102',ix_);
        pb.xnames{iv} = 'Q0102';
        ig = ig_('obj0102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0102',ix_);
        pb.xnames{iv} = 'Q0102';
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0203',ix_);
        pb.xnames{iv} = 'Q0203';
        ig = ig_('obj0203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0203',ix_);
        pb.xnames{iv} = 'Q0203';
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0204',ix_);
        pb.xnames{iv} = 'Q0204';
        ig = ig_('obj0204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0204',ix_);
        pb.xnames{iv} = 'Q0204';
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0305',ix_);
        pb.xnames{iv} = 'Q0305';
        ig = ig_('obj0305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0305',ix_);
        pb.xnames{iv} = 'Q0305';
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0405',ix_);
        pb.xnames{iv} = 'Q0405';
        ig = ig_('obj0405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0405',ix_);
        pb.xnames{iv} = 'Q0405';
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0406',ix_);
        pb.xnames{iv} = 'Q0406';
        ig = ig_('obj0406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0406',ix_);
        pb.xnames{iv} = 'Q0406';
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0607',ix_);
        pb.xnames{iv} = 'Q0607';
        ig = ig_('obj0607');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0607',ix_);
        pb.xnames{iv} = 'Q0607';
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0705',ix_);
        pb.xnames{iv} = 'Q0705';
        ig = ig_('obj0705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0705',ix_);
        pb.xnames{iv} = 'Q0705';
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q01u0',ix_);
        pb.xnames{iv} = 'Q01u0';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q01u0',ix_);
        pb.xnames{iv} = 'Q01u0';
        ig = ig_('c8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y02up',ix_);
        pb.xnames{iv} = 'y02up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y02up',ix_);
        pb.xnames{iv} = 'y02up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y03up',ix_);
        pb.xnames{iv} = 'y03up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y03up',ix_);
        pb.xnames{iv} = 'y03up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y04up',ix_);
        pb.xnames{iv} = 'y04up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y04up',ix_);
        pb.xnames{iv} = 'y04up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y05up',ix_);
        pb.xnames{iv} = 'y05up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y05up',ix_);
        pb.xnames{iv} = 'y05up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y06up',ix_);
        pb.xnames{iv} = 'y06up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y06up',ix_);
        pb.xnames{iv} = 'y06up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','y07up',ix_);
        pb.xnames{iv} = 'y07up';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 210+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 210;
        end
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','y07up',ix_);
        pb.xnames{iv} = 'y07up';
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu02',ix_);
        pb.xnames{iv} = 'yqu02';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -175;
        end
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu02',ix_);
        pb.xnames{iv} = 'yqu02';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu03',ix_);
        pb.xnames{iv} = 'yqu03';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -190+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -190;
        end
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu03',ix_);
        pb.xnames{iv} = 'yqu03';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu04',ix_);
        pb.xnames{iv} = 'yqu04';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -185;
        end
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu04',ix_);
        pb.xnames{iv} = 'yqu04';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu05',ix_);
        pb.xnames{iv} = 'yqu05';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -180+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -180;
        end
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu05',ix_);
        pb.xnames{iv} = 'yqu05';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu06',ix_);
        pb.xnames{iv} = 'yqu06';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -195+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -195;
        end
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu06',ix_);
        pb.xnames{iv} = 'yqu06';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu07',ix_);
        pb.xnames{iv} = 'yqu07';
        ig = ig_('obj');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -190+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -190;
        end
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','yqu07',ix_);
        pb.xnames{iv} = 'yqu07';
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0201',ix_);
        pb.xnames{iv} = 'Q0201';
        ig = ig_('c1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0302',ix_);
        pb.xnames{iv} = 'Q0302';
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0402',ix_);
        pb.xnames{iv} = 'Q0402';
        ig = ig_('c2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0503',ix_);
        pb.xnames{iv} = 'Q0503';
        ig = ig_('c3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0504',ix_);
        pb.xnames{iv} = 'Q0504';
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0604',ix_);
        pb.xnames{iv} = 'Q0604';
        ig = ig_('c4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0507',ix_);
        pb.xnames{iv} = 'Q0507';
        ig = ig_('c5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [iv,ix_] = s2mpjlib('ii','Q0706',ix_);
        pb.xnames{iv} = 'Q0706';
        ig = ig_('c6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yupu0',ix_);
        pb.xnames{iv} = 'yupu0';
        ig = ig_('c8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        ig = ig_('c9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [iv,ix_] = s2mpjlib('ii','yu0uq',ix_);
        pb.xnames{iv} = 'yu0uq';
        ig = ig_('c8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('c10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
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
        pbm.gconst(ig_('c1')) = 1120;
        pbm.gconst(ig_('c2')) = -100;
        pbm.gconst(ig_('c3')) = -100;
        pbm.gconst(ig_('c4')) = -120;
        pbm.gconst(ig_('c5')) = -270;
        pbm.gconst(ig_('c6')) = -330;
        pbm.gconst(ig_('c7')) = -200;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('Q0102')) = 1200;
        pb.xupper(ix_('Q0203')) = 1200;
        pb.xupper(ix_('Q0204')) = 1200;
        pb.xupper(ix_('Q0305')) = 1200;
        pb.xupper(ix_('Q0405')) = 1200;
        pb.xupper(ix_('Q0406')) = 1200;
        pb.xupper(ix_('Q0607')) = 1200;
        pb.xupper(ix_('Q0705')) = 1200;
        pb.xupper(ix_('Q01u0')) = 1200;
        pb.xupper(ix_('y02up')) = 1200;
        pb.xupper(ix_('y03up')) = 1200;
        pb.xupper(ix_('y04up')) = 1200;
        pb.xupper(ix_('y05up')) = 1200;
        pb.xupper(ix_('y06up')) = 1200;
        pb.xupper(ix_('y07up')) = 1200;
        pb.xupper(ix_('yqu02')) = 1200;
        pb.xupper(ix_('yqu03')) = 1200;
        pb.xupper(ix_('yqu04')) = 1200;
        pb.xupper(ix_('yqu05')) = 1200;
        pb.xupper(ix_('yqu06')) = 1200;
        pb.xupper(ix_('yqu07')) = 1200;
        pb.xupper(ix_('Q0201')) = 1200;
        pb.xupper(ix_('Q0302')) = 1200;
        pb.xupper(ix_('Q0402')) = 1200;
        pb.xupper(ix_('Q0503')) = 1200;
        pb.xupper(ix_('Q0504')) = 1200;
        pb.xupper(ix_('Q0604')) = 1200;
        pb.xupper(ix_('Q0507')) = 1200;
        pb.xupper(ix_('Q0706')) = 1200;
        pb.xupper(ix_('yupu0')) = 1200;
        pb.xupper(ix_('yu0uq')) = 1200;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gPOWER',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('obj0102');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0203');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0204');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0305');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0405');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0406');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0607');
        pbm.grftype{ig} = 'gPOWER';
        ig = ig_('obj0705');
        pbm.grftype{ig} = 'gPOWER';
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION           1.054938D+04
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'ONR2-MN-31-10';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gPOWER'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^2.852;
        if(nargout>1)
            g_ = 2.852*GVAR_^1.852;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 5.282*GVAR_^.852;
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

