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
%    classification = 'C-CONR2-MN-31-10'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'WATER';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        irA  = [];
        icA  = [];
        valA = [];
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
        ngrp   = ig_.Count;
        [iv,ix_] = s2mpjlib('ii','Q0102',ix_);
        pb.xnames{iv} = 'Q0102';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0102');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c1');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0102',ix_);
        pb.xnames{iv} = 'Q0102';
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0203',ix_);
        pb.xnames{iv} = 'Q0203';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0203');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0203',ix_);
        pb.xnames{iv} = 'Q0203';
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0204',ix_);
        pb.xnames{iv} = 'Q0204';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0204');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0204',ix_);
        pb.xnames{iv} = 'Q0204';
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0305',ix_);
        pb.xnames{iv} = 'Q0305';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0305');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0305',ix_);
        pb.xnames{iv} = 'Q0305';
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0405',ix_);
        pb.xnames{iv} = 'Q0405';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0405');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0405',ix_);
        pb.xnames{iv} = 'Q0405';
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0406',ix_);
        pb.xnames{iv} = 'Q0406';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0406');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0406',ix_);
        pb.xnames{iv} = 'Q0406';
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0607',ix_);
        pb.xnames{iv} = 'Q0607';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0607');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0607',ix_);
        pb.xnames{iv} = 'Q0607';
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0705',ix_);
        pb.xnames{iv} = 'Q0705';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj0705');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0705',ix_);
        pb.xnames{iv} = 'Q0705';
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q01u0',ix_);
        pb.xnames{iv} = 'Q01u0';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c1');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q01u0',ix_);
        pb.xnames{iv} = 'Q01u0';
        icA(end+1) = iv;
        irA(end+1) = ig_('c8');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y02up',ix_);
        pb.xnames{iv} = 'y02up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y02up',ix_);
        pb.xnames{iv} = 'y02up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y03up',ix_);
        pb.xnames{iv} = 'y03up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y03up',ix_);
        pb.xnames{iv} = 'y03up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y04up',ix_);
        pb.xnames{iv} = 'y04up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y04up',ix_);
        pb.xnames{iv} = 'y04up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y05up',ix_);
        pb.xnames{iv} = 'y05up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y05up',ix_);
        pb.xnames{iv} = 'y05up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y06up',ix_);
        pb.xnames{iv} = 'y06up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y06up',ix_);
        pb.xnames{iv} = 'y06up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','y07up',ix_);
        pb.xnames{iv} = 'y07up';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = 210;
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','y07up',ix_);
        pb.xnames{iv} = 'y07up';
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu02',ix_);
        pb.xnames{iv} = 'yqu02';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -175;
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu02',ix_);
        pb.xnames{iv} = 'yqu02';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yqu03',ix_);
        pb.xnames{iv} = 'yqu03';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -190;
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu03',ix_);
        pb.xnames{iv} = 'yqu03';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yqu04',ix_);
        pb.xnames{iv} = 'yqu04';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -185;
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu04',ix_);
        pb.xnames{iv} = 'yqu04';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yqu05',ix_);
        pb.xnames{iv} = 'yqu05';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -180;
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu05',ix_);
        pb.xnames{iv} = 'yqu05';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yqu06',ix_);
        pb.xnames{iv} = 'yqu06';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -195;
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu06',ix_);
        pb.xnames{iv} = 'yqu06';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yqu07',ix_);
        pb.xnames{iv} = 'yqu07';
        icA(end+1) = iv;
        irA(end+1) = ig_('obj');
        valA(end+1) = -190;
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','yqu07',ix_);
        pb.xnames{iv} = 'yqu07';
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0201',ix_);
        pb.xnames{iv} = 'Q0201';
        icA(end+1) = iv;
        irA(end+1) = ig_('c1');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0302',ix_);
        pb.xnames{iv} = 'Q0302';
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0402',ix_);
        pb.xnames{iv} = 'Q0402';
        icA(end+1) = iv;
        irA(end+1) = ig_('c2');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0503',ix_);
        pb.xnames{iv} = 'Q0503';
        icA(end+1) = iv;
        irA(end+1) = ig_('c3');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0504',ix_);
        pb.xnames{iv} = 'Q0504';
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0604',ix_);
        pb.xnames{iv} = 'Q0604';
        icA(end+1) = iv;
        irA(end+1) = ig_('c4');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','Q0507',ix_);
        pb.xnames{iv} = 'Q0507';
        icA(end+1) = iv;
        irA(end+1) = ig_('c5');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = -1;
        [iv,ix_] = s2mpjlib('ii','Q0706',ix_);
        pb.xnames{iv} = 'Q0706';
        icA(end+1) = iv;
        irA(end+1) = ig_('c6');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c7');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yupu0',ix_);
        pb.xnames{iv} = 'yupu0';
        icA(end+1) = iv;
        irA(end+1) = ig_('c8');
        valA(end+1) = -1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c9');
        valA(end+1) = 1;
        [iv,ix_] = s2mpjlib('ii','yu0uq',ix_);
        pb.xnames{iv} = 'yu0uq';
        icA(end+1) = iv;
        irA(end+1) = ig_('c8');
        valA(end+1) = 1;
        icA(end+1) = iv;
        irA(end+1) = ig_('c10');
        valA(end+1) = -1;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
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
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gPOWER',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
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
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CONR2-MN-31-10';
        pb.x0          = zeros(pb.n,1);
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

