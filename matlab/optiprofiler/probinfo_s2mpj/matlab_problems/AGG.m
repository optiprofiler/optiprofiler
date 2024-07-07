function varargout = AGG(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    An LP with degeneracy
% 
%    Source:
%    The NETLIB collection of test problems.
% 
%    SIF input: (already in MPS format)
% 
%    classification = 'LLR2-AN-163-488'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'AGG';

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
        [ig,ig_] = s2mpjlib('ii','CAP00101',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00101';
        [ig,ig_] = s2mpjlib('ii','CAP00201',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00201';
        [ig,ig_] = s2mpjlib('ii','CAP00301',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00301';
        [ig,ig_] = s2mpjlib('ii','CAP00401',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00401';
        [ig,ig_] = s2mpjlib('ii','CAP00501',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00501';
        [ig,ig_] = s2mpjlib('ii','CAP00601',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00601';
        [ig,ig_] = s2mpjlib('ii','CAP00701',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00701';
        [ig,ig_] = s2mpjlib('ii','CAP00801',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00801';
        [ig,ig_] = s2mpjlib('ii','CAP00901',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00901';
        [ig,ig_] = s2mpjlib('ii','CAP01001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01001';
        [ig,ig_] = s2mpjlib('ii','CAP01101',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01101';
        [ig,ig_] = s2mpjlib('ii','CAP01201',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01201';
        [ig,ig_] = s2mpjlib('ii','CAP01301',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01301';
        [ig,ig_] = s2mpjlib('ii','CAP01401',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01401';
        [ig,ig_] = s2mpjlib('ii','CAP01501',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01501';
        [ig,ig_] = s2mpjlib('ii','CAP01601',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01601';
        [ig,ig_] = s2mpjlib('ii','CAP01701',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01701';
        [ig,ig_] = s2mpjlib('ii','CAP01801',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01801';
        [ig,ig_] = s2mpjlib('ii','CAP02001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02001';
        [ig,ig_] = s2mpjlib('ii','CAP02201',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02201';
        [ig,ig_] = s2mpjlib('ii','CAP02501',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02501';
        [ig,ig_] = s2mpjlib('ii','CAP02701',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02701';
        [ig,ig_] = s2mpjlib('ii','CAP03101',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03101';
        [ig,ig_] = s2mpjlib('ii','CAP03301',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03301';
        [ig,ig_] = s2mpjlib('ii','CAP03701',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03701';
        [ig,ig_] = s2mpjlib('ii','CAP04101',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04101';
        [ig,ig_] = s2mpjlib('ii','CAP04301',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04301';
        [ig,ig_] = s2mpjlib('ii','CAP04601',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04601';
        [ig,ig_] = s2mpjlib('ii','CAP05501',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05501';
        [ig,ig_] = s2mpjlib('ii','CAP06101',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06101';
        [ig,ig_] = s2mpjlib('ii','CAP06201',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06201';
        [ig,ig_] = s2mpjlib('ii','CAP06301',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06301';
        [ig,ig_] = s2mpjlib('ii','CAP06501',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06501';
        [ig,ig_] = s2mpjlib('ii','CAP00102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00102';
        [ig,ig_] = s2mpjlib('ii','CAP00202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00202';
        [ig,ig_] = s2mpjlib('ii','CAP00302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00302';
        [ig,ig_] = s2mpjlib('ii','CAP00402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00402';
        [ig,ig_] = s2mpjlib('ii','CAP00502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00502';
        [ig,ig_] = s2mpjlib('ii','CAP00602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00602';
        [ig,ig_] = s2mpjlib('ii','CAP00702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00702';
        [ig,ig_] = s2mpjlib('ii','CAP00802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00802';
        [ig,ig_] = s2mpjlib('ii','CAP00902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00902';
        [ig,ig_] = s2mpjlib('ii','CAP01002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01002';
        [ig,ig_] = s2mpjlib('ii','CAP01102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01102';
        [ig,ig_] = s2mpjlib('ii','CAP01202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01202';
        [ig,ig_] = s2mpjlib('ii','CAP01302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01302';
        [ig,ig_] = s2mpjlib('ii','CAP01402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01402';
        [ig,ig_] = s2mpjlib('ii','CAP01502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01502';
        [ig,ig_] = s2mpjlib('ii','CAP01602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01602';
        [ig,ig_] = s2mpjlib('ii','CAP01702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01702';
        [ig,ig_] = s2mpjlib('ii','CAP01802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01802';
        [ig,ig_] = s2mpjlib('ii','CAP01902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01902';
        [ig,ig_] = s2mpjlib('ii','CAP02002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02002';
        [ig,ig_] = s2mpjlib('ii','CAP02102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02102';
        [ig,ig_] = s2mpjlib('ii','CAP02202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02202';
        [ig,ig_] = s2mpjlib('ii','CAP02302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02302';
        [ig,ig_] = s2mpjlib('ii','CAP02402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02402';
        [ig,ig_] = s2mpjlib('ii','CAP02502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02502';
        [ig,ig_] = s2mpjlib('ii','CAP02602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02602';
        [ig,ig_] = s2mpjlib('ii','CAP02702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02702';
        [ig,ig_] = s2mpjlib('ii','CAP02802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02802';
        [ig,ig_] = s2mpjlib('ii','CAP02902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02902';
        [ig,ig_] = s2mpjlib('ii','CAP03002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03002';
        [ig,ig_] = s2mpjlib('ii','CAP03102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03102';
        [ig,ig_] = s2mpjlib('ii','CAP03202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03202';
        [ig,ig_] = s2mpjlib('ii','CAP03302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03302';
        [ig,ig_] = s2mpjlib('ii','CAP03402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03402';
        [ig,ig_] = s2mpjlib('ii','CAP03502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03502';
        [ig,ig_] = s2mpjlib('ii','CAP03602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03602';
        [ig,ig_] = s2mpjlib('ii','CAP03702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03702';
        [ig,ig_] = s2mpjlib('ii','CAP03802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03802';
        [ig,ig_] = s2mpjlib('ii','CAP03902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03902';
        [ig,ig_] = s2mpjlib('ii','CAP04002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04002';
        [ig,ig_] = s2mpjlib('ii','CAP04102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04102';
        [ig,ig_] = s2mpjlib('ii','CAP04202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04202';
        [ig,ig_] = s2mpjlib('ii','CAP04302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04302';
        [ig,ig_] = s2mpjlib('ii','CAP04402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04402';
        [ig,ig_] = s2mpjlib('ii','CAP04502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04502';
        [ig,ig_] = s2mpjlib('ii','CAP04602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04602';
        [ig,ig_] = s2mpjlib('ii','CAP04702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04702';
        [ig,ig_] = s2mpjlib('ii','CAP04802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04802';
        [ig,ig_] = s2mpjlib('ii','CAP04902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04902';
        [ig,ig_] = s2mpjlib('ii','CAP05002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05002';
        [ig,ig_] = s2mpjlib('ii','CAP05102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05102';
        [ig,ig_] = s2mpjlib('ii','CAP05202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05202';
        [ig,ig_] = s2mpjlib('ii','CAP05302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05302';
        [ig,ig_] = s2mpjlib('ii','CAP05402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05402';
        [ig,ig_] = s2mpjlib('ii','CAP05502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05502';
        [ig,ig_] = s2mpjlib('ii','CAP05602',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05602';
        [ig,ig_] = s2mpjlib('ii','CAP05702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05702';
        [ig,ig_] = s2mpjlib('ii','CAP05802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05802';
        [ig,ig_] = s2mpjlib('ii','CAP05902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05902';
        [ig,ig_] = s2mpjlib('ii','CAP06002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06002';
        [ig,ig_] = s2mpjlib('ii','CAP06102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06102';
        [ig,ig_] = s2mpjlib('ii','CAP06202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06202';
        [ig,ig_] = s2mpjlib('ii','CAP06302',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06302';
        [ig,ig_] = s2mpjlib('ii','CAP06402',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06402';
        [ig,ig_] = s2mpjlib('ii','CAP06502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06502';
        [ig,ig_] = s2mpjlib('ii','CAP00103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00103';
        [ig,ig_] = s2mpjlib('ii','CAP00203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00203';
        [ig,ig_] = s2mpjlib('ii','CAP00303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00303';
        [ig,ig_] = s2mpjlib('ii','CAP00403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00403';
        [ig,ig_] = s2mpjlib('ii','CAP00503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00503';
        [ig,ig_] = s2mpjlib('ii','CAP00603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00603';
        [ig,ig_] = s2mpjlib('ii','CAP00703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00703';
        [ig,ig_] = s2mpjlib('ii','CAP00803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00803';
        [ig,ig_] = s2mpjlib('ii','CAP00903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00903';
        [ig,ig_] = s2mpjlib('ii','CAP01003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01003';
        [ig,ig_] = s2mpjlib('ii','CAP01103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01103';
        [ig,ig_] = s2mpjlib('ii','CAP01203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01203';
        [ig,ig_] = s2mpjlib('ii','CAP01303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01303';
        [ig,ig_] = s2mpjlib('ii','CAP01403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01403';
        [ig,ig_] = s2mpjlib('ii','CAP01503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01503';
        [ig,ig_] = s2mpjlib('ii','CAP01603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01603';
        [ig,ig_] = s2mpjlib('ii','CAP01703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01703';
        [ig,ig_] = s2mpjlib('ii','CAP01803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01803';
        [ig,ig_] = s2mpjlib('ii','CAP01903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01903';
        [ig,ig_] = s2mpjlib('ii','CAP02003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02003';
        [ig,ig_] = s2mpjlib('ii','CAP02103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02103';
        [ig,ig_] = s2mpjlib('ii','CAP02203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02203';
        [ig,ig_] = s2mpjlib('ii','CAP02303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02303';
        [ig,ig_] = s2mpjlib('ii','CAP02403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02403';
        [ig,ig_] = s2mpjlib('ii','CAP02503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02503';
        [ig,ig_] = s2mpjlib('ii','CAP02603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02603';
        [ig,ig_] = s2mpjlib('ii','CAP02703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02703';
        [ig,ig_] = s2mpjlib('ii','CAP02803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02803';
        [ig,ig_] = s2mpjlib('ii','CAP02903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02903';
        [ig,ig_] = s2mpjlib('ii','CAP03003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03003';
        [ig,ig_] = s2mpjlib('ii','CAP03103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03103';
        [ig,ig_] = s2mpjlib('ii','CAP03203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03203';
        [ig,ig_] = s2mpjlib('ii','CAP03303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03303';
        [ig,ig_] = s2mpjlib('ii','CAP03403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03403';
        [ig,ig_] = s2mpjlib('ii','CAP03503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03503';
        [ig,ig_] = s2mpjlib('ii','CAP03603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03603';
        [ig,ig_] = s2mpjlib('ii','CAP03703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03703';
        [ig,ig_] = s2mpjlib('ii','CAP03803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03803';
        [ig,ig_] = s2mpjlib('ii','CAP03903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03903';
        [ig,ig_] = s2mpjlib('ii','CAP04003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04003';
        [ig,ig_] = s2mpjlib('ii','CAP04103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04103';
        [ig,ig_] = s2mpjlib('ii','CAP04203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04203';
        [ig,ig_] = s2mpjlib('ii','CAP04303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04303';
        [ig,ig_] = s2mpjlib('ii','CAP04403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04403';
        [ig,ig_] = s2mpjlib('ii','CAP04503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04503';
        [ig,ig_] = s2mpjlib('ii','CAP04603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04603';
        [ig,ig_] = s2mpjlib('ii','CAP04703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04703';
        [ig,ig_] = s2mpjlib('ii','CAP04803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04803';
        [ig,ig_] = s2mpjlib('ii','CAP04903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04903';
        [ig,ig_] = s2mpjlib('ii','CAP05003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05003';
        [ig,ig_] = s2mpjlib('ii','CAP05103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05103';
        [ig,ig_] = s2mpjlib('ii','CAP05203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05203';
        [ig,ig_] = s2mpjlib('ii','CAP05303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05303';
        [ig,ig_] = s2mpjlib('ii','CAP05403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05403';
        [ig,ig_] = s2mpjlib('ii','CAP05503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05503';
        [ig,ig_] = s2mpjlib('ii','CAP05603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05603';
        [ig,ig_] = s2mpjlib('ii','CAP05703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05703';
        [ig,ig_] = s2mpjlib('ii','CAP05803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05803';
        [ig,ig_] = s2mpjlib('ii','CAP05903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05903';
        [ig,ig_] = s2mpjlib('ii','CAP06003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06003';
        [ig,ig_] = s2mpjlib('ii','CAP06103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06103';
        [ig,ig_] = s2mpjlib('ii','CAP06203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06203';
        [ig,ig_] = s2mpjlib('ii','CAP06303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06303';
        [ig,ig_] = s2mpjlib('ii','CAP06403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06403';
        [ig,ig_] = s2mpjlib('ii','CAP06503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06503';
        [ig,ig_] = s2mpjlib('ii','CAP00104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00104';
        [ig,ig_] = s2mpjlib('ii','CAP00204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00204';
        [ig,ig_] = s2mpjlib('ii','CAP00304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00304';
        [ig,ig_] = s2mpjlib('ii','CAP00404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00404';
        [ig,ig_] = s2mpjlib('ii','CAP00504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00504';
        [ig,ig_] = s2mpjlib('ii','CAP00604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00604';
        [ig,ig_] = s2mpjlib('ii','CAP00704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00704';
        [ig,ig_] = s2mpjlib('ii','CAP00804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00804';
        [ig,ig_] = s2mpjlib('ii','CAP00904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00904';
        [ig,ig_] = s2mpjlib('ii','CAP01004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01004';
        [ig,ig_] = s2mpjlib('ii','CAP01104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01104';
        [ig,ig_] = s2mpjlib('ii','CAP01204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01204';
        [ig,ig_] = s2mpjlib('ii','CAP01304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01304';
        [ig,ig_] = s2mpjlib('ii','CAP01404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01404';
        [ig,ig_] = s2mpjlib('ii','CAP01504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01504';
        [ig,ig_] = s2mpjlib('ii','CAP01604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01604';
        [ig,ig_] = s2mpjlib('ii','CAP01704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01704';
        [ig,ig_] = s2mpjlib('ii','CAP01804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01804';
        [ig,ig_] = s2mpjlib('ii','CAP01904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01904';
        [ig,ig_] = s2mpjlib('ii','CAP02004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02004';
        [ig,ig_] = s2mpjlib('ii','CAP02104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02104';
        [ig,ig_] = s2mpjlib('ii','CAP02204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02204';
        [ig,ig_] = s2mpjlib('ii','CAP02304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02304';
        [ig,ig_] = s2mpjlib('ii','CAP02404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02404';
        [ig,ig_] = s2mpjlib('ii','CAP02504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02504';
        [ig,ig_] = s2mpjlib('ii','CAP02604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02604';
        [ig,ig_] = s2mpjlib('ii','CAP02704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02704';
        [ig,ig_] = s2mpjlib('ii','CAP02804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02804';
        [ig,ig_] = s2mpjlib('ii','CAP02904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02904';
        [ig,ig_] = s2mpjlib('ii','CAP03004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03004';
        [ig,ig_] = s2mpjlib('ii','CAP03104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03104';
        [ig,ig_] = s2mpjlib('ii','CAP03204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03204';
        [ig,ig_] = s2mpjlib('ii','CAP03304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03304';
        [ig,ig_] = s2mpjlib('ii','CAP03404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03404';
        [ig,ig_] = s2mpjlib('ii','CAP03504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03504';
        [ig,ig_] = s2mpjlib('ii','CAP03604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03604';
        [ig,ig_] = s2mpjlib('ii','CAP03704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03704';
        [ig,ig_] = s2mpjlib('ii','CAP03804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03804';
        [ig,ig_] = s2mpjlib('ii','CAP03904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03904';
        [ig,ig_] = s2mpjlib('ii','CAP04004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04004';
        [ig,ig_] = s2mpjlib('ii','CAP04104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04104';
        [ig,ig_] = s2mpjlib('ii','CAP04204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04204';
        [ig,ig_] = s2mpjlib('ii','CAP04304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04304';
        [ig,ig_] = s2mpjlib('ii','CAP04404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04404';
        [ig,ig_] = s2mpjlib('ii','CAP04504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04504';
        [ig,ig_] = s2mpjlib('ii','CAP04604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04604';
        [ig,ig_] = s2mpjlib('ii','CAP04704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04704';
        [ig,ig_] = s2mpjlib('ii','CAP04804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04804';
        [ig,ig_] = s2mpjlib('ii','CAP04904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04904';
        [ig,ig_] = s2mpjlib('ii','CAP05004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05004';
        [ig,ig_] = s2mpjlib('ii','CAP05104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05104';
        [ig,ig_] = s2mpjlib('ii','CAP05204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05204';
        [ig,ig_] = s2mpjlib('ii','CAP05304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05304';
        [ig,ig_] = s2mpjlib('ii','CAP05404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05404';
        [ig,ig_] = s2mpjlib('ii','CAP05504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05504';
        [ig,ig_] = s2mpjlib('ii','CAP05604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05604';
        [ig,ig_] = s2mpjlib('ii','CAP05704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05704';
        [ig,ig_] = s2mpjlib('ii','CAP05804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05804';
        [ig,ig_] = s2mpjlib('ii','CAP05904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05904';
        [ig,ig_] = s2mpjlib('ii','CAP06004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06004';
        [ig,ig_] = s2mpjlib('ii','CAP06104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06104';
        [ig,ig_] = s2mpjlib('ii','CAP06204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06204';
        [ig,ig_] = s2mpjlib('ii','CAP06304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06304';
        [ig,ig_] = s2mpjlib('ii','CAP06404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06404';
        [ig,ig_] = s2mpjlib('ii','CAP06504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06504';
        [ig,ig_] = s2mpjlib('ii','CAP00105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00105';
        [ig,ig_] = s2mpjlib('ii','CAP00205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00205';
        [ig,ig_] = s2mpjlib('ii','CAP00305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00305';
        [ig,ig_] = s2mpjlib('ii','CAP00405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00405';
        [ig,ig_] = s2mpjlib('ii','CAP00505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00505';
        [ig,ig_] = s2mpjlib('ii','CAP00605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00605';
        [ig,ig_] = s2mpjlib('ii','CAP00705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00705';
        [ig,ig_] = s2mpjlib('ii','CAP00805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00805';
        [ig,ig_] = s2mpjlib('ii','CAP00905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00905';
        [ig,ig_] = s2mpjlib('ii','CAP01005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01005';
        [ig,ig_] = s2mpjlib('ii','CAP01105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01105';
        [ig,ig_] = s2mpjlib('ii','CAP01205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01205';
        [ig,ig_] = s2mpjlib('ii','CAP01305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01305';
        [ig,ig_] = s2mpjlib('ii','CAP01405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01405';
        [ig,ig_] = s2mpjlib('ii','CAP01505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01505';
        [ig,ig_] = s2mpjlib('ii','CAP01605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01605';
        [ig,ig_] = s2mpjlib('ii','CAP01705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01705';
        [ig,ig_] = s2mpjlib('ii','CAP01805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01805';
        [ig,ig_] = s2mpjlib('ii','CAP01905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01905';
        [ig,ig_] = s2mpjlib('ii','CAP02005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02005';
        [ig,ig_] = s2mpjlib('ii','CAP02105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02105';
        [ig,ig_] = s2mpjlib('ii','CAP02205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02205';
        [ig,ig_] = s2mpjlib('ii','CAP02305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02305';
        [ig,ig_] = s2mpjlib('ii','CAP02405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02405';
        [ig,ig_] = s2mpjlib('ii','CAP02505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02505';
        [ig,ig_] = s2mpjlib('ii','CAP02605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02605';
        [ig,ig_] = s2mpjlib('ii','CAP02705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02705';
        [ig,ig_] = s2mpjlib('ii','CAP02805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02805';
        [ig,ig_] = s2mpjlib('ii','CAP02905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02905';
        [ig,ig_] = s2mpjlib('ii','CAP03005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03005';
        [ig,ig_] = s2mpjlib('ii','CAP03105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03105';
        [ig,ig_] = s2mpjlib('ii','CAP03205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03205';
        [ig,ig_] = s2mpjlib('ii','CAP03305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03305';
        [ig,ig_] = s2mpjlib('ii','CAP03405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03405';
        [ig,ig_] = s2mpjlib('ii','CAP03505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03505';
        [ig,ig_] = s2mpjlib('ii','CAP03605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03605';
        [ig,ig_] = s2mpjlib('ii','CAP03705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03705';
        [ig,ig_] = s2mpjlib('ii','CAP03805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03805';
        [ig,ig_] = s2mpjlib('ii','CAP03905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03905';
        [ig,ig_] = s2mpjlib('ii','CAP04005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04005';
        [ig,ig_] = s2mpjlib('ii','CAP04105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04105';
        [ig,ig_] = s2mpjlib('ii','CAP04205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04205';
        [ig,ig_] = s2mpjlib('ii','CAP04305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04305';
        [ig,ig_] = s2mpjlib('ii','CAP04405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04405';
        [ig,ig_] = s2mpjlib('ii','CAP04505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04505';
        [ig,ig_] = s2mpjlib('ii','CAP04605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04605';
        [ig,ig_] = s2mpjlib('ii','CAP04705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04705';
        [ig,ig_] = s2mpjlib('ii','CAP04805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04805';
        [ig,ig_] = s2mpjlib('ii','CAP04905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04905';
        [ig,ig_] = s2mpjlib('ii','CAP05005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05005';
        [ig,ig_] = s2mpjlib('ii','CAP05105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05105';
        [ig,ig_] = s2mpjlib('ii','CAP05205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05205';
        [ig,ig_] = s2mpjlib('ii','CAP05305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05305';
        [ig,ig_] = s2mpjlib('ii','CAP05405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05405';
        [ig,ig_] = s2mpjlib('ii','CAP05505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05505';
        [ig,ig_] = s2mpjlib('ii','CAP05605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05605';
        [ig,ig_] = s2mpjlib('ii','CAP05705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05705';
        [ig,ig_] = s2mpjlib('ii','CAP05805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05805';
        [ig,ig_] = s2mpjlib('ii','CAP05905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05905';
        [ig,ig_] = s2mpjlib('ii','CAP06005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06005';
        [ig,ig_] = s2mpjlib('ii','CAP06105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06105';
        [ig,ig_] = s2mpjlib('ii','CAP06205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06205';
        [ig,ig_] = s2mpjlib('ii','CAP06305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06305';
        [ig,ig_] = s2mpjlib('ii','CAP06405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06405';
        [ig,ig_] = s2mpjlib('ii','CAP06505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06505';
        [ig,ig_] = s2mpjlib('ii','CAP00106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00106';
        [ig,ig_] = s2mpjlib('ii','CAP00206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00206';
        [ig,ig_] = s2mpjlib('ii','CAP00306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00306';
        [ig,ig_] = s2mpjlib('ii','CAP00406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00406';
        [ig,ig_] = s2mpjlib('ii','CAP00506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00506';
        [ig,ig_] = s2mpjlib('ii','CAP00606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00606';
        [ig,ig_] = s2mpjlib('ii','CAP00706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00706';
        [ig,ig_] = s2mpjlib('ii','CAP00806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00806';
        [ig,ig_] = s2mpjlib('ii','CAP00906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP00906';
        [ig,ig_] = s2mpjlib('ii','CAP01006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01006';
        [ig,ig_] = s2mpjlib('ii','CAP01106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01106';
        [ig,ig_] = s2mpjlib('ii','CAP01206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01206';
        [ig,ig_] = s2mpjlib('ii','CAP01306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01306';
        [ig,ig_] = s2mpjlib('ii','CAP01406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01406';
        [ig,ig_] = s2mpjlib('ii','CAP01506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01506';
        [ig,ig_] = s2mpjlib('ii','CAP01606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01606';
        [ig,ig_] = s2mpjlib('ii','CAP01706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01706';
        [ig,ig_] = s2mpjlib('ii','CAP01806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01806';
        [ig,ig_] = s2mpjlib('ii','CAP01906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP01906';
        [ig,ig_] = s2mpjlib('ii','CAP02006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02006';
        [ig,ig_] = s2mpjlib('ii','CAP02106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02106';
        [ig,ig_] = s2mpjlib('ii','CAP02206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02206';
        [ig,ig_] = s2mpjlib('ii','CAP02306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02306';
        [ig,ig_] = s2mpjlib('ii','CAP02406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02406';
        [ig,ig_] = s2mpjlib('ii','CAP02506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02506';
        [ig,ig_] = s2mpjlib('ii','CAP02606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02606';
        [ig,ig_] = s2mpjlib('ii','CAP02706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02706';
        [ig,ig_] = s2mpjlib('ii','CAP02806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02806';
        [ig,ig_] = s2mpjlib('ii','CAP02906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP02906';
        [ig,ig_] = s2mpjlib('ii','CAP03006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03006';
        [ig,ig_] = s2mpjlib('ii','CAP03106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03106';
        [ig,ig_] = s2mpjlib('ii','CAP03206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03206';
        [ig,ig_] = s2mpjlib('ii','CAP03306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03306';
        [ig,ig_] = s2mpjlib('ii','CAP03406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03406';
        [ig,ig_] = s2mpjlib('ii','CAP03506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03506';
        [ig,ig_] = s2mpjlib('ii','CAP03606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03606';
        [ig,ig_] = s2mpjlib('ii','CAP03706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03706';
        [ig,ig_] = s2mpjlib('ii','CAP03806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03806';
        [ig,ig_] = s2mpjlib('ii','CAP03906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP03906';
        [ig,ig_] = s2mpjlib('ii','CAP04006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04006';
        [ig,ig_] = s2mpjlib('ii','CAP04106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04106';
        [ig,ig_] = s2mpjlib('ii','CAP04206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04206';
        [ig,ig_] = s2mpjlib('ii','CAP04306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04306';
        [ig,ig_] = s2mpjlib('ii','CAP04406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04406';
        [ig,ig_] = s2mpjlib('ii','CAP04506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04506';
        [ig,ig_] = s2mpjlib('ii','CAP04606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04606';
        [ig,ig_] = s2mpjlib('ii','CAP04706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04706';
        [ig,ig_] = s2mpjlib('ii','CAP04806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04806';
        [ig,ig_] = s2mpjlib('ii','CAP04906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP04906';
        [ig,ig_] = s2mpjlib('ii','CAP05006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05006';
        [ig,ig_] = s2mpjlib('ii','CAP05106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05106';
        [ig,ig_] = s2mpjlib('ii','CAP05206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05206';
        [ig,ig_] = s2mpjlib('ii','CAP05306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05306';
        [ig,ig_] = s2mpjlib('ii','CAP05406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05406';
        [ig,ig_] = s2mpjlib('ii','CAP05506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05506';
        [ig,ig_] = s2mpjlib('ii','CAP05606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05606';
        [ig,ig_] = s2mpjlib('ii','CAP05706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05706';
        [ig,ig_] = s2mpjlib('ii','CAP05806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05806';
        [ig,ig_] = s2mpjlib('ii','CAP05906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP05906';
        [ig,ig_] = s2mpjlib('ii','CAP06006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06006';
        [ig,ig_] = s2mpjlib('ii','CAP06106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06106';
        [ig,ig_] = s2mpjlib('ii','CAP06206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06206';
        [ig,ig_] = s2mpjlib('ii','CAP06306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06306';
        [ig,ig_] = s2mpjlib('ii','CAP06406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06406';
        [ig,ig_] = s2mpjlib('ii','CAP06506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CAP06506';
        [ig,ig_] = s2mpjlib('ii','MND00102',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00102';
        [ig,ig_] = s2mpjlib('ii','MXD00102',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00102';
        [ig,ig_] = s2mpjlib('ii','MND00202',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00202';
        [ig,ig_] = s2mpjlib('ii','MXD00202',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00202';
        [ig,ig_] = s2mpjlib('ii','MND00502',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00502';
        [ig,ig_] = s2mpjlib('ii','MXD00502',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00502';
        [ig,ig_] = s2mpjlib('ii','MND00702',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00702';
        [ig,ig_] = s2mpjlib('ii','MXD00702',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00702';
        [ig,ig_] = s2mpjlib('ii','MND00802',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00802';
        [ig,ig_] = s2mpjlib('ii','MXD00802',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00802';
        [ig,ig_] = s2mpjlib('ii','MND00902',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00902';
        [ig,ig_] = s2mpjlib('ii','MXD00902',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00902';
        [ig,ig_] = s2mpjlib('ii','MND01002',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND01002';
        [ig,ig_] = s2mpjlib('ii','MXD01002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD01002';
        [ig,ig_] = s2mpjlib('ii','MND00103',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00103';
        [ig,ig_] = s2mpjlib('ii','MXD00103',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00103';
        [ig,ig_] = s2mpjlib('ii','MND00203',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00203';
        [ig,ig_] = s2mpjlib('ii','MXD00203',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00203';
        [ig,ig_] = s2mpjlib('ii','MND00303',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00303';
        [ig,ig_] = s2mpjlib('ii','MXD00303',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00303';
        [ig,ig_] = s2mpjlib('ii','MND00403',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00403';
        [ig,ig_] = s2mpjlib('ii','MXD00403',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00403';
        [ig,ig_] = s2mpjlib('ii','MND00503',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00503';
        [ig,ig_] = s2mpjlib('ii','MXD00503',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00503';
        [ig,ig_] = s2mpjlib('ii','MND00603',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00603';
        [ig,ig_] = s2mpjlib('ii','MXD00603',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00603';
        [ig,ig_] = s2mpjlib('ii','MND00703',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00703';
        [ig,ig_] = s2mpjlib('ii','MXD00703',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00703';
        [ig,ig_] = s2mpjlib('ii','MND00803',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00803';
        [ig,ig_] = s2mpjlib('ii','MXD00803',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00803';
        [ig,ig_] = s2mpjlib('ii','MND00903',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00903';
        [ig,ig_] = s2mpjlib('ii','MXD00903',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00903';
        [ig,ig_] = s2mpjlib('ii','MND01003',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND01003';
        [ig,ig_] = s2mpjlib('ii','MXD01003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD01003';
        [ig,ig_] = s2mpjlib('ii','MND00104',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00104';
        [ig,ig_] = s2mpjlib('ii','MXD00104',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00104';
        [ig,ig_] = s2mpjlib('ii','MND00204',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00204';
        [ig,ig_] = s2mpjlib('ii','MXD00204',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00204';
        [ig,ig_] = s2mpjlib('ii','MND00304',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00304';
        [ig,ig_] = s2mpjlib('ii','MXD00304',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00304';
        [ig,ig_] = s2mpjlib('ii','MND00404',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00404';
        [ig,ig_] = s2mpjlib('ii','MXD00404',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00404';
        [ig,ig_] = s2mpjlib('ii','MND00504',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00504';
        [ig,ig_] = s2mpjlib('ii','MXD00504',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00504';
        [ig,ig_] = s2mpjlib('ii','MND00604',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00604';
        [ig,ig_] = s2mpjlib('ii','MXD00604',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00604';
        [ig,ig_] = s2mpjlib('ii','MND00704',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00704';
        [ig,ig_] = s2mpjlib('ii','MXD00704',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00704';
        [ig,ig_] = s2mpjlib('ii','MND00804',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00804';
        [ig,ig_] = s2mpjlib('ii','MXD00804',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00804';
        [ig,ig_] = s2mpjlib('ii','MND00904',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00904';
        [ig,ig_] = s2mpjlib('ii','MXD00904',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00904';
        [ig,ig_] = s2mpjlib('ii','MND01004',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND01004';
        [ig,ig_] = s2mpjlib('ii','MXD01004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD01004';
        [ig,ig_] = s2mpjlib('ii','MND00105',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00105';
        [ig,ig_] = s2mpjlib('ii','MXD00105',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00105';
        [ig,ig_] = s2mpjlib('ii','MND00205',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00205';
        [ig,ig_] = s2mpjlib('ii','MXD00205',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00205';
        [ig,ig_] = s2mpjlib('ii','MND00305',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00305';
        [ig,ig_] = s2mpjlib('ii','MXD00305',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00305';
        [ig,ig_] = s2mpjlib('ii','MND00405',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00405';
        [ig,ig_] = s2mpjlib('ii','MXD00405',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00405';
        [ig,ig_] = s2mpjlib('ii','MND00505',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00505';
        [ig,ig_] = s2mpjlib('ii','MXD00505',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00505';
        [ig,ig_] = s2mpjlib('ii','MND00605',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00605';
        [ig,ig_] = s2mpjlib('ii','MXD00605',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00605';
        [ig,ig_] = s2mpjlib('ii','MND00705',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00705';
        [ig,ig_] = s2mpjlib('ii','MXD00705',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00705';
        [ig,ig_] = s2mpjlib('ii','MND00805',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00805';
        [ig,ig_] = s2mpjlib('ii','MXD00805',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00805';
        [ig,ig_] = s2mpjlib('ii','MND00905',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00905';
        [ig,ig_] = s2mpjlib('ii','MXD00905',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00905';
        [ig,ig_] = s2mpjlib('ii','MND01005',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND01005';
        [ig,ig_] = s2mpjlib('ii','MXD01005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD01005';
        [ig,ig_] = s2mpjlib('ii','MND00106',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00106';
        [ig,ig_] = s2mpjlib('ii','MXD00106',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00106';
        [ig,ig_] = s2mpjlib('ii','MND00206',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00206';
        [ig,ig_] = s2mpjlib('ii','MXD00206',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00206';
        [ig,ig_] = s2mpjlib('ii','MND00306',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00306';
        [ig,ig_] = s2mpjlib('ii','MXD00306',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00306';
        [ig,ig_] = s2mpjlib('ii','MND00406',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00406';
        [ig,ig_] = s2mpjlib('ii','MXD00406',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00406';
        [ig,ig_] = s2mpjlib('ii','MND00506',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00506';
        [ig,ig_] = s2mpjlib('ii','MXD00506',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00506';
        [ig,ig_] = s2mpjlib('ii','MND00606',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00606';
        [ig,ig_] = s2mpjlib('ii','MXD00606',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00606';
        [ig,ig_] = s2mpjlib('ii','MND00706',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00706';
        [ig,ig_] = s2mpjlib('ii','MXD00706',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00706';
        [ig,ig_] = s2mpjlib('ii','MND00806',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00806';
        [ig,ig_] = s2mpjlib('ii','MXD00806',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00806';
        [ig,ig_] = s2mpjlib('ii','MND00906',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND00906';
        [ig,ig_] = s2mpjlib('ii','MXD00906',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD00906';
        [ig,ig_] = s2mpjlib('ii','MND01006',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'MND01006';
        [ig,ig_] = s2mpjlib('ii','MXD01006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MXD01006';
        [ig,ig_] = s2mpjlib('ii','INV00101',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00101';
        [ig,ig_] = s2mpjlib('ii','INV00201',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00201';
        [ig,ig_] = s2mpjlib('ii','INV00301',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00301';
        [ig,ig_] = s2mpjlib('ii','INV00401',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00401';
        [ig,ig_] = s2mpjlib('ii','INV00501',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00501';
        [ig,ig_] = s2mpjlib('ii','INV00601',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00601';
        [ig,ig_] = s2mpjlib('ii','INV00102',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00102';
        [ig,ig_] = s2mpjlib('ii','INV00202',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00202';
        [ig,ig_] = s2mpjlib('ii','INV00302',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00302';
        [ig,ig_] = s2mpjlib('ii','INV00402',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00402';
        [ig,ig_] = s2mpjlib('ii','INV00502',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00502';
        [ig,ig_] = s2mpjlib('ii','INV00602',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00602';
        [ig,ig_] = s2mpjlib('ii','INV00103',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00103';
        [ig,ig_] = s2mpjlib('ii','INV00203',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00203';
        [ig,ig_] = s2mpjlib('ii','INV00303',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00303';
        [ig,ig_] = s2mpjlib('ii','INV00403',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00403';
        [ig,ig_] = s2mpjlib('ii','INV00503',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00503';
        [ig,ig_] = s2mpjlib('ii','INV00603',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00603';
        [ig,ig_] = s2mpjlib('ii','INV00104',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00104';
        [ig,ig_] = s2mpjlib('ii','INV00204',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00204';
        [ig,ig_] = s2mpjlib('ii','INV00304',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00304';
        [ig,ig_] = s2mpjlib('ii','INV00404',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00404';
        [ig,ig_] = s2mpjlib('ii','INV00504',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00504';
        [ig,ig_] = s2mpjlib('ii','INV00604',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00604';
        [ig,ig_] = s2mpjlib('ii','INV00105',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00105';
        [ig,ig_] = s2mpjlib('ii','INV00205',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00205';
        [ig,ig_] = s2mpjlib('ii','INV00305',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00305';
        [ig,ig_] = s2mpjlib('ii','INV00405',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00405';
        [ig,ig_] = s2mpjlib('ii','INV00505',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00505';
        [ig,ig_] = s2mpjlib('ii','INV00605',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00605';
        [ig,ig_] = s2mpjlib('ii','INV00106',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00106';
        [ig,ig_] = s2mpjlib('ii','INV00206',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00206';
        [ig,ig_] = s2mpjlib('ii','INV00306',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00306';
        [ig,ig_] = s2mpjlib('ii','INV00406',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00406';
        [ig,ig_] = s2mpjlib('ii','INV00506',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00506';
        [ig,ig_] = s2mpjlib('ii','INV00606',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INV00606';
        [ig,ig_] = s2mpjlib('ii','OBJECTIV',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = numEntries(ig_);
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01773+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01773;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('CAP01501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01775+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01775;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00132;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MND00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.76829+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.76829;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.1405+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.1405;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MXD00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01248+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01248;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00141+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00141;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MXD00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -31.09+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -31.09;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.05277+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.05277;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -29.52+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -29.52;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01385+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01385;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00143;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.85602+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.85602;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01636+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01636;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MXD00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01438+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01438;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MND00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('CAP01502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01639+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01639;
        end
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01385+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01385;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01438+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01438;
        end
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01636+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01636;
        end
        ig = ig_('CAP01503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01639+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01639;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.03;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00143;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.85602+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.85602;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.05277+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.05277;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('CAP01504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01578+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01578;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -26.61+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -26.61;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01576;
        end
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89501+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89501;
        end
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.01378+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.01378;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00147+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00147;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00125;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01498;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('CAP01505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01446+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01446;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('CAP01506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03024;
        end
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00272+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00272;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.26+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.26;
        end
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.90879+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.90879;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.09488+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.09488;
        end
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01702+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01702;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('CAP01505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01704+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01704;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('CAP01406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03074;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00135+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00135;
        end
        [iv,ix_] = s2mpjlib('ii','Y00106',ix_);
        pb.xnames{iv} = 'Y00106';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00132;
        end
        ig = ig_('CAP00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02059+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02059;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -31.15+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -31.15;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01176+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01176;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00141+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00141;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.52969+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.52969;
        end
        ig = ig_('CAP00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00965+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00965;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.3791+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.3791;
        end
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('MXD00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('CAP00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02083;
        end
        ig = ig_('MXD00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('MND00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('MXD00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00143;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01336+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01336;
        end
        ig = ig_('CAP00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01123+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01123;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .019+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .019;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.27302+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.27302;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.63578+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.63578;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01923+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01923;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -29.57+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -29.57;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('MND00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01923+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01923;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.08;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00143;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.27302+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.27302;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.63578+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.63578;
        end
        ig = ig_('CAP00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01123+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01123;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01336+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01336;
        end
        ig = ig_('CAP00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .019+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .019;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('MXD00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.22587+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.22587;
        end
        ig = ig_('MXD00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('CAP00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01194+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01194;
        end
        ig = ig_('CAP00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01408+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01408;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00125;
        end
        ig = ig_('CAP00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0183;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.68292+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.68292;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00147+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00147;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01852+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01852;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -26.66+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -26.66;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('MND00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('MND00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.31+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.31;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.32394+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.32394;
        end
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.90879+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.90879;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00272+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00272;
        end
        ig = ig_('CAP00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('CAP00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03024;
        end
        ig = ig_('CAP00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01976+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01976;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('CAP00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03259+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03259;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00135+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00135;
        end
        [iv,ix_] = s2mpjlib('ii','Y00206',ix_);
        pb.xnames{iv} = 'Y00206';
        ig = ig_('MXD00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('CAP00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01455+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01455;
        end
        ig = ig_('MXD00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.46;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('CAP00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01426+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01426;
        end
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06023+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06023;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('MND00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00165+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00165;
        end
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.67787+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.67787;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        ig = ig_('CAP00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04611+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04611;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.76124+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.76124;
        end
        ig = ig_('MND00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('MXD00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06135+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06135;
        end
        ig = ig_('CAP00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0472+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0472;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -30.04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -30.04;
        end
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00166+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00166;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('MXD00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.81339+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.81339;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.62573+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.62573;
        end
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01343+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01343;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('CAP00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01316+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01316;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01343+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01343;
        end
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00166+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00166;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.62573+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.62573;
        end
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.81339+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.81339;
        end
        ig = ig_('CAP00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0472+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0472;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.52+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.52;
        end
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06135+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06135;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        ig = ig_('CAP00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01316+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01316;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00167;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.83656+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.83656;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.60255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.60255;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06184+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06184;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('CAP00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01268+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01268;
        end
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -27.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -27.08;
        end
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01293+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01293;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00013;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('CAP00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04769+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04769;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0018+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0018;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.43911+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.43911;
        end
        ig = ig_('CAP00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06037;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('CAP00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .07478+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .07478;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01397+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01397;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.71+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.71;
        end
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01369+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01369;
        end
        [iv,ix_] = s2mpjlib('ii','Y00306',ix_);
        pb.xnames{iv} = 'Y00306';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.65075+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.65075;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00165+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00165;
        end
        ig = ig_('CAP01701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00928+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00928;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.42+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.42;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('CAP01602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05177+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05177;
        end
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('CAP01601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0086+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0086;
        end
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.37298+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.37298;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('MXD00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.06613+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.06613;
        end
        ig = ig_('MND00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('CAP01702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .061+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .061;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.09482+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.09482;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00166+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00166;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('CAP01702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00857+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00857;
        end
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -30.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -30.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('CAP01603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05243;
        end
        ig = ig_('CAP01602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00794+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00794;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('CAP01703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06171;
        end
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.34429+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.34429;
        end
        ig = ig_('MND00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('MXD00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('MXD00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00166+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00166;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.34429+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.34429;
        end
        ig = ig_('MND00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('CAP01604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05243;
        end
        ig = ig_('CAP01704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06171;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.09482+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.09482;
        end
        ig = ig_('CAP01603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00794+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00794;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('CAP01703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00857+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00857;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.48+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.48;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('MND00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00765+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00765;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00825+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00825;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -27.04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -27.04;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00013;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.10758+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.10758;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('CAP01605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05272+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05272;
        end
        ig = ig_('CAP01705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06203+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06203;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00167;
        end
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33154+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33154;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('MXD00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('CAP01706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .07028+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .07028;
        end
        ig = ig_('CAP01605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00826;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.35806+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.35806;
        end
        ig = ig_('MND00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('MXD00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00014;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.43911+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.43911;
        end
        ig = ig_('CAP01606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06037;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('CAP01705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00891+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00891;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.67;
        end
        [iv,ix_] = s2mpjlib('ii','Y00406',ix_);
        pb.xnames{iv} = 'Y00406';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0018+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0018;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('CAP00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00401;
        end
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.216;
        end
        ig = ig_('MND00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('CAP06002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00225;
        end
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00124+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00124;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03973+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03973;
        end
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -33.21+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -33.21;
        end
        ig = ig_('MXD00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('CAP00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00331+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00331;
        end
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.18801+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.18801;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('MXD00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03074;
        end
        [iv,ix_] = s2mpjlib('ii','Y00503',ix_);
        pb.xnames{iv} = 'Y00503';
        ig = ig_('MND00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00124+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00124;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('MND00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03341+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03341;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('CAP00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00134+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00134;
        end
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00225;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.29601+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.29601;
        end
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.108+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.108;
        end
        [iv,ix_] = s2mpjlib('ii','Y00504',ix_);
        pb.xnames{iv} = 'Y00504';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.74+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.74;
        end
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04304;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00133+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00133;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -27.28+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -27.28;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.104;
        end
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00216;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00159+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00159;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('CAP00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00129+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00129;
        end
        ig = ig_('CAP00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03346+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03346;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04145+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04145;
        end
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.30001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.30001;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00505',ix_);
        pb.xnames{iv} = 'Y00505';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.51633+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.51633;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03753+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03753;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.40401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.40401;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.9;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00234+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00234;
        end
        ig = ig_('CAP06006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00349+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00349;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.05616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.05616;
        end
        ig = ig_('CAP00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03475;
        end
        [iv,ix_] = s2mpjlib('ii','Y00506',ix_);
        pb.xnames{iv} = 'Y00506';
        ig = ig_('CAP00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04304;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MXD00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.29601+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.29601;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03208+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03208;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MND00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00225;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03734+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03734;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.054;
        end
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00124+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00124;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.054;
        end
        ig = ig_('CAP01703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00233;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('MND00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00167;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('CAP01701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00078;
        end
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00603',ix_);
        pb.xnames{iv} = 'Y00603';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -33.21+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -33.21;
        end
        ig = ig_('CAP01601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .001;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.74+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.74;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.054;
        end
        ig = ig_('MXD00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00225;
        end
        ig = ig_('MND00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('CAP01603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03308+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03308;
        end
        ig = ig_('CAP01703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03812+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03812;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00167;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.35001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.35001;
        end
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00124+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00124;
        end
        ig = ig_('CAP01704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00233;
        end
        [iv,ix_] = s2mpjlib('ii','Y00604',ix_);
        pb.xnames{iv} = 'Y00604';
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00133+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00133;
        end
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.40401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.40401;
        end
        ig = ig_('CAP01705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00375+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00375;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('MXD00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03671+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03671;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('MND00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -27.28+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -27.28;
        end
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.30001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.30001;
        end
        [iv,ix_] = s2mpjlib('ii','Y00605',ix_);
        pb.xnames{iv} = 'Y00605';
        ig = ig_('CAP01604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03185;
        end
        ig = ig_('CAP01605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0029+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0029;
        end
        [iv,ix_] = s2mpjlib('ii','Y00606',ix_);
        pb.xnames{iv} = 'Y00606';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00234+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00234;
        end
        ig = ig_('MXD00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00606',ix_);
        pb.xnames{iv} = 'Y00606';
        ig = ig_('CAP06006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00349+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00349;
        end
        ig = ig_('MND00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00606',ix_);
        pb.xnames{iv} = 'Y00606';
        ig = ig_('CAP01706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04046+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04046;
        end
        ig = ig_('CAP01606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03475;
        end
        [iv,ix_] = s2mpjlib('ii','Y00606',ix_);
        pb.xnames{iv} = 'Y00606';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -25.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -25.9;
        end
        ig = ig_('CAP01705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03965+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03965;
        end
        [iv,ix_] = s2mpjlib('ii','Y00606',ix_);
        pb.xnames{iv} = 'Y00606';
        ig = ig_('CAP01605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0344+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0344;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('CAP00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00217;
        end
        ig = ig_('CAP06402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .001;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('MND00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('MND00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('CAP06403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0009;
        end
        ig = ig_('MXD00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('MXD00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.24142+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.24142;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('MND00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0111+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0111;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -87.07+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -87.07;
        end
        ig = ig_('MXD00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('INV00601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.37243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.37243;
        end
        ig = ig_('CAP00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00264+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00264;
        end
        [iv,ix_] = s2mpjlib('ii','Y00703',ix_);
        pb.xnames{iv} = 'Y00703';
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01192;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('CAP00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00159+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00159;
        end
        ig = ig_('MND00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('MXD00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01215+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01215;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('MND00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('CAP06403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .001;
        end
        ig = ig_('INV00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.24828+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.24828;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00108+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00108;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -82.66+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -82.66;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('MND00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('CAP06404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0009;
        end
        ig = ig_('MXD00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00704',ix_);
        pb.xnames{iv} = 'Y00704';
        ig = ig_('INV00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.36556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.36556;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('INV00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.37476+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.37476;
        end
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('MND00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -78.48+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -78.48;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('INV00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.23909+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.23909;
        end
        ig = ig_('CAP06405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00094+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00094;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('CAP06404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00097+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00097;
        end
        ig = ig_('MND00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('MXD00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01305+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01305;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('CAP00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01221+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01221;
        end
        ig = ig_('CAP00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00153+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00153;
        end
        [iv,ix_] = s2mpjlib('ii','Y00705',ix_);
        pb.xnames{iv} = 'Y00705';
        ig = ig_('MXD00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('CAP00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01409;
        end
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00056+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00056;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('MND00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -74.51+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -74.51;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('CAP00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01522+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01522;
        end
        ig = ig_('CAP06406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00191+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00191;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('CAP06405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        ig = ig_('CAP00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01374+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01374;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('INV00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.19366+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.19366;
        end
        ig = ig_('CAP00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01484+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01484;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('INV00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.61385+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.61385;
        end
        ig = ig_('CAP00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0011+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0011;
        end
        [iv,ix_] = s2mpjlib('ii','Y00706',ix_);
        pb.xnames{iv} = 'Y00706';
        ig = ig_('MXD00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.74296+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.74296;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MXD00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01774+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01774;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('CAP06501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00033+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00033;
        end
        ig = ig_('CAP00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00241;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.64004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.64004;
        end
        ig = ig_('MND00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('CAP00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00401;
        end
        ig = ig_('CAP00601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01267+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01267;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.54157+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.54157;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.06;
        end
        ig = ig_('MXD00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0023+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0023;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MND00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00802',ix_);
        pb.xnames{iv} = 'Y00802';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('CAP06502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0003;
        end
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.68169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.68169;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('CAP00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00222+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00222;
        end
        ig = ig_('MND00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.49991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.49991;
        end
        ig = ig_('CAP00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.86+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.86;
        end
        ig = ig_('MXD00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('CAP06503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00232+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00232;
        end
        ig = ig_('CAP00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01872+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01872;
        end
        [iv,ix_] = s2mpjlib('ii','Y00803',ix_);
        pb.xnames{iv} = 'Y00803';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('CAP00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.66+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.66;
        end
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.68169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.68169;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('CAP00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01872+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01872;
        end
        ig = ig_('CAP06503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0003;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('CAP06504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00232+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00232;
        end
        ig = ig_('CAP00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.49991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.49991;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00804',ix_);
        pb.xnames{iv} = 'Y00804';
        ig = ig_('CAP00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00222+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00222;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('CAP00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        ig = ig_('CAP06504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00029+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00029;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00233;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.48139+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.48139;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('CAP00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01915;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.48+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.48;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('CAP00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00214+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00214;
        end
        ig = ig_('CAP00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00428+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00428;
        end
        [iv,ix_] = s2mpjlib('ii','Y00805',ix_);
        pb.xnames{iv} = 'Y00805';
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.70021+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.70021;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('CAP00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00642+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00642;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('CAP06505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00031+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00031;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('INV00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.1816+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.1816;
        end
        ig = ig_('CAP06506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00262+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00262;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('CAP00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03041;
        end
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.51991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.51991;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('CAP00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00231+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00231;
        end
        ig = ig_('CAP00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01217;
        end
        [iv,ix_] = s2mpjlib('ii','Y00806',ix_);
        pb.xnames{iv} = 'Y00806';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.3;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MND00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.06;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MXD00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01267+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01267;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.54157+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.54157;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('CAP06501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00033+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00033;
        end
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.64004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.64004;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00401;
        end
        ig = ig_('CAP01102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01774+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01774;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00241;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('CAP06502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0023+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0023;
        end
        ig = ig_('MXD00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00902',ix_);
        pb.xnames{iv} = 'Y00902';
        ig = ig_('MND00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('CAP06503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00232+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00232;
        end
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.68169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.68169;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('CAP06502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0003;
        end
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.49991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.49991;
        end
        ig = ig_('CAP01103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01872+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01872;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.86+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.86;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('MND00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00222+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00222;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('MXD00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00903',ix_);
        pb.xnames{iv} = 'Y00903';
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('MND00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.49991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.49991;
        end
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('CAP06504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00232+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00232;
        end
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.66+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.66;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        ig = ig_('CAP01104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01872+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01872;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00222+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00222;
        end
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.68169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.68169;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('CAP06503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0003;
        end
        ig = ig_('CAP01103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        [iv,ix_] = s2mpjlib('ii','Y00904',ix_);
        pb.xnames{iv} = 'Y00904';
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.48139+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.48139;
        end
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.70021+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.70021;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('MND00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('CAP06505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00233;
        end
        ig = ig_('CAP01104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00214+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00214;
        end
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00428+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00428;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('CAP06504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00029+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00029;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.48+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.48;
        end
        [iv,ix_] = s2mpjlib('ii','Y00905',ix_);
        pb.xnames{iv} = 'Y00905';
        ig = ig_('CAP01105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01915;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('CAP01106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03041;
        end
        ig = ig_('CAP06505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00031+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00031;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('INV00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.1816+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.1816;
        end
        ig = ig_('CAP01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00642+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00642;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('MND00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.3;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('CAP06506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00262+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00262;
        end
        ig = ig_('CAP01105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01217;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('MXD00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.51991+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.51991;
        end
        [iv,ix_] = s2mpjlib('ii','Y00906',ix_);
        pb.xnames{iv} = 'Y00906';
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00231+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00231;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('CAP06002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00699;
        end
        ig = ig_('MND00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01383+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01383;
        end
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.2326+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.2326;
        end
        ig = ig_('INV00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.22411+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.22411;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('CAP01201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        ig = ig_('MND00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('MXD00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00115+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00115;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -33.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -33.1;
        end
        ig = ig_('MXD00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00995+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00995;
        end
        [iv,ix_] = s2mpjlib('ii','Y01003',ix_);
        pb.xnames{iv} = 'Y01003';
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00384+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00384;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00043+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00043;
        end
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01498;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00699;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -31.42+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -31.42;
        end
        ig = ig_('MXD00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.11205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.11205;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01082+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01082;
        end
        ig = ig_('MND00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01004',ix_);
        pb.xnames{iv} = 'Y01004';
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00384+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00384;
        end
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.34466+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.34466;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00674+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00674;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0041;
        end
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.34881+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.34881;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -29.83+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -29.83;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1079;
        end
        ig = ig_('CAP01204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01083;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00055+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00055;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01005',ix_);
        pb.xnames{iv} = 'Y01005';
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01443+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01443;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('INV00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.45671+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.45671;
        end
        ig = ig_('INV00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.57325+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.57325;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00727+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00727;
        end
        ig = ig_('CAP01205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01215+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01215;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01558;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.32+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.32;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.05827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.05827;
        end
        ig = ig_('CAP01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01498;
        end
        [iv,ix_] = s2mpjlib('ii','Y01006',ix_);
        pb.xnames{iv} = 'Y01006';
        ig = ig_('CAP01206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01125;
        end
        ig = ig_('CAP06006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01084+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01084;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('MXD00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('CAP06002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00699;
        end
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00384+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00384;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.2326+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.2326;
        end
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01383+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01383;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('CAP00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00115+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00115;
        end
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00995+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00995;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('MND00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('CAP00701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -33.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -33.1;
        end
        ig = ig_('MXD00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01103',ix_);
        pb.xnames{iv} = 'Y01103';
        ig = ig_('INV00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.22411+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.22411;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00043+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00043;
        end
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01498;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('CAP06003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00699;
        end
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.34466+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.34466;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00384+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00384;
        end
        ig = ig_('MXD00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('MND00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -31.42+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -31.42;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01082+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01082;
        end
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01104',ix_);
        pb.xnames{iv} = 'Y01104';
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.11205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.11205;
        end
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.34881+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.34881;
        end
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -29.83+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -29.83;
        end
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01443+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01443;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        ig = ig_('MND00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('CAP00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01083;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00055+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00055;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('CAP06004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00674+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00674;
        end
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0041;
        end
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1079;
        end
        [iv,ix_] = s2mpjlib('ii','Y01105',ix_);
        pb.xnames{iv} = 'Y01105';
        ig = ig_('MXD00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('CAP00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01215+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01215;
        end
        ig = ig_('MXD00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('MND00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.32+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.32;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('CAP06006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01084+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01084;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01558;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('INV00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.45671+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.45671;
        end
        ig = ig_('INV00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.57325+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.57325;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('CAP00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01498;
        end
        ig = ig_('CAP06005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00727+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00727;
        end
        [iv,ix_] = s2mpjlib('ii','Y01106',ix_);
        pb.xnames{iv} = 'Y01106';
        ig = ig_('CAP00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01125;
        end
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.05827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.05827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.86606+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.86606;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.6+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.6;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MXD00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MND00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MND00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00186+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00186;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00706;
        end
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0099+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0099;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.28869+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.28869;
        end
        [iv,ix_] = s2mpjlib('ii','Y01202',ix_);
        pb.xnames{iv} = 'Y01202';
        ig = ig_('CAP01001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00198;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0072;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.32+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.32;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('MND00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01005;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00183;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01203',ix_);
        pb.xnames{iv} = 'Y01203';
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00171;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('CAP01204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0072;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01005;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.05;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00171;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00183;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01204',ix_);
        pb.xnames{iv} = 'Y01204';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01012+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01012;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('CAP01205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00727+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00727;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89814+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89814;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('CAP01204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00165+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00165;
        end
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00176+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00176;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.79+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.79;
        end
        [iv,ix_] = s2mpjlib('ii','Y01205',ix_);
        pb.xnames{iv} = 'Y01205';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.25661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.25661;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('CAP01205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00178+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00178;
        end
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0019+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0019;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.27714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.27714;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP01206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00892+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00892;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('CAP01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01188+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01188;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.55+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.55;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.15475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.15475;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01206',ix_);
        pb.xnames{iv} = 'Y01206';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00706;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00198;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0099+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0099;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.28869+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.28869;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.86606+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.86606;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('MXD00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.62+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.62;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('MXD00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('CAP00701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00186+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00186;
        end
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01302',ix_);
        pb.xnames{iv} = 'Y01302';
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01005;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0072;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00183;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.34+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.34;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('MXD00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01303',ix_);
        pb.xnames{iv} = 'Y01303';
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00171;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.07+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.07;
        end
        ig = ig_('MND00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01005;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        ig = ig_('MXD00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('CAP00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0072;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00171;
        end
        [iv,ix_] = s2mpjlib('ii','Y01304',ix_);
        pb.xnames{iv} = 'Y01304';
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00183;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('CAP00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00165+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00165;
        end
        ig = ig_('MXD00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00176+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00176;
        end
        ig = ig_('CAP00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00727+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00727;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('MND00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01012+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01012;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.81+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.81;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89814+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89814;
        end
        [iv,ix_] = s2mpjlib('ii','Y01305',ix_);
        pb.xnames{iv} = 'Y01305';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.25661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.25661;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('CAP00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00178+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00178;
        end
        ig = ig_('MND00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('MXD00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.57+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.57;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('CAP00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01188+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01188;
        end
        ig = ig_('CAP00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00892+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00892;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0019+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0019;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.15475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.15475;
        end
        [iv,ix_] = s2mpjlib('ii','Y01306',ix_);
        pb.xnames{iv} = 'Y01306';
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.27714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.27714;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MND00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MND00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('CAP01001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00412+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00412;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.46;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.7014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.7014;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('CAP06201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00442+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00442;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00066;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('INV00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.501+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.501;
        end
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00824+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00824;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0058+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0058;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('CAP01201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00348+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00348;
        end
        ig = ig_('CAP06101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01402',ix_);
        pb.xnames{iv} = 'Y01402';
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00856+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00856;
        end
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.18+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.18;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00607+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00607;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('MND00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('MXD00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('CAP01202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00321;
        end
        ig = ig_('CAP01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0038+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0038;
        end
        [iv,ix_] = s2mpjlib('ii','Y01403',ix_);
        pb.xnames{iv} = 'Y01403';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00321;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('CAP01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0038+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0038;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.92+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.92;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00856+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00856;
        end
        [iv,ix_] = s2mpjlib('ii','Y01404',ix_);
        pb.xnames{iv} = 'Y01404';
        ig = ig_('CAP01204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00607+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00607;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.44534+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.44534;
        end
        ig = ig_('CAP01205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00619+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00619;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0087+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0087;
        end
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.75707+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.75707;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00074;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('CAP01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00366+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00366;
        end
        ig = ig_('CAP01204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0031+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0031;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00005;
        end
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00452+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00452;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01405',ix_);
        pb.xnames{iv} = 'Y01405';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.67;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('CAP01206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00929+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00929;
        end
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('CAP01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01237+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01237;
        end
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.48096+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.48096;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.44+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.44;
        end
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.20241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.20241;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('CAP06106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('CAP06206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00525+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00525;
        end
        ig = ig_('CAP01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00396+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00396;
        end
        [iv,ix_] = s2mpjlib('ii','Y01406',ix_);
        pb.xnames{iv} = 'Y01406';
        ig = ig_('CAP01205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00334+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00334;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00066;
        end
        ig = ig_('CAP06101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.46;
        end
        ig = ig_('INV00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.501+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.501;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00348+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00348;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('CAP06201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        ig = ig_('CAP00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00412+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00412;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00442+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00442;
        end
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00824+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00824;
        end
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.7014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.7014;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('MXD00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('MND00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01502',ix_);
        pb.xnames{iv} = 'Y01502';
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0058+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0058;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00607+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00607;
        end
        ig = ig_('MXD00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('CAP00702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00321;
        end
        ig = ig_('MND00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0038+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0038;
        end
        [iv,ix_] = s2mpjlib('ii','Y01503',ix_);
        pb.xnames{iv} = 'Y01503';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.18+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.18;
        end
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00856+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00856;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('CAP00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00607+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00607;
        end
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00856+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00856;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('CAP00703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00321;
        end
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('MND00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        ig = ig_('CAP00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0038+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0038;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01504',ix_);
        pb.xnames{iv} = 'Y01504';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.92+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.92;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('MND00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00005;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('MXD00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00074;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.67;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('CAP00704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0031+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0031;
        end
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00452+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00452;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('CAP00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00366+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00366;
        end
        ig = ig_('CAP00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00619+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00619;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.44534+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.44534;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0087+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0087;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01505',ix_);
        pb.xnames{iv} = 'Y01505';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.75707+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.75707;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.20241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.20241;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.48096+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.48096;
        end
        ig = ig_('MXD00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('CAP00706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00929+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00929;
        end
        ig = ig_('CAP00705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00334+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00334;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        ig = ig_('CAP00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00396+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00396;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('CAP00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01237+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01237;
        end
        ig = ig_('CAP06106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('MND00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00525+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00525;
        end
        [iv,ix_] = s2mpjlib('ii','Y01506',ix_);
        pb.xnames{iv} = 'Y01506';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.44+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.44;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.99+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.99;
        end
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('CAP00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00063+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00063;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.28869+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.28869;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.86606+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.86606;
        end
        ig = ig_('MND00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('MXD00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00131;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00654;
        end
        [iv,ix_] = s2mpjlib('ii','Y01602',ix_);
        pb.xnames{iv} = 'Y01602';
        ig = ig_('CAP00801');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00017+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00017;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('CAP00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00121+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00121;
        end
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.79+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.79;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00664;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01603',ix_);
        pb.xnames{iv} = 'Y01603';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('CAP00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00121+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00121;
        end
        ig = ig_('CAP00803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00664;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.59+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.59;
        end
        ig = ig_('CAP00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01604',ix_);
        pb.xnames{iv} = 'Y01604';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.41+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.41;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.25661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.25661;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('CAP00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00116+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00116;
        end
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('CAP00804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89814+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89814;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('CAP00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00668+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00668;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01605',ix_);
        pb.xnames{iv} = 'Y01605';
        ig = ig_('CAP00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.27714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.27714;
        end
        ig = ig_('CAP00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00784+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00784;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.24+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.24;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('CAP00805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00016+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00016;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('CAP00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00125;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('CAP00806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.15475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.15475;
        end
        [iv,ix_] = s2mpjlib('ii','Y01606',ix_);
        pb.xnames{iv} = 'Y01606';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('CAP00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00654;
        end
        ig = ig_('CAP01301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00017+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00017;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('MND00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.28869+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.28869;
        end
        ig = ig_('MXD00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('CAP01302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00063+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00063;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.97+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.97;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.86606+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.86606;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('CAP00901');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00131;
        end
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01702',ix_);
        pb.xnames{iv} = 'Y01702';
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.77+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.77;
        end
        ig = ig_('CAP01302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00121+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00121;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00664;
        end
        [iv,ix_] = s2mpjlib('ii','Y01703',ix_);
        pb.xnames{iv} = 'Y01703';
        ig = ig_('CAP01303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.88827+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.88827;
        end
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('CAP01303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.26648+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.26648;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.58+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.58;
        end
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('CAP00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00121+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00121;
        end
        ig = ig_('CAP01304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('CAP00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00664;
        end
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01704',ix_);
        pb.xnames{iv} = 'Y01704';
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('CAP01305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00065;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.25661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.25661;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89814+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89814;
        end
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.4;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('CAP00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00116+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00116;
        end
        ig = ig_('CAP01304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00015;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00668+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00668;
        end
        [iv,ix_] = s2mpjlib('ii','Y01705',ix_);
        pb.xnames{iv} = 'Y01705';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.15475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.15475;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.22+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.22;
        end
        ig = ig_('CAP01305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00016+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00016;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('CAP01306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.27714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.27714;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('CAP00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00125;
        end
        ig = ig_('CAP00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00784+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00784;
        end
        [iv,ix_] = s2mpjlib('ii','Y01706',ix_);
        pb.xnames{iv} = 'Y01706';
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.93+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.93;
        end
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.14434+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.14434;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('CAP06301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP01802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00071+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00071;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00662+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00662;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0104;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('CAP01801');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00009;
        end
        ig = ig_('MXD00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01802',ix_);
        pb.xnames{iv} = 'Y01802';
        ig = ig_('MND00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.02151+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.02151;
        end
        ig = ig_('MXD00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('MND00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00071+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00071;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('CAP01802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00008;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('CAP06302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP01803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00668+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00668;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.73+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.73;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.13324+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.13324;
        end
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01803',ix_);
        pb.xnames{iv} = 'Y01803';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00668+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00668;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.13324+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.13324;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.02151+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.02151;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('MND00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.54+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.54;
        end
        [iv,ix_] = s2mpjlib('ii','Y01804',ix_);
        pb.xnames{iv} = 'Y01804';
        ig = ig_('MXD00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00071+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00071;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('MXD00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('CAP01805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0004;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00068+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00068;
        end
        ig = ig_('CAP06304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00671+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00671;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.02644+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.02644;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('CAP01804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00008;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.12831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.12831;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.36+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.36;
        end
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01805',ix_);
        pb.xnames{iv} = 'Y01805';
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.13857+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.13857;
        end
        ig = ig_('MXD00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('CAP01406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00739+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00739;
        end
        ig = ig_('CAP01806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('MND00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.19+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.19;
        end
        ig = ig_('CAP06305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00074;
        end
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.15475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.15475;
        end
        [iv,ix_] = s2mpjlib('ii','Y01806',ix_);
        pb.xnames{iv} = 'Y01806';
        ig = ig_('CAP01805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00009;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.67;
        end
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('CAP00901');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00272+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00272;
        end
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00066;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('CAP06101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('CAP01302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00052+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00052;
        end
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00031+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00031;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('MND01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('MXD01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.7014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.7014;
        end
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('INV00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.501+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.501;
        end
        ig = ig_('CAP00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00544+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00544;
        end
        [iv,ix_] = s2mpjlib('ii','Y01902',ix_);
        pb.xnames{iv} = 'Y01902';
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00442+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00442;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.49+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.49;
        end
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('CAP00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00565+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00565;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('CAP01302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00029+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00029;
        end
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('MND01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        ig = ig_('CAP01303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00054;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('CAP00902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00251+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00251;
        end
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01903',ix_);
        pb.xnames{iv} = 'Y01903';
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('CAP01303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00029+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00029;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('CAP01304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00054;
        end
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.73994+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.73994;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('CAP00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00565+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00565;
        end
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.46246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.46246;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.31+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.31;
        end
        ig = ig_('CAP00903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00251+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00251;
        end
        [iv,ix_] = s2mpjlib('ii','Y01904',ix_);
        pb.xnames{iv} = 'Y01904';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('CAP00904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00242+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00242;
        end
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('CAP01304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00028+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00028;
        end
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00452+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00452;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.14+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.14;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.44534+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.44534;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('CAP01305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00056+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00056;
        end
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00005;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00575+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00575;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00074;
        end
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01905',ix_);
        pb.xnames{iv} = 'Y01905';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.75707+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.75707;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('CAP01305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0003;
        end
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.48096+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.48096;
        end
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('CAP01306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.99+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.99;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('CAP00905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00261+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00261;
        end
        ig = ig_('CAP06106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        [iv,ix_] = s2mpjlib('ii','Y01906',ix_);
        pb.xnames{iv} = 'Y01906';
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.20241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.20241;
        end
        ig = ig_('CAP00906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00817+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00817;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.65+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.65;
        end
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00066;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP01802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0006;
        end
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('MND01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.85171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.85171;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00561+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00561;
        end
        ig = ig_('MXD01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('MND01002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP01801');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00023+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00023;
        end
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP01401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00209+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00209;
        end
        ig = ig_('CAP06201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP06101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('INV00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.3507+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.3507;
        end
        [iv,ix_] = s2mpjlib('ii','Y02002',ix_);
        pb.xnames{iv} = 'Y02002';
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00442+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00442;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.46;
        end
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00062+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00062;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('CAP01402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00192;
        end
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('MND01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00577+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00577;
        end
        ig = ig_('CAP01802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00022+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00022;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.32373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.32373;
        end
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.87868+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.87868;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('MXD01003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP06102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('CAP06202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02003',ix_);
        pb.xnames{iv} = 'Y02003';
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('CAP01804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00062+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00062;
        end
        ig = ig_('CAP01803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00022+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00022;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('CAP06103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('CAP06203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00577+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00577;
        end
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00449+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00449;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        ig = ig_('MND01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.32373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.32373;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('MXD01004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('CAP01403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00192;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.87868+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.87868;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02004',ix_);
        pb.xnames{iv} = 'Y02004';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.29+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.29;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('CAP06104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00005;
        end
        ig = ig_('CAP01804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00021+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00021;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.89067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.89067;
        end
        ig = ig_('CAP06204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00074+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00074;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.31174+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.31174;
        end
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00585+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00585;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('CAP01404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00185;
        end
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00067;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.12+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.12;
        end
        ig = ig_('CAP01805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00062+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00062;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00452+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00452;
        end
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MXD01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02005',ix_);
        pb.xnames{iv} = 'Y02005';
        ig = ig_('MND01005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.96+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.96;
        end
        ig = ig_('CAP06206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00525+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00525;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('CAP06106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33667+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33667;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('MXD01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('MND01006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.20241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.20241;
        end
        ig = ig_('CAP01806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00083;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('CAP01405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .002;
        end
        ig = ig_('CAP01406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0077;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('CAP06205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0008;
        end
        ig = ig_('CAP06105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','Y02006',ix_);
        pb.xnames{iv} = 'Y02006';
        ig = ig_('CAP01805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00022+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00022;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00677+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00677;
        end
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02109+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02109;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP03101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0033+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0033;
        end
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0013;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00072+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00072;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00936+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00936;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00546+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00546;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0395+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0395;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP01902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01285;
        end
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01142+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01142;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04825+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04825;
        end
        ig = ig_('CAP02501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00294+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00294;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00356+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00356;
        end
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .10013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .10013;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('CAP02802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06479+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06479;
        end
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01241;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 280.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 280.0;
        end
        ig = ig_('CAP03102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08414+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08414;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00947+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00947;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00042;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00226+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00226;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01324+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01324;
        end
        ig = ig_('CAP03102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0021+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0021;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 280.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 280.0;
        end
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01226+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01226;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00211+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00211;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08241;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .09784+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .09784;
        end
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02167;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00702+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00702;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05183;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01437+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01437;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04242+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04242;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00149+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00149;
        end
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06185;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .09897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .09897;
        end
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00113+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00113;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01226+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01226;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05183+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05183;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00141+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00141;
        end
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01394+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01394;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00105+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00105;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04242+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04242;
        end
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06185;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01437+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01437;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08346+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08346;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02167+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02167;
        end
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00016+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00016;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00149+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00149;
        end
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00702+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00702;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00973+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00973;
        end
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 280.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 280.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08854+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08854;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .10487+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .10487;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('INV00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 280.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 280.0;
        end
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06479+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06479;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01608;
        end
        ig = ig_('CAP01906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03394+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03394;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01139+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01139;
        end
        ig = ig_('CAP03106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .12694+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .12694;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01223+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01223;
        end
        ig = ig_('CAP02506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01535;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .07621+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .07621;
        end
        ig = ig_('CAP02706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .15193+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .15193;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00546+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00546;
        end
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01036;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00016+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00016;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0011+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0011;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00148+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00148;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00119+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00119;
        end
        [iv,ix_] = s2mpjlib('ii','X00106',ix_);
        pb.xnames{iv} = 'X00106';
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01285;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00766+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00766;
        end
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01492+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01492;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00051;
        end
        ig = ig_('CAP02701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00333;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00296+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00296;
        end
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01161+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01161;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01487+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01487;
        end
        ig = ig_('CAP03101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00008+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00008;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0098+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0098;
        end
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00459+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00459;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0206+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0206;
        end
        ig = ig_('CAP01902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01768+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01768;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02157+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02157;
        end
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02094+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02094;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00342+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00342;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03082+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03082;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00079;
        end
        ig = ig_('CAP03102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06879+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06879;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 75.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 75.0;
        end
        ig = ig_('CAP03202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02505+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02505;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01722+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01722;
        end
        ig = ig_('CAP02602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02373;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08608;
        end
        ig = ig_('CAP02302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01551+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01551;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02391+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02391;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06574+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06574;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01601+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01601;
        end
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01644+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01644;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08418+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08418;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03395+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03395;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02265;
        end
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01241;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01481+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01481;
        end
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00731+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00731;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00987+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00987;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00217;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02265;
        end
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00216;
        end
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 75.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 75.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01572+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01572;
        end
        ig = ig_('CAP02304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02131;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00494+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00494;
        end
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00073+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00073;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .024;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0015+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0015;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .024;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00731+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00731;
        end
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00117;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('INV00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 75.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 75.0;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01601+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01601;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01012+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01012;
        end
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02265;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03395+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03395;
        end
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00025+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00025;
        end
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01644+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01644;
        end
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00073+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00073;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00494+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00494;
        end
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01572+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01572;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .08519+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .08519;
        end
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02391+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02391;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02265;
        end
        ig = ig_('CAP02305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02131;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06574+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06574;
        end
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01307+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01307;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('CAP02304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01481+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01481;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP03206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03992+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03992;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00122;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01722+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01722;
        end
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01768+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01768;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01526+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01526;
        end
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01086+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01086;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .09047+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .09047;
        end
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00026+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00026;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00766+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00766;
        end
        ig = ig_('INV00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 75.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 75.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02373;
        end
        ig = ig_('CAP02305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01551+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01551;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00157+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00157;
        end
        ig = ig_('CAP02806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02064+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02064;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06887;
        end
        ig = ig_('CAP02506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01457+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01457;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP01906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0326+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0326;
        end
        ig = ig_('CAP02306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03611+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03611;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11036;
        end
        ig = ig_('CAP02606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04529+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04529;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP02006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0111+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0111;
        end
        ig = ig_('CAP02106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01225;
        end
        [iv,ix_] = s2mpjlib('ii','X00206',ix_);
        pb.xnames{iv} = 'X00206';
        ig = ig_('CAP03106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .09969+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .09969;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02505+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02505;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP03002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00998+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00998;
        end
        ig = ig_('CAP02501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00285;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .12349+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .12349;
        end
        ig = ig_('CAP02403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00837+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00837;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        ig = ig_('CAP02903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00719+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00719;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00392+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00392;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01187+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01187;
        end
        ig = ig_('CAP02902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .028+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .028;
        end
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11034+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11034;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02005;
        end
        ig = ig_('CAP02001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00079;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01664;
        end
        ig = ig_('CAP02203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02513+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02513;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02041;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00977+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00977;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00054;
        end
        ig = ig_('CAP03003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01077;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP03202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01831;
        end
        ig = ig_('CAP02701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP03102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00594+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00594;
        end
        ig = ig_('CAP02402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01363+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01363;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP02602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01608;
        end
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02973+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02973;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02077;
        end
        ig = ig_('CAP01902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01255;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01979+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01979;
        end
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00205;
        end
        ig = ig_('CAP02204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03078;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1083;
        end
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03066;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0152+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0152;
        end
        ig = ig_('CAP02404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00899+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00899;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00376+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00376;
        end
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0216;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('CAP03003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00953+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00953;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01589+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01589;
        end
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01748+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01748;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00567+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00567;
        end
        ig = ig_('CAP02403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00791+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00791;
        end
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00467+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00467;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00983+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00983;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01225;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01948+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01948;
        end
        ig = ig_('CAP02203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11861+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11861;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02409;
        end
        ig = ig_('CAP03004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01122;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02036;
        end
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02873;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01266;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01535;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02036;
        end
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01333;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03078;
        end
        ig = ig_('CAP02904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0152+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0152;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11861+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11861;
        end
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01589+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01589;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01535;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .10954+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .10954;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03066;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00567+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00567;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP03004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00953+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00953;
        end
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00138+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00138;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02873;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0216;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01225;
        end
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01948+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01948;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00253+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00253;
        end
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01748+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01748;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00025+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00025;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02409;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01009;
        end
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP02905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00791+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00791;
        end
        ig = ig_('CAP02405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00899+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00899;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('CAP03005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01122;
        end
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00467+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00467;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP03206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03908+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03908;
        end
        ig = ig_('CAP03006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02075+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02075;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05014;
        end
        ig = ig_('CAP02206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .14939+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .14939;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP03106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01793+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01793;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1174+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1174;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00594+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00594;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00265;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02056+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02056;
        end
        ig = ig_('CAP02006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0111+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0111;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP03005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00998+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00998;
        end
        ig = ig_('CAP02905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01664;
        end
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02041;
        end
        ig = ig_('CAP02405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01363+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01363;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00026+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00026;
        end
        ig = ig_('CAP02906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02312+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02312;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04408+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04408;
        end
        ig = ig_('CAP02505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01542;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01083;
        end
        ig = ig_('CAP02506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01472+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01472;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02201+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02201;
        end
        ig = ig_('CAP02205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .12425+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .12425;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .13616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .13616;
        end
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00145+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00145;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01255;
        end
        ig = ig_('CAP01906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03234+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03234;
        end
        [iv,ix_] = s2mpjlib('ii','X00306',ix_);
        pb.xnames{iv} = 'X00306';
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01608;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01831;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01583+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01583;
        end
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01286;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        ig = ig_('CAP04503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        ig = ig_('CAP03301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        ig = ig_('CAP03702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03304;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP03402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02558;
        end
        ig = ig_('CAP03802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02733+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02733;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01621+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01621;
        end
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01089+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01089;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00591+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00591;
        end
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00431+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00431;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00439;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00685+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00685;
        end
        ig = ig_('CAP04302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02743+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02743;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01706;
        end
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00362;
        end
        ig = ig_('CAP03302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01548+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01548;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02252+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02252;
        end
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02393;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02556;
        end
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03192;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01566+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01566;
        end
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00255;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01004;
        end
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP05402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00786+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00786;
        end
        ig = ig_('CAP05302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01148+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01148;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00027+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00027;
        end
        ig = ig_('CAP05102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0132;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP05002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01633+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01633;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 153.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 153.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('CAP04902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01393;
        end
        ig = ig_('CAP04101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00047+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00047;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00685+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00685;
        end
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0132;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02733+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02733;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00439;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00431+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00431;
        end
        ig = ig_('CAP04704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01583+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01583;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01633+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01633;
        end
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02558;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00591+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00591;
        end
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01089+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01089;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01148+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01148;
        end
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03304;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00786+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00786;
        end
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01554+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01554;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01706;
        end
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02743+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02743;
        end
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01286;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 153.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 153.0;
        end
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03192;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02393;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01621+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01621;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00255;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        ig = ig_('CAP04503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02252+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02252;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01004;
        end
        ig = ig_('CAP04504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01217;
        end
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00362;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01393;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('CAP04004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02556;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00313+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00313;
        end
        ig = ig_('CAP04505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01772+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01772;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01496+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01496;
        end
        ig = ig_('CAP04004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01561+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01561;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02632;
        end
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00678+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00678;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01306+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01306;
        end
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00095+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00095;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03287+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03287;
        end
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00485+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00485;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01807;
        end
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02464+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02464;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02616;
        end
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02435+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02435;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00488+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00488;
        end
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01534+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01534;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00967+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00967;
        end
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01341+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01341;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01573+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01573;
        end
        ig = ig_('CAP04504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02169;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 153.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 153.0;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01271+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01271;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00502+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00502;
        end
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01106+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01106;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00757+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00757;
        end
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00163+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00163;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02642+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02642;
        end
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01118+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01118;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00634+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00634;
        end
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01084+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01084;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00745+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00745;
        end
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00483+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00483;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01172+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01172;
        end
        ig = ig_('CAP04705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0162+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0162;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03181+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03181;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00542;
        end
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01194+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01194;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01373;
        end
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01699;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01448;
        end
        ig = ig_('CAP04705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01044+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01044;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01656+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01656;
        end
        ig = ig_('CAP04505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02342+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02342;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02853+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02853;
        end
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00817+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00817;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01171;
        end
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01266;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01686+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01686;
        end
        ig = ig_('CAP04906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01824+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01824;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03436+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03436;
        end
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02661;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP03306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01809+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01809;
        end
        ig = ig_('CAP04406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01808+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01808;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02586+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02586;
        end
        ig = ig_('CAP04606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01696+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01696;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0394+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0394;
        end
        ig = ig_('CAP05006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02318+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02318;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP05106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01759;
        end
        ig = ig_('CAP05306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0174+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0174;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP05406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01875+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01875;
        end
        ig = ig_('CAP04306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03319+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03319;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01616;
        end
        ig = ig_('CAP04206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03519+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03519;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP04106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01267+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01267;
        end
        ig = ig_('CAP04006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04177+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04177;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP03806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04439;
        end
        ig = ig_('CAP03706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03666+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03666;
        end
        [iv,ix_] = s2mpjlib('ii','X00406',ix_);
        pb.xnames{iv} = 'X00406';
        ig = ig_('CAP03406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05751+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05751;
        end
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02843+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02843;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP03002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00998+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00998;
        end
        ig = ig_('CAP02802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01664;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00392+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00392;
        end
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02005;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        ig = ig_('CAP03102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00594+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00594;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11034+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11034;
        end
        ig = ig_('CAP03202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01831;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        ig = ig_('CAP02902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .028+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .028;
        end
        ig = ig_('CAP02903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00719+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00719;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP03003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01077;
        end
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02077;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00285;
        end
        ig = ig_('CAP02201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00079;
        end
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01979+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01979;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        ig = ig_('CAP02203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02513+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02513;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00054+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00054;
        end
        ig = ig_('CAP02202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .12349+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .12349;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP01902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01255;
        end
        ig = ig_('CAP02402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01363+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01363;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02041;
        end
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01187+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01187;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02973+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02973;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00977+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00977;
        end
        [iv,ix_] = s2mpjlib('ii','X00503',ix_);
        pb.xnames{iv} = 'X00503';
        ig = ig_('CAP02602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01608;
        end
        ig = ig_('CAP02403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00837+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00837;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP03103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00567+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00567;
        end
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0216;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP01903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00983+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00983;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01266;
        end
        ig = ig_('CAP03203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01748+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01748;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0152+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0152;
        end
        ig = ig_('CAP02103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01948+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01948;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01589+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01589;
        end
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1083;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11861+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11861;
        end
        ig = ig_('CAP02403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01535;
        end
        ig = ig_('CAP03003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00953+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00953;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP03004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01122;
        end
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00467+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00467;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        ig = ig_('CAP02702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00376+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00376;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02409;
        end
        ig = ig_('CAP02904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00791+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00791;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02873;
        end
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03066;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00205;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01225;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        ig = ig_('CAP02002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02036;
        end
        ig = ig_('CAP02204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03078;
        end
        [iv,ix_] = s2mpjlib('ii','X00504',ix_);
        pb.xnames{iv} = 'X00504';
        ig = ig_('CAP02404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00899+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00899;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00253+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00253;
        end
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02036+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02036;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00899+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00899;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00467+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00467;
        end
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .10954+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .10954;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00025+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00025;
        end
        ig = ig_('CAP02404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01301+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01301;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00791+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00791;
        end
        ig = ig_('CAP02104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01948+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01948;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00138+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00138;
        end
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02409;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02873;
        end
        ig = ig_('CAP02204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .11861+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .11861;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01009+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01009;
        end
        ig = ig_('CAP02604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01535;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP03004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00953+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00953;
        end
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01333;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0152+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0152;
        end
        ig = ig_('CAP02205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03078;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP01904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01198+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01198;
        end
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00077+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00077;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01589+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01589;
        end
        ig = ig_('CAP03204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01748+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01748;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03066+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03066;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0216+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0216;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01225+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01225;
        end
        ig = ig_('CAP03104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00567+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00567;
        end
        [iv,ix_] = s2mpjlib('ii','X00505',ix_);
        pb.xnames{iv} = 'X00505';
        ig = ig_('CAP03005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01122;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00265+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00265;
        end
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 43.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 43.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        ig = ig_('CAP02504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00145+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00145;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1174+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1174;
        end
        ig = ig_('CAP03205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01831;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP03105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00594+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00594;
        end
        ig = ig_('CAP02805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01664+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01664;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP03005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00998+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00998;
        end
        ig = ig_('CAP02106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05014;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP03006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02075+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02075;
        end
        ig = ig_('CAP02706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .13616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .13616;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP01905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01255;
        end
        ig = ig_('CAP02806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02056+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02056;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02312+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02312;
        end
        ig = ig_('CAP02606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04408+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04408;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00026+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00026;
        end
        ig = ig_('CAP03106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01793+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01793;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01083;
        end
        ig = ig_('CAP02205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .12425+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .12425;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP03206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03908+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03908;
        end
        ig = ig_('CAP01906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03234+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03234;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01472+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01472;
        end
        ig = ig_('CAP02406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02201+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02201;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01363+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01363;
        end
        ig = ig_('CAP02505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01542;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02041+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02041;
        end
        ig = ig_('CAP02206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .14939+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .14939;
        end
        [iv,ix_] = s2mpjlib('ii','X00506',ix_);
        pb.xnames{iv} = 'X00506';
        ig = ig_('CAP02605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01608+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01608;
        end
        ig = ig_('CAP02006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0111+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0111;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00591+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00591;
        end
        ig = ig_('CAP05102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0132;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00255;
        end
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01706;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00027+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00027;
        end
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01089+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01089;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP05302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01148+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01148;
        end
        ig = ig_('CAP05402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00786+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00786;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03192;
        end
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00362;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP03301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00006;
        end
        ig = ig_('CAP04101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00047+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00047;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        ig = ig_('CAP04702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01004;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00431+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00431;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02556;
        end
        ig = ig_('CAP04402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02393;
        end
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP03302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01548+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01548;
        end
        ig = ig_('CAP04902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01393;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02252+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02252;
        end
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00439;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01286;
        end
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00685+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00685;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01566+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01566;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 143.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 143.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02743+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02743;
        end
        ig = ig_('CAP03802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02733+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02733;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01621+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01621;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP05002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01633+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01633;
        end
        ig = ig_('CAP03702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03304;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        ig = ig_('CAP03402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02558;
        end
        [iv,ix_] = s2mpjlib('ii','X00603',ix_);
        pb.xnames{iv} = 'X00603';
        ig = ig_('CAP04102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0117;
        end
        ig = ig_('CAP04703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01583+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01583;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01393;
        end
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01706+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01706;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00576+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00576;
        end
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01633+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01633;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00362;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02393+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02393;
        end
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02743+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02743;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01126;
        end
        ig = ig_('CAP04004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02556+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02556;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02558+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02558;
        end
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00685+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00685;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01004+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01004;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00439;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 143.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 143.0;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00255;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01148+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01148;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01593+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01593;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03192+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03192;
        end
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02733+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02733;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01583+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01583;
        end
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00431+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00431;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02252+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02252;
        end
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00104+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00104;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0132;
        end
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00786+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00786;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01286;
        end
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01089+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01089;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03304+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03304;
        end
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00521+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00521;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00591+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00591;
        end
        ig = ig_('CAP04003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01621+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01621;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01217;
        end
        ig = ig_('CAP04504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01688+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01688;
        end
        [iv,ix_] = s2mpjlib('ii','X00604',ix_);
        pb.xnames{iv} = 'X00604';
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01554+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01554;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00095+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00095;
        end
        ig = ig_('CAP04505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01772+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01772;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02616;
        end
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01807;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02435+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02435;
        end
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00485+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00485;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03287+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03287;
        end
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00678+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00678;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01306+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01306;
        end
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00163+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00163;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01118+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01118;
        end
        ig = ig_('CAP04705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0162+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0162;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00483+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00483;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 143.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 143.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00745+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00745;
        end
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00488+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00488;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00634+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00634;
        end
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00313+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00313;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02642+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02642;
        end
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01084+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01084;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01271+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01271;
        end
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01106+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01106;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03181+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03181;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01496+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01496;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00967+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00967;
        end
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02632;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00502+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00502;
        end
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00757+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00757;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02464+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02464;
        end
        ig = ig_('CAP04504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02169;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01534+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01534;
        end
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01341+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01341;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01172+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01172;
        end
        ig = ig_('CAP04004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01561+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01561;
        end
        [iv,ix_] = s2mpjlib('ii','X00605',ix_);
        pb.xnames{iv} = 'X00605';
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01573+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01573;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01696+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01696;
        end
        ig = ig_('CAP04206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03519+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03519;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02853+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02853;
        end
        ig = ig_('CAP04106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01267+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01267;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04177+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04177;
        end
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01171+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01171;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01824+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01824;
        end
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00542;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP03806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04439+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04439;
        end
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03436+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03436;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03319+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03319;
        end
        ig = ig_('CAP04506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0394+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0394;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01808+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01808;
        end
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02843+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02843;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01616;
        end
        ig = ig_('CAP05106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01759;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP05406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01875+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01875;
        end
        ig = ig_('CAP05006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02318+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02318;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01686+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01686;
        end
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02661+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02661;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01266;
        end
        ig = ig_('CAP05306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0174+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0174;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02586+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02586;
        end
        ig = ig_('CAP03306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01809+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01809;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02342+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02342;
        end
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01373+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01373;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00817+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00817;
        end
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01699+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01699;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01194+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01194;
        end
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01448;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP03406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05751+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05751;
        end
        ig = ig_('CAP04705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01044+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01044;
        end
        [iv,ix_] = s2mpjlib('ii','X00606',ix_);
        pb.xnames{iv} = 'X00606';
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01656+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01656;
        end
        ig = ig_('CAP03706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03666+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03666;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00376+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00376;
        end
        ig = ig_('CAP03902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00969+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00969;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00093+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00093;
        end
        ig = ig_('CAP03802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01945+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01945;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00862+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00862;
        end
        ig = ig_('CAP04102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01116+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01116;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0244+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0244;
        end
        ig = ig_('CAP04402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00496+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00496;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00978+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00978;
        end
        ig = ig_('CAP05703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02434+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02434;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01463+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01463;
        end
        ig = ig_('CAP03302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01962+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01962;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02055+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02055;
        end
        ig = ig_('CAP03903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0079;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01478+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01478;
        end
        ig = ig_('CAP03701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0026+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0026;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0313+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0313;
        end
        ig = ig_('CAP05803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01677+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01677;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0227+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0227;
        end
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00396+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00396;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01299;
        end
        ig = ig_('CAP03702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0251+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0251;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02908+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02908;
        end
        ig = ig_('CAP05302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01822+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01822;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP04301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00012+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00012;
        end
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01977+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01977;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP03503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03477+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03477;
        end
        ig = ig_('CAP05902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01455+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01455;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 424.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 424.0;
        end
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00539+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00539;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00542;
        end
        ig = ig_('CAP04601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00258+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00258;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01259+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01259;
        end
        ig = ig_('CAP05501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00229+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00229;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02266;
        end
        ig = ig_('CAP04101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00203+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00203;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01024;
        end
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04854+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04854;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01247+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01247;
        end
        ig = ig_('CAP05802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01936+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01936;
        end
        [iv,ix_] = s2mpjlib('ii','X00703',ix_);
        pb.xnames{iv} = 'X00703';
        ig = ig_('CAP05603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01208+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01208;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01977+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01977;
        end
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0244+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0244;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01455+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01455;
        end
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0292+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0292;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02723+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02723;
        end
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00542;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01822+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01822;
        end
        ig = ig_('CAP05803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01936+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01936;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00862+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00862;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 424.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 424.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02434+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02434;
        end
        ig = ig_('CAP03504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03477+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03477;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00045+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00045;
        end
        ig = ig_('CAP05604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01208+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01208;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01218+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01218;
        end
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01247+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01247;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP04602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00125+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00125;
        end
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04854+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04854;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02266+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02266;
        end
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01299;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01677+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01677;
        end
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00396+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00396;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP04203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00978+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00978;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00376+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00376;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0079+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0079;
        end
        ig = ig_('CAP03803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01945+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01945;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01463+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01463;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00539+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00539;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0313+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0313;
        end
        ig = ig_('CAP05503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01374+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01374;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP04102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00101+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00101;
        end
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00496+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00496;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01612+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01612;
        end
        ig = ig_('CAP05904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02055+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02055;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0201+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0201;
        end
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0227+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0227;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00114+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00114;
        end
        ig = ig_('CAP03702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00047+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00047;
        end
        [iv,ix_] = s2mpjlib('ii','X00704',ix_);
        pb.xnames{iv} = 'X00704';
        ig = ig_('CAP05703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01024;
        end
        ig = ig_('CAP03903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00969+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00969;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01221+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01221;
        end
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00043+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00043;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01292+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01292;
        end
        ig = ig_('CAP05704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00986+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00986;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 424.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 424.0;
        end
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00571+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00571;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01401;
        end
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01267+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01267;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02109+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02109;
        end
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02668+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02668;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02472+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02472;
        end
        ig = ig_('CAP03904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00933+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00933;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01864+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01864;
        end
        ig = ig_('CAP05805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01749+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01749;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP04204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00942+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00942;
        end
        ig = ig_('CAP05503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0011+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0011;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01873;
        end
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02044+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02044;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00522+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00522;
        end
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00429+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00429;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01616+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01616;
        end
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00504+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00504;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01958+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01958;
        end
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02476+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02476;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00045+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00045;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0083;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00826;
        end
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01317+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01317;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02342+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02342;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0012+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0012;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0497+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0497;
        end
        ig = ig_('CAP03505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03531+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03531;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02812+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02812;
        end
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00057+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00057;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02182+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02182;
        end
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00098+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00098;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00477+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00477;
        end
        ig = ig_('CAP03404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03014+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03014;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP05504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01378+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01378;
        end
        ig = ig_('CAP05304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01755+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01755;
        end
        [iv,ix_] = s2mpjlib('ii','X00705',ix_);
        pb.xnames{iv} = 'X00705';
        ig = ig_('CAP03504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01409+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01409;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01789+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01789;
        end
        ig = ig_('CAP05906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0351+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0351;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01795+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01795;
        end
        ig = ig_('CAP03306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0243;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03799+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03799;
        end
        ig = ig_('CAP05506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01488+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01488;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01401+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01401;
        end
        ig = ig_('CAP03706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0277+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0277;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01319+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01319;
        end
        ig = ig_('CAP03806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04215+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04215;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03613+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03613;
        end
        ig = ig_('CAP03906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01759;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00022+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00022;
        end
        ig = ig_('CAP05706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03458+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03458;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0494+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0494;
        end
        ig = ig_('CAP04306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03316+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03316;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .07984+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .07984;
        end
        ig = ig_('CAP04606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01736+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01736;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03418+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03418;
        end
        ig = ig_('CAP05606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03474+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03474;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 424.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 424.0;
        end
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0293+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0293;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00564+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00564;
        end
        ig = ig_('CAP05305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01895+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01895;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02162+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02162;
        end
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00897;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03255;
        end
        ig = ig_('CAP03505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01522+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01522;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02023+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02023;
        end
        ig = ig_('CAP05505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01607+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01607;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP03905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01007+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01007;
        end
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01425+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01425;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01017+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01017;
        end
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01875+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01875;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03037;
        end
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00515+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00515;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0006;
        end
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00053+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00053;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02356+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02356;
        end
        ig = ig_('CAP05705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01065+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01065;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0006+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0006;
        end
        ig = ig_('CAP05905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01513+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01513;
        end
        [iv,ix_] = s2mpjlib('ii','X00706',ix_);
        pb.xnames{iv} = 'X00706';
        ig = ig_('CAP05805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02013+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02013;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00692+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00692;
        end
        ig = ig_('CAP05102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02245+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02245;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01289+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01289;
        end
        ig = ig_('CAP03701');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00081+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00081;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0323+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0323;
        end
        ig = ig_('CAP05002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01034+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01034;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01151+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01151;
        end
        ig = ig_('CAP04902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01364+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01364;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00097+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00097;
        end
        ig = ig_('CAP03603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01934+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01934;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02063+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02063;
        end
        ig = ig_('CAP05902');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01835+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01835;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00453+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00453;
        end
        ig = ig_('CAP04803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00494+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00494;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02611+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02611;
        end
        ig = ig_('CAP05702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01309;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03388+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03388;
        end
        ig = ig_('CAP05402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00679+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00679;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00098+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00098;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00037;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP05603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00898+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00898;
        end
        ig = ig_('CAP05202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01958+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01958;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00397+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00397;
        end
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00678+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00678;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('INV00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 26.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 26.0;
        end
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00271+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00271;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP03302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01518+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01518;
        end
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP03602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03149+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03149;
        end
        ig = ig_('CAP04602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01598+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01598;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00529+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00529;
        end
        ig = ig_('CAP03702');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05747+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05747;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP03301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00033+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00033;
        end
        ig = ig_('CAP05803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01413+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01413;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01196+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01196;
        end
        ig = ig_('CAP04302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03242+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03242;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04802');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01307+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01307;
        end
        ig = ig_('CAP05703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02476+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02476;
        end
        [iv,ix_] = s2mpjlib('ii','X00803',ix_);
        pb.xnames{iv} = 'X00803';
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00489+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00489;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01151+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01151;
        end
        ig = ig_('CAP04103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01294+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01294;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00397+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00397;
        end
        ig = ig_('CAP05204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0323+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0323;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00692+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00692;
        end
        ig = ig_('CAP04603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01695+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01695;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03388+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03388;
        end
        ig = ig_('CAP03703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05828+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05828;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('INV00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 26.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 26.0;
        end
        ig = ig_('CAP03303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0155;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP03603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03149+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03149;
        end
        ig = ig_('CAP05904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02063+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02063;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00679+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00679;
        end
        ig = ig_('CAP05804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01413+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01413;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02476+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02476;
        end
        ig = ig_('CAP05203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01958+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01958;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01307+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01307;
        end
        ig = ig_('CAP05103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02245+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02245;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00898+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00898;
        end
        ig = ig_('CAP05003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01034+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01034;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05703');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01309;
        end
        ig = ig_('CAP04903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01364+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01364;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00002;
        end
        ig = ig_('CAP03604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01934+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01934;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00678+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00678;
        end
        ig = ig_('CAP05903');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01835+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01835;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01289+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01289;
        end
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00489+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00489;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00453+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00453;
        end
        ig = ig_('CAP04804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00494+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00494;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00529+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00529;
        end
        ig = ig_('CAP05803');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02611+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02611;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00037;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00271+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00271;
        end
        [iv,ix_] = s2mpjlib('ii','X00804',ix_);
        pb.xnames{iv} = 'X00804';
        ig = ig_('CAP04303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03242+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03242;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0005+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0005;
        end
        ig = ig_('CAP04805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00542+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00542;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02131;
        end
        ig = ig_('CAP05604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03262+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03262;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00099+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00099;
        end
        ig = ig_('CAP04104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01246+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01246;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0126;
        end
        ig = ig_('CAP04604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01632;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02515+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02515;
        end
        ig = ig_('CAP04804');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01258+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01258;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00609+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00609;
        end
        ig = ig_('INV00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 26.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 26.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01767+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01767;
        end
        ig = ig_('CAP04404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0051;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01309;
        end
        ig = ig_('CAP03304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01493+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01493;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00776+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00776;
        end
        ig = ig_('CAP05805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0151+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0151;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP03604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03033+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03033;
        end
        ig = ig_('CAP05205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03303;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP03605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0205;
        end
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00328+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00328;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00717+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00717;
        end
        ig = ig_('CAP05104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02161+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02161;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP03704');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05612+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05612;
        end
        ig = ig_('CAP05004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00996+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00996;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01886+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01886;
        end
        ig = ig_('CAP04304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03122+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03122;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00613+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00613;
        end
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01176+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01176;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP04904');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01313+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01313;
        end
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00503+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00503;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01024+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01024;
        end
        ig = ig_('CAP05404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00654;
        end
        [iv,ix_] = s2mpjlib('ii','X00805',ix_);
        pb.xnames{iv} = 'X00805';
        ig = ig_('CAP05705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02525+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02525;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01361+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01361;
        end
        ig = ig_('CAP04805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01359+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01359;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05805');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02716+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02716;
        end
        ig = ig_('CAP05605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03523+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03523;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00707+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00707;
        end
        ig = ig_('CAP05205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02037;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02334+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02334;
        end
        ig = ig_('CAP05005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01076+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01076;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01418+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01418;
        end
        ig = ig_('CAP04605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01763+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01763;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05905');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01908+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01908;
        end
        ig = ig_('CAP04405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .00551+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .00551;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03372+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03372;
        end
        ig = ig_('CAP04105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01346+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01346;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP03705');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06061+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06061;
        end
        ig = ig_('CAP03605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03275+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03275;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP03305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01612+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01612;
        end
        ig = ig_('CAP03306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01821+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01821;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01817+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01817;
        end
        ig = ig_('CAP05806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04025+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04025;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03785+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03785;
        end
        ig = ig_('CAP05606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04286+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04286;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01831+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01831;
        end
        ig = ig_('CAP05206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05189+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05189;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP05106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .02937+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .02937;
        end
        ig = ig_('CAP05006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01713+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01713;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04806');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .018+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .018;
        end
        ig = ig_('CAP05906');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03898+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03898;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01732+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01732;
        end
        ig = ig_('CAP04406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01818+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01818;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP04306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .03731+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .03731;
        end
        ig = ig_('CAP04106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01296+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01296;
        end
        [iv,ix_] = s2mpjlib('ii','X00806',ix_);
        pb.xnames{iv} = 'X00806';
        ig = ig_('CAP03706');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .06226+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .06226;
        end
        ig = ig_('CAP03606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .05083+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .05083;
        end
        [iv,ix_] = s2mpjlib('ii','I00101',ix_);
        pb.xnames{iv} = 'I00101';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00101',ix_);
        pb.xnames{iv} = 'I00101';
        ig = ig_('INV00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00102',ix_);
        pb.xnames{iv} = 'I00102';
        ig = ig_('INV00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00102',ix_);
        pb.xnames{iv} = 'I00102';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00103',ix_);
        pb.xnames{iv} = 'I00103';
        ig = ig_('INV00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00103',ix_);
        pb.xnames{iv} = 'I00103';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00104',ix_);
        pb.xnames{iv} = 'I00104';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00104',ix_);
        pb.xnames{iv} = 'I00104';
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00105',ix_);
        pb.xnames{iv} = 'I00105';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00105',ix_);
        pb.xnames{iv} = 'I00105';
        ig = ig_('INV00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00106',ix_);
        pb.xnames{iv} = 'I00106';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00106');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00201',ix_);
        pb.xnames{iv} = 'I00201';
        ig = ig_('INV00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00201',ix_);
        pb.xnames{iv} = 'I00201';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00202',ix_);
        pb.xnames{iv} = 'I00202';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00202',ix_);
        pb.xnames{iv} = 'I00202';
        ig = ig_('INV00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00203',ix_);
        pb.xnames{iv} = 'I00203';
        ig = ig_('INV00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00203',ix_);
        pb.xnames{iv} = 'I00203';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00204',ix_);
        pb.xnames{iv} = 'I00204';
        ig = ig_('INV00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00204',ix_);
        pb.xnames{iv} = 'I00204';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00205',ix_);
        pb.xnames{iv} = 'I00205';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00205',ix_);
        pb.xnames{iv} = 'I00205';
        ig = ig_('INV00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00206',ix_);
        pb.xnames{iv} = 'I00206';
        ig = ig_('INV00206');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00301',ix_);
        pb.xnames{iv} = 'I00301';
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00301',ix_);
        pb.xnames{iv} = 'I00301';
        ig = ig_('INV00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00302',ix_);
        pb.xnames{iv} = 'I00302';
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00302',ix_);
        pb.xnames{iv} = 'I00302';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00303',ix_);
        pb.xnames{iv} = 'I00303';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00303',ix_);
        pb.xnames{iv} = 'I00303';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00304',ix_);
        pb.xnames{iv} = 'I00304';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00304',ix_);
        pb.xnames{iv} = 'I00304';
        ig = ig_('INV00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00305',ix_);
        pb.xnames{iv} = 'I00305';
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00305',ix_);
        pb.xnames{iv} = 'I00305';
        ig = ig_('INV00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00306',ix_);
        pb.xnames{iv} = 'I00306';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00306');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00401',ix_);
        pb.xnames{iv} = 'I00401';
        ig = ig_('INV00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00401',ix_);
        pb.xnames{iv} = 'I00401';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00402',ix_);
        pb.xnames{iv} = 'I00402';
        ig = ig_('INV00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00402',ix_);
        pb.xnames{iv} = 'I00402';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00403',ix_);
        pb.xnames{iv} = 'I00403';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00403',ix_);
        pb.xnames{iv} = 'I00403';
        ig = ig_('INV00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00404',ix_);
        pb.xnames{iv} = 'I00404';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00404',ix_);
        pb.xnames{iv} = 'I00404';
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00405',ix_);
        pb.xnames{iv} = 'I00405';
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00405',ix_);
        pb.xnames{iv} = 'I00405';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00406',ix_);
        pb.xnames{iv} = 'I00406';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00406');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00501',ix_);
        pb.xnames{iv} = 'I00501';
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00501',ix_);
        pb.xnames{iv} = 'I00501';
        ig = ig_('INV00501');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00502',ix_);
        pb.xnames{iv} = 'I00502';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00502',ix_);
        pb.xnames{iv} = 'I00502';
        ig = ig_('INV00502');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00503',ix_);
        pb.xnames{iv} = 'I00503';
        ig = ig_('INV00503');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00503',ix_);
        pb.xnames{iv} = 'I00503';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00504',ix_);
        pb.xnames{iv} = 'I00504';
        ig = ig_('INV00504');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00504',ix_);
        pb.xnames{iv} = 'I00504';
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00505',ix_);
        pb.xnames{iv} = 'I00505';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00505');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00505',ix_);
        pb.xnames{iv} = 'I00505';
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00506',ix_);
        pb.xnames{iv} = 'I00506';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        ig = ig_('INV00506');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00601',ix_);
        pb.xnames{iv} = 'I00601';
        ig = ig_('INV00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00601',ix_);
        pb.xnames{iv} = 'I00601';
        ig = ig_('INV00601');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00602',ix_);
        pb.xnames{iv} = 'I00602';
        ig = ig_('INV00602');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00602',ix_);
        pb.xnames{iv} = 'I00602';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00603',ix_);
        pb.xnames{iv} = 'I00603';
        ig = ig_('INV00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('INV00603');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00603',ix_);
        pb.xnames{iv} = 'I00603';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00604',ix_);
        pb.xnames{iv} = 'I00604';
        ig = ig_('INV00604');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('INV00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00604',ix_);
        pb.xnames{iv} = 'I00604';
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00605',ix_);
        pb.xnames{iv} = 'I00605';
        ig = ig_('INV00605');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
        end
        [iv,ix_] = s2mpjlib('ii','I00605',ix_);
        pb.xnames{iv} = 'I00605';
        ig = ig_('INV00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00606',ix_);
        pb.xnames{iv} = 'I00606';
        ig = ig_('INV00606');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('OBJECTIV');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.08;
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
        pbm.gconst(ig_('CAP00101')) = 23995.8;
        pbm.gconst(ig_('CAP00201')) = 2133.0;
        pbm.gconst(ig_('CAP00301')) = 6398.9;
        pbm.gconst(ig_('CAP00401')) = 7998.6;
        pbm.gconst(ig_('CAP00501')) = 29328.1;
        pbm.gconst(ig_('CAP00601')) = 7465.3;
        pbm.gconst(ig_('CAP00701')) = 11731.3;
        pbm.gconst(ig_('CAP00801')) = 1599.7;
        pbm.gconst(ig_('CAP00901')) = 40824.0;
        pbm.gconst(ig_('CAP01001')) = 27216.0;
        pbm.gconst(ig_('CAP01101')) = 6237.0;
        pbm.gconst(ig_('CAP01201')) = 2835.0;
        pbm.gconst(ig_('CAP01301')) = 1701.0;
        pbm.gconst(ig_('CAP01401')) = 6479.7;
        pbm.gconst(ig_('CAP01501')) = 540.0;
        pbm.gconst(ig_('CAP01601')) = 540.0;
        pbm.gconst(ig_('CAP01701')) = 540.0;
        pbm.gconst(ig_('CAP01801')) = 1080.0;
        pbm.gconst(ig_('CAP02001')) = 993.6;
        pbm.gconst(ig_('CAP02201')) = 441.6;
        pbm.gconst(ig_('CAP02501')) = 552.0;
        pbm.gconst(ig_('CAP02701')) = 2550.2;
        pbm.gconst(ig_('CAP03101')) = 728.6;
        pbm.gconst(ig_('CAP03301')) = 1296.0;
        pbm.gconst(ig_('CAP03701')) = 3013.2;
        pbm.gconst(ig_('CAP04101')) = 751.7;
        pbm.gconst(ig_('CAP04301')) = 1224.7;
        pbm.gconst(ig_('CAP04601')) = 1503.4;
        pbm.gconst(ig_('CAP05501')) = 1769.0;
        pbm.gconst(ig_('CAP06101')) = 3864.0;
        pbm.gconst(ig_('CAP06201')) = 1449.0;
        pbm.gconst(ig_('CAP06301')) = 3732.7;
        pbm.gconst(ig_('CAP06501')) = 2649.6;
        pbm.gconst(ig_('CAP00102')) = 21329.6;
        pbm.gconst(ig_('CAP00202')) = 1896.0;
        pbm.gconst(ig_('CAP00302')) = 5687.9;
        pbm.gconst(ig_('CAP00402')) = 7109.9;
        pbm.gconst(ig_('CAP00502')) = 26069.5;
        pbm.gconst(ig_('CAP00602')) = 6635.9;
        pbm.gconst(ig_('CAP00702')) = 10427.8;
        pbm.gconst(ig_('CAP00802')) = 1422.0;
        pbm.gconst(ig_('CAP00902')) = 36288.0;
        pbm.gconst(ig_('CAP01002')) = 24192.0;
        pbm.gconst(ig_('CAP01102')) = 5544.0;
        pbm.gconst(ig_('CAP01202')) = 2520.0;
        pbm.gconst(ig_('CAP01302')) = 1512.0;
        pbm.gconst(ig_('CAP01402')) = 5759.8;
        pbm.gconst(ig_('CAP01502')) = 480.0;
        pbm.gconst(ig_('CAP01602')) = 480.0;
        pbm.gconst(ig_('CAP01702')) = 480.0;
        pbm.gconst(ig_('CAP01802')) = 960.0;
        pbm.gconst(ig_('CAP01902')) = 652.8;
        pbm.gconst(ig_('CAP02002')) = 864.0;
        pbm.gconst(ig_('CAP02102')) = 859.7;
        pbm.gconst(ig_('CAP02202')) = 384.0;
        pbm.gconst(ig_('CAP02302')) = 573.1;
        pbm.gconst(ig_('CAP02402')) = 576.0;
        pbm.gconst(ig_('CAP02502')) = 480.0;
        pbm.gconst(ig_('CAP02602')) = 1488.0;
        pbm.gconst(ig_('CAP02702')) = 2217.6;
        pbm.gconst(ig_('CAP02802')) = 1392.0;
        pbm.gconst(ig_('CAP02902')) = 432.0;
        pbm.gconst(ig_('CAP03002')) = 480.0;
        pbm.gconst(ig_('CAP03102')) = 633.6;
        pbm.gconst(ig_('CAP03202')) = 1219.2;
        pbm.gconst(ig_('CAP03302')) = 1152.0;
        pbm.gconst(ig_('CAP03402')) = 2102.4;
        pbm.gconst(ig_('CAP03502')) = 115.2;
        pbm.gconst(ig_('CAP03602')) = 1209.6;
        pbm.gconst(ig_('CAP03702')) = 2678.4;
        pbm.gconst(ig_('CAP03802')) = 1157.8;
        pbm.gconst(ig_('CAP03902')) = 299.5;
        pbm.gconst(ig_('CAP04002')) = 328.3;
        pbm.gconst(ig_('CAP04102')) = 668.2;
        pbm.gconst(ig_('CAP04202')) = 1036.8;
        pbm.gconst(ig_('CAP04302')) = 1088.6;
        pbm.gconst(ig_('CAP04402')) = 362.9;
        pbm.gconst(ig_('CAP04502')) = 1002.2;
        pbm.gconst(ig_('CAP04602')) = 1336.3;
        pbm.gconst(ig_('CAP04702')) = 1336.3;
        pbm.gconst(ig_('CAP04802')) = 668.2;
        pbm.gconst(ig_('CAP04902')) = 668.2;
        pbm.gconst(ig_('CAP05002')) = 668.2;
        pbm.gconst(ig_('CAP05102')) = 679.7;
        pbm.gconst(ig_('CAP05202')) = 1670.4;
        pbm.gconst(ig_('CAP05302')) = 576.0;
        pbm.gconst(ig_('CAP05402')) = 725.8;
        pbm.gconst(ig_('CAP05502')) = 1572.5;
        pbm.gconst(ig_('CAP05602')) = 985.0;
        pbm.gconst(ig_('CAP05702')) = 898.6;
        pbm.gconst(ig_('CAP05802')) = 985.0;
        pbm.gconst(ig_('CAP05902')) = 673.9;
        pbm.gconst(ig_('CAP06002')) = 3960.0;
        pbm.gconst(ig_('CAP06102')) = 3360.0;
        pbm.gconst(ig_('CAP06202')) = 1260.0;
        pbm.gconst(ig_('CAP06302')) = 3317.9;
        pbm.gconst(ig_('CAP06402')) = 1728.0;
        pbm.gconst(ig_('CAP06502')) = 2304.0;
        pbm.gconst(ig_('CAP00103')) = 23107.0;
        pbm.gconst(ig_('CAP00203')) = 2054.0;
        pbm.gconst(ig_('CAP00303')) = 6161.9;
        pbm.gconst(ig_('CAP00403')) = 7702.3;
        pbm.gconst(ig_('CAP00503')) = 28241.9;
        pbm.gconst(ig_('CAP00603')) = 7188.9;
        pbm.gconst(ig_('CAP00703')) = 11296.8;
        pbm.gconst(ig_('CAP00803')) = 1540.5;
        pbm.gconst(ig_('CAP00903')) = 39312.0;
        pbm.gconst(ig_('CAP01003')) = 26208.0;
        pbm.gconst(ig_('CAP01103')) = 6006.0;
        pbm.gconst(ig_('CAP01203')) = 2730.0;
        pbm.gconst(ig_('CAP01303')) = 1638.0;
        pbm.gconst(ig_('CAP01403')) = 6239.7;
        pbm.gconst(ig_('CAP01503')) = 520.0;
        pbm.gconst(ig_('CAP01603')) = 520.0;
        pbm.gconst(ig_('CAP01703')) = 520.0;
        pbm.gconst(ig_('CAP01803')) = 1040.0;
        pbm.gconst(ig_('CAP01903')) = 685.4;
        pbm.gconst(ig_('CAP02003')) = 907.2;
        pbm.gconst(ig_('CAP02103')) = 902.7;
        pbm.gconst(ig_('CAP02203')) = 403.2;
        pbm.gconst(ig_('CAP02303')) = 601.8;
        pbm.gconst(ig_('CAP02403')) = 604.8;
        pbm.gconst(ig_('CAP02503')) = 504.0;
        pbm.gconst(ig_('CAP02603')) = 1562.4;
        pbm.gconst(ig_('CAP02703')) = 2328.5;
        pbm.gconst(ig_('CAP02803')) = 1461.6;
        pbm.gconst(ig_('CAP02903')) = 453.6;
        pbm.gconst(ig_('CAP03003')) = 504.0;
        pbm.gconst(ig_('CAP03103')) = 665.3;
        pbm.gconst(ig_('CAP03203')) = 1280.2;
        pbm.gconst(ig_('CAP03303')) = 1248.0;
        pbm.gconst(ig_('CAP03403')) = 2277.6;
        pbm.gconst(ig_('CAP03503')) = 124.8;
        pbm.gconst(ig_('CAP03603')) = 1310.4;
        pbm.gconst(ig_('CAP03703')) = 2901.6;
        pbm.gconst(ig_('CAP03803')) = 1254.2;
        pbm.gconst(ig_('CAP03903')) = 324.5;
        pbm.gconst(ig_('CAP04003')) = 355.7;
        pbm.gconst(ig_('CAP04103')) = 723.8;
        pbm.gconst(ig_('CAP04203')) = 1123.2;
        pbm.gconst(ig_('CAP04303')) = 1179.4;
        pbm.gconst(ig_('CAP04403')) = 393.1;
        pbm.gconst(ig_('CAP04503')) = 1085.8;
        pbm.gconst(ig_('CAP04603')) = 1447.7;
        pbm.gconst(ig_('CAP04703')) = 1447.7;
        pbm.gconst(ig_('CAP04803')) = 723.8;
        pbm.gconst(ig_('CAP04903')) = 723.8;
        pbm.gconst(ig_('CAP05003')) = 723.8;
        pbm.gconst(ig_('CAP05103')) = 736.3;
        pbm.gconst(ig_('CAP05203')) = 1809.6;
        pbm.gconst(ig_('CAP05303')) = 624.0;
        pbm.gconst(ig_('CAP05403')) = 786.2;
        pbm.gconst(ig_('CAP05503')) = 1703.5;
        pbm.gconst(ig_('CAP05603')) = 1067.0;
        pbm.gconst(ig_('CAP05703')) = 973.4;
        pbm.gconst(ig_('CAP05803')) = 1067.0;
        pbm.gconst(ig_('CAP05903')) = 730.1;
        pbm.gconst(ig_('CAP06003')) = 4158.0;
        pbm.gconst(ig_('CAP06103')) = 3528.0;
        pbm.gconst(ig_('CAP06203')) = 1323.0;
        pbm.gconst(ig_('CAP06303')) = 3594.4;
        pbm.gconst(ig_('CAP06403')) = 1814.4;
        pbm.gconst(ig_('CAP06503')) = 2419.2;
        pbm.gconst(ig_('CAP00104')) = 23107.0;
        pbm.gconst(ig_('CAP00204')) = 2054.0;
        pbm.gconst(ig_('CAP00304')) = 6161.9;
        pbm.gconst(ig_('CAP00404')) = 7702.3;
        pbm.gconst(ig_('CAP00504')) = 28241.9;
        pbm.gconst(ig_('CAP00604')) = 7188.9;
        pbm.gconst(ig_('CAP00704')) = 11296.8;
        pbm.gconst(ig_('CAP00804')) = 1540.5;
        pbm.gconst(ig_('CAP00904')) = 39312.0;
        pbm.gconst(ig_('CAP01004')) = 26208.0;
        pbm.gconst(ig_('CAP01104')) = 6006.0;
        pbm.gconst(ig_('CAP01204')) = 2730.0;
        pbm.gconst(ig_('CAP01304')) = 1638.0;
        pbm.gconst(ig_('CAP01404')) = 6239.7;
        pbm.gconst(ig_('CAP01504')) = 520.0;
        pbm.gconst(ig_('CAP01604')) = 520.0;
        pbm.gconst(ig_('CAP01704')) = 520.0;
        pbm.gconst(ig_('CAP01804')) = 1040.0;
        pbm.gconst(ig_('CAP01904')) = 718.1;
        pbm.gconst(ig_('CAP02004')) = 950.4;
        pbm.gconst(ig_('CAP02104')) = 945.6;
        pbm.gconst(ig_('CAP02204')) = 422.4;
        pbm.gconst(ig_('CAP02304')) = 630.4;
        pbm.gconst(ig_('CAP02404')) = 633.6;
        pbm.gconst(ig_('CAP02504')) = 528.0;
        pbm.gconst(ig_('CAP02604')) = 1636.8;
        pbm.gconst(ig_('CAP02704')) = 2439.4;
        pbm.gconst(ig_('CAP02804')) = 1531.2;
        pbm.gconst(ig_('CAP02904')) = 475.2;
        pbm.gconst(ig_('CAP03004')) = 528.0;
        pbm.gconst(ig_('CAP03104')) = 697.0;
        pbm.gconst(ig_('CAP03204')) = 1341.1;
        pbm.gconst(ig_('CAP03304')) = 1248.0;
        pbm.gconst(ig_('CAP03404')) = 2277.6;
        pbm.gconst(ig_('CAP03504')) = 124.8;
        pbm.gconst(ig_('CAP03604')) = 1310.4;
        pbm.gconst(ig_('CAP03704')) = 2901.6;
        pbm.gconst(ig_('CAP03804')) = 1254.2;
        pbm.gconst(ig_('CAP03904')) = 324.5;
        pbm.gconst(ig_('CAP04004')) = 355.7;
        pbm.gconst(ig_('CAP04104')) = 723.8;
        pbm.gconst(ig_('CAP04204')) = 1123.2;
        pbm.gconst(ig_('CAP04304')) = 1179.4;
        pbm.gconst(ig_('CAP04404')) = 393.1;
        pbm.gconst(ig_('CAP04504')) = 1085.8;
        pbm.gconst(ig_('CAP04604')) = 1447.7;
        pbm.gconst(ig_('CAP04704')) = 1447.7;
        pbm.gconst(ig_('CAP04804')) = 723.8;
        pbm.gconst(ig_('CAP04904')) = 723.8;
        pbm.gconst(ig_('CAP05004')) = 723.8;
        pbm.gconst(ig_('CAP05104')) = 736.3;
        pbm.gconst(ig_('CAP05204')) = 1809.6;
        pbm.gconst(ig_('CAP05304')) = 624.0;
        pbm.gconst(ig_('CAP05404')) = 786.2;
        pbm.gconst(ig_('CAP05504')) = 1703.5;
        pbm.gconst(ig_('CAP05604')) = 1067.0;
        pbm.gconst(ig_('CAP05704')) = 973.4;
        pbm.gconst(ig_('CAP05804')) = 1067.0;
        pbm.gconst(ig_('CAP05904')) = 730.1;
        pbm.gconst(ig_('CAP06004')) = 4356.0;
        pbm.gconst(ig_('CAP06104')) = 3696.0;
        pbm.gconst(ig_('CAP06204')) = 1386.0;
        pbm.gconst(ig_('CAP06304')) = 3594.4;
        pbm.gconst(ig_('CAP06404')) = 1900.8;
        pbm.gconst(ig_('CAP06504')) = 2534.4;
        pbm.gconst(ig_('CAP00105')) = 23995.8;
        pbm.gconst(ig_('CAP00205')) = 2133.0;
        pbm.gconst(ig_('CAP00305')) = 6398.9;
        pbm.gconst(ig_('CAP00405')) = 7998.6;
        pbm.gconst(ig_('CAP00505')) = 29328.1;
        pbm.gconst(ig_('CAP00605')) = 7465.3;
        pbm.gconst(ig_('CAP00705')) = 11731.3;
        pbm.gconst(ig_('CAP00805')) = 1599.7;
        pbm.gconst(ig_('CAP00905')) = 40824.0;
        pbm.gconst(ig_('CAP01005')) = 27216.0;
        pbm.gconst(ig_('CAP01105')) = 6237.0;
        pbm.gconst(ig_('CAP01205')) = 2835.0;
        pbm.gconst(ig_('CAP01305')) = 1701.0;
        pbm.gconst(ig_('CAP01405')) = 6479.7;
        pbm.gconst(ig_('CAP01505')) = 540.0;
        pbm.gconst(ig_('CAP01605')) = 540.0;
        pbm.gconst(ig_('CAP01705')) = 540.0;
        pbm.gconst(ig_('CAP01805')) = 1080.0;
        pbm.gconst(ig_('CAP01905')) = 718.1;
        pbm.gconst(ig_('CAP02005')) = 950.4;
        pbm.gconst(ig_('CAP02105')) = 945.6;
        pbm.gconst(ig_('CAP02205')) = 422.4;
        pbm.gconst(ig_('CAP02305')) = 630.4;
        pbm.gconst(ig_('CAP02405')) = 633.6;
        pbm.gconst(ig_('CAP02505')) = 528.0;
        pbm.gconst(ig_('CAP02605')) = 1636.8;
        pbm.gconst(ig_('CAP02705')) = 2439.4;
        pbm.gconst(ig_('CAP02805')) = 1531.2;
        pbm.gconst(ig_('CAP02905')) = 475.2;
        pbm.gconst(ig_('CAP03005')) = 528.0;
        pbm.gconst(ig_('CAP03105')) = 697.0;
        pbm.gconst(ig_('CAP03205')) = 1341.1;
        pbm.gconst(ig_('CAP03305')) = 1296.0;
        pbm.gconst(ig_('CAP03405')) = 2365.2;
        pbm.gconst(ig_('CAP03505')) = 129.6;
        pbm.gconst(ig_('CAP03605')) = 1360.8;
        pbm.gconst(ig_('CAP03705')) = 3013.2;
        pbm.gconst(ig_('CAP03805')) = 1302.5;
        pbm.gconst(ig_('CAP03905')) = 337.0;
        pbm.gconst(ig_('CAP04005')) = 369.4;
        pbm.gconst(ig_('CAP04105')) = 751.7;
        pbm.gconst(ig_('CAP04205')) = 1166.4;
        pbm.gconst(ig_('CAP04305')) = 1224.7;
        pbm.gconst(ig_('CAP04405')) = 408.2;
        pbm.gconst(ig_('CAP04505')) = 1127.5;
        pbm.gconst(ig_('CAP04605')) = 1503.4;
        pbm.gconst(ig_('CAP04705')) = 1503.4;
        pbm.gconst(ig_('CAP04805')) = 751.7;
        pbm.gconst(ig_('CAP04905')) = 751.7;
        pbm.gconst(ig_('CAP05005')) = 751.7;
        pbm.gconst(ig_('CAP05105')) = 764.6;
        pbm.gconst(ig_('CAP05205')) = 1879.2;
        pbm.gconst(ig_('CAP05305')) = 648.0;
        pbm.gconst(ig_('CAP05405')) = 816.5;
        pbm.gconst(ig_('CAP05505')) = 1769.0;
        pbm.gconst(ig_('CAP05605')) = 1108.1;
        pbm.gconst(ig_('CAP05705')) = 1010.9;
        pbm.gconst(ig_('CAP05805')) = 1108.1;
        pbm.gconst(ig_('CAP05905')) = 758.2;
        pbm.gconst(ig_('CAP06005')) = 4356.0;
        pbm.gconst(ig_('CAP06105')) = 3696.0;
        pbm.gconst(ig_('CAP06205')) = 1386.0;
        pbm.gconst(ig_('CAP06305')) = 3732.7;
        pbm.gconst(ig_('CAP06405')) = 1900.8;
        pbm.gconst(ig_('CAP06505')) = 2534.4;
        pbm.gconst(ig_('CAP00106')) = 22218.3;
        pbm.gconst(ig_('CAP00206')) = 1975.0;
        pbm.gconst(ig_('CAP00306')) = 5924.9;
        pbm.gconst(ig_('CAP00406')) = 7406.1;
        pbm.gconst(ig_('CAP00506')) = 27155.7;
        pbm.gconst(ig_('CAP00606')) = 6912.4;
        pbm.gconst(ig_('CAP00706')) = 10862.3;
        pbm.gconst(ig_('CAP00806')) = 1481.2;
        pbm.gconst(ig_('CAP00906')) = 37800.0;
        pbm.gconst(ig_('CAP01006')) = 25200.0;
        pbm.gconst(ig_('CAP01106')) = 5775.0;
        pbm.gconst(ig_('CAP01206')) = 2625.0;
        pbm.gconst(ig_('CAP01306')) = 1575.0;
        pbm.gconst(ig_('CAP01406')) = 1999.9;
        pbm.gconst(ig_('CAP01506')) = 166.7;
        pbm.gconst(ig_('CAP01606')) = 166.7;
        pbm.gconst(ig_('CAP01706')) = 166.7;
        pbm.gconst(ig_('CAP01806')) = 333.3;
        pbm.gconst(ig_('CAP01906')) = 685.4;
        pbm.gconst(ig_('CAP02006')) = 907.2;
        pbm.gconst(ig_('CAP02106')) = 902.7;
        pbm.gconst(ig_('CAP02206')) = 403.2;
        pbm.gconst(ig_('CAP02306')) = 601.8;
        pbm.gconst(ig_('CAP02406')) = 604.8;
        pbm.gconst(ig_('CAP02506')) = 504.0;
        pbm.gconst(ig_('CAP02606')) = 1562.4;
        pbm.gconst(ig_('CAP02706')) = 2328.5;
        pbm.gconst(ig_('CAP02806')) = 1461.6;
        pbm.gconst(ig_('CAP02906')) = 453.6;
        pbm.gconst(ig_('CAP03006')) = 504.0;
        pbm.gconst(ig_('CAP03106')) = 665.3;
        pbm.gconst(ig_('CAP03206')) = 1280.2;
        pbm.gconst(ig_('CAP03306')) = 1200.0;
        pbm.gconst(ig_('CAP03406')) = 2190.0;
        pbm.gconst(ig_('CAP03506')) = 120.0;
        pbm.gconst(ig_('CAP03606')) = 1260.0;
        pbm.gconst(ig_('CAP03706')) = 2790.0;
        pbm.gconst(ig_('CAP03806')) = 1206.0;
        pbm.gconst(ig_('CAP03906')) = 312.0;
        pbm.gconst(ig_('CAP04006')) = 342.0;
        pbm.gconst(ig_('CAP04106')) = 696.0;
        pbm.gconst(ig_('CAP04206')) = 1080.0;
        pbm.gconst(ig_('CAP04306')) = 1134.0;
        pbm.gconst(ig_('CAP04406')) = 378.0;
        pbm.gconst(ig_('CAP04506')) = 1044.0;
        pbm.gconst(ig_('CAP04606')) = 1392.0;
        pbm.gconst(ig_('CAP04706')) = 1392.0;
        pbm.gconst(ig_('CAP04806')) = 696.0;
        pbm.gconst(ig_('CAP04906')) = 696.0;
        pbm.gconst(ig_('CAP05006')) = 696.0;
        pbm.gconst(ig_('CAP05106')) = 708.0;
        pbm.gconst(ig_('CAP05206')) = 1740.0;
        pbm.gconst(ig_('CAP05306')) = 600.0;
        pbm.gconst(ig_('CAP05406')) = 756.0;
        pbm.gconst(ig_('CAP05506')) = 1638.0;
        pbm.gconst(ig_('CAP05606')) = 1026.0;
        pbm.gconst(ig_('CAP05706')) = 936.0;
        pbm.gconst(ig_('CAP05806')) = 1026.0;
        pbm.gconst(ig_('CAP05906')) = 702.0;
        pbm.gconst(ig_('CAP06006')) = 4158.0;
        pbm.gconst(ig_('CAP06106')) = 3528.0;
        pbm.gconst(ig_('CAP06206')) = 1323.0;
        pbm.gconst(ig_('CAP06306')) = 3456.2;
        pbm.gconst(ig_('CAP06406')) = 1814.4;
        pbm.gconst(ig_('CAP06506')) = 2419.2;
        pbm.gconst(ig_('MXD00102')) = 56200.0;
        pbm.gconst(ig_('MXD00202')) = 245714.0;
        pbm.gconst(ig_('MXD00502')) = 1201600.0;
        pbm.gconst(ig_('MXD00702')) = 78000.0;
        pbm.gconst(ig_('MXD00802')) = 180000.0;
        pbm.gconst(ig_('MXD00902')) = 1391199.0;
        pbm.gconst(ig_('MXD01002')) = 450000.0;
        pbm.gconst(ig_('MXD00103')) = 98200.0;
        pbm.gconst(ig_('MXD00203')) = 377143.0;
        pbm.gconst(ig_('MXD00303')) = 142824.0;
        pbm.gconst(ig_('MXD00403')) = 16000.0;
        pbm.gconst(ig_('MXD00503')) = 1984799.0;
        pbm.gconst(ig_('MXD00603')) = 34200.0;
        pbm.gconst(ig_('MXD00703')) = 155200.0;
        pbm.gconst(ig_('MXD00803')) = 258000.0;
        pbm.gconst(ig_('MXD00903')) = 2414999.0;
        pbm.gconst(ig_('MXD01003')) = 706000.0;
        pbm.gconst(ig_('MND00104')) = 89091.0;
        pbm.gconst(ig_('MXD00104')) = 232200.0;
        pbm.gconst(ig_('MXD00204')) = 491428.0;
        pbm.gconst(ig_('MND00304')) = 13640.0;
        pbm.gconst(ig_('MXD00304')) = 288706.0;
        pbm.gconst(ig_('MXD00404')) = 20400.0;
        pbm.gconst(ig_('MND00504')) = 433908.0;
        pbm.gconst(ig_('MXD00504')) = 2669599.0;
        pbm.gconst(ig_('MXD00604')) = 54200.0;
        pbm.gconst(ig_('MND00704')) = 99890.0;
        pbm.gconst(ig_('MXD00704')) = 214400.0;
        pbm.gconst(ig_('MND00804')) = 149985.0;
        pbm.gconst(ig_('MXD00804')) = 336000.0;
        pbm.gconst(ig_('MND00904')) = 68537.0;
        pbm.gconst(ig_('MXD00904')) = 3635799.0;
        pbm.gconst(ig_('MND01004')) = 118071.0;
        pbm.gconst(ig_('MXD01004')) = 966000.0;
        pbm.gconst(ig_('MND00105')) = 129087.0;
        pbm.gconst(ig_('MXD00105')) = 366200.0;
        pbm.gconst(ig_('MND00205')) = 150791.0;
        pbm.gconst(ig_('MXD00205')) = 492571.0;
        pbm.gconst(ig_('MND00305')) = 200788.0;
        pbm.gconst(ig_('MXD00305')) = 486353.0;
        pbm.gconst(ig_('MND00405')) = 11699.0;
        pbm.gconst(ig_('MXD00405')) = 25400.0;
        pbm.gconst(ig_('MND00505')) = 884246.0;
        pbm.gconst(ig_('MXD00505')) = 3241599.0;
        pbm.gconst(ig_('MND00605')) = 24096.0;
        pbm.gconst(ig_('MXD00605')) = 68200.0;
        pbm.gconst(ig_('MND00705')) = 125187.0;
        pbm.gconst(ig_('MXD00705')) = 273600.0;
        pbm.gconst(ig_('MND00805')) = 179982.0;
        pbm.gconst(ig_('MXD00805')) = 414000.0;
        pbm.gconst(ig_('MND00905')) = 1849407.0;
        pbm.gconst(ig_('MXD00905')) = 4798598.0;
        pbm.gconst(ig_('MND01005')) = 334408.0;
        pbm.gconst(ig_('MXD01005')) = 1226000.0;
        pbm.gconst(ig_('MND00106')) = 169083.0;
        pbm.gconst(ig_('MXD00106')) = 500200.0;
        pbm.gconst(ig_('MND00206')) = 150791.0;
        pbm.gconst(ig_('MXD00206')) = 493714.0;
        pbm.gconst(ig_('MND00306')) = 233184.0;
        pbm.gconst(ig_('MXD00306')) = 702823.0;
        pbm.gconst(ig_('MND00406')) = 11699.0;
        pbm.gconst(ig_('MXD00406')) = 30800.0;
        pbm.gconst(ig_('MND00506')) = 1301415.0;
        pbm.gconst(ig_('MXD00506')) = 3817599.0;
        pbm.gconst(ig_('MND00606')) = 30096.0;
        pbm.gconst(ig_('MXD00606')) = 80200.0;
        pbm.gconst(ig_('MND00706')) = 125187.0;
        pbm.gconst(ig_('MXD00706')) = 345200.0;
        pbm.gconst(ig_('MND00806')) = 209979.0;
        pbm.gconst(ig_('MXD00806')) = 492000.0;
        pbm.gconst(ig_('MND00906')) = 1849407.0;
        pbm.gconst(ig_('MXD00906')) = 6141396.0;
        pbm.gconst(ig_('MND01006')) = 644932.0;
        pbm.gconst(ig_('MXD01006')) = 1486000.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                   ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-AN-163-488';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

