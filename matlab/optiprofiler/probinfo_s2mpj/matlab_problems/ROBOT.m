function varargout = ROBOT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    This program solves the displacement optimization problem in 
%    REDUNDANT ROBOTS. A redundant robot is one which has more links than 
%    the dimensions it moves in.  Because this redundancy allows almost
%    infinite combinations of joint angles for a particular orientation of
%    end-effector of a robot, choosing an optimum combination has always been a
%    problem of research. 
%    The ROBOT considered here is a 7 link robot moving in 2 dimensional space.
% 
%    Source: an exercize for L. Watson course on LANCELOT in the Spring 1993.
%    B.Benhabib, R.G.Fenton and A.A.Goldberg, 
%    "Analytical trajectory optimization of seven degrees of freedom redundant
%    robot",  
%    Transactions of the Canadian Society for Mechanical Engineering,
%    vol.11(4), 1987, pp 197-200.
% 
%    SIF input: Manish Sabu at Virginia Tech., Spring 1993.
%               Minor modifications by Ph. L. Toint, April 1993.
% 
%    classification = 'C-CQOR2-MY-14-2'
% 
%  This segment describes the initial values of angles (by THnIN)
%   and final position of the end effector (by XPOS and YPOS)
%   these values can be changed here according to the needs of the user.
%  The segment also defines the upper and lower bounds of the various joint 
%   angles (by HIGH and DOWN)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ROBOT';

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
        v_('TH1IN') = 0.0;
        v_('TH2IN') = 0.0;
        v_('TH3IN') = 0.0;
        v_('TH4IN') = 0.0;
        v_('TH5IN') = 0.0;
        v_('TH6IN') = 0.0;
        v_('TH7IN') = 0.0;
        v_('XPOS') = 4.0;
        v_('YPOS') = 4.0;
        v_('HIGH') = 2.356194;
        v_('DOWN') = -2.356194;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','TH1',ix_);
        pb.xnames{iv} = 'TH1';
        [iv,ix_] = s2mpjlib('ii','TH2',ix_);
        pb.xnames{iv} = 'TH2';
        [iv,ix_] = s2mpjlib('ii','TH3',ix_);
        pb.xnames{iv} = 'TH3';
        [iv,ix_] = s2mpjlib('ii','TH4',ix_);
        pb.xnames{iv} = 'TH4';
        [iv,ix_] = s2mpjlib('ii','TH5',ix_);
        pb.xnames{iv} = 'TH5';
        [iv,ix_] = s2mpjlib('ii','TH6',ix_);
        pb.xnames{iv} = 'TH6';
        [iv,ix_] = s2mpjlib('ii','TH7',ix_);
        pb.xnames{iv} = 'TH7';
        [iv,ix_] = s2mpjlib('ii','TH1I',ix_);
        pb.xnames{iv} = 'TH1I';
        [iv,ix_] = s2mpjlib('ii','TH2I',ix_);
        pb.xnames{iv} = 'TH2I';
        [iv,ix_] = s2mpjlib('ii','TH3I',ix_);
        pb.xnames{iv} = 'TH3I';
        [iv,ix_] = s2mpjlib('ii','TH4I',ix_);
        pb.xnames{iv} = 'TH4I';
        [iv,ix_] = s2mpjlib('ii','TH5I',ix_);
        pb.xnames{iv} = 'TH5I';
        [iv,ix_] = s2mpjlib('ii','TH6I',ix_);
        pb.xnames{iv} = 'TH6I';
        [iv,ix_] = s2mpjlib('ii','TH7I',ix_);
        pb.xnames{iv} = 'TH7I';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CONSTR1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CONSTR1';
        [ig,ig_] = s2mpjlib('ii','CONSTR2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CONSTR2';
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
        pbm.gconst(ig_('CONSTR1')) = v_('XPOS');
        pbm.gconst(ig_('CONSTR2')) = v_('YPOS');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('TH1'),1) = v_('DOWN');
        pb.xupper(ix_('TH1')) = v_('HIGH');
        pb.xlower(ix_('TH2'),1) = v_('DOWN');
        pb.xupper(ix_('TH2')) = v_('HIGH');
        pb.xlower(ix_('TH3'),1) = v_('DOWN');
        pb.xupper(ix_('TH3')) = v_('HIGH');
        pb.xlower(ix_('TH4'),1) = v_('DOWN');
        pb.xupper(ix_('TH4')) = v_('HIGH');
        pb.xlower(ix_('TH5'),1) = v_('DOWN');
        pb.xupper(ix_('TH5')) = v_('HIGH');
        pb.xlower(ix_('TH6'),1) = v_('DOWN');
        pb.xupper(ix_('TH6')) = v_('HIGH');
        pb.xlower(ix_('TH7'),1) = v_('DOWN');
        pb.xupper(ix_('TH7')) = v_('HIGH');
        pb.xlower(ix_('TH1I'),1) = v_('TH1IN');
        pb.xupper(ix_('TH1I'),1) = v_('TH1IN');
        pb.xlower(ix_('TH2I'),1) = v_('TH2IN');
        pb.xupper(ix_('TH2I'),1) = v_('TH2IN');
        pb.xlower(ix_('TH3I'),1) = v_('TH3IN');
        pb.xupper(ix_('TH3I'),1) = v_('TH3IN');
        pb.xlower(ix_('TH4I'),1) = v_('TH4IN');
        pb.xupper(ix_('TH4I'),1) = v_('TH4IN');
        pb.xlower(ix_('TH5I'),1) = v_('TH5IN');
        pb.xupper(ix_('TH5I'),1) = v_('TH5IN');
        pb.xlower(ix_('TH6I'),1) = v_('TH6IN');
        pb.xupper(ix_('TH6I'),1) = v_('TH6IN');
        pb.xlower(ix_('TH7I'),1) = v_('TH7IN');
        pb.xupper(ix_('TH7I'),1) = v_('TH7IN');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('TH1'),1) = 0.0;
        pb.x0(ix_('TH2'),1) = 0.0;
        pb.x0(ix_('TH3'),1) = 0.0;
        pb.x0(ix_('TH4'),1) = 0.0;
        pb.x0(ix_('TH5'),1) = 0.0;
        pb.x0(ix_('TH6'),1) = 0.0;
        pb.x0(ix_('TH7'),1) = 0.0;
        pb.x0(ix_('TH1I'),1) = 0.0;
        pb.x0(ix_('TH2I'),1) = 0.0;
        pb.x0(ix_('TH3I'),1) = 0.0;
        pb.x0(ix_('TH4I'),1) = 0.0;
        pb.x0(ix_('TH5I'),1) = 0.0;
        pb.x0(ix_('TH6I'),1) = 0.0;
        pb.x0(ix_('TH7I'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eCOSTH',iet_);
        elftv{it}{1} = 'THETAC';
        [it,iet_] = s2mpjlib( 'ii', 'eSINTH',iet_);
        elftv{it}{1} = 'THETAS';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'TH1SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH1I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH2SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH2I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH3SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH3I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH4SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH4I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH5SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH5I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH6SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH6I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'TH7SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eISQ';
        ielftype(ie) = iet_('eISQ');
        vname = 'TH7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TH7I';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C1TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C2TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C3TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C4TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C5TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C6TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C7TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOSTH';
        ielftype(ie) = iet_('eCOSTH');
        vname = 'TH7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAC',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S1TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S2TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S3TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S4TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S5TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S6TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'S7TH';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINTH';
        ielftype(ie) = iet_('eSINTH');
        vname = 'TH7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('THETAS',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TH1SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('TH2SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TH3SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('TH4SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TH5SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('TH6SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TH7SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('CONSTR1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C2TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C3TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C4TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C5TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C6TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C7TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.5;
        ig = ig_('CONSTR2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('S1TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('S2TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('S3TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('S4TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('S5TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('S6TH');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('S7TH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.5;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            5.46283877
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MY-14-2';
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

    case 'eISQ'

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

    case 'eCOSTH'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP = cos(EV_(1));
        varargout{1} = TMP;
        if(nargout>1)
            g_(1,1) = -sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -TMP;
                varargout{3} = H_;
            end
        end

    case 'eSINTH'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP = sin(EV_(1));
        varargout{1} = TMP;
        if(nargout>1)
            g_(1,1) = cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -TMP;
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

