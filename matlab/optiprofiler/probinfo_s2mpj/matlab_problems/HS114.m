function varargout = HS114(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS114
%    *********
% 
%    An alkylation process problem.
% 
%    Source: problem 114 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: J.M. Collin, Jan 1990.
% 
%    classification = 'C-CQOR2-MY-10-11'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS114';

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
        v_('N') = 10;
        v_('1') = 1;
        v_('A') = 0.99;
        v_('B') = 0.9;
        v_('-A') = -1.0*v_('A');
        v_('-B') = -1.0*v_('B');
        v_('INVA') = 1.0/v_('A');
        v_('INVB') = 1.0/v_('B');
        v_('-INVA') = -1.0*v_('INVA');
        v_('-INVB') = -1.0*v_('INVB');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.035;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 10.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 3.36;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 5.04;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = -0.222;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = v_('-B');
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 3.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = v_('-A');
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 0.222;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = v_('INVB');
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = v_('INVA');
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -3.0;
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.12;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = v_('-A');
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.098;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 0.325;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = v_('-A');
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.12;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = v_('INVA');
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -1.098;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -0.325;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = v_('INVA');
        [ig,ig_] = s2mpjlib('ii','C9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C9';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.22;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C10';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -1.0;
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
        pbm.gconst(ig_('C1')) = -35.82;
        pbm.gconst(ig_('C2')) = 133.0;
        pbm.gconst(ig_('C3')) = 35.82;
        pbm.gconst(ig_('C4')) = -133.0;
        pbm.gconst(ig_('C6')) = -57.425;
        pbm.gconst(ig_('C8')) = 57.425;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.00001;
        pb.xupper(ix_('X1')) = 2000.0;
        pb.xlower(ix_('X2'),1) = 0.00001;
        pb.xupper(ix_('X2')) = 16000.0;
        pb.xlower(ix_('X3'),1) = 0.00001;
        pb.xupper(ix_('X3')) = 120.0;
        pb.xlower(ix_('X4'),1) = 0.00001;
        pb.xupper(ix_('X4')) = 5000.0;
        pb.xlower(ix_('X5'),1) = 0.00001;
        pb.xupper(ix_('X5')) = 2000.0;
        pb.xlower(ix_('X6'),1) = 85.0;
        pb.xupper(ix_('X6')) = 93.0;
        pb.xlower(ix_('X7'),1) = 90.0;
        pb.xupper(ix_('X7')) = 95.0;
        pb.xlower(ix_('X8'),1) = 3.0;
        pb.xupper(ix_('X8')) = 12.0;
        pb.xlower(ix_('X9'),1) = 1.2;
        pb.xupper(ix_('X9')) = 4.0;
        pb.xlower(ix_('X10'),1) = 145.0;
        pb.xupper(ix_('X10')) = 162.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 1745.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 1745.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 12000.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 12000.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 110.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 110.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 3048.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 3048.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 1974.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 1974.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 89.2;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 89.2;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 92.8;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 92.8;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 8.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8')),1) = 8.0;
        end
        if(isKey(ix_,'X9'))
            pb.x0(ix_('X9'),1) = 3.6;
        else
            pb.y0(find(pbm.congrps==ig_('X9')),1) = 3.6;
        end
        if(isKey(ix_,'X10'))
            pb.x0(ix_('X10'),1) = 145.0;
        else
            pb.y0(find(pbm.congrps==ig_('X10')),1) = 145.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eWSQ',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eRAP1',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'S';
        elftv{it}{3} = 'T';
        elftv{it}{4} = 'U';
        elftp{it}{1} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eRAP2',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OE1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.063;
        ename = 'CE1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.13167;
        ename = 'CE2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.00667;
        ename = 'CE3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eWSQ';
        ielftype(ie) = iet_('eWSQ');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.038;
        ename = 'CE4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.13167;
        ename = 'CE5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD2';
        ielftype(ie) = iet_('ePROD2');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00667;
        ename = 'CE6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eWSQ';
        ielftype(ie) = iet_('eWSQ');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.038;
        ename = 'CE7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eRAP1';
        ielftype(ie) = iet_('eRAP1');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('S',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('T',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('W',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 98000.0;
        ename = 'CE8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eRAP2';
        ielftype(ie) = iet_('eRAP2');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('CE2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C7');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('CE5');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C10');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C11');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('CE8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -1768.80696
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MY-10-11';
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

    case 'eWSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)^2;
        if(nargout>1)
            g_(1,1) = 2.0*pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*pbm.elpar{iel_}(1);
                varargout{3} = H_;
            end
        end

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2);
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2)^2;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2)^2;
            g_(2,1) = 2.0*pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 2.0*pbm.elpar{iel_}(1)*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*pbm.elpar{iel_}(1)*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eRAP1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        WX = pbm.elpar{iel_}(1)*EV_(1);
        DENOM = EV_(2)*EV_(3)+1000.0*EV_(4);
        varargout{1} = WX/DENOM;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)/DENOM;
            g_(2,1) = -WX*EV_(3)/DENOM^2;
            g_(3,1) = -WX*EV_(2)/DENOM^2;
            g_(4,1) = -1000.0*WX/DENOM^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = -pbm.elpar{iel_}(1)*EV_(3)/DENOM^2;
                H_(2,1) = H_(1,2);
                H_(1,3) = -pbm.elpar{iel_}(1)*EV_(2)/DENOM^2;
                H_(3,1) = H_(1,3);
                H_(1,4) = -1000.0*pbm.elpar{iel_}(1)/DENOM^2;
                H_(4,1) = H_(1,4);
                H_(2,2) = 2.0*WX*EV_(3)*EV_(3)/DENOM^3;
                H_(2,3) = 2.0*WX*EV_(2)*EV_(3)/DENOM^3-WX/DENOM^2;
                H_(3,2) = H_(2,3);
                H_(2,4) = 2000.0*WX*EV_(3)/DENOM^3;
                H_(4,2) = H_(2,4);
                H_(3,3) = 2.0*WX*EV_(2)*EV_(2)/DENOM^3;
                H_(3,4) = 2000.0*WX*EV_(2)/DENOM^3;
                H_(4,3) = H_(3,4);
                H_(4,4) = 2000000.0*WX/DENOM^3;
                varargout{3} = H_;
            end
        end

    case 'eRAP2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(2,3) = U_(2,3)+1;
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)/IV_(2);
        if(nargout>1)
            g_(1,1) = 1.0/IV_(2);
            g_(2,1) = -IV_(1)/IV_(2)^2;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0/IV_(2)^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*IV_(1)/IV_(2)^3;
                varargout{3} = U_.'*H_*U_;
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

