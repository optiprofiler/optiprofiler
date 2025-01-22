function varargout = HS59(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS59
%    *********
% 
%    Source: problem 59 (corrected) in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn and Ph.L. Toint, November 1990.
% 
%    classification = 'C-COQR2-AN-2-3'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS59';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 3.8112;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 6.8306;
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON1';
        [ig,ig_] = s2mpjlib('ii','CON2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','CON3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -5.0;
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
        pbm.gconst(ig_('OBJ')) = 75.196;
        pbm.gconst(ig_('CON1')) = 700.0;
        pbm.gconst(ig_('CON3')) = -275.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('X1')) = 75.0;
        pb.xupper(ix_('X2')) = 65.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 90.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 90.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 10.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 10.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePOWER',iet_);
        elftv{it}{1} = 'V1';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'ePOWPR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        [it,iet_] = s2mpjlib( 'ii', 'eEXPR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'eSHPOW',iet_);
        elftv{it}{1} = 'V1';
        elftp{it}{1} = 'P';
        elftp{it}{2} = 'S';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.0;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.0;
        ename = 'E9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSHPOW';
        ielftype(ie) = iet_('eSHPOW');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('S',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        ename = 'E10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        ename = 'E13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E14';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        ename = 'E15';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEXPR';
        ielftype(ie) = iet_('eEXPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E16';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWPR';
        ielftype(ie) = iet_('ePOWPR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePOWER';
        ielftype(ie) = iet_('ePOWER');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSHPOW';
        ielftype(ie) = iet_('eSHPOW');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('S',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 50.0;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0020567;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = -1.0345e-5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.030234;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.28134e-3;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.266e-7;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        pbm.grelw{ig}(posel) = -0.25645;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0034604;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        pbm.grelw{ig}(posel) = -1.3514e-5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 28.106;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        pbm.grelw{ig}(posel) = 5.2375e-6;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.3e-8;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E12');
        pbm.grelw{ig}(posel) = -7.0e-10;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -3.405e-4;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E14');
        pbm.grelw{ig}(posel) = 1.6638e-6;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.8673;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E16');
        pbm.grelw{ig}(posel) = -3.5256e-5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E17');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.12694;
        ig = ig_('CON1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('CON2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E17');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.008;
        ig = ig_('CON3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -7.8027894
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COQR2-AN-2-3';
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

    case 'ePOWER'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^pbm.elpar{iel_}(1);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(1)^(pbm.elpar{iel_}(1)-1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) =...
                      pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0)*EV_(1)^(pbm.elpar{iel_}(1)-2.0);
                varargout{3} = H_;
            end
        end

    case 'ePOWPR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V1P1 = EV_(1)^pbm.elpar{iel_}(1);
        V2P2 = EV_(2)^pbm.elpar{iel_}(2);
        P1M1 = pbm.elpar{iel_}(1)-1.0;
        P2M1 = pbm.elpar{iel_}(2)-1.0;
        EH1 = max(0.0e0,pbm.elpar{iel_}(1)-2.0);
        EH2 = max(0.0e0,pbm.elpar{iel_}(2)-2.0);
        varargout{1} = V1P1*V2P2;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(1)^P1M1*V2P2;
            g_(2,1) = V1P1*pbm.elpar{iel_}(2)*EV_(2)^P2M1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.elpar{iel_}(1)*P1M1*EV_(1)^EH1*V2P2;
                H_(1,2) = pbm.elpar{iel_}(1)*EV_(1)^P1M1*pbm.elpar{iel_}(2)*EV_(2)^P2M1;
                H_(2,1) = H_(1,2);
                H_(2,2) = V1P1*pbm.elpar{iel_}(2)*P2M1*EV_(2)^EH2;
                varargout{3} = H_;
            end
        end

    case 'eEXPR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E = exp(0.0005*EV_(1)*EV_(2));
        varargout{1} = E;
        if(nargout>1)
            g_(1,1) = 0.0005*EV_(2)*E;
            g_(2,1) = 0.0005*EV_(1)*E;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.5e-7*EV_(2)*EV_(2)*E;
                H_(1,2) = 2.5e-7*EV_(1)*EV_(2)*E+0.0005*E;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.5e-7*EV_(1)*EV_(1)*E;
                varargout{3} = H_;
            end
        end

    case 'eSHPOW'

        EV_  = varargin{1};
        iel_ = varargin{2};
        VMS = EV_(1)-pbm.elpar{iel_}(2);
        varargout{1} = VMS^pbm.elpar{iel_}(1);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*VMS^(pbm.elpar{iel_}(1)-1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) =...
                      pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0)*VMS^(pbm.elpar{iel_}(1)-2.0);
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

