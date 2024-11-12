function varargout = HS104(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS104
%    *********
% 
%    Source: problem 104 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'C-COOR2-AN-8-5'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS104';

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
        v_('N') = 8;
        v_('1') = 1;
        v_('4') = 4;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C1';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0e-1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0e-1;
        end
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C2';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0e-1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0e-1;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0e-1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0e-1;
        end
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C3';
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C4';
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C5';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
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
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0e+0*ones(ngrp,1);
        pbm.gconst(ig_('OBJ')) = -1.0e+1;
        pbm.gconst(ig_('C5')) = -9.0e+0;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        grange(ig_('C5')) = 3.2e+0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 1.0e-1*ones(pb.n,1);
        pb.xupper = 1.0e+1*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 6.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 6.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 3.0;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 3.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.4;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.4;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 0.2;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 0.2;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 6.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 6.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 6.0;
        else
            pb.y0(find(pbm.congrps==ig('X6')),1) = 6.0;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 1.0;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig('X8')),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
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
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.67;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.67;
        ename = 'OE2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.67;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.67;
        ename = 'C1E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'C2E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'C3E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        ename = 'C3E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.71;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        ename = 'C3E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.3;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'C4E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        ename = 'C4E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.71;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        ename = 'C4E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-1,1.0e+1,[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.3;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0e-1;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OE2');
        pbm.grelw{ig}(posel) = 4.0e-1;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.88e-2;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C2E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.88e-2;
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C3E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C3E2');
        pbm.grelw{ig}(posel) = 2.0e+0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C3E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.88e-2;
        ig = ig_('C4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C4E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C4E2');
        pbm.grelw{ig}(posel) = 2.0e+0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C4E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.88e-2;
        ig = ig_('C5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0e-1;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OE2');
        pbm.grelw{ig}(posel) = 4.0e-1;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               3.9511634396
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-8-5';
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
        varargout{1} = (EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^pbm.elpar{iel_}(2));
        if(nargout>1)
            g_(1,1) =...
                  pbm.elpar{iel_}(1)*(EV_(1)^(pbm.elpar{iel_}(1)-1.0))*(EV_(2)^pbm.elpar{iel_}(2));
            g_(2,1) =...
                  pbm.elpar{iel_}(2)*(EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^(pbm.elpar{iel_}(2)-1.0));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) =...
                      pbm.elpar{iel_}(1)*(EV_(1)^(pbm.elpar{iel_}(1)-2.0))*(pbm.elpar{iel_}(1)-1.0)*(EV_(2)^pbm.elpar{iel_}(2));
                H_(1,2) =...
                      pbm.elpar{iel_}(1)*(EV_(1)^(pbm.elpar{iel_}(1)-1.0))*pbm.elpar{iel_}(2)*(EV_(2)^(pbm.elpar{iel_}(2)-1.0));
                H_(2,1) = H_(1,2);
                H_(2,2) =...
                      pbm.elpar{iel_}(2)*(pbm.elpar{iel_}(2)-1.0)*(EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^(pbm.elpar{iel_}(2)-2.0));
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

