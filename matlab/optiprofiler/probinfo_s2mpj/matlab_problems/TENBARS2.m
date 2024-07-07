function varargout = TENBARS2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TENBARS2
%    *********
%    The "ten bar truss" structural optimization problem,
%    version P2.
% 
%    The problem is to minimize the cross section areas of the bars
%    in the structure
% 
%      /|
%      /|>o------o------o
%      /|  \    /|\    /|
%           \  / | \  / |
%            \/  |  \/  |
%            /\  |  /\  |
%           /  \ | /  \ |
%      /|  /    \|/    \|
%      /|>o------o------o
%      /|
% 
%    submitted to vertical forces of equal magnitude (P0) applied at
%    the two free lower nodes, subject to limits of nodal displacements.
% 
%    Source:
%    K. Svanberg,
%    private communication,  August 1990.
%    See also
%    K. Svanberg,
%    "On local and global minima in structural optimization",
%    in "New directions in optimum structural design" (Atrek, Ragsdell
%    and Zienkiwewicz, eds.), Wiley, 1984.
% 
%    SIF input: Ph. Toint, August 1990.
% 
%    classification = 'LOR2-MY-18-8'
% 
%    Constants
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TENBARS2';

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
        v_('1') = 1;
        v_('2') = 2;
        v_('8') = 8;
        v_('10') = 10;
        v_('SQ2') = sqrt(2.0);
        v_('SQ8') = sqrt(8.0);
        v_('1/SQ8') = 1.0/v_('SQ8');
        v_('-1/SQ8') = -1.0*v_('1/SQ8');
        v_('C0') = 2.53106;
        v_('-P0') = -589.884;
        v_('C0SQ2') = v_('C0')*v_('SQ2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('8')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        for I=v_('1'):v_('10')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0SQ2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0SQ2');
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0SQ2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0SQ2');
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0SQ2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0SQ2');
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0SQ2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0SQ2');
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('C0')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('C0');
        end
        for I=v_('1'):v_('8')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(I)];
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        pbm.gconst(ig_('C4')) = v_('-P0');
        pbm.gconst(ig_('C8')) = v_('-P0');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for I=v_('1'):v_('10')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.645;
        end
        for I=v_('2'):v_('2'):v_('8')
            pb.xlower(ix_(['U',int2str(I)]),1) = -50.8;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXU',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'eXUPV',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'U';
        elftv{it}{3} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eXUMV',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'U';
        elftv{it}{3} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eXBIG',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'U';
        elftv{it}{3} = 'V';
        elftv{it}{4} = 'W';
        elftv{it}{5} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'EA';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXU';
        ielftype(ie) = iet_('eXU');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXUPV';
        ielftype(ie) = iet_('eXUPV');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EC';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eXUMV';
            ielftype(ie) = iet_('eXUMV');
        end
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'ED';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXBIG';
        ielftype(ie) = iet_('eXBIG');
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EE';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eXUMV';
            ielftype(ie) = iet_('eXUMV');
        end
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EF';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eXUMV';
            ielftype(ie) = iet_('eXUMV');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EG';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXU';
        ielftype(ie) = iet_('eXU');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EH';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXBIG';
        ielftype(ie) = iet_('eXBIG');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EI';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eXUMV';
            ielftype(ie) = iet_('eXUMV');
        end
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJ';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eXUMV';
            ielftype(ie) = iet_('eXUMV');
        end
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('U',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EA');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EC');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EB');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ED');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EB');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EE');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ED');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/SQ8');
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EI');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EG');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EF');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EE');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EC');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EH');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C7');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ED');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EI');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ED');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/SQ8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN-A             2302.55
% LO SOLTN-B             2277.9458
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MY-18-8';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eXU'

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

    case 'eXUPV'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,2) = U_(1,2)+1;
        U_(1,3) = U_(1,3)+1;
        U_(2,1) = U_(2,1)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(2)*IV_(1);
        if(nargout>1)
            g_(2,1) = IV_(1);
            g_(1,1) = IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = 1.0;
                H_(1,2) = H_(2,1);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXUMV'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,2) = U_(1,2)+1;
        U_(1,3) = U_(1,3)-1;
        U_(2,1) = U_(2,1)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(2)*IV_(1);
        if(nargout>1)
            g_(2,1) = IV_(1);
            g_(1,1) = IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = 1.0;
                H_(1,2) = H_(2,1);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eXBIG'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,5);
        U_(1,2) = U_(1,2)+1;
        U_(1,3) = U_(1,3)+1;
        U_(1,4) = U_(1,4)-1;
        U_(1,5) = U_(1,5)-1;
        U_(2,1) = U_(2,1)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(2)*IV_(1);
        if(nargout>1)
            g_(2,1) = IV_(1);
            g_(1,1) = IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = 1.0;
                H_(1,2) = H_(2,1);
                varargout{3} = U_.'*H_*U_;
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

