function varargout = HS91(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS91
%    *********
% 
%    A time-optimal heat conduction problem.
% 
%    Source: problem 91 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, September 1991.
%      SAVEs removed December 3rd 2014
% 
%    classification = 'QOR2-MN-5-1'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS91';

switch(action)

    case 'setup'

    pb.name      = 'HS91';
    pb.sifpbname = 'HS91';
    pbm.name     = 'HS91';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 5;
        v_('EPS') = 0.01;
        v_('EPSSQR') = v_('EPS')*v_('EPS');
        v_('-EPSSQR') = -1.0*v_('EPSSQR');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','CON',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON';
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('CON')) = v_('-EPSSQR');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('2'):v_('N')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = 0.5;
            end
        end
        for I=v_('2'):v_('2'):v_('N')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = -0.5;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = -0.5;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'SQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'H',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['O',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQR';
            ielftype(ie) = iet_('SQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'H';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'H';
        ielftype(ie) = iet_('H');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['O',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
        end
        ig = ig_('CON');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('H');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QOR2-MN-5-1';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
