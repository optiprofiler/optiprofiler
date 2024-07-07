function varargout = HS85(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS85
%    *********
% 
%    The problem is to optimize the net profit of an hypothetical
%    wood-pulp plant. The constraints include the usual material
%    and energy balances as well as several empirical equations.
% 
%    Source: problem 85 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, September 1991.
%      SAVEs removed December 3rd 2014
% 
%    classification = 'OOI2-MN-5-21'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS85';

switch(action)

    case 'setup'

    pb.name      = 'HS85';
    pb.sifpbname = 'HS85';
    pbm.name     = 'HS85';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 5;
        v_('A2') = 17.505;
        v_('A3') = 11.275;
        v_('A4') = 214.228;
        v_('A5') = 7.458;
        v_('A6') = 0.961;
        v_('A7') = 1.612;
        v_('A8') = 0.146;
        v_('A9') = 107.99;
        v_('A10') = 922.693;
        v_('A11') = 926.832;
        v_('A12') = 18.766;
        v_('A13') = 1072.163;
        v_('A14') = 8961.448;
        v_('A15') = 0.063;
        v_('A16') = 71084.33;
        v_('A17') = 2802713.0;
        v_('B2') = 1053.6667;
        v_('B3') = 35.03;
        v_('B4') = 665.585;
        v_('B5') = 584.463;
        v_('B6') = 265.916;
        v_('B7') = 7.046;
        v_('B8') = 0.222;
        v_('B9') = 273.366;
        v_('B10') = 1286.105;
        v_('B11') = 1444.046;
        v_('B12') = 537.141;
        v_('B13') = 3247.039;
        v_('B14') = 26844.086;
        v_('B15') = 0.386;
        v_('B16') = 140000.0;
        v_('B17') = 12146108.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('17') = 17;
        v_('19') = 19;
        v_('20') = 20;
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
        [ig,ig_] = s2xlib('ii','CON0',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON0';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.5;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('1'):v_('20')
            [ig,ig_] = s2xlib('ii',['CON',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CON',int2str(I)];
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('OBJ')) = 0.1365;
        pbm.gconst(ig_('CON1')) = 213.1;
        for I=v_('2'):v_('17')
            pbm.gconst(ig_(['CON',int2str(I)])) = v_(['A',int2str(I)]);
        end
        pbm.gconst(ig_('CON19')) = -21.0;
        pbm.gconst(ig_('CON20')) = 110.6;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        for I=v_('2'):v_('17')
            v_('DIF') = v_(['B',int2str(I)])-v_(['A',int2str(I)]);
            grange(ig_(['CON',int2str(I)])) = v_('DIF');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 7.044148e+2;
        pb.xupper(ix_('X1')) = 9.063855e+2;
        pb.xlower(ix_('X2'),1) = 6.86e+1;
        pb.xupper(ix_('X2')) = 2.8888e+2;
        pb.xupper(ix_('X3')) = 1.3475e+2;
        pb.xlower(ix_('X4'),1) = 1.930e+2;
        pb.xupper(ix_('X4')) = 2.870966e+2;
        pb.xlower(ix_('X5'),1) = 2.50e+1;
        pb.xupper(ix_('X5')) = 8.41988e+1;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 9.0e+2;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 9.0e+2;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 8.0e+1;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 8.0e+1;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 1.15e+2;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 1.15e+2;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 2.67e+2;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 2.67e+2;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 2.7e+1;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 2.7e+1;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'Y',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftp{it}{1} = 'PI';
        [it,iet_] = s2xlib( 'ii', 'C',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftp{it}{1} = 'PI';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('19')
            v_('PI') = I;
            v_('PI') = 0.01+v_('PI');
            ename = ['C',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'C';
            ielftype(ie) = iet_('C');
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
            [~,posep] = ismember('PI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('PI');
        end
        for I=v_('1'):v_('20')
            v_('PI') = I;
            v_('PI') = 0.01+v_('PI');
            ename = ['Y',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'Y';
            ielftype(ie) = iet_('Y');
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
            [~,posep] = ismember('PI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('PI');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Y17');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -5.843e-7;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Y14');
        pbm.grelw{ig}(posel) = 1.17e-4;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Y13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.358e-5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Y16');
        pbm.grelw{ig}(posel) = 1.502e-6;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('Y12');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0321;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('Y5');
        pbm.grelw{ig}(posel) = 0.00423;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0e-4;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C19');
        pbm.grelw{ig}(posel) = 37.48;
        for I=v_('1'):v_('20')
            ig = ig_(['CON',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Y',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOI2-MN-5-21';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
