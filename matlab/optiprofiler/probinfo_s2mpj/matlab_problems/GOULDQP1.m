function varargout = GOULDQP1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : GOULDQP1
%    *********
% 
%    Source: problem 118 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981, as modified by N.I.M. Gould in "An algorithm 
%    for large-scale quadratic programming", IMA J. Num. Anal (1991),
%    11, 299-324, problem class 1.
% 
%    SIF input: B Baudson, Jan 1990 modified by Nick Gould, Jan, 2011
% 
%    classification = 'C-CQLR2-AN-32-17'
% 
%    Other useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GOULDQP1';

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
        v_('0') = 0;
        v_('1') = 1;
        v_('4') = 4;
        v_('5') = 5;
        v_('8') = 8;
        v_('12') = 12;
        v_('15') = 15;
        v_('17') = 17;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('15')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for K=v_('1'):v_('4')
            [iv,ix_] = s2mpjlib('ii',['AS',int2str(K)],ix_);
            pb.xnames{iv} = ['AS',int2str(K)];
        end
        for K=v_('1'):v_('4')
            [iv,ix_] = s2mpjlib('ii',['CS',int2str(K)],ix_);
            pb.xnames{iv} = ['CS',int2str(K)];
        end
        for K=v_('1'):v_('4')
            [iv,ix_] = s2mpjlib('ii',['BS',int2str(K)],ix_);
            pb.xnames{iv} = ['BS',int2str(K)];
        end
        for K=v_('1'):v_('5')
            [iv,ix_] = s2mpjlib('ii',['DS',int2str(K)],ix_);
            pb.xnames{iv} = ['DS',int2str(K)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for K=v_('0'):v_('4')
            v_('3K') = 3*K;
            v_('3K+1') = 1+v_('3K');
            v_('3K+2') = 2+v_('3K');
            v_('3K+3') = 3+v_('3K');
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+1')))]);
            valA(end+1) = 2.3;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+2')))]);
            valA(end+1) = 1.7;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+3')))]);
            valA(end+1) = 2.2;
        end
        for K=v_('1'):v_('4')
            v_('3K') = 3*K;
            v_('3K+1') = 1+v_('3K');
            v_('3K+2') = 2+v_('3K');
            v_('3K+3') = 3+v_('3K');
            v_('3K-2') = -2+v_('3K');
            v_('3K-1') = -1+v_('3K');
            [ig,ig_] = s2mpjlib('ii',['A',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(K)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K-2')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['AS',int2str(K)]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['B',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['B',int2str(K)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+3')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['BS',int2str(K)]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(K)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K+2')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('3K-1')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['CS',int2str(K)]);
            valA(end+1) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','D1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DS1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','D2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DS2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','D3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DS3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','D4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DS4');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','D5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DS5');
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
        for K=v_('1'):v_('4')
            pbm.gconst(ig_(['A',int2str(K)])) = -7.0;
            pbm.gconst(ig_(['B',int2str(K)])) = -7.0;
            pbm.gconst(ig_(['C',int2str(K)])) = -7.0;
        end
        pbm.gconst(ig_('D1')) = 60.0;
        pbm.gconst(ig_('D2')) = 50.0;
        pbm.gconst(ig_('D3')) = 70.0;
        pbm.gconst(ig_('D4')) = 85.0;
        pbm.gconst(ig_('D5')) = 100.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 8.0;
        pb.xupper(ix_('X1')) = 21.0;
        pb.xlower(ix_('X2'),1) = 43.0;
        pb.xupper(ix_('X2')) = 57.0;
        pb.xlower(ix_('X3'),1) = 3.0;
        pb.xupper(ix_('X3')) = 16.0;
        for K=v_('1'):v_('4')
            v_('3K') = 3*K;
            v_('3K+1') = 1+v_('3K');
            v_('3K+2') = 2+v_('3K');
            v_('3K+3') = 3+v_('3K');
            pb.xupper(ix_(['X',int2str(round(v_('3K+1')))])) = 90.0;
            pb.xupper(ix_(['X',int2str(round(v_('3K+2')))])) = 120.0;
            pb.xupper(ix_(['X',int2str(round(v_('3K+3')))])) = 60.0;
        end
        for K=v_('1'):v_('4')
            pb.xlower(ix_(['AS',int2str(K)]),1) = 0.0;
            pb.xupper(ix_(['AS',int2str(K)])) = 13.0;
            pb.xlower(ix_(['BS',int2str(K)]),1) = 0.0;
            pb.xupper(ix_(['BS',int2str(K)])) = 13.0;
            pb.xlower(ix_(['CS',int2str(K)]),1) = 0.0;
            pb.xupper(ix_(['CS',int2str(K)])) = 14.0;
        end
        pb.xlower(ix_('DS1'),1) = 0.0;
        pb.xlower(ix_('DS2'),1) = 0.0;
        pb.xlower(ix_('DS3'),1) = 0.0;
        pb.xlower(ix_('DS4'),1) = 0.0;
        pb.xlower(ix_('DS5'),1) = 0.0;
        pb.xupper(ix_('DS1')) = 60.0;
        pb.xupper(ix_('DS2')) = 50.0;
        pb.xupper(ix_('DS3')) = 70.0;
        pb.xupper(ix_('DS4')) = 85.0;
        pb.xupper(ix_('DS5')) = 100.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 20.0*ones(pb.n,1);
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 55.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 55.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 15.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 15.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 60.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 60.0;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 60.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8')),1) = 60.0;
        end
        if(isKey(ix_,'X11'))
            pb.x0(ix_('X11'),1) = 60.0;
        else
            pb.y0(find(pbm.congrps==ig_('X11')),1) = 60.0;
        end
        if(isKey(ix_,'X14'))
            pb.x0(ix_('X14'),1) = 60.0;
        else
            pb.y0(find(pbm.congrps==ig_('X14')),1) = 60.0;
        end
        if(isKey(ix_,'AS1'))
            pb.x0(ix_('AS1'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('AS1')),1) = 7.0;
        end
        if(isKey(ix_,'BS1'))
            pb.x0(ix_('BS1'),1) = 12.0;
        else
            pb.y0(find(pbm.congrps==ig_('BS1')),1) = 12.0;
        end
        if(isKey(ix_,'CS1'))
            pb.x0(ix_('CS1'),1) = 12.0;
        else
            pb.y0(find(pbm.congrps==ig_('CS1')),1) = 12.0;
        end
        if(isKey(ix_,'AS2'))
            pb.x0(ix_('AS2'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('AS2')),1) = 7.0;
        end
        if(isKey(ix_,'BS2'))
            pb.x0(ix_('BS2'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('BS2')),1) = 7.0;
        end
        if(isKey(ix_,'CS2'))
            pb.x0(ix_('CS2'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('CS2')),1) = 7.0;
        end
        if(isKey(ix_,'AS3'))
            pb.x0(ix_('AS3'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('AS3')),1) = 7.0;
        end
        if(isKey(ix_,'BS3'))
            pb.x0(ix_('BS3'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('BS3')),1) = 7.0;
        end
        if(isKey(ix_,'CS3'))
            pb.x0(ix_('CS3'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('CS3')),1) = 7.0;
        end
        if(isKey(ix_,'AS4'))
            pb.x0(ix_('AS4'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('AS4')),1) = 7.0;
        end
        if(isKey(ix_,'BS4'))
            pb.x0(ix_('BS4'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('BS4')),1) = 7.0;
        end
        if(isKey(ix_,'CS4'))
            pb.x0(ix_('CS4'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('CS4')),1) = 7.0;
        end
        if(isKey(ix_,'DS1'))
            pb.x0(ix_('DS1'),1) = 30.0;
        else
            pb.y0(find(pbm.congrps==ig_('DS1')),1) = 30.0;
        end
        if(isKey(ix_,'DS2'))
            pb.x0(ix_('DS2'),1) = 50.0;
        else
            pb.y0(find(pbm.congrps==ig_('DS2')),1) = 50.0;
        end
        if(isKey(ix_,'DS3'))
            pb.x0(ix_('DS3'),1) = 30.0;
        else
            pb.y0(find(pbm.congrps==ig_('DS3')),1) = 30.0;
        end
        if(isKey(ix_,'DS4'))
            pb.x0(ix_('DS4'),1) = 15.0;
        else
            pb.y0(find(pbm.congrps==ig_('DS4')),1) = 15.0;
        end
        if(isKey(ix_,'DS5'))
            pb.x0(ix_('DS5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('DS5')),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('15')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],20.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.00015;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 10.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 25.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E12');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.00015;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.0001;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.00015;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -3.485333E+3
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-AN-32-17';
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

