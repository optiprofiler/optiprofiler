function varargout = HS84(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    *********
%    Problem : HS84
%    *********
% 
%    Source: problem 84 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn, March 1991.
% 
%    classification = 'C-CQQR2-AN-5-3'
% 
%    Set useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS84';

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
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_('15') = 15;
        v_('16') = 16;
        v_('17') = 17;
        v_('18') = 18;
        v_('19') = 19;
        v_('20') = 20;
        v_('21') = 21;
        v_(['A',int2str(round(v_('1')))]) = -24345.0;
        v_(['A',int2str(round(v_('2')))]) = -8720288.849;
        v_(['MA',int2str(round(v_('2')))]) = -1.0*v_(['A',int2str(round(v_('2')))]);
        v_(['A',int2str(round(v_('3')))]) = 150512.5253;
        v_(['MA',int2str(round(v_('3')))]) = -1.0*v_(['A',int2str(round(v_('3')))]);
        v_(['A',int2str(round(v_('4')))]) = -156.6950325;
        v_(['MA',int2str(round(v_('4')))]) = -1.0*v_(['A',int2str(round(v_('4')))]);
        v_(['A',int2str(round(v_('5')))]) = 476470.3222;
        v_(['MA',int2str(round(v_('5')))]) = -1.0*v_(['A',int2str(round(v_('5')))]);
        v_(['A',int2str(round(v_('6')))]) = 729482.8271;
        v_(['MA',int2str(round(v_('6')))]) = -1.0*v_(['A',int2str(round(v_('6')))]);
        v_(['A',int2str(round(v_('7')))]) = -145421.402;
        v_(['A',int2str(round(v_('8')))]) = 2931.1506;
        v_(['A',int2str(round(v_('9')))]) = -40.427932;
        v_(['A',int2str(round(v_('10')))]) = 5106.192;
        v_(['A',int2str(round(v_('11')))]) = 15711.36;
        v_(['A',int2str(round(v_('12')))]) = -155011.1084;
        v_(['A',int2str(round(v_('13')))]) = 4360.53352;
        v_(['A',int2str(round(v_('14')))]) = 12.9492344;
        v_(['A',int2str(round(v_('15')))]) = 10236.884;
        v_(['A',int2str(round(v_('16')))]) = 13176.786;
        v_(['A',int2str(round(v_('17')))]) = -326669.5104;
        v_(['A',int2str(round(v_('18')))]) = 7390.68412;
        v_(['A',int2str(round(v_('19')))]) = -27.8986976;
        v_(['A',int2str(round(v_('20')))]) = 16643.076;
        v_(['A',int2str(round(v_('21')))]) = 30988.146;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('5')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = v_(['MA',int2str(round(v_('2')))]);
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = v_(['A',int2str(round(v_('7')))]);
        [ig,ig_] = s2mpjlib('ii','CON2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = v_(['A',int2str(round(v_('12')))]);
        [ig,ig_] = s2mpjlib('ii','CON3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = v_(['A',int2str(round(v_('17')))]);
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
        pbm.gconst(ig_('OBJ')) = v_(['A',int2str(round(v_('1')))]);
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        grange(ig_('CON1')) = 294000.0;
        grange(ig_('CON2')) = 294000.0;
        grange(ig_('CON3')) = 277200.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('2')))]),1) = 1.2;
        pb.xlower(ix_(['X',int2str(round(v_('3')))]),1) = 20.0;
        pb.xlower(ix_(['X',int2str(round(v_('4')))]),1) = 9.0;
        pb.xlower(ix_(['X',int2str(round(v_('5')))]),1) = 6.5;
        pb.xupper(ix_(['X',int2str(round(v_('1')))])) = 1000.0;
        pb.xupper(ix_(['X',int2str(round(v_('2')))])) = 2.4;
        pb.xupper(ix_(['X',int2str(round(v_('3')))])) = 60.0;
        pb.xupper(ix_(['X',int2str(round(v_('4')))])) = 9.3;
        pb.xupper(ix_(['X',int2str(round(v_('5')))])) = 7.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('1')))]),1) = 2.52;
        pb.x0(ix_(['X',int2str(round(v_('2')))]),1) = 2.0;
        pb.x0(ix_(['X',int2str(round(v_('3')))]),1) = 37.5;
        pb.x0(ix_(['X',int2str(round(v_('4')))]),1) = 9.25;
        pb.x0(ix_(['X',int2str(round(v_('5')))]),1) = 6.8;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('4')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            v_('IP1') = 1+I;
            vname = ['X',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('IP1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['MA',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['MA',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['MA',int2str(round(v_('5')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['MA',int2str(round(v_('6')))]);
        ig = ig_('CON1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('8')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('9')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('10')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('11')))]);
        ig = ig_('CON2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('13')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('14')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('15')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('16')))]);
        ig = ig_('CON3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('18')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('19')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('20')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('21')))]);
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQQR2-AN-5-3';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

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

