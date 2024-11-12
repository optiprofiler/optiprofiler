function varargout = HS105(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS105
%    *********
% 
%    Source: problem 105 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
%    bug correction (line 351) Ph. Toint, May 2024
% 
%    classification = 'C-COLR2-AY-8-1'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS105';

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
        v_('235') = 235;
        v_('Y1') = 95.0;
        v_('Y2') = 105.0;
        v_('LOW') = 3;
        v_('UP') = 6;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 110.0;
        end
        v_('LOW') = 7;
        v_('UP') = 10;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 115.0;
        end
        v_('LOW') = 11;
        v_('UP') = 25;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 120.0;
        end
        v_('LOW') = 26;
        v_('UP') = 40;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 125.0;
        end
        v_('LOW') = 41;
        v_('UP') = 55;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 130.0;
        end
        v_('LOW') = 56;
        v_('UP') = 68;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 135.0;
        end
        v_('LOW') = 69;
        v_('UP') = 89;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 140.0;
        end
        v_('LOW') = 90;
        v_('UP') = 101;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 145.0;
        end
        v_('LOW') = 102;
        v_('UP') = 118;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 150.0;
        end
        v_('LOW') = 119;
        v_('UP') = 122;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 155.0;
        end
        v_('LOW') = 123;
        v_('UP') = 142;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 160.0;
        end
        v_('LOW') = 143;
        v_('UP') = 150;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 165.0;
        end
        v_('LOW') = 151;
        v_('UP') = 167;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 170.0;
        end
        v_('LOW') = 168;
        v_('UP') = 175;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 175.0;
        end
        v_('LOW') = 176;
        v_('UP') = 181;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 180.0;
        end
        v_('LOW') = 182;
        v_('UP') = 187;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 185.0;
        end
        v_('LOW') = 188;
        v_('UP') = 194;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 190.0;
        end
        v_('LOW') = 195;
        v_('UP') = 198;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 195.0;
        end
        v_('LOW') = 199;
        v_('UP') = 201;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 200.0;
        end
        v_('LOW') = 202;
        v_('UP') = 204;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 205.0;
        end
        v_('LOW') = 205;
        v_('UP') = 212;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 210.0;
        end
        v_('Y213') = 215.0;
        v_('LOW') = 214;
        v_('UP') = 219;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 220.0;
        end
        v_('LOW') = 220;
        v_('UP') = 224;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 230.0;
        end
        v_('Y225') = 235.0;
        v_('LOW') = 226;
        v_('UP') = 232;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 240.0;
        end
        v_('Y233') = 245.0;
        v_('LOW') = 234;
        v_('UP') = 235;
        for I=v_('LOW'):v_('UP')
            v_(['Y',int2str(I)]) = 250.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('235')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('C1')) = -1.0e+0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.001;
        pb.xupper(ix_('X1')) = 0.499;
        pb.xlower(ix_('X2'),1) = 0.001;
        pb.xupper(ix_('X2')) = 0.499;
        pb.xlower(ix_('X3'),1) = 100.0;
        pb.xupper(ix_('X3')) = 180.0;
        pb.xlower(ix_('X4'),1) = 130.0;
        pb.xupper(ix_('X4')) = 210.0;
        pb.xlower(ix_('X5'),1) = 170.0;
        pb.xupper(ix_('X5')) = 240.0;
        pb.xlower(ix_('X6'),1) = 5.0;
        pb.xupper(ix_('X6')) = 25.0;
        pb.xlower(ix_('X7'),1) = 5.0;
        pb.xupper(ix_('X7')) = 25.0;
        pb.xlower(ix_('X8'),1) = 5.0;
        pb.xupper(ix_('X8')) = 25.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 0.1;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 0.1;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 0.2;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 0.2;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 100.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 100.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 125.0;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 125.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 175.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 175.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 11.2;
        else
            pb.y0(find(pbm.congrps==ig('X6')),1) = 11.2;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 13.2;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 13.2;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 15.8;
        else
            pb.y0(find(pbm.congrps==ig('X8')),1) = 15.8;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eABI',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'YI';
        [it,iet_] = s2mpjlib( 'ii', 'eCI',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X5';
        elftv{it}{4} = 'X8';
        elftp{it}{1} = 'YI';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('235')
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eABI';
            ielftype(ie) = iet_('eABI');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('YI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eABI';
            ielftype(ie) = iet_('eABI');
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('YI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCI';
            ielftype(ie) = iet_('eCI');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X8';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('YI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('235')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gLOG';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               1138.416240
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AY-8-1';
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

    case 'eABI'

        EV_  = varargin{1};
        iel_ = varargin{2};
        R = EV_(1)/EV_(2);
        D = (pbm.elpar{iel_}(1)-EV_(3))/EV_(2);
        E = exp(-5.0e-1*D*D);
        DDV2 = -D/EV_(2);
        DEV2 = E*(-D)*DDV2;
        DDV3 = -1.0e+0/EV_(2);
        DEV3 = E*(-D)*DDV3;
        varargout{1} = R*E;
        if(nargout>1)
            g_(1,1) = E/EV_(2);
            g_(2,1) = (D*D-1.0e+0)*R*E/EV_(2);
            g_(3,1) = D*R*E/EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = (DEV2-E/EV_(2))/EV_(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = DEV3/EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/EV_(2)))*R/EV_(2);
                H_(2,3) = (DDV2*E+D*DEV2-2.0e+0*D*E/EV_(2))*R/EV_(2);
                H_(3,2) = H_(2,3);
                H_(3,3) = (DDV3*E+D*DEV3)*R/EV_(2);
                varargout{3} = H_;
            end
        end

    case 'eCI'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(1,1) = U_(1,1)-1;
        U_(1,2) = U_(1,2)-1;
        U_(2,4) = U_(2,4)+1;
        U_(3,3) = U_(3,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        R = (1.0e+0+IV_(1))/IV_(2);
        D = (pbm.elpar{iel_}(1)-IV_(3))/IV_(2);
        E = exp(-5.0e-1*D*D);
        DDV2 = -D/IV_(2);
        DEV2 = E*(-D)*DDV2;
        DDV3 = -1.0e+0/IV_(2);
        DEV3 = E*(-D)*DDV3;
        varargout{1} = R*E;
        if(nargout>1)
            g_(1,1) = E/IV_(2);
            g_(2,1) = (D*D-1.0e+0)*R*E/IV_(2);
            g_(3,1) = D*R*E/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = (DEV2-E/IV_(2))/IV_(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = DEV3/IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/IV_(2)))*R/IV_(2);
                H_(2,3) = (DDV2*E+D*DEV2-2.0e+0*D*E/IV_(2))*R/IV_(2);
                H_(3,2) = H_(2,3);
                H_(3,3) = (DDV3*E+D*DEV3)*R/IV_(2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'g_globs'

        pbm = varargin{1};
        pbm.gfpar(1) = 3.9894228040143270e-01;    % this is P
        varargout{1} = pbm;

    case 'gLOG'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = -log(pbm.gfpar(1)*GVAR_);
        if(nargout>1)
            g_ = -1/GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 1/GVAR_^2;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,1];
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

