function varargout = HS57(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS57
%    *********
% 
%    Source: problem 57 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn, April 1990
% 
%    classification = 'C-CSQR2-AN-2-1'
% 
%    Problem parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS57';

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
        v_('A1') = 8.0;
        v_('A2') = 8.0;
        v_('A3') = 10.0;
        v_('A4') = 10.0;
        v_('A5') = 10.0;
        v_('A6') = 10.0;
        v_('A7') = 12.0;
        v_('A8') = 12.0;
        v_('A9') = 12.0;
        v_('A10') = 12.0;
        v_('A11') = 14.0;
        v_('A12') = 14.0;
        v_('A13') = 14.0;
        v_('A14') = 16.0;
        v_('A15') = 16.0;
        v_('A16') = 16.0;
        v_('A17') = 18.0;
        v_('A18') = 18.0;
        v_('A19') = 20.0;
        v_('A20') = 20.0;
        v_('A21') = 20.0;
        v_('A22') = 22.0;
        v_('A23') = 22.0;
        v_('A24') = 22.0;
        v_('A25') = 24.0;
        v_('A26') = 24.0;
        v_('A27') = 24.0;
        v_('A28') = 26.0;
        v_('A29') = 26.0;
        v_('A30') = 26.0;
        v_('A31') = 28.0;
        v_('A32') = 28.0;
        v_('A33') = 30.0;
        v_('A34') = 30.0;
        v_('A35') = 30.0;
        v_('A36') = 32.0;
        v_('A37') = 32.0;
        v_('A38') = 34.0;
        v_('A39') = 36.0;
        v_('A40') = 36.0;
        v_('A41') = 38.0;
        v_('A42') = 38.0;
        v_('A43') = 40.0;
        v_('A44') = 42.0;
        v_('B1') = 0.49;
        v_('B2') = 0.49;
        v_('B3') = 0.48;
        v_('B4') = 0.47;
        v_('B5') = 0.48;
        v_('B6') = 0.47;
        v_('B7') = 0.46;
        v_('B8') = 0.46;
        v_('B9') = 0.45;
        v_('B10') = 0.43;
        v_('B11') = 0.45;
        v_('B12') = 0.43;
        v_('B13') = 0.43;
        v_('B14') = 0.44;
        v_('B15') = 0.43;
        v_('B16') = 0.43;
        v_('B17') = 0.46;
        v_('B18') = 0.45;
        v_('B19') = 0.42;
        v_('B20') = 0.42;
        v_('B21') = 0.43;
        v_('B22') = 0.41;
        v_('B23') = 0.41;
        v_('B24') = 0.40;
        v_('B25') = 0.42;
        v_('B26') = 0.40;
        v_('B27') = 0.40;
        v_('B28') = 0.41;
        v_('B29') = 0.40;
        v_('B30') = 0.41;
        v_('B31') = 0.41;
        v_('B32') = 0.40;
        v_('B33') = 0.40;
        v_('B34') = 0.40;
        v_('B35') = 0.38;
        v_('B36') = 0.41;
        v_('B37') = 0.40;
        v_('B38') = 0.40;
        v_('B39') = 0.41;
        v_('B40') = 0.38;
        v_('B41') = 0.40;
        v_('B42') = 0.40;
        v_('B43') = 0.39;
        v_('B44') = 0.39;
        v_('1') = 1;
        v_('44') = 44;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON1';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.49+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.49;
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
        pbm.gconst(ig_('CON1')) = 0.09;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.4;
        pb.xlower(ix_('X2'),1) = -4.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 0.42;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 0.42;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 5.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 5.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOBSQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'AA';
        elftp{it}{2} = 'BB';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('44')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eOBSQ';
            ielftype(ie) = iet_('eOBSQ');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('AA',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(I)]);
            [~,posep] = ismember('BB',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['B',int2str(I)]);
        end
        ename = 'PR';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('44')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_('CON1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PR');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                 0.02845966
% LO SOLTN                 0.03063791
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CSQR2-AN-2-1';
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

    case 'eOBSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        AM8 = pbm.elpar{iel_}(1)-8.0;
        CMV1 = 0.49-EV_(1);
        E = exp(-EV_(2)*AM8);
        DED2 = -AM8*E;
        R = pbm.elpar{iel_}(2)-EV_(1)-CMV1*E;
        DRD1 = E-1.0;
        DRD2 = -CMV1*DED2;
        D2RD22 = -CMV1*AM8*AM8*E;
        varargout{1} = R*R;
        if(nargout>1)
            g_(1,1) = 2.0*R*DRD1;
            g_(2,1) = 2.0*R*DRD2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*DRD1*DRD1;
                H_(1,2) = 2.0*(DRD2*DRD1+R*DED2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*(DRD2*DRD2+R*D2RD22);
                varargout{3} = H_;
            end
        end

    case 'en2PR'

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

