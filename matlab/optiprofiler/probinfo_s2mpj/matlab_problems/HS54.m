function varargout = HS54(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS54
%    *********
% 
%    Source: problem 54, incorrectly stated in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    Betts problem 11.7, JOTA 21, 1977, pp.137-174.
%    SIF input: A.R. Conn, April 1990 and Nick Gould, October 1990
% 
%    classification = 'OLR2-AN-6-1'
% 
%    some useful parameters, including N, the number of variables.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS54';

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
        v_('N') = 6;
        v_('1') = 1;
        v_('6') = 6;
        v_('RHO') = 2.0e-1;
        v_('RHOSQR') = v_('RHO')*v_('RHO');
        v_('1-RHOSQR') = 1.0-v_('RHOSQR');
        v_('FACTOR') = 1.0/v_('1-RHOSQR');
        v_('MU1') = 1.0e+4;
        v_('MU2') = 1.0e+0;
        v_('MU3') = 2.0e+6;
        v_('MU4') = 1.0e+1;
        v_('MU5') = 1.0e-3;
        v_('MU6') = 1.0e+8;
        v_('SIGMA1') = 8.0e+3;
        v_('SIGMA2') = 1.0e+0;
        v_('SIGMA3') = 7.0e+6;
        v_('SIGMA4') = 5.0e+1;
        v_('SIGMA5') = 5.0e-2;
        v_('SIGMA6') = 5.0e+8;
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
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON1';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.0e+3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.0e+3;
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
        v_('0.2SI1') = 2.0e-1*v_('SIGMA1');
        v_('2000SI2') = 2.0e+3*v_('SIGMA2');
        v_('4000MU2') = 4.0e+3*v_('MU2');
        v_('RHS') = v_('MU1')+v_('4000MU2');
        v_('RHS') = v_('RHS')+v_('0.2SI1');
        v_('RHS') = v_('RHS')+v_('2000SI2');
        pbm.gconst(ig_('CON1')) = v_('RHS');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('X1')) = 2.0e+4;
        pb.xlower(ix_('X2'),1) = - 1.0e+1;
        pb.xupper(ix_('X2')) = 1.0e+1;
        pb.xupper(ix_('X3')) = 1.0e+7;
        pb.xupper(ix_('X4')) = 2.0e+1;
        pb.xlower(ix_('X5'),1) = - 1.0e+0;
        pb.xupper(ix_('X5')) = 1.0e+0;
        pb.xupper(ix_('X6')) = 2.0e+8;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 6.0e+3;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 6.0e+3;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 1.5e+0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 1.5e+0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 4.0e+6;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 4.0e+6;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 2.0e+0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 2.0e+0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 3.0e-3;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 3.0e-3;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 5.0e+7;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 5.0e+7;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V1';
        elftp{it}{1} = 'MU';
        elftp{it}{2} = 'SIGMA';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'MU1';
        elftp{it}{2} = 'MU2';
        elftp{it}{3} = 'SIGMA1';
        elftp{it}{4} = 'SIGMA2';
        elftp{it}{5} = 'RHO';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('6')
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eSQR';
                ielftype(ie) = iet_('eSQR');
            end
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('MU',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['MU',int2str(I)]);
            [~,posep] = ismember('SIGMA',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['SIGMA',int2str(I)]);
        end
        ename = 'F1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RHO',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('RHO');
        [~,posep] = ismember('MU1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('MU1');
        [~,posep] = ismember('MU2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('MU2');
        [~,posep] = ismember('SIGMA1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('SIGMA1');
        [~,posep] = ismember('SIGMA2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('SIGMA2');
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gNORMAL',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        pbm.grftype{ig} = 'gNORMAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('FACTOR');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('FACTOR');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.0e+0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        pbm.grelw{ig}(posel) = 1.0e+0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('FACTOR');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.90807482
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OLR2-AN-6-1';
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

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V1MP = (EV_(1)-pbm.elpar{iel_}(1))/pbm.elpar{iel_}(2);
        varargout{1} = V1MP^2;
        if(nargout>1)
            g_(1,1) = 2.0*V1MP/pbm.elpar{iel_}(2);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0/pbm.elpar{iel_}(2)^2;
                varargout{3} = H_;
            end
        end

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TERM1 = (EV_(1)-pbm.elpar{iel_}(1))/pbm.elpar{iel_}(3);
        TERM2 = (EV_(2)-pbm.elpar{iel_}(2))/pbm.elpar{iel_}(4);
        RHO2 = pbm.elpar{iel_}(5)+pbm.elpar{iel_}(5);
        varargout{1} = RHO2*TERM1*TERM2;
        if(nargout>1)
            g_(1,1) = RHO2*TERM2/pbm.elpar{iel_}(3);
            g_(2,1) = RHO2*TERM1/pbm.elpar{iel_}(4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = RHO2/(pbm.elpar{iel_}(3)*pbm.elpar{iel_}(4));
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gNORMAL'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        EXPHV = exp(-0.5*GVAR_);
        varargout{1} = -EXPHV;
        if(nargout>1)
            g_ = 5.0e-1*EXPHV;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -2.5e-1*EXPHV;
                varargout{3} = H_;
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

