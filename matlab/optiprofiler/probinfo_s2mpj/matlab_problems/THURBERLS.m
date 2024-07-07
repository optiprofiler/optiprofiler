function varargout = THURBERLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : THURBERLS
%    *********
% 
%    NIST Data fitting problem THURBERLS.
% 
%    Fit: y = (b1 + b2*x + b3*x**2 + b4*x**3) / 
%             (1 + b5*x + b6*x**2 + b7*x**3) + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Thurber, R., NIST (197?).  
%      Semiconductor electron mobility modeling.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'SUR2-MN-7-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'THURBERLS';

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
        v_('M') = 37;
        v_('N') = 7;
        v_('1') = 1;
        v_('X1') = -3.067;
        v_('X2') = -2.981;
        v_('X3') = -2.921;
        v_('X4') = -2.912;
        v_('X5') = -2.840;
        v_('X6') = -2.797;
        v_('X7') = -2.702;
        v_('X8') = -2.699;
        v_('X9') = -2.633;
        v_('X10') = -2.481;
        v_('X11') = -2.363;
        v_('X12') = -2.322;
        v_('X13') = -1.501;
        v_('X14') = -1.460;
        v_('X15') = -1.274;
        v_('X16') = -1.212;
        v_('X17') = -1.100;
        v_('X18') = -1.046;
        v_('X19') = -0.915;
        v_('X20') = -0.714;
        v_('X21') = -0.566;
        v_('X22') = -0.545;
        v_('X23') = -0.400;
        v_('X24') = -0.309;
        v_('X25') = -0.109;
        v_('X26') = -0.103;
        v_('X27') = 0.010;
        v_('X28') = 0.119;
        v_('X29') = 0.377;
        v_('X30') = 0.790;
        v_('X31') = 0.963;
        v_('X32') = 1.006;
        v_('X33') = 1.115;
        v_('X34') = 1.572;
        v_('X35') = 1.841;
        v_('X36') = 2.047;
        v_('X37') = 2.200;
        v_('Y1') = 80.574;
        v_('Y2') = 84.248;
        v_('Y3') = 87.264;
        v_('Y4') = 87.195;
        v_('Y5') = 89.076;
        v_('Y6') = 89.608;
        v_('Y7') = 89.868;
        v_('Y8') = 90.101;
        v_('Y9') = 92.405;
        v_('Y10') = 95.854;
        v_('Y11') = 100.696;
        v_('Y12') = 101.060;
        v_('Y13') = 401.672;
        v_('Y14') = 390.724;
        v_('Y15') = 567.534;
        v_('Y16') = 635.316;
        v_('Y17') = 733.054;
        v_('Y18') = 759.087;
        v_('Y19') = 894.206;
        v_('Y20') = 990.785;
        v_('Y21') = 1090.109;
        v_('Y22') = 1080.914;
        v_('Y23') = 1122.643;
        v_('Y24') = 1178.351;
        v_('Y25') = 1260.531;
        v_('Y26') = 1273.514;
        v_('Y27') = 1288.339;
        v_('Y28') = 1327.543;
        v_('Y29') = 1353.863;
        v_('Y30') = 1414.509;
        v_('Y31') = 1425.208;
        v_('Y32') = 1421.384;
        v_('Y33') = 1442.962;
        v_('Y34') = 1464.350;
        v_('Y35') = 1468.705;
        v_('Y36') = 1447.894;
        v_('Y37') = 1457.628;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('B1'),1) = 1000.0;
        pb.x0(ix_('B2'),1) = 1000.0;
        pb.x0(ix_('B3'),1) = 400.0;
        pb.x0(ix_('B4'),1) = 40.0;
        pb.x0(ix_('B5'),1) = 0.7;
        pb.x0(ix_('B6'),1) = 0.3;
        pb.x0(ix_('B7'),1) = 0.03;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE19',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE19';
            ielftype(ie) = iet_('eE19');
            vname = 'B1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-7-0';
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

    case 'eE19'

        EV_  = varargin{1};
        iel_ = varargin{2};
        X2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        X3 = X2*pbm.elpar{iel_}(1);
        X4 = X3*pbm.elpar{iel_}(1);
        X5 = X4*pbm.elpar{iel_}(1);
        X6 = X5*pbm.elpar{iel_}(1);
        T = EV_(1)+EV_(2)*pbm.elpar{iel_}(1)+EV_(3)*X2+EV_(4)*X3;
        D = 1.0e0+EV_(5)*pbm.elpar{iel_}(1)+EV_(6)*X2+EV_(7)*X3;
        D2 = D*D;
        TD3 = 0.5e0*D2*D;
        varargout{1} = T/D;
        if(nargout>1)
            g_(1,1) = 1.0e0/D;
            g_(2,1) = pbm.elpar{iel_}(1)/D;
            g_(3,1) = X2/D;
            g_(4,1) = X3/D;
            g_(5,1) = -pbm.elpar{iel_}(1)*T/D2;
            g_(6,1) = -X2*T/D2;
            g_(7,1) = -X3*T/D2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(7,7);
                H_(1,5) = -pbm.elpar{iel_}(1)/D2;
                H_(5,1) = H_(1,5);
                H_(1,6) = -X2/D2;
                H_(6,1) = H_(1,6);
                H_(1,7) = -X3/D2;
                H_(7,1) = H_(1,7);
                H_(2,5) = -X2/D2;
                H_(5,2) = H_(2,5);
                H_(2,6) = -X3/D2;
                H_(6,2) = H_(2,6);
                H_(2,7) = -X4/D2;
                H_(7,2) = H_(2,7);
                H_(3,5) = -X3/D2;
                H_(5,3) = H_(3,5);
                H_(3,6) = -X4/D2;
                H_(6,3) = H_(3,6);
                H_(3,7) = -X5/D2;
                H_(7,3) = H_(3,7);
                H_(4,5) = -X4/D2;
                H_(5,4) = H_(4,5);
                H_(4,6) = -X5/D2;
                H_(6,4) = H_(4,6);
                H_(4,7) = -X6/D2;
                H_(7,4) = H_(4,7);
                H_(5,5) = X2*T/TD3;
                H_(5,6) = X3*T/TD3;
                H_(6,5) = H_(5,6);
                H_(5,7) = X4*T/TD3;
                H_(7,5) = H_(5,7);
                H_(6,6) = X4*T/TD3;
                H_(6,7) = X5*T/TD3;
                H_(7,6) = H_(6,7);
                H_(7,7) = X6*T/TD3;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
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

