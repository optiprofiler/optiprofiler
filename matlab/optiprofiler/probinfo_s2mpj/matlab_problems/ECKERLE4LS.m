function varargout = ECKERLE4LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ECKERLE4LS
%    *********
% 
%    NIST Data fitting problem ECKERLE4.
% 
%    Fit: y = (b1/b2) * exp[-0.5*((x-b3)/b2)**2] + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Eckerle, K., NIST (197?).  
%      Circular Interference Transmittance Study.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'C-CSUR2-MN-3-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ECKERLE4LS';

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
        v_('M') = 35;
        v_('N') = 3;
        v_('1') = 1;
        v_('X1') = 400.000000;
        v_('X2') = 405.000000;
        v_('X3') = 410.000000;
        v_('X4') = 415.000000;
        v_('X5') = 420.000000;
        v_('X6') = 425.000000;
        v_('X7') = 430.000000;
        v_('X8') = 435.000000;
        v_('X9') = 436.500000;
        v_('X10') = 438.000000;
        v_('X11') = 439.500000;
        v_('X12') = 441.000000;
        v_('X13') = 442.500000;
        v_('X14') = 444.000000;
        v_('X15') = 445.500000;
        v_('X16') = 447.000000;
        v_('X17') = 448.500000;
        v_('X18') = 450.000000;
        v_('X19') = 451.500000;
        v_('X20') = 453.000000;
        v_('X21') = 454.500000;
        v_('X22') = 456.000000;
        v_('X23') = 457.500000;
        v_('X24') = 459.000000;
        v_('X25') = 460.500000;
        v_('X26') = 462.000000;
        v_('X27') = 463.500000;
        v_('X28') = 465.000000;
        v_('X29') = 470.000000;
        v_('X30') = 475.000000;
        v_('X31') = 480.000000;
        v_('X32') = 485.000000;
        v_('X33') = 490.000000;
        v_('X34') = 495.000000;
        v_('X35') = 500.000000;
        v_('Y1') = 0.0001575;
        v_('Y2') = 0.0001699;
        v_('Y3') = 0.0002350;
        v_('Y4') = 0.0003102;
        v_('Y5') = 0.0004917;
        v_('Y6') = 0.0008710;
        v_('Y7') = 0.0017418;
        v_('Y8') = 0.0046400;
        v_('Y9') = 0.0065895;
        v_('Y10') = 0.0097302;
        v_('Y11') = 0.0149002;
        v_('Y12') = 0.0237310;
        v_('Y13') = 0.0401683;
        v_('Y14') = 0.0712559;
        v_('Y15') = 0.1264458;
        v_('Y16') = 0.2073413;
        v_('Y17') = 0.2902366;
        v_('Y18') = 0.3445623;
        v_('Y19') = 0.3698049;
        v_('Y20') = 0.3668534;
        v_('Y21') = 0.3106727;
        v_('Y22') = 0.2078154;
        v_('Y23') = 0.1164354;
        v_('Y24') = 0.0616764;
        v_('Y25') = 0.0337200;
        v_('Y26') = 0.0194023;
        v_('Y27') = 0.0117831;
        v_('Y28') = 0.0074357;
        v_('Y29') = 0.0022732;
        v_('Y30') = 0.0008800;
        v_('Y31') = 0.0004579;
        v_('Y32') = 0.0002345;
        v_('Y33') = 0.0001586;
        v_('Y34') = 0.0001143;
        v_('Y35') = 0.0000710;
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
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('B1'),1) = 1.0;
        pb.x0(ix_('B2'),1) = 10.0;
        pb.x0(ix_('B3'),1) = 500.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE';
            ielftype(ie) = iet_('eE');
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
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
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
        pb.pbclass = 'C-CSUR2-MN-3-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

    case 'eE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V3MX = EV_(3)-pbm.elpar{iel_}(1);
        TV3MX = 2.0*V3MX;
        V3MX2 = V3MX^2;
        V22 = EV_(2)^2;
        V23 = EV_(2)*V22;
        V24 = V22*V22;
        V25 = V22*V23;
        V26 = V23*V23;
        V27 = V23*V24;
        E = exp(-0.5*V3MX2/V22);
        V1E = EV_(1)*E;
        DIFF = V3MX2/V24-1.0/V22;
        varargout{1} = EV_(1)*E/EV_(2);
        if(nargout>1)
            g_(1,1) = E/EV_(2);
            g_(2,1) = V1E*DIFF;
            g_(3,1) = -V1E*V3MX/V23;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = E*DIFF;
                H_(2,1) = H_(1,2);
                H_(1,3) = -0.5*E*TV3MX/V23;
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*V1E/V23-5.0*V1E*V3MX2/V25+V1E*V3MX^4/V27;
                H_(2,3) = 1.5*V1E*TV3MX/V24-0.5*V1E*TV3MX*V3MX2/V26;
                H_(3,2) = H_(2,3);
                H_(3,3) = 0.5*V1E*V3MX2/V25-V1E/V23;
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

