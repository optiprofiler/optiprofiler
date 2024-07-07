function varargout = CHWIRUT2LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHWIRUT2LS
%    *********
% 
%    NIST Data fitting problem CHWIRUT2.
% 
%    Fit: y = exp[-b1*x]/(b2+b3*x) + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Chwirut, D., NIST (197?).  
%      Ultrasonic Reference Block Study. 
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'SUR2-MN-3-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHWIRUT2LS';

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
        v_('M') = 54;
        v_('N') = 3;
        v_('1') = 1;
        v_('X1') = 0.5;
        v_('X2') = 1.0;
        v_('X3') = 1.75;
        v_('X4') = 3.75;
        v_('X5') = 5.75;
        v_('X6') = 0.875;
        v_('X7') = 2.25;
        v_('X8') = 3.25;
        v_('X9') = 5.25;
        v_('X10') = 0.75;
        v_('X11') = 1.75;
        v_('X12') = 2.75;
        v_('X13') = 4.75;
        v_('X14') = 0.625;
        v_('X15') = 1.25;
        v_('X16') = 2.25;
        v_('X17') = 4.25;
        v_('X18') = 0.5;
        v_('X19') = 3.0;
        v_('X20') = 0.75;
        v_('X21') = 3.0;
        v_('X22') = 1.5;
        v_('X23') = 6.0;
        v_('X24') = 3.0;
        v_('X25') = 6.0;
        v_('X26') = 1.5;
        v_('X27') = 3.0;
        v_('X28') = 0.5;
        v_('X29') = 2.0;
        v_('X30') = 4.0;
        v_('X31') = 0.75;
        v_('X32') = 2.0;
        v_('X33') = 5.0;
        v_('X34') = 0.75;
        v_('X35') = 2.25;
        v_('X36') = 3.75;
        v_('X37') = 5.75;
        v_('X38') = 3.0;
        v_('X39') = 0.75;
        v_('X40') = 2.5;
        v_('X41') = 4.0;
        v_('X42') = 0.75;
        v_('X43') = 2.5;
        v_('X44') = 4.0;
        v_('X45') = 0.75;
        v_('X46') = 2.5;
        v_('X47') = 4.0;
        v_('X48') = 0.5;
        v_('X49') = 6.0;
        v_('X50') = 3.0;
        v_('X51') = 0.5;
        v_('X52') = 2.75;
        v_('X53') = 0.5;
        v_('X54') = 1.75;
        v_('Y1') = 92.9;
        v_('Y2') = 57.1;
        v_('Y3') = 31.05;
        v_('Y4') = 11.5875;
        v_('Y5') = 8.025;
        v_('Y6') = 63.6;
        v_('Y7') = 21.4;
        v_('Y8') = 14.25;
        v_('Y9') = 8.475;
        v_('Y10') = 63.8;
        v_('Y11') = 26.8;
        v_('Y12') = 16.4625;
        v_('Y13') = 7.125;
        v_('Y14') = 67.3;
        v_('Y15') = 41.0;
        v_('Y16') = 21.15;
        v_('Y17') = 8.175;
        v_('Y18') = 81.50;
        v_('Y19') = 13.12;
        v_('Y20') = 59.9;
        v_('Y21') = 14.62;
        v_('Y22') = 32.9;
        v_('Y23') = 5.44;
        v_('Y24') = 12.56;
        v_('Y25') = 5.44;
        v_('Y26') = 32.0;
        v_('Y27') = 13.95;
        v_('Y28') = 75.8;
        v_('Y29') = 20.0;
        v_('Y30') = 10.42;
        v_('Y31') = 59.5;
        v_('Y32') = 21.67;
        v_('Y33') = 8.55;
        v_('Y34') = 62.0;
        v_('Y35') = 20.2;
        v_('Y36') = 7.76;
        v_('Y37') = 3.75;
        v_('Y38') = 11.81;
        v_('Y39') = 54.7;
        v_('Y40') = 23.7;
        v_('Y41') = 11.55;
        v_('Y42') = 61.3;
        v_('Y43') = 17.7;
        v_('Y44') = 8.74;
        v_('Y45') = 59.2;
        v_('Y46') = 16.3;
        v_('Y47') = 8.62;
        v_('Y48') = 81.0;
        v_('Y49') = 4.87;
        v_('Y50') = 14.62;
        v_('Y51') = 81.7;
        v_('Y52') = 17.17;
        v_('Y53') = 81.3;
        v_('Y54') = 28.9;
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
        pb.x0(ix_('B1'),1) = 0.1;
        pb.x0(ix_('B2'),1) = 0.01;
        pb.x0(ix_('B3'),1) = 0.02;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE16',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
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
            pbm.elftype{ie} = 'eE16';
            ielftype(ie) = iet_('eE16');
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
        pb.pbclass = 'SUR2-MN-3-0';
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

    case 'eE16'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E = exp(-EV_(1)*pbm.elpar{iel_}(1));
        EX = E*pbm.elpar{iel_}(1);
        EX2 = EX*pbm.elpar{iel_}(1);
        V2PV3X = EV_(2)+EV_(3)*pbm.elpar{iel_}(1);
        V2PV32 = V2PV3X*V2PV3X;
        V2PV33 = V2PV3X*V2PV32;
        varargout{1} = E/V2PV3X;
        if(nargout>1)
            g_(1,1) = -EX/V2PV3X;
            g_(2,1) = -E/V2PV32;
            g_(3,1) = -EX/V2PV32;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = EX2/V2PV3X;
                H_(1,2) = EX/V2PV32;
                H_(2,1) = H_(1,2);
                H_(1,3) = EX2/V2PV32;
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*E/V2PV33;
                H_(2,3) = 2.0*EX/V2PV33;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*EX2/V2PV33;
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

