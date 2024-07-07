function varargout = MGH10SLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MGH10SLS
%    *********
% 
%    NIST Data fitting problem MGH10 given as an inconsistent set of
%    nonlinear equations (scaled version)
% 
%    Fit: y = 0.01 * b1 * exp[ 1000 * b2 / (x + 100 * b3 ) ] + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Meyer, R. R. (1970).
%      Theoretical and computational aspects of nonlinear
%      regression.  In Nonlinear Programming, Rosen,
%      Mangasarian and Ritter (Eds).
%      New York, NY: Academic Press, pp. 465-486.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'SUR2-MN-3-0'
%    Least-squares version of MGH10S.SIF, Nick Gould, Jan 2020
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MGH10SLS';

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
        v_('M') = 16;
        v_('N') = 3;
        v_('1') = 1;
        v_('X1') = 5.00E+01;
        v_('X2') = 5.50E+01;
        v_('X3') = 6.00E+01;
        v_('X4') = 6.50E+01;
        v_('X5') = 7.00E+01;
        v_('X6') = 7.50E+01;
        v_('X7') = 8.00E+01;
        v_('X8') = 8.50E+01;
        v_('X9') = 9.00E+01;
        v_('X10') = 9.50E+01;
        v_('X11') = 1.00E+02;
        v_('X12') = 1.05E+02;
        v_('X13') = 1.10E+02;
        v_('X14') = 1.15E+02;
        v_('X15') = 1.20E+02;
        v_('X16') = 1.25E+02;
        v_('Y1') = 3.478e+04;
        v_('Y2') = 2.861e+04;
        v_('Y3') = 2.365e+04;
        v_('Y4') = 1.963e+04;
        v_('Y5') = 1.637e+04;
        v_('Y6') = 1.372e+04;
        v_('Y7') = 1.154e+04;
        v_('Y8') = 9.744e+03;
        v_('Y9') = 8.261e+03;
        v_('Y10') = 7.030e+03;
        v_('Y11') = 6.005e+03;
        v_('Y12') = 5.147e+03;
        v_('Y13') = 4.427e+03;
        v_('Y14') = 3.820e+03;
        v_('Y15') = 3.307e+03;
        v_('Y16') = 2.872e+03;
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
        pb.x0(ix_('B1'),1) = 200.0;
        pb.x0(ix_('B2'),1) = 400.0;
        pb.x0(ix_('B3'),1) = 250.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE12',iet_);
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
            pbm.elftype{ie} = 'eE12';
            ielftype(ie) = iet_('eE12');
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

    case 'eE12'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S1 = 0.01;
        S2 = 1000.0;
        S3 = 100.0;
        XPV3 = pbm.elpar{iel_}(1)+S3*EV_(3);
        XPV32 = XPV3*XPV3;
        XPV33 = XPV3*XPV32;
        E = S1*exp(S2*EV_(2)/XPV3);
        F = EV_(1)*E;
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = E;
            g_(2,1) = S2*F/XPV3;
            g_(3,1) = -S2*S3*F*EV_(2)/XPV32;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = S2*E/XPV3;
                H_(2,1) = H_(1,2);
                H_(1,3) = -S2*S3*EV_(2)*E/XPV32;
                H_(3,1) = H_(1,3);
                H_(2,2) = S2*S2*F/XPV32;
                H_(2,3) = -S2*S3*F*(1.0/XPV32+S2*EV_(2)/XPV33);
                H_(3,2) = H_(2,3);
                H_(3,3) = S2*S3*S3*F*EV_(2)*(S2*EV_(2)/XPV3^4+2.0/XPV33);
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

