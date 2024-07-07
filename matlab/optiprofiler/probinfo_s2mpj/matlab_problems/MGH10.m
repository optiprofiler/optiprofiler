function varargout = MGH10(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MGH10
%    *********
% 
%    NIST Data fitting problem MGH10 given as an inconsistent set of
%    nonlinear equations.
% 
%    Fit: y = b1 * exp[b2/(x+b3)] + e
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
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'NOR2-MN-3-16'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MGH10';

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
        v_('Y1') = 3.478E+04;
        v_('Y2') = 2.861E+04;
        v_('Y3') = 2.365E+04;
        v_('Y4') = 1.963E+04;
        v_('Y5') = 1.637E+04;
        v_('Y6') = 1.372E+04;
        v_('Y7') = 1.154E+04;
        v_('Y8') = 9.744E+03;
        v_('Y9') = 8.261E+03;
        v_('Y10') = 7.030E+03;
        v_('Y11') = 6.005E+03;
        v_('Y12') = 5.147E+03;
        v_('Y13') = 4.427E+03;
        v_('Y14') = 3.820E+03;
        v_('Y15') = 3.307E+03;
        v_('Y16') = 2.872E+03;
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
            gtype{ig}  = '==';
            cnames{ig} = ['F',int2str(I)];
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'B1'))
            pb.x0(ix_('B1'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 2.0;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = 400000.0;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = 400000.0;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = 25000.0;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = 25000.0;
        end
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
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-3-16';
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
        V3PX = EV_(3)+pbm.elpar{iel_}(1);
        V3PX2 = V3PX*V3PX;
        V3PX3 = V3PX*V3PX2;
        E = exp(EV_(2)/V3PX);
        F = EV_(1)*E;
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = E;
            g_(2,1) = F/V3PX;
            g_(3,1) = -F*EV_(2)/V3PX2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = E/V3PX;
                H_(2,1) = H_(1,2);
                H_(1,3) = -EV_(2)*E/V3PX2;
                H_(3,1) = H_(1,3);
                H_(2,2) = F/V3PX2;
                H_(2,3) = -F*(1.0/V3PX2+EV_(2)/V3PX3);
                H_(3,2) = H_(2,3);
                H_(3,3) = F*EV_(2)*(EV_(2)/V3PX^4+2.0/V3PX3);
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

