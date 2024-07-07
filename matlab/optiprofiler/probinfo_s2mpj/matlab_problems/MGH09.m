function varargout = MGH09(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MGH09
%    *********
% 
%    NIST Data fitting problem MGH09 given as an inconsistent set of
%    nonlinear equations.
% 
%    Fit: y = b1*(x**2+x*b2) / (x**2+x*b3+b4) + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Kowalik, J.S., and M. R. Osborne, (1978).  
%      Methods for Unconstrained Optimization Problems.  
%      New York, NY:  Elsevier North-Holland.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'NOR2-MN-4-11'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MGH09';

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
        v_('M') = 11;
        v_('N') = 4;
        v_('1') = 1;
        v_('X1') = 4.0E+0;
        v_('X2') = 2.0E+0;
        v_('X3') = 1.0E+0;
        v_('X4') = 5.00E-1;
        v_('X5') = 2.50E-1;
        v_('X6') = 1.67E-1;
        v_('X7') = 1.25E-1;
        v_('X8') = 1.00E-1;
        v_('X9') = 8.33E-2;
        v_('X10') = 7.14E-2;
        v_('X11') = 6.25E-2;
        v_('Y1') = 1.957E-1;
        v_('Y2') = 1.947E-1;
        v_('Y3') = 1.735E-1;
        v_('Y4') = 1.60E-1;
        v_('Y5') = 8.44E-2;
        v_('Y6') = 6.27E-2;
        v_('Y7') = 4.56E-2;
        v_('Y8') = 3.42E-2;
        v_('Y9') = 3.23E-2;
        v_('Y10') = 2.35E-2;
        v_('Y11') = 2.46E-2;
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
            pb.x0(ix_('B1'),1) = 25.0;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 25.0;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = 39.0;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = 39.0;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = 41.5;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = 41.5;
        end
        if(isKey(ix_,'B4'))
            pb.x0(ix_('B4'),1) = 39.0;
        else
            pb.y0(find(pbm.congrps==ig_('B4')),1) = 39.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE10',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
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
            pbm.elftype{ie} = 'eE10';
            ielftype(ie) = iet_('eE10');
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
        pb.pbclass = 'NOR2-MN-4-11';
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

    case 'eE10'

        EV_  = varargin{1};
        iel_ = varargin{2};
        X2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        T = EV_(2)*pbm.elpar{iel_}(1)+X2;
        B = EV_(4)+EV_(3)*pbm.elpar{iel_}(1)+X2;
        B2 = B*B;
        B3 = B*B2;
        V1X = EV_(1)*pbm.elpar{iel_}(1);
        V1X2 = EV_(1)*X2;
        V1T = EV_(1)*T;
        V1XT = V1X*T;
        varargout{1} = V1T/B;
        if(nargout>1)
            g_(1,1) = T/B;
            g_(2,1) = V1X/B;
            g_(3,1) = -V1XT/B2;
            g_(4,1) = -V1T/B2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = pbm.elpar{iel_}(1)/B;
                H_(2,1) = H_(1,2);
                H_(1,3) = -pbm.elpar{iel_}(1)*T/B2;
                H_(3,1) = H_(1,3);
                H_(1,4) = -T/B2;
                H_(4,1) = H_(1,4);
                H_(2,3) = -V1X2/B2;
                H_(3,2) = H_(2,3);
                H_(2,4) = -V1X/B2;
                H_(4,2) = H_(2,4);
                H_(3,3) = 2.0*V1X2*T/B3;
                H_(3,4) = 2.0*V1XT/B3;
                H_(4,3) = H_(3,4);
                H_(4,4) = 2.0*V1T/B3;
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

