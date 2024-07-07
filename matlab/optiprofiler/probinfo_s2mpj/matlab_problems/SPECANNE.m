function varargout = SPECANNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPECANNE
%    *********
% 
%    Source: a problem in spectral analysis suggested
%    by J. Eriksson and P. Lindstrom in "A Parallel Algorithm
%    for Bound Constrained Nonlinear Least Squares", UMEA TR S-901 87
% 
%    SIF input: Michael Ferris, July 1993
%    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
% 
%    classification = 'NOR2-AN-V-V'
% 
%    Number of Gaussians
% 
%       Alternative values for the SIF file parameters:
% IE K                   1              $-PARAMETER
% IE K                   2              $-PARAMETER
% IE K                   3              $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SPECANNE';

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
        if(nargs<1)
            v_('K') = 3;  %  SIF file default value
        else
            v_('K') = varargin{1};
        end
        v_('N') = 3;
        v_('M') = 5000;
        v_('RealM') = v_('M');
        v_('H') = 25.0/v_('RealM');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('ONE') = 1.0;
        v_('ROOTP5') = sqrt(0.5);
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for p=v_('1'):v_('K')
            for j=v_('1'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(p),',',int2str(j)],ix_);
                pb.xnames{iv} = ['X',int2str(p),',',int2str(j)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for p=v_('1'):v_('K')
            for I=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(p),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['OBJ',int2str(p),',',int2str(I)];
                pbm.gscale(ig,1) = v_('ROOTP5');
            end
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
        v_('SOLN1,1') = 19.0;
        v_('SOLN1,2') = 4.2;
        v_('SOLN1,3') = 1.2;
        v_('SOLN2,1') = 8.0;
        v_('SOLN2,2') = 2.5;
        v_('SOLN2,3') = 4.6;
        v_('SOLN3,1') = 10.0;
        v_('SOLN3,2') = 2.0;
        v_('SOLN3,3') = 2.6;
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('IH') = v_('H')*v_('RI');
            v_('TI') = v_('ONE')+v_('IH');
            v_('Differ') = v_('TI')-v_('SOLN1,2');
            v_('Numer') = v_('Differ')*v_('Differ');
            v_('Denom') = v_('SOLN1,3')*v_('SOLN1,3');
            v_('Differ') = v_('Numer')/v_('Denom');
            v_('Ratio') = 0.0-v_('Differ');
            v_('ERat') = exp(v_('Ratio'));
            v_('Yi1') = v_('SOLN1,1')*v_('ERat');
            v_('Differ') = v_('TI')-v_('SOLN2,2');
            v_('Numer') = v_('Differ')*v_('Differ');
            v_('Denom') = v_('SOLN2,3')*v_('SOLN2,3');
            v_('Differ') = v_('Numer')/v_('Denom');
            v_('Ratio') = 0.0-v_('Differ');
            v_('ERat') = exp(v_('Ratio'));
            v_('Yi2') = v_('SOLN2,1')*v_('ERat');
            v_('Differ') = v_('TI')-v_('SOLN3,2');
            v_('Numer') = v_('Differ')*v_('Differ');
            v_('Denom') = v_('SOLN3,3')*v_('SOLN3,3');
            v_('Differ') = v_('Numer')/v_('Denom');
            v_('Ratio') = 0.0-v_('Differ');
            v_('ERat') = exp(v_('Ratio'));
            v_('Yi3') = v_('SOLN3,1')*v_('ERat');
            for p=v_('1'):v_('K')
                pbm.gconst(ig_(['OBJ',int2str(p),',',int2str(I)])) = v_(['Yi',int2str(p)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        v_('LOWER1,1') = 15.0;
        v_('LOWER1,2') = 3.5;
        v_('LOWER1,3') = 0.3;
        v_('LOWER2,1') = 5.0;
        v_('LOWER2,2') = 2.2;
        v_('LOWER2,3') = 2.6;
        v_('LOWER3,1') = 5.0;
        v_('LOWER3,2') = 1.2;
        v_('LOWER3,3') = 1.3;
        v_('UPPER1,1') = 31.0;
        v_('UPPER1,2') = 6.3;
        v_('UPPER1,3') = 3.7;
        v_('UPPER2,1') = 15.0;
        v_('UPPER2,2') = 5.3;
        v_('UPPER2,3') = 6.2;
        v_('UPPER3,1') = 14.0;
        v_('UPPER3,2') = 3.3;
        v_('UPPER3,3') = 2.8;
        for p=v_('1'):v_('K')
            for j=v_('1'):v_('N')
                pb.xlower(ix_(['X',int2str(p),',',int2str(j)]),1) =...
                      v_(['LOWER',int2str(p),',',int2str(j)]);
                pb.xupper(ix_(['X',int2str(p),',',int2str(j)])) =...
                      v_(['UPPER',int2str(p),',',int2str(j)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('START1,1') = 25.0;
        v_('START1,2') = 5.2;
        v_('START1,3') = 3.2;
        v_('START2,1') = 7.0;
        v_('START2,2') = 4.1;
        v_('START2,3') = 3.6;
        v_('START3,1') = 11.6;
        v_('START3,2') = 1.9;
        v_('START3,3') = 2.2;
        for p=v_('1'):v_('K')
            for j=v_('1'):v_('N')
                pb.x0(ix_(['X',int2str(p),',',int2str(j)]),1) =...
                      v_(['START',int2str(p),',',int2str(j)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXPSQ',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        elftv{it}{3} = 'W';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for p=v_('1'):v_('K')
            for I=v_('1'):v_('M')
                ename = ['E',int2str(p),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eEXPSQ';
                ielftype(ie) = iet_('eEXPSQ');
                vname = ['X',int2str(p),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(p),',',int2str(round(v_('2')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(p),',',int2str(round(v_('3')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('W',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                v_('RI') = I;
                v_('IH') = v_('H')*v_('RI');
                v_('TI') = v_('ONE')+v_('IH');
                [~,posep] = ismember('T',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('TI');
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for p=v_('1'):v_('K')
            for I=v_('1'):v_('M')
                ig = ig_(['OBJ',int2str(p),',',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(p),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-V-V';
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

    case 'eEXPSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        R = (pbm.elpar{iel_}(1)-EV_(2))^2;
        S = EV_(3)^2;
        E = exp(-R/S);
        varargout{1} = EV_(1)*E;
        if(nargout>1)
            g_(1,1) = E;
            g_(2,1) = 2.0*(pbm.elpar{iel_}(1)-EV_(2))*EV_(1)*E/S;
            g_(3,1) = 2.0*R*EV_(1)*E/(S*EV_(3));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = 2.0*(pbm.elpar{iel_}(1)-EV_(2))*E/S;
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0*R*E/(S*EV_(3));
                H_(3,1) = H_(1,3);
                H_(2,2) = (2.0*EV_(1)*E/S)*(2.0*R/S-1.0);
                H_(2,3) = 4.0*(pbm.elpar{iel_}(1)-EV_(2))*EV_(1)*E/(S*EV_(3))*(R/S-1.0);
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*R*EV_(1)*E/(S^3)*(2.0*R-3.0*S);
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

