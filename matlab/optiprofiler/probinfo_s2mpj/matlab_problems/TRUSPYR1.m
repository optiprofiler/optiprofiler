function varargout = TRUSPYR1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    This is a structural optimization problem.
%    The problem is to minimize the weight of a given
%    8-bar truss structure formed as a pyramid for a given external load.
%    There are upper bounds on the strain energy and lower bounds
%    on the cross-sectional areas of the bars.
% 
%    Source:
%    K. Svanberg, 
%    "Local and global optima",
%    Proceedings of the NATO/DFG ASI on Optimization of large structural
%    systems, 
%    G. I. N. Rozvany, ed., Kluwer, 1993, pp. 579-588.
% 
%    SIF input: A. Forsgren, Royal Institute of Technology, December 1993.
% 
%    classification = 'LQR2-MN-11-4'
% 
%    Number of bars
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TRUSPYR1';

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
        v_('NBAR') = 8;
        v_('NDIM') = 3;
        v_('1') = 1;
        v_('2') = 2;
        v_('NBAR/2') = fix(v_('NBAR')/v_('2'));
        v_('8.0') = 8.0;
        v_('SQRT17') = sqrt(17.0);
        v_('SQRT18') = sqrt(18.0);
        v_('SQRT105') = sqrt(105.0);
        v_('P1') = 40.0;
        v_('P2') = 20.0;
        v_('P3') = 200.0;
        v_('Q1') = 2.0/v_('SQRT105');
        v_('Q2') = 1.0/v_('SQRT105');
        v_('Q3') = 10.0/v_('SQRT105');
        v_('ALPHA') = 0.3291726437;
        for J=v_('1'):v_('NBAR/2')
            v_(['L',int2str(J)]) = v_('SQRT17')/v_('8.0');
            v_('J+4') = J+v_('NBAR/2');
            v_(['L',int2str(round(v_('J+4')))]) = v_('SQRT18')/v_('8.0');
        end
        v_('E') = 21.0;
        v_('R1,1') = 0.250;
        v_('R2,1') = 0.250;
        v_('R3,1') = 0.375;
        v_('R1,2') = 0.250;
        v_('R2,2') = -0.250;
        v_('R3,2') = 0.375;
        v_('R1,3') = -0.250;
        v_('R2,3') = -0.250;
        v_('R3,3') = 0.375;
        v_('R1,4') = -0.250;
        v_('R2,4') = 0.250;
        v_('R3,4') = 0.375;
        v_('R1,5') = 0.375;
        v_('R2,5') = 0.000;
        v_('R3,5') = 0.375;
        v_('R1,6') = 0.000;
        v_('R2,6') = -0.375;
        v_('R3,6') = 0.375;
        v_('R1,7') = -0.375;
        v_('R2,7') = 0.000;
        v_('R3,7') = 0.375;
        v_('R1,8') = 0.000;
        v_('R2,8') = 0.375;
        v_('R3,8') = 0.375;
        for J=v_('1'):v_('NBAR')
            v_(['L2',int2str(J)]) = v_(['L',int2str(J)])*v_(['L',int2str(J)]);
            v_(['L3',int2str(J)]) = v_(['L2',int2str(J)])*v_(['L',int2str(J)]);
            v_(['GAMMA',int2str(J)]) = v_('E')/v_(['L3',int2str(J)]);
            v_(['DL2',int2str(J)]) = v_(['L2',int2str(J)])/v_('E');
            v_(['W',int2str(J)]) = 0.78*v_(['L',int2str(J)]);
            v_(['STRUP',int2str(J)]) = 10.0*v_(['DL2',int2str(J)]);
            for I=v_('1'):v_('NDIM')
                v_(['RG',int2str(I),',',int2str(J)]) =...
                      v_(['GAMMA',int2str(J)])*v_(['R',int2str(I),',',int2str(J)]);
                for K=v_('1'):v_('NDIM')
                    v_(['RR',int2str(I),',',int2str(J),',',int2str(K)]) =...
                          v_(['RG',int2str(I),',',int2str(J)])*v_(['R',int2str(K),',',int2str(J)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('NBAR')
            [iv,ix_] = s2mpjlib('ii',['XAREA',int2str(J)],ix_);
            pb.xnames{iv} = ['XAREA',int2str(J)];
        end
        for I=v_('1'):v_('NDIM')
            [iv,ix_] = s2mpjlib('ii',['DISPL',int2str(I)],ix_);
            pb.xnames{iv} = ['DISPL',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('NBAR')
            [ig,ig_] = s2mpjlib('ii','WEIGHT',ig_);
            gtype{ig} = '<>';
            iv = ix_(['XAREA',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['W',int2str(J)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['W',int2str(J)]);
            end
        end
        for K=v_('1'):v_('NDIM')
            [ig,ig_] = s2mpjlib('ii',['EQUIL',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EQUIL',int2str(K)];
        end
        for I=v_('1'):v_('NDIM')
            [ig,ig_] = s2mpjlib('ii','STREN',ig_);
            gtype{ig}  = '<=';
            cnames{ig} = 'STREN';
            iv = ix_(['DISPL',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['Q',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['Q',int2str(I)]);
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
        for K=v_('1'):v_('NDIM')
            pbm.gconst(ig_(['EQUIL',int2str(K)])) = v_(['P',int2str(K)]);
        end
        pbm.gconst(ig_('STREN')) = v_('ALPHA');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for J=v_('1'):v_('NBAR')
            pb.xlower(ix_(['XAREA',int2str(J)]),1) = 1.0;
        end
        for I=v_('1'):v_('NDIM')
            pb.xlower(ix_(['DISPL',int2str(I)])) = -Inf;
            pb.xupper(ix_(['DISPL',int2str(I)]),1) = +Inf;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NDIM')
            for J=v_('1'):v_('NBAR')
                ename = ['UX',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['DISPL',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['XAREA',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NDIM')
            for J=v_('1'):v_('NBAR')
                for K=v_('1'):v_('NDIM')
                    ig = ig_(['EQUIL',int2str(K)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['UX',int2str(I),',',int2str(J)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['RR',int2str(I),',',int2str(J),',',int2str(K)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Objective function value corresponding to the global minimizer above
        pb.objlower = 1.2287408808;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-MN-11-4';
        pb.x0          = zeros(pb.n,1);
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

