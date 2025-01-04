function varargout = KISSING2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem: A second formulation of the KISSING NUMBER PROBLEM
%                                                                    
%    Source: This problem is associated to the family of Hard-Spheres 
%    problem. It belongs to the family of sphere packing problems, a 
%    class of challenging problems dating from the beginning of the 
%    17th century which is related to practical problems in Chemistry, 
%    Biology and Physics. Given a fixed unit sphere at the origin in R^n, 
%    the problem consists of arranging a further m unit spheres so that 
%    sum of the distances to these spheres is as small as possible.
%    This problem may be reduced to a nonconvex nonlinear optimization 
%    problem with a potentially large number of (nonoptimal) points 
%    satisfying optimality conditions. We have, thus, a class of problems 
%    indexed by the parameters m and n, that provides a suitable 
%    set of test problems for evaluating nonlinear programming codes.
%    After some algebric manipulations, we can formulate this problem as
%               m
%     Minimize sum <p_i,p_i> - m n
%              i=1
%     subject to
%        
%      <p_i - p_j, p_i - p_j> >= 4 for all different pair of indices i, j
%      and  
%      <p_i, p_i> >= 4 for all indices i
%   
%      as well as n(n-1)/2 normalisation constraints fixing components.
%      The goal is to find an objective value equal to 0.
%      [1]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
%            N. J. C. Sloane, Springer-Verlag, NY, 1988.
%    SIF input: Nick Gould, September 2000
% 
%    classification = 'C-CQQR2-RN-V-V'
% 
% **********************************************************************
% 
%    Number of points: m
% 
%       Alternative values for the SIF file parameters:
% IE m                   24             $-PARAMETER  number of points
% IE m                   25             $-PARAMETER  number of points
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'KISSING2';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
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
        if(nargs<1)
            v_('m') = 25;  %  SIF file default value
        else
            v_('m') = varargin{1};
        end
% IE m                   100            $-PARAMETER  number of points
% IE n                    4             $-PARAMETER  dimension of sphere
        if(nargs<2)
            v_('n') = 4;  %  SIF file default value
        else
            v_('n') = varargin{2};
        end
% IE n                    8             $-PARAMETER  dimension of sphere
        v_('1') = 1;
        v_('2') = 2;
        v_('n-1') = v_('n')-v_('1');
        v_('rm') = v_('m');
        v_('rn') = v_('n');
        v_('RM+N') = v_('rm')+v_('rn');
        v_('mn') = v_('rm')*v_('rn');
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('PI/m') = v_('PI')/v_('rm');
        v_('2PI/m') = 2.0*v_('PI/m');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('m')
            for J=v_('1'):v_('n')
                [iv,ix_] = s2mpjlib('ii',['P',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['P',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('m')
            for J=v_('1'):v_('m')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(I),',',int2str(J)];
            end
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
        pbm.gconst(ig_('OBJ')) = v_('mn');
        for I=v_('1'):v_('m')
            for J=v_('1'):v_('m')
                pbm.gconst(ig_(['C',int2str(I),',',int2str(J)])) = 4.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('m')
            for J=v_('1'):v_('n')
                pb.xlower(ix_(['P',int2str(I),',',int2str(J)])) = -Inf;
                pb.xupper(ix_(['P',int2str(I),',',int2str(J)]),1) = +Inf;
            end
        end
        for I=v_('2'):v_('n')
            for J=I:v_('n')
                pb.xlower(ix_(['P',int2str(I),',',int2str(J)]),1) = 0.0;
                pb.xupper(ix_(['P',int2str(I),',',int2str(J)]),1) = 0.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('m')
            v_('RI') = I;
            v_('2PIi/m') = v_('2PI/m')*v_('RI');
            v_('cos') = cos(v_('2PIi/m'));
            v_('sin') = sin(v_('2PIi/m'));
            v_('cos') = v_('cos')+v_('cos');
            v_('sin') = v_('sin')+v_('sin');
            if(isKey(ix_,['P',int2str(I),',',int2str(round(v_('1')))]))
                pb.x0(ix_(['P',int2str(I),',',int2str(round(v_('1')))]),1) = v_('cos');
            else
                pb.y0(find(pbm.congrps==ig_(['P',int2str(I),',',int2str(round(v_('1')))])),1) = v_('cos');
            end
            for J=v_('2'):v_('n-1')
                if(isKey(ix_,['P',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['P',int2str(I),',',int2str(J)]),1) = v_('sin');
                else
                    pb.y0(find(pbm.congrps==ig_(['P',int2str(I),',',int2str(J)])),1) =...
                          v_('sin');
                end
            end
        end
        if(isKey(ix_,['P',int2str(round(v_('m'))),',',int2str(round(v_('n')))]))
            pb.x0(ix_(['P',int2str(round(v_('m'))),',',int2str(round(v_('n')))]),1) =...
                  v_('cos');
        else
            pb.y0(find(pbm.congrps==ig_(['P',int2str(round(v_('m'))),',',int2str(round(v_('n')))])),1) = v_('cos');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD1',iet_);
        elftv{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD2',iet_);
        elftv{it}{1} = 'Q';
        elftv{it}{2} = 'R';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('m')
            v_('I-') = -1+I;
            v_('I+') = 1+I;
            for J=v_('1'):v_('I-')
                for K=v_('1'):v_('n')
                    ename = ['E',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'ePROD2';
                    ielftype(ie) = iet_('ePROD2');
                    vname = ['P',int2str(I),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('Q',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['P',int2str(J),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('R',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
            for K=v_('1'):v_('n')
                ename = ['E',int2str(I),',',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD1';
                ielftype(ie) = iet_('ePROD1');
                vname = ['P',int2str(I),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('P',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for J=v_('I+'):v_('m')
                for K=v_('1'):v_('n')
                    ename = ['E',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'ePROD2';
                    ielftype(ie) = iet_('ePROD2');
                    vname = ['P',int2str(I),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('Q',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['P',int2str(J),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('R',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('m')
            for K=v_('1'):v_('n')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
            for J=v_('1'):v_('m')
                for K=v_('1'):v_('n')
                    ig = ig_(['C',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J),',',int2str(K)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION             0.00000D+00   $ n=4, m = 24
% XL SOLUTION             6.48030D+00   $ n=4, m = 25 one of many local solutions
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQQR2-RN-V-V';
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

    case 'ePROD1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'ePROD2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
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

