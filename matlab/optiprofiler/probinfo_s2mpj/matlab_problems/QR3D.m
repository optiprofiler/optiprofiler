function varargout = QR3D(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : QR3D
%    *********
% 
%    Find the QR factorization of a tridiagonal matrix A.
%    The problem is formulated as a system of quadratic equations
%    whose unknowns are the elements of the orthogonal matrix Q and of
%    the upper triangular matrix R.  In this version of the problem,
%    the banded structure of R is not imposed as a constraint. See problem
%    QR3DBD for the case where this structure is explicitly used.
% 
%    The problem is non-convex.
% 
%    Source:
%    Ph. Toint, private communication.
% 
%    SIF input: Ph. Toint, Nov 1993
% 
%    classification = 'C-CNQR2-AN-V-V'
% 
%    Define the matrix order M  ( M >= 3 ).
%    There are M * ( 3M + 1) / 2 variables and equations.
% 
%       Alternative values for the SIF file parameters:
% IE M                   5              $-PARAMETER  n =  40
% IE M                   10             $-PARAMETER  n = 155  original value
% IE M                   20             $-PARAMETER  n = 610
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'QR3D';

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
        if(nargs<1)
            v_('M') = 5;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('M-1') = -1+v_('M');
        v_('RM') = v_('M');
        v_('2/M') = 2.0/v_('RM');
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('1')))]) = v_('2/M');
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 0.0;
        for I=v_('2'):v_('M-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('1-I') = -1*v_('I-1');
            v_('R1-I') = v_('1-I');
            v_('1-I/M') = v_('R1-I')/v_('RM');
            v_('2I') = 2*I;
            v_('R2I') = v_('2I');
            v_('2I/M') = v_('R2I')/v_('RM');
            v_(['A',int2str(I),',',int2str(round(v_('I-1')))]) = v_('1-I/M');
            v_(['A',int2str(I),',',int2str(I)]) = v_('2I/M');
            v_(['A',int2str(I),',',int2str(round(v_('I+1')))]) = v_('1-I/M');
        end
        v_('RM-1') = v_('M-1');
        v_('1-M') = -1.0*v_('RM-1');
        v_('1-M/M') = v_('1-M')/v_('RM');
        v_('2M') = 2.0*v_('RM');
        v_(['A',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]) =...
              v_('1-M/M');
        v_(['A',int2str(round(v_('M'))),',',int2str(round(v_('M')))]) = v_('2M');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('M')
                [iv,ix_] = s2mpjlib('ii',['Q',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Q',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('M')
            for J=I:v_('M')
                [iv,ix_] = s2mpjlib('ii',['R',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['R',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            for J=I:v_('M')
                [ig,ig_] = s2mpjlib('ii',['O',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['O',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['F',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['F',int2str(I),',',int2str(J)];
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['O',int2str(I),',',int2str(I)])) = 1.0;
        end
        pbm.gconst(ig_(['F',int2str(round(v_('1'))),',',int2str(round(v_('1')))])) =...
              v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
        pbm.gconst(ig_(['F',int2str(round(v_('1'))),',',int2str(round(v_('2')))])) =...
              v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
        for I=v_('2'):v_('M-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            pbm.gconst(ig_(['F',int2str(I),',',int2str(round(v_('I-1')))])) =...
                  v_(['A',int2str(I),',',int2str(round(v_('I-1')))]);
            pbm.gconst(ig_(['F',int2str(I),',',int2str(I)])) =...
                  v_(['A',int2str(I),',',int2str(I)]);
            pbm.gconst(ig_(['F',int2str(I),',',int2str(round(v_('I+1')))])) =...
                  v_(['A',int2str(I),',',int2str(round(v_('I+1')))]);
        end
        pbm.gconst(ig_(['F',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))])) = v_(['A',int2str(round(v_('M'))),',',int2str(round(v_('M-1')))]);
        pbm.gconst(ig_(['F',int2str(round(v_('M'))),',',int2str(round(v_('M')))])) =...
              v_(['A',int2str(round(v_('M'))),',',int2str(round(v_('M')))]);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('M')
            pb.xlower(ix_(['R',int2str(I),',',int2str(I)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('M')
            pb.x0(ix_(['Q',int2str(I),',',int2str(I)]),1) = 1.0;
        end
        for I=v_('1'):v_('M-1')
            v_('I+1') = 1+I;
            pb.x0(ix_(['R',int2str(I),',',int2str(I)]),1) =...
                  v_(['A',int2str(I),',',int2str(I)]);
            pb.x0(ix_(['R',int2str(I),',',int2str(round(v_('I+1')))]),1) =...
                  v_(['A',int2str(I),',',int2str(round(v_('I+1')))]);
        end
        pb.x0(ix_(['R',int2str(round(v_('M'))),',',int2str(round(v_('M')))]),1) =...
              v_(['A',int2str(round(v_('M'))),',',int2str(round(v_('M')))]);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('M')
            for J=I:v_('M')
                for K=v_('1'):v_('M')
                    ename = ['C',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                    if(newelt)
                        pbm.elftype{ie} = 'en2PR';
                        ielftype(ie) = iet_('en2PR');
                    end
                    vname = ['Q',int2str(I),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['Q',int2str(J),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('M')
                for K=v_('1'):J
                    ename = ['B',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                    if(newelt)
                        pbm.elftype{ie} = 'en2PR';
                        ielftype(ie) = iet_('en2PR');
                    end
                    vname = ['Q',int2str(I),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['R',int2str(K),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('V2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('M')
            for J=I:v_('M')
                for K=v_('1'):v_('M')
                    ig = ig_(['O',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J),',',int2str(K)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('M')
                for K=v_('1'):J
                    ig = ig_(['F',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J),',',int2str(K)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNQR2-AN-V-V';
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

