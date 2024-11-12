function varargout = OPTMASS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OPTMASS
%    *********
% 
%    A constrained optimal control problem
%    adapted from Gawande and Dunn
% 
%    The problem is that of a particle of unit mass moving on a
%    frictionless plane under the action of a controlling force whose
%    magnitude may not exceed unity. At time=0, the particle moves through
%    the origin of the plane in the direction of the positive x-axis with
%    speed SPEED.  The cost function incorporates two conflicting control
%    objectives, namely: maximization of the particle's final (at time=1)
%    distance from the origin and minimization of its final speed.  By
%    increasing the  value of the penalty constant PEN, more stress can be
%    placed on the latter objective.
% 
%    Gawande and Dunn originally use a starting point (in the control
%    only) that is much closer to the solution than the one chosen
%    here.
% 
%    Source:
%    M. Gawande and J. Dunn,
%    "A Projected Newton Method in a Cartesian Product of Balls",
%    JOTA 59(1): 59-69, 1988.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'C-CQQR2-AN-V-V'
% 
%    Number of discretization steps in the time interval
%    The number of variables is 6 * (N + 2) -2 , 4 of which are fixed.
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER n = 70    original value
% IE N                   100            $-PARAMETER n = 610
% IE N                   200            $-PARAMETER n = 1210
% IE N                   500            $-PARAMETER n = 3010
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OPTMASS';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   1000           $-PARAMETER n = 6010
% IE N                   5000           $-PARAMETER n = 30010
        v_('SPEED') = 0.01;
        v_('PEN') = 0.335;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N+1') = 1+v_('N');
        v_('RN') = v_('N');
        v_('1/N') = 1.0/v_('RN');
        v_('-1/N') = -1.0*v_('1/N');
        v_('1/N2') = v_('1/N')*v_('1/N');
        v_('-1/2N2') = -0.5*v_('1/N2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            for J=v_('1'):v_('2')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['X',int2str(J),',',int2str(I)];
                [iv,ix_] = s2mpjlib('ii',['V',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['V',int2str(J),',',int2str(I)];
                [iv,ix_] = s2mpjlib('ii',['F',int2str(J),',',int2str(I)],ix_);
                pb.xnames{iv} = ['F',int2str(J),',',int2str(I)];
            end
        end
        for J=v_('1'):v_('2')
            [iv,ix_] =...
                  s2mpjlib('ii',['X',int2str(J),',',int2str(round(v_('N+1')))],ix_);
            pb.xnames{iv} = ['X',int2str(J),',',int2str(round(v_('N+1')))];
            [iv,ix_] =...
                  s2mpjlib('ii',['V',int2str(J),',',int2str(round(v_('N+1')))],ix_);
            pb.xnames{iv} = ['V',int2str(J),',',int2str(round(v_('N+1')))];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','F',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N+1')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('2')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(J),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['A',int2str(J),',',int2str(I)];
                iv = ix_(['X',int2str(J),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(J),',',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['V',int2str(J),',',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-1/N')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-1/N');
                end
                iv = ix_(['F',int2str(J),',',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-1/2N2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-1/2N2');
                end
                [ig,ig_] = s2mpjlib('ii',['B',int2str(J),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['B',int2str(J),',',int2str(I)];
                iv = ix_(['V',int2str(J),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['V',int2str(J),',',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['F',int2str(J),',',int2str(round(v_('I-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-1/N')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-1/N');
                end
            end
        end
        for I=v_('0'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(I)];
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
        for I=v_('0'):v_('N')
            pbm.gconst(ig_(['C',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['V',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = v_('SPEED');
        pb.xupper(ix_(['V',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = v_('SPEED');
        pb.xlower(ix_(['V',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['V',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        if(isKey(ix_,['V',int2str(round(v_('1'))),',',int2str(round(v_('0')))]))
            pb.x0(ix_(['V',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) =...
                  v_('SPEED');
        else
            pb.y0(find(pbm.congrps==ig_(['V',int2str(round(v_('1'))),',',int2str(round(v_('0')))])),1) = v_('SPEED');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'O1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = ['X',int2str(round(v_('1'))),',',int2str(round(v_('N+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = ['X',int2str(round(v_('2'))),',',int2str(round(v_('N+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = ['V',int2str(round(v_('1'))),',',int2str(round(v_('N+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = ['V',int2str(round(v_('2'))),',',int2str(round(v_('N+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('0'):v_('N')
            for J=v_('1'):v_('2')
                ename = ['D',int2str(J),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = ['F',int2str(J),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('F');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('O2');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('PEN');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('PEN');
        for I=v_('0'):v_('N')
            ig = ig_(['C',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('1'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('2'))),',',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           -0.04647
% LO SOLTN(100)          ???
% LO SOLTN(200)          ???
% LO SOLTN(500)          ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQQR2-AN-V-V';
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

    case 'eSQ'

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

