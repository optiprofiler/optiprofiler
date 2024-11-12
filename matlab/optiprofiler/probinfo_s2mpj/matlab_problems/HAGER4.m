function varargout = HAGER4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HAGER4
%    *********
% 
%    A nonlinear optimal control problem, by W. Hager.
% 
%    NOTE: The solution for x given in the article below by Hager has
%    a typo. On the interval [1/2, 1], x(t) = (exp(2t) + exp(t))/d. In
%    other words, the minus sign in the article should be a plus sign.
% 
%    Source: problem P4 in
%    W.W. Hager,
%    "Multiplier Methods for Nonlinear Optimal Control",
%    SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.
% 
%    SIF input: Ph. Toint, April 1991.
% 
%    classification = 'C-COLR2-AN-V-V'
% 
%    Number of discretized points in [0,1]
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   2500           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HAGER4';

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
% IE N                   5000           $-PARAMETER
        v_('1/H') = v_('N');
        v_('H') = 1.0/v_('1/H');
        v_('H/2') = 0.5*v_('H');
        v_('1/H-1') = -1.0+v_('1/H');
        v_('-1/H') = -1.0*v_('1/H');
        v_('1/HSQ') = v_('1/H')*v_('1/H');
        v_('1/2HSQ') = 0.5*v_('1/HSQ');
        v_('0') = 0;
        v_('1') = 1;
        for I=v_('0'):v_('N')
            v_('RI') = I;
            v_(['T',int2str(I)]) = v_('RI')*v_('H');
            v_('-2TI') = -2.0*v_(['T',int2str(I)]);
            v_(['Z',int2str(I)]) = exp(v_('-2TI'));
        end
        for I=v_('0'):v_('1')
            v_(['A',int2str(I)]) = -0.5*v_(['Z',int2str(I)]);
            v_('TI+1/2') = 0.5+v_(['T',int2str(I)]);
            v_(['B',int2str(I)]) = v_(['A',int2str(I)])*v_('TI+1/2');
            v_('TISQ') = v_(['T',int2str(I)])*v_(['T',int2str(I)]);
            v_('TIETC') = v_('TISQ')+v_('TI+1/2');
            v_(['C',int2str(I)]) = v_(['A',int2str(I)])*v_('TIETC');
        end
        v_('DA') =...
              v_(['A',int2str(round(v_('1')))])-v_(['A',int2str(round(v_('0')))]);
        v_('SCDA') = 0.5*v_('DA');
        v_('DB') =...
              v_(['B',int2str(round(v_('1')))])-v_(['B',int2str(round(v_('0')))]);
        v_('SCDB') = v_('DB')*v_('1/H');
        v_('DC') =...
              v_(['C',int2str(round(v_('1')))])-v_(['C',int2str(round(v_('0')))]);
        v_('SCDC') = v_('DC')*v_('1/2HSQ');
        v_('E') = exp(1.0);
        v_('3E') = 3.0*v_('E');
        v_('1+3E') = 1.0+v_('3E');
        v_('1-E') = 1.0-v_('E');
        v_('2-2E') = 2.0*v_('1-E');
        v_('XX0') = v_('1+3E')/v_('2-2E');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['S',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['S',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H-1');
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            v_('ETI') = exp(v_(['T',int2str(I)]));
            v_('-ETI') = -1.0*v_('ETI');
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-ETI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-ETI');
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = v_('XX0');
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = v_('XX0');
        for I=v_('1'):v_('N')
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('0')))]),1) = v_('XX0');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eELT',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'D';
        elftp{it}{2} = 'E';
        elftp{it}{3} = 'F';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            ename = ['EL',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eELT';
            ielftype(ie) = iet_('eELT');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('DD') = v_('SCDA')*v_(['Z',int2str(round(v_('I-1')))]);
            v_('EE') = v_('SCDB')*v_(['Z',int2str(round(v_('I-1')))]);
            v_('FF') = v_('SCDC')*v_(['Z',int2str(round(v_('I-1')))]);
            [~,posep] = ismember('D',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('DD');
            [~,posep] = ismember('E',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('EE');
            [~,posep] = ismember('F',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('FF');
            ename = ['U',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EL',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['U',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           2.833914199
% LO SOLTN(50)           2.799810928
% LO SOLTN(100)          2.796761851
% LO SOLTN(500)          2.794513229
% LO SOLTN(1000)         2.794244187
% LO SOLTN(5000)         ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AN-V-V';
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

    case 'eELT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(1)+pbm.elpar{iel_}(2)*...
             EV_(1)*(EV_(2)-EV_(1))+pbm.elpar{iel_}(3)*(EV_(2)-EV_(1))^2;
        if(nargout>1)
            g_(1,1) = 2.0*pbm.elpar{iel_}(1)*EV_(1)+pbm.elpar{iel_}(2)*(EV_(2)-2.0*EV_(1))-...
                 2.0*pbm.elpar{iel_}(3)*(EV_(2)-EV_(1));
            g_(2,1) = pbm.elpar{iel_}(2)*EV_(1)+2.0*pbm.elpar{iel_}(3)*(EV_(2)-EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*(pbm.elpar{iel_}(1)-pbm.elpar{iel_}(2)+pbm.elpar{iel_}(3));
                H_(1,2) = pbm.elpar{iel_}(2)-2.0*pbm.elpar{iel_}(3);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*pbm.elpar{iel_}(3);
                varargout{3} = H_;
            end
        end

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

