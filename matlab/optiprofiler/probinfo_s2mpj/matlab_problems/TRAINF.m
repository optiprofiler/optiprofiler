function varargout = TRAINF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TRAINF
%    *********
% 
%    This is an optimal control problem.
%    The problem is to minimize the energy spent to move a train 
%    from the beginning of a flat track to its end in a given time.  The train
%    is slowed down by some drag (assumed to be quadratic in the the velocity).
%    The control variables are the acceleration force (UA) and the braking
%    force (UB) applied on the train.
% 
%    Source:
%    J. Kautsky and N. K. Nichols,
%    "OTEP-2: Optimal Train Energy Programme, mark 2",
%    Numerical Analysis Report NA/4/83,
%    Department of Mathematics, University of Reading, 1983.
% 
%    SIF input: N. Nichols and Ph. Toint, April 1993
% 
%    classification = 'C-CQQR2-MN-V-V'
% 
%    Problem variants
% 
%       Alternative values for the SIF file parameters:
% RE TIME                4.8            $-PARAMETER  travel time
% RE LENGTH              6.0            $-PARAMETER  length of track
% 
% RE TIME                2.0            $-PARAMETER  travel time
% RE LENGTH              2.0            $-PARAMETER  length of track
% 
% RE TIME                1.5            $-PARAMETER  travel time
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TRAINF';

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
            v_('TIME') = 1.5;  %  SIF file default value
        else
            v_('TIME') = varargin{1};
        end
% RE LENGTH              2.0            $-PARAMETER  length of track
        if(nargs<2)
            v_('LENGTH') = 2;  %  SIF file default value
        else
            v_('LENGTH') = varargin{2};
        end
% IE N                   11             $-PARAMETER
% IE N                   51             $-PARAMETER
% IE N                   101            $-PARAMETER     original value
% IE N                   201            $-PARAMETER
% IE N                   501            $-PARAMETER
% IE N                   1001           $-PARAMETER
        if(nargs<3)
            v_('N') = 11;  %  SIF file default value
        else
            v_('N') = varargin{3};
        end
% IE N                   5001           $-PARAMETER
% IE N                   10001          $-PARAMETER
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = v_('TIME')/v_('RN');
        v_('H/2') = 0.5*v_('H');
        v_('-H') = -1.0*v_('H');
        v_('-H/2') = -1.0*v_('H/2');
        v_('UAMAX') = 10.0;
        v_('UBMIN') = -2.0;
        v_('VMAX') = 10.0;
        v_('A') = 0.3;
        v_('B') = 0.14;
        v_('C') = 0.16;
        v_('0') = 0;
        v_('1') = 1;
        v_('BH/2') = v_('B')*v_('H/2');
        v_('1+BH/2') = 1.0+v_('BH/2');
        v_('BH/2-1') = -1.0+v_('BH/2');
        v_('-AH') = v_('A')*v_('-H');
        v_('LENGTH/N') = v_('LENGTH')/v_('RN');
        v_('CH/2') = v_('C')*v_('H/2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['V',int2str(I)],ix_);
            pb.xnames{iv} = ['V',int2str(I)];
        end
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['UA',int2str(I)],ix_);
            pb.xnames{iv} = ['UA',int2str(I)];
        end
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['UB',int2str(I)],ix_);
            pb.xnames{iv} = ['UB',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','ENERGY',ig_);
        gtype{ig} = '<>';
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['XEQ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['XEQ',int2str(I)];
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['V',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
            end
            iv = ix_(['V',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
            end
            [ig,ig_] = s2mpjlib('ii',['VEQ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['VEQ',int2str(I)];
            iv = ix_(['V',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1+BH/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1+BH/2');
            end
            iv = ix_(['V',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('BH/2-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('BH/2-1');
            end
            iv = ix_(['UA',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
            end
            iv = ix_(['UA',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
            end
            iv = ix_(['UB',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
            end
            iv = ix_(['UB',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-H/2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-H/2');
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
        for I=v_('0'):v_('N-1')
            pbm.gconst(ig_(['VEQ',int2str(I)])) = v_('-AH');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['V',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['V',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['UA',int2str(round(v_('0')))]),1) = v_('UAMAX');
        pb.xupper(ix_(['UA',int2str(round(v_('0')))]),1) = v_('UAMAX');
        pb.xlower(ix_(['UB',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['UB',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('N-1')
            pb.xlower(ix_(['X',int2str(I)])) = -Inf;
            pb.xupper(ix_(['X',int2str(I)]),1) = +Inf;
            pb.xlower(ix_(['V',int2str(I)])) = -Inf;
            pb.xupper(ix_(['V',int2str(I)]),1) = +Inf;
            pb.xlower(ix_(['UA',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['UA',int2str(I)])) = v_('UAMAX');
            pb.xlower(ix_(['UB',int2str(I)]),1) = v_('UBMIN');
            pb.xupper(ix_(['UB',int2str(I)])) = 0.0;
        end
        pb.xlower(ix_(['X',int2str(round(v_('N')))]),1) = v_('LENGTH');
        pb.xupper(ix_(['X',int2str(round(v_('N')))]),1) = v_('LENGTH');
        pb.xlower(ix_(['V',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['V',int2str(round(v_('N')))]),1) = 0.0;
        pb.xlower(ix_(['UA',int2str(round(v_('N')))]),1) = 0.0;
        pb.xupper(ix_(['UA',int2str(round(v_('N')))]),1) = 0.0;
        pb.xlower(ix_(['UB',int2str(round(v_('N')))]),1) = v_('UBMIN');
        pb.xupper(ix_(['UB',int2str(round(v_('N')))]),1) = v_('UBMIN');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['V',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['UA',int2str(round(v_('0')))]),1) = v_('UAMAX');
        pb.x0(ix_(['UB',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('N-1')
            v_('RI') = I;
            v_('PI') = v_('LENGTH/N')*v_('RI');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('PI');
            pb.x0(ix_(['V',int2str(I)]),1) = v_('LENGTH/N');
            pb.x0(ix_(['UA',int2str(I)]),1) = 0.0;
            pb.x0(ix_(['UB',int2str(I)]),1) = 0.0;
        end
        pb.x0(ix_(['X',int2str(round(v_('N')))]),1) = v_('LENGTH');
        pb.x0(ix_(['V',int2str(round(v_('N')))]),1) = 0.0;
        pb.x0(ix_(['UA',int2str(round(v_('N')))]),1) = 0.0;
        pb.x0(ix_(['UB',int2str(round(v_('N')))]),1) = v_('UBMIN');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'UU';
        elftv{it}{2} = 'VV';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'VVV';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('N')
            v_('I+1') = 1+I;
            ename = ['VISQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['V',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VVV',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('N-1')
            ename = ['UV',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = ['UA',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('UU',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['V',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VV',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            ig = ig_(['VEQ',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['VISQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('CH/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['VISQ',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('CH/2');
        end
        for I=v_('1'):v_('N-1')
            ig = ig_('ENERGY');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['UV',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION            3.09751881012
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQQR2-MN-V-V';
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

% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

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

    case 'ePROD'

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

