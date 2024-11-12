function varargout = HANGING(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HANGING
%    *********
% 
%    A catenary problem in 3 dimensions.  A rectangular grid is hung from its
%    4 corners under gravity.  The problem is to determine the resulting shape.
% 
%    Source:  
%    an example in a talk by Nesterova and Vial, LLN, 1994.
% 
%    SIF input: Ph. Toint, November 1994.
% 
%    classification = 'C-CLQR2-AY-V-V'
% 
%    dimension of the grid
% 
%       Alternative values for the SIF file parameters:
% IE NX                  3              $-PARAMETER n = 27
% IE NY                  3              $-PARAMETER
% 
% IE NX                  5              $-PARAMETER n = 90
% IE NY                  6              $-PARAMETER
% 
% IE NX                  10             $-PARAMETER n = 300  original value
% IE NY                  10             $-PARAMETER
% 
% IE NX                  20             $-PARAMETER n = 1800
% IE NY                  30             $-PARAMETER
% 
% IE NX                  40             $-PARAMETER n = 3600
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HANGING';

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
            v_('NX') = 3;  %  SIF file default value
        else
            v_('NX') = varargin{1};
        end
% IE NY                  30             $-PARAMETER
        if(nargs<2)
            v_('NY') = 3;  %  SIF file default value
        else
            v_('NY') = varargin{2};
        end
        v_('LX') = 1.8;
        v_('LY') = 1.8;
        v_('1') = 1;
        v_('NX-1') = -1+v_('NX');
        v_('NY-1') = -1+v_('NY');
        v_('LX2') = v_('LX')*v_('LX');
        v_('LY2') = v_('LY')*v_('LY');
        v_('RNX') = v_('NX');
        v_('RNY') = v_('NY');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['Z',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Z',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
                gtype{ig} = '<>';
                iv = ix_(['Z',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY-1')
                [ig,ig_] = s2mpjlib('ii',['RC',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['RC',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('NX-1')
            for J=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['DC',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['DC',int2str(I),',',int2str(J)];
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
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY-1')
                pbm.gconst(ig_(['RC',int2str(I),',',int2str(J)])) = v_('LX2');
            end
        end
        for I=v_('1'):v_('NX-1')
            for J=v_('1'):v_('NY')
                pbm.gconst(ig_(['DC',int2str(I),',',int2str(J)])) = v_('LY2');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Z',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Z',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = v_('RNX');
        pb.xupper(ix_(['X',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = v_('RNX');
        pb.xlower(ix_(['Y',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Z',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Z',int2str(round(v_('NX'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = v_('RNY');
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = v_('RNY');
        pb.xlower(ix_(['Z',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        pb.xupper(ix_(['Z',int2str(round(v_('1'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = v_('RNX');
        pb.xupper(ix_(['X',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = v_('RNX');
        pb.xlower(ix_(['Y',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = v_('RNY');
        pb.xupper(ix_(['Y',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = v_('RNY');
        pb.xlower(ix_(['Z',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        pb.xupper(ix_(['Z',int2str(round(v_('NX'))),',',int2str(round(v_('NY')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('NX')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            for J=v_('1'):v_('NY')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('RI-1');
                pb.x0(ix_(['Y',int2str(I),',',int2str(J)]),1) = v_('RJ-1');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('NY-1')
            v_('J+1') = 1+J;
            for I=v_('1'):v_('NX')
                ename = ['RX',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['RY',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['RZ',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['Z',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('1'):v_('NX-1')
            v_('I+1') = 1+I;
            for J=v_('1'):v_('NY')
                ename = ['DX',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DY',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['DZ',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eISQ';
                    ielftype(ie) = iet_('eISQ');
                end
                vname = ['Z',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YY',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY-1')
                ig = ig_(['RC',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RX',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['RY',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RZ',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for I=v_('1'):v_('NX-1')
            for J=v_('1'):v_('NY')
                ig = ig_(['DC',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['DX',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['DY',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['DZ',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(3,3)          -6.1184107487
% LO SOLTN(5,6)          -77.260229515
% LO SOLTN(10,10)        -620.17603242
% LO SOLTN(20,30)        -1025.4292887
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLQR2-AY-V-V';
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

    case 'eISQ'

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

