function varargout = BLOWEYC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BLOWEYC
%    *********
% 
%    A nonconvex quadratic program proposed by 
%    James Blowey (University of Durham)
% 
%    Given function v(s) and f(s) = v(s) + A(inv) v(s), s in [0,1],
%    minimize 
% 
%         (u(s) - v(s))(trans) ( A + A(inv) ) u(s) - (u(s) - v(s))(trans)f(s)
% 
%    where 
% 
%       u(s) in [-1,1] and int[0,1] u(s) ds = int[0,1] v(s) ds
% 
%    and A is the 
% 
%       "- Laplacian with Neumann boundary conditions on a uniform mesh"
% 
%    The troublesome term A(inv) u(s) is replaced by the additional 
%    variable w(s) and the constraint A w(s) = u(s)
% 
%    The function v(s) is chosen to be 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BLOWEYC';

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
%    classification = 'C-CQLR2-MN-V-V'
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER  n = 22, m = 12
% IE N                   100            $-PARAMETER  n = 202, m = 102
% IE N                   1000           $-PARAMETER  n = 2002, m = 1002
% IE N                   2000           $-PARAMETER  n = 4002, m = 2002
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   4000           $-PARAMETER  n = 8002, m = 4002
% IE N                   8000           $-PARAMETER  n = 16002, m = 8002
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('ONE') = 1.0;
        v_('-ONE') = -1.0;
        v_('TWO') = 2.0;
        v_('-TWO') = -2.0;
        v_('RN') = v_('N');
        v_('N**2') = v_('RN')*v_('RN');
        v_('N-1') = -1+v_('N');
        v_('1/N**2') = v_('ONE')/v_('N**2');
        v_('-1/N**2') = v_('-ONE')/v_('N**2');
        v_('-2/N**2') = 2.0*v_('-1/N**2');
        v_('2N**2') = 2.0*v_('N**2');
        v_('-2N**2') = -2.0*v_('N**2');
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N/5') = fix(v_('N')/v_('5'));
        v_('NA') = v_('N/5');
        v_('A') = v_('NA');
        v_('A') = v_('A')/v_('RN');
        v_('NA+1') = 1+v_('NA');
        v_('NB') = v_('N/2');
        v_('B') = v_('NB');
        v_('B') = v_('B')/v_('RN');
        v_('NB+1') = 1+v_('NB');
        v_('NC') = v_('N/2');
        v_('C') = v_('NC');
        v_('C') = v_('C')/v_('RN');
        v_('NC+1') = 1+v_('NC');
        v_('ND') = v_('N/5')*v_('4');
        v_('D') = v_('ND');
        v_('D') = v_('D')/v_('RN');
        v_('ND+1') = 1+v_('ND');
        v_('INT') = v_('ONE');
        v_('INT') = v_('INT')+v_('A');
        v_('INT') = v_('INT')+v_('B');
        v_('INT') = v_('INT')-v_('C');
        v_('INT') = v_('INT')-v_('D');
        v_('INT') = v_('INT')*v_('RN');
        for I=v_('0'):v_('NA')
            v_(['V',int2str(I)]) = 1.0;
        end
        v_('STEP') = v_('B')-v_('A');
        v_('STEP') = v_('STEP')*v_('RN');
        v_('STEP') = v_('TWO')/v_('STEP');
        for I=v_('NA+1'):v_('NB')
            v_('J') = I-v_('NA');
            v_('RJ') = v_('J');
            v_('VAL') = v_('RJ')*v_('STEP');
            v_('VAL') = v_('ONE')-v_('VAL');
            v_(['V',int2str(I)]) = v_('VAL');
        end
        for I=v_('NB+1'):v_('NC')
            v_(['V',int2str(I)]) = -1.0;
        end
        v_('STEP') = v_('D')-v_('C');
        v_('STEP') = v_('STEP')*v_('RN');
        v_('STEP') = v_('TWO')/v_('STEP');
        for I=v_('NC+1'):v_('ND')
            v_('J') = I-v_('NC');
            v_('RJ') = v_('J');
            v_('VAL') = v_('RJ')*v_('STEP');
            v_('VAL') = v_('-ONE')+v_('VAL');
            v_(['V',int2str(I)]) = v_('VAL');
        end
        for I=v_('ND'):v_('N')
            v_(['V',int2str(I)]) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I)],ix_);
            pb.xnames{iv} = ['W',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('0'):v_('N')
            v_('VAL') = v_(['V',int2str(I)])*v_('-1/N**2');
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = v_('VAL');
            v_('VAL') = v_(['V',int2str(I)])*v_('-2/N**2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = v_('VAL');
        end
        v_('VAL') =...
              v_(['V',int2str(round(v_('1')))])-v_(['V',int2str(round(v_('0')))]);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('0')))]);
        valA(end+1) = v_('VAL');
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('VAL') = -2.0*v_(['V',int2str(I)]);
            v_('VAL') = v_('VAL')+v_(['V',int2str(round(v_('I-1')))]);
            v_('VAL') = v_('VAL')+v_(['V',int2str(round(v_('I+1')))]);
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = v_('VAL');
        end
        v_('VAL') =...
              v_(['V',int2str(round(v_('N-1')))])-v_(['V',int2str(round(v_('N')))]);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('N')))]);
        valA(end+1) = v_('VAL');
        [ig,ig_] = s2mpjlib('ii','INT',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INT';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('0')))]);
        valA(end+1) = 0.5;
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('0')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('0')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('0')))]);
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('1')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('0')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('0')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['W',int2str(round(v_('0')))]);
        valA(end+1) = v_('-1/N**2');
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['CON',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CON',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = 2.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(round(v_('I+1')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(round(v_('I-1')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = v_('-1/N**2');
            [ig,ig_] = s2mpjlib('ii','INT',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'INT';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('N')))]);
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('N-1')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('N')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['W',int2str(round(v_('N')))]);
        valA(end+1) = v_('-1/N**2');
        [ig,ig_] = s2mpjlib('ii','INT',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'INT';
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['U',int2str(round(v_('N')))]);
        valA(end+1) = 0.5;
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
        pbm.gconst(ig_('INT')) = v_('INT');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('0'):v_('N')
            pb.xlower(ix_(['U',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['U',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        for I=v_('0'):v_('N')
            pb.x0(ix_(['U',int2str(I)]),1) = v_(['V',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('N')
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            ename = ['D',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['O',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['D',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['D',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['U',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['O',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-TWO');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('ONE');
        for I=v_('1'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['O',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-TWO');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('TWO');
        end
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['D',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('ONE');
        for I=v_('0'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/N**2');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION            -2.38816D+02   $ N = 10 
% XL SOLUTION            -2.65340D+03   $ N = 100
% XL SOLUTION            -2.67211D+04   $ N = 1000
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-MN-V-V';
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

