function varargout = HS86(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS86
%    *********
% 
%    Source: problem 86 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'C-COLR2-AN-5-10'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS86';

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
        v_('N') = 5;
        v_('M') = 10;
        v_('1') = 1;
        v_('E1') = -15.0;
        v_('E2') = -27.0;
        v_('E3') = -36.0;
        v_('E4') = -18.0;
        v_('E5') = -12.0;
        v_('C1,1') = 30.0;
        v_('C2,1') = -20.0;
        v_('C3,1') = -10.0;
        v_('C4,1') = 32.0;
        v_('C5,1') = -10.0;
        v_('C1,2') = -20.0;
        v_('C2,2') = 39.0;
        v_('C3,2') = -6.0;
        v_('C4,2') = -31.0;
        v_('C5,2') = 32.0;
        v_('C1,3') = -10.0;
        v_('C2,3') = -6.0;
        v_('C3,3') = 10.0;
        v_('C4,3') = -6.0;
        v_('C5,3') = -10.0;
        v_('C1,4') = 32.0;
        v_('C2,4') = -31.0;
        v_('C3,4') = -6.0;
        v_('C4,4') = 39.0;
        v_('C5,4') = -20.0;
        v_('C1,5') = -10.0;
        v_('C2,5') = 32.0;
        v_('C3,5') = -10.0;
        v_('C4,5') = -20.0;
        v_('C5,5') = 30.0;
        v_('D1') = 4.0;
        v_('D2') = 8.0;
        v_('D3') = 10.0;
        v_('D4') = 6.0;
        v_('D5') = 2.0;
        v_('A1,1') = -16.0;
        v_('A2,1') = 0.0;
        v_('A3,1') = -3.5;
        v_('A4,1') = 0.0;
        v_('A5,1') = 0.0;
        v_('A6,1') = 2.0;
        v_('A7,1') = -1.0;
        v_('A8,1') = -1.0;
        v_('A9,1') = 1.0;
        v_('A10,1') = 1.0;
        v_('A1,2') = 2.0;
        v_('A2,2') = -2.0;
        v_('A3,2') = 0.0;
        v_('A4,2') = -2.0;
        v_('A5,2') = -9.0;
        v_('A6,2') = 0.0;
        v_('A7,2') = -1.0;
        v_('A8,2') = -2.0;
        v_('A9,2') = 2.0;
        v_('A10,2') = 1.0;
        v_('A1,3') = 0.0;
        v_('A2,3') = 0.0;
        v_('A3,3') = 2.0;
        v_('A4,3') = 0.0;
        v_('A5,3') = -2.0;
        v_('A6,3') = -4.0;
        v_('A7,3') = -1.0;
        v_('A8,3') = -3.0;
        v_('A9,3') = 3.0;
        v_('A10,3') = 1.0;
        v_('A1,4') = 1.0;
        v_('A2,4') = 4.0;
        v_('A3,4') = 0.0;
        v_('A4,4') = -4.0;
        v_('A5,4') = 1.0;
        v_('A6,4') = 0.0;
        v_('A7,4') = -1.0;
        v_('A8,4') = -2.0;
        v_('A9,4') = 4.0;
        v_('A10,4') = 1.0;
        v_('A1,5') = 0.0;
        v_('A2,5') = 2.0;
        v_('A3,5') = 0.0;
        v_('A4,5') = -1.0;
        v_('A5,5') = -2.8;
        v_('A6,5') = 0.0;
        v_('A7,5') = -1.0;
        v_('A8,5') = -1.0;
        v_('A9,5') = 5.0;
        v_('A10,5') = 1.0;
        v_('B1') = -40.0;
        v_('B2') = -2.0;
        v_('B3') = -0.25;
        v_('B4') = -4.0;
        v_('B5') = -4.0;
        v_('B6') = -1.0;
        v_('B7') = -40.0;
        v_('B8') = -60.0;
        v_('B9') = 5.0;
        v_('B10') = 1.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(J)]);
            valA(end+1) = v_(['E',int2str(J)]);
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_(['A',int2str(I),',',int2str(J)]);
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
            pbm.gconst(ig_(['C',int2str(I)])) = v_(['B',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        pb.y0 = 0.0*ones(pb.m,1);
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCUBE',iet_);
        elftv{it}{1} = 'XJ';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'XI';
        elftv{it}{2} = 'XJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            ename = ['D',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCUBE';
            ielftype(ie) = iet_('eCUBE');
            vname = ['X',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('XJ',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            for I=v_('1'):v_('N')
                ename = ['C',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('XI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('XJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['D',int2str(J)]);
            for I=v_('1'):v_('N')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_(['C',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -32.34867897
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AN-5-10';
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

    case 'eCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 3.0e+0*EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0e+0*EV_(1);
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
                H_(1,2) = 1.0e+0;
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

