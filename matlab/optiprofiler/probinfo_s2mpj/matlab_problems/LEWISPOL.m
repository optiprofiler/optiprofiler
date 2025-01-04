function varargout = LEWISPOL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Adrian Lewis Polynomial Problem,
%    The problem is a transformation of a number theory integer
%    programming problem.
% 
%    Source:
%    A. Lewis, private communication.
% 
%    SIF input: A.R. Conn and Ph. Toint, March 1990.
% 
%    classification = 'C-CQOR2-AN-6-9'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LEWISPOL';

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
        v_('N') = 6;
        v_('DEG') = 3;
        v_('PEN') = 1.0e4;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('DEG-1') = -1+v_('DEG');
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for J=v_('0'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['A',int2str(J)],ix_);
            pb.xnames{iv} = ['A',int2str(J)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for J=v_('0'):v_('N-1')
            v_(['C',int2str(round(v_('0'))),',',int2str(J)]) = 1.0;
            [ig,ig_] = s2mpjlib('ii','D0',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'D0';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['A',int2str(J)]);
            valA(end+1) = v_(['C',int2str(round(v_('0'))),',',int2str(J)]);
        end
        for I=v_('1'):v_('DEG-1')
            v_('I-1') = -1+I;
            for J=I:v_('N-1')
                v_('RJ') = J;
                v_(['C',int2str(I),',',int2str(J)]) =...
                      v_(['C',int2str(round(v_('I-1'))),',',int2str(J)])*v_('RJ');
                [ig,ig_] = s2mpjlib('ii',['D',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['D',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['A',int2str(J)]);
                valA(end+1) = v_(['C',int2str(I),',',int2str(J)]);
            end
        end
        for J=v_('0'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['INT',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['INT',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['A',int2str(J)]);
            valA(end+1) = -1.0;
            pbm.gscale(ig,1) = v_('PEN');
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
        v_(['CT',int2str(round(v_('0')))]) = -1.0;
        pbm.gconst(ig_('D0')) = v_(['CT',int2str(round(v_('0')))]);
        for I=v_('1'):v_('DEG-1')
            v_('I-1') = -1+I;
            v_('-I') = -1*I;
            v_('N+1-I') = v_('N+1')+v_('-I');
            v_('VAL') = v_('N+1-I');
            v_(['CT',int2str(I)]) = v_(['CT',int2str(round(v_('I-1')))])*v_('VAL');
            pbm.gconst(ig_(['D',int2str(I)])) = v_(['CT',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -10.0*ones(pb.n,1);
        pb.xupper = 10.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'A0'))
            pb.x0(ix_('A0'),1) = -1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A0')),1) = -1.0;
        end
        if(isKey(ix_,'A1'))
            pb.x0(ix_('A1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A1')),1) = 1.0;
        end
        if(isKey(ix_,'A2'))
            pb.x0(ix_('A2'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A2')),1) = 1.0;
        end
        if(isKey(ix_,'A3'))
            pb.x0(ix_('A3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('A3')),1) = 0.0;
        end
        if(isKey(ix_,'A4'))
            pb.x0(ix_('A4'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A4')),1) = 1.0;
        end
        if(isKey(ix_,'A5'))
            pb.x0(ix_('A5'),1) = -1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A5')),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eCB',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('0'):v_('N-1')
            ename = ['O',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['A',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCB';
            ielftype(ie) = iet_('eCB');
            vname = ['A',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('0'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['O',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['INT',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-AN-6-9';
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

    case 'eCB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3;
        if(nargout>1)
            g_(1,1) = 3.0*EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*EV_(1);
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

