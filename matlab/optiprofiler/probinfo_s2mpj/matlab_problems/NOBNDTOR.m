function varargout = NOBNDTOR(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NOBNDTOR
%    *********
% 
%    The quadratic elastic torsion problem
% 
%    The problem comes from the obstacle problem on a square.
% 
%    The square is discretized into (px-1)(py-1) little squares. The
%    heights of the considered surface above the corners of these little
%    squares are the problem variables,  There are px**2 of them.
% 
%    The dimension of the problem is specified by Q, which is half the
%    number discretization points along one of the coordinate
%    direction.
%    Since the number of variables is P**2, it is given by 4Q**2
% 
%    Source: problem 1 (c=5, starting point U = upper bound) in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    A variant of TORSION1 in which some of the variables are
%    unconstrained.
% 
%    classification = 'C-CQBR2-AY-V-0'
% 
%    Q is half the number of discretized points along the X axis
% 
%       Alternative values for the SIF file parameters:
% IE Q                   2              $-PARAMETER  n= 16
% IE Q                   5              $-PARAMETER  n= 100   original value
% IE Q                   11             $-PARAMETER  n= 484
% IE Q                   16             $-PARAMETER  n= 1024
% IE Q                   37             $-PARAMETER  n= 5476
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NOBNDTOR';

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
            v_('Q') = 3;  %  SIF file default value
        else
            v_('Q') = varargin{1};
        end
% IE Q                   50             $-PARAMETER  n= 10000
% IE Q                   61             $-PARAMETER  n= 14884
% IE Q                   100            $-PARAMETER  n= 40000
        v_('C') = 5.0;
        v_('Q+1') = 1+v_('Q');
        v_('P') = v_('Q')+v_('Q');
        v_('P-1') = -1+v_('P');
        v_('1/H') = v_('P-1');
        v_('H') = 1.0/v_('1/H');
        v_('H2') = v_('H')*v_('H');
        v_('C0') = v_('H2')*v_('C');
        v_('LC') = -1.0*v_('C0');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for J=v_('1'):v_('P')
            for I=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('LC');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for J=v_('1'):v_('P')
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('P-1')
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        for I=v_('2'):v_('Q')
            for J=v_('2'):I
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('UPPL') = v_('RJ-1')*v_('H');
                v_('UPPL') = 1.0e+21;
                v_('LOWL') = -1.0*v_('UPPL');
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWL');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPL');
            end
            v_('MI') = -1*I;
            v_('P-I') = v_('P')+v_('MI');
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('UPPM') = v_('RI-1')*v_('H');
            v_('UPPM') = 1.0e+21;
            v_('LOWM') = -1.0*v_('UPPM');
            v_('P-I+1') = 1+v_('P-I');
            for J=I:v_('P-I+1')
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWM');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPM');
            end
            for J=v_('P-I+1'):v_('P-1')
                v_('MJ') = -1*J;
                v_('P-J') = v_('P')+v_('MJ');
                v_('RP-J') = v_('P-J');
                v_('UPPR') = v_('RP-J')*v_('H');
                v_('UPPR') = 1.0e+21;
                v_('LOWR') = -1.0*v_('UPPR');
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWR');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPR');
            end
        end
        for I=v_('Q+1'):v_('P-1')
            v_('MI') = -1*I;
            v_('P-I') = v_('P')+v_('MI');
            v_('P-I+1') = 1+v_('P-I');
            for J=v_('2'):v_('P-I+1')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('UPPL') = v_('RJ-1')*v_('H');
                v_('LOWL') = -1.0*v_('UPPL');
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWL');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPL');
            end
            v_('RP-I') = v_('P-I');
            v_('UPPM') = v_('RP-I')*v_('H');
            v_('LOWM') = -1.0*v_('UPPM');
            for J=v_('P-I+1'):I
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWM');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPM');
            end
            for J=I:v_('P-1')
                v_('MJ') = -1*J;
                v_('P-J') = v_('P')+v_('MJ');
                v_('RP-J') = v_('P-J');
                v_('UPPR') = v_('RP-J')*v_('H');
                v_('LOWR') = -1.0*v_('UPPR');
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('LOWR');
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = v_('UPPR');
            end
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for J=v_('1'):v_('P')
            pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.x0(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('P-1')
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        for I=v_('2'):v_('Q')
            for J=v_('2'):I
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('UPPL') = v_('RJ-1')*v_('H');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPL');
            end
            v_('MI') = -1*I;
            v_('P-I') = v_('P')+v_('MI');
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('UPPM') = v_('RI-1')*v_('H');
            v_('P-I+1') = 1+v_('P-I');
            for J=I:v_('P-I+1')
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPM');
            end
            for J=v_('P-I+1'):v_('P-1')
                v_('MJ') = -1*J;
                v_('P-J') = v_('P')+v_('MJ');
                v_('RP-J') = v_('P-J');
                v_('UPPR') = v_('RP-J')*v_('H');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPR');
            end
        end
        for I=v_('Q+1'):v_('P-1')
            v_('MI') = -1*I;
            v_('P-I') = v_('P')+v_('MI');
            v_('P-I+1') = 1+v_('P-I');
            for J=v_('2'):v_('P-I+1')
                v_('J-1') = -1+J;
                v_('RJ-1') = v_('J-1');
                v_('UPPL') = v_('RJ-1')*v_('H');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPL');
            end
            v_('RP-I') = v_('P-I');
            v_('UPPM') = v_('RP-I')*v_('H');
            for J=v_('P-I+1'):I
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPM');
            end
            for J=I:v_('P-1')
                v_('MJ') = -1*J;
                v_('P-J') = v_('P')+v_('MJ');
                v_('RP-J') = v_('P-J');
                v_('UPPR') = v_('RP-J')*v_('H');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('UPPR');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('P-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('2'):v_('P-1')
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['C',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['D',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(I),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                ig = ig_(['G',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = .25;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 0.25;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = .25;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['D',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 0.25;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(2)            -5.1851852D-1
% LO SOLTN(5)            -4.9234185D-1
% LO SOLTN(11)           -4.5608771D-1
% LO SOLTN(16)           ???
% LO SOLTN(37)           ???
% LO SOLTN(50)           ???
% LO SOLTN(61)           ???
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CQBR2-AY-V-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

