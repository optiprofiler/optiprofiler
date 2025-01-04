function varargout = JNLBRNG2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : JNLBRNG2
%    *********
% 
%    The quadratic journal bearing problem (with excentricity = 0.5)
%    This is a variant of the problem stated in the report quoted below.
%    It corresponds to the problem as distributed in MINPACK-2.
% 
%    Source:
%    J. More' and G. Toraldo,
%    "On the Solution of Large Quadratic-Programming Problems with Bound
%    Constraints", 
%    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
% 
%    SIF input: Ph. Toint, Dec 1989.
%    modified by Peihuang Chen, according to MINPACK-2, Apr 1992
% 
%    classification = 'C-CQBR2-AY-V-0'
% 
%    The rectangle is discretized into (pt-1)(py-1) little rectangles. The
%    heights of the considered surface above the corners of these little
%    rectangles are the problem variables,  There are px*py of them.
% 
%    PT is the number of points along the T (\theta) side of the rectangle
%    PY is the number of points along the Y side of the rectangle
% 
%       Alternative values for the SIF file parameters:
% IE PT                  4              $-PARAMETER  n=16
% IE PY                  4              $-PARAMETER
% 
% IE PT                  10             $-PARAMETER  n=100
% IE PY                  10             $-PARAMETER
% 
% IE PT                  23             $-PARAMETER  n=529
% IE PY                  23             $-PARAMETER
% 
% IE PT                  32             $-PARAMETER  n=1024
% IE PY                  32             $-PARAMETER
% 
% IE PT                  34             $-PARAMETER  n=1156
% IE PY                  34             $-PARAMETER
% 
% IE PT                  75             $-PARAMETER  n=5625   original value
% IE PY                  75             $-PARAMETER           original value
% 
% IE PT                  100            $-PARAMETER  n=10000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'JNLBRNG2';

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
            v_('PT') = 5;  %  SIF file default value
        else
            v_('PT') = varargin{1};
        end
% IE PY                  100            $-PARAMETER
        if(nargs<2)
            v_('PY') = 5;  %  SIF file default value
        else
            v_('PY') = varargin{2};
        end
% IE PT                  125            $-PARAMETER  n=15625
% IE PY                  125            $-PARAMETER
        if(nargs<3)
            v_('EX') = 0.5;  %  SIF file default value
        else
            v_('EX') = varargin{3};
        end
        v_('PI/4') = atan(1.0);
        v_('LT') = 8.0*v_('PI/4');
        v_('LY') = 20.0;
        v_('SIX') = 6.0;
        v_('PT-1') = -1+v_('PT');
        v_('RPT-1') = v_('PT-1');
        v_('HT1') = 1.0/v_('RPT-1');
        v_('HT') = v_('HT1')*v_('LT');
        v_('1/HT') = 1.0/v_('HT');
        v_('PY-1') = -1+v_('PY');
        v_('RPY-1') = v_('PY-1');
        v_('HY1') = 1.0/v_('RPY-1');
        v_('HY') = v_('HY1')*v_('LY');
        v_('1/HY') = 1.0/v_('HY');
        v_('HTHY') = v_('HT')*v_('HY');
        v_('HT/HY') = v_('HT')*v_('1/HY');
        v_('HY/HT') = v_('HY')*v_('1/HT');
        v_('EXHTHY') = v_('HTHY')*v_('EX');
        v_('CLINC') = -1.0*v_('EXHTHY');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('PT')
            for J=v_('1'):v_('PY')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('2'):v_('PT-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HT');
            v_('SXI1') = sin(v_('XI1'));
            v_('COEFF') = v_('SXI1')*v_('CLINC');
            for J=v_('2'):v_('PY-1')
                [ig,ig_] = s2mpjlib('ii','G',ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = v_('COEFF');
            end
        end
        for I=v_('1'):v_('PT-1')
            for J=v_('1'):v_('PY-1')
                [ig,ig_] = s2mpjlib('ii',['GR',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = 2.0;
            end
        end
        for I=v_('2'):v_('PT')
            for J=v_('2'):v_('PY')
                [ig,ig_] = s2mpjlib('ii',['GL',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = 2.0;
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
        for J=v_('1'):v_('PY')
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(round(v_('PT'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('PT'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('PT-1')
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('PY')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('PY')))]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('2'):v_('PT-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HT');
            v_('SXI1') = sin(v_('XI1'));
            for J=v_('2'):v_('PY-1')
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('SXI1');
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
        for I=v_('1'):v_('PT-1')
            v_('I+1') = 1+I;
            for J=v_('1'):v_('PY-1')
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
            end
        end
        for I=v_('2'):v_('PT')
            v_('I-1') = -1+I;
            for J=v_('2'):v_('PY')
                v_('J-1') = -1+J;
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
        for I=v_('1'):v_('PT-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HT');
            v_('CXI1') = cos(v_('XI1'));
            v_('ECX') = v_('CXI1')*v_('EX');
            v_('ECX1') = 1.0+v_('ECX');
            v_('E12') = v_('ECX1')*v_('ECX1');
            v_('WI') = v_('ECX1')*v_('E12');
            v_('2WI') = v_('WI')+v_('WI');
            v_('XI+1') = v_('XI1')+v_('HT');
            v_('CXI+1') = cos(v_('XI+1'));
            v_('E+CX0') = v_('CXI+1')*v_('EX');
            v_('E+CX1') = 1.0+v_('E+CX0');
            v_('E22') = v_('E+CX1')*v_('E+CX1');
            v_('WI+1') = v_('E+CX1')*v_('E22');
            v_('PM0') = v_('2WI')+v_('WI+1');
            v_('PM1') = v_('PM0')/v_('SIX');
            v_('LA/HY2') = v_('PM1')*v_('HT/HY');
            v_('LA/HT2') = v_('PM1')*v_('HY/HT');
            for J=v_('1'):v_('PY-1')
                ig = ig_(['GR',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('LA/HT2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('LA/HY2');
            end
        end
        for I=v_('2'):v_('PT')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('XI1') = v_('RI-1')*v_('HT');
            v_('CXI1') = cos(v_('XI1'));
            v_('ECX') = v_('CXI1')*v_('EX');
            v_('ECX1') = 1.0+v_('ECX');
            v_('E12') = v_('ECX1')*v_('ECX1');
            v_('WI') = v_('ECX1')*v_('E12');
            v_('2WI') = v_('WI')+v_('WI');
            v_('XI-1') = v_('XI1')-v_('HT');
            v_('CXI-1') = cos(v_('XI-1'));
            v_('E-CX0') = v_('CXI-1')*v_('EX');
            v_('E-CX1') = 1.0+v_('E-CX0');
            v_('E32') = v_('E-CX1')*v_('E-CX1');
            v_('WI-1') = v_('E-CX1')*v_('E32');
            v_('PL0') = v_('2WI')+v_('WI-1');
            v_('PL1') = v_('PL0')/v_('SIX');
            v_('MU/HY2') = v_('PL1')*v_('HT/HY');
            v_('MU/HT2') = v_('PL1')*v_('HY/HT');
            for J=v_('2'):v_('PY')
                ig = ig_(['GL',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('MU/HT2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['D',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('MU/HY2');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(4)            -0.4764000 
% LO SOLTN(10)           -0.3952800
% LO SOLTN(23)           -0.4102400
% LO SOLTN(32)           -0.4124900
% LO SOLTN(75)           -0.4146600
% LO SOLTN(100)          -0.4148700
% LO SOLTN(125)          -0.4149600
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

