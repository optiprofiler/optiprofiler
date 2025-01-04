function varargout = TRAINH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TRAINH
%    *********
% 
%    The problem is to minimize the energy spent to move a train 
%    from the beginning of a track to its end in a given time.  The train
%    is slowed down by some drag (assumed to be quadratic in the the velocity).
%    The track follows the slope of a hill.  The track geometry is given
%    by the equation
% 
%           1                    1  ns-1                          x - z_i
%    g(x) = - ( a_1 + a_{ns} ) + -- SUM ( s_{i+1} - s_i ) arctan( ------- )
%           2                    pi  1                              eps
% 
%    where the z_i are the breakpoints between sections of the track, and where 
%    the s_i are the "slopes" on these sections (eps is a regularization
%    parameter). Here we have a track of the overall shape
% 
%                      ______
%                     /      \      z0 = 0, z1 = 2, z2 = 4, z3 = 6
%                    /        \     s1 = 2, s2 = 0, s3 = -2
%                   /          \    eps = 0.05
% 
%    The control variables are the acceleration force (UA) and the braking
%    force (UB) applied on the train.
% 
%    Source: adapted from
%    J. Kautsky and N. K. Nichols,
%    "OTEP-2: Optimal Train Energy Programme, mark 2",
%    Numerical Analysis Report NA/4/83,
%    Department of Mathematics, University of Reading, 1983.
% 
%    SIF input: N. Nichols and Ph. Toint, April 1993
% 
%    classification = 'C-CQOR2-MN-V-V'
% 
%    Number of discretized points in the interval
% 
%       Alternative values for the SIF file parameters:
% IE N                   11             $-PARAMETER n=48, m=22
% IE N                   51             $-PARAMETER n=208, m=102
% IE N                   101            $-PARAMETER n=408, m=202  original value
% IE N                   201            $-PARAMETER n=808, m=402
% IE N                   501            $-PARAMETER n=2008, m=1002 
% IE N                   1001           $-PARAMETER n=4008, m=2002
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TRAINH';

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
            v_('N') = 11;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5001           $-PARAMETER n=20008, m=10002
        if(nargs<2)
            v_('TIME') = 4.8;  %  SIF file default value
        else
            v_('TIME') = varargin{2};
        end
        if(nargs<3)
            v_('LENGTH') = 6.0;  %  SIF file default value
        else
            v_('LENGTH') = varargin{3};
        end
        if(nargs<4)
            v_('NS') = 3;  %  SIF file default value
        else
            v_('NS') = varargin{4};
        end
        if(nargs<5)
            v_('Z1') = 2.0;  %  SIF file default value
        else
            v_('Z1') = varargin{5};
        end
        if(nargs<6)
            v_('Z2') = 4.0;  %  SIF file default value
        else
            v_('Z2') = varargin{6};
        end
        if(nargs<7)
            v_('S1') = 2.0;  %  SIF file default value
        else
            v_('S1') = varargin{7};
        end
        if(nargs<8)
            v_('S2') = 0.0;  %  SIF file default value
        else
            v_('S2') = varargin{8};
        end
        if(nargs<9)
            v_('S3') = -2.0;  %  SIF file default value
        else
            v_('S3') = varargin{9};
        end
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
        v_('EPS') = 0.05;
        v_('0') = 0;
        v_('1') = 1;
        v_('PI') = 3.1415926535;
        v_('NS-1') = -1+v_('NS');
        v_('BH/2') = v_('B')*v_('H/2');
        v_('1+BH/2') = 1.0+v_('BH/2');
        v_('BH/2-1') = -1.0+v_('BH/2');
        v_('-AH') = v_('A')*v_('-H');
        v_('LENGTH/N') = v_('LENGTH')/v_('RN');
        v_('CH/2') = v_('C')*v_('H/2');
        v_('SUMS') =...
              v_(['S',int2str(round(v_('1')))])+v_(['S',int2str(round(v_('NS')))]);
        v_('-AVS') = -0.5*v_('SUMS');
        v_('-AVSH') = v_('-AVS')*v_('H');
        v_('CNST') = v_('-AH')+v_('-AVSH');
        v_('H/2PI') = v_('H/2')/v_('PI');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
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
        [ig,ig_] = s2mpjlib('ii','ENERGY',ig_);
        gtype{ig} = '<>';
        for I=v_('0'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['XEQ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['XEQ',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('-H/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I)]);
            valA(end+1) = v_('-H/2');
            [ig,ig_] = s2mpjlib('ii',['VEQ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['VEQ',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('1+BH/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['V',int2str(I)]);
            valA(end+1) = v_('BH/2-1');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['UA',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('-H/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['UA',int2str(I)]);
            valA(end+1) = v_('-H/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['UB',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('-H/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['UB',int2str(I)]);
            valA(end+1) = v_('-H/2');
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
            pbm.gconst(ig_(['VEQ',int2str(I)])) = v_('CNST');
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
        [it,iet_] = s2mpjlib( 'ii', 'eATAN',iet_);
        elftv{it}{1} = 'XX';
        elftp{it}{1} = 'ZZ';
        elftp{it}{2} = 'E';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
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
            for J=v_('1'):v_('NS-1')
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eATAN';
                ielftype(ie) = iet_('eATAN');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XX',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('ZZ',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['Z',int2str(J)]);
                [~,posep] = ismember('E',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('EPS');
            end
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
            for J=v_('1'):v_('NS-1')
                v_('J+1') = 1+J;
                v_('DS') = v_(['S',int2str(round(v_('J+1')))])-v_(['S',int2str(J)]);
                v_('WJ') = v_('DS')*v_('H/2PI');
                ig = ig_(['VEQ',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('WJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('I+1'))),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('WJ');
            end
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
% LO SOLUTION(11)        12.423025536
% LO SOLUTION(51)        12.295777964
% LO SOLUTION(101)       12.306399739
% LO SOLUTION(201)       12.309848614
% LO SOLUTION(1001)      12.307801327
% LO SOLUTION(5001)      12.221148056
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MN-V-V';
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

    case 'eATAN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DX = EV_(1)-pbm.elpar{iel_}(1);
        E2 = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(2);
        varargout{1} = atan(DX/pbm.elpar{iel_}(2));
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(2)/(E2+DX*DX);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -2.0*DX*pbm.elpar{iel_}(2)/(E2+DX*DX)^2;
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

