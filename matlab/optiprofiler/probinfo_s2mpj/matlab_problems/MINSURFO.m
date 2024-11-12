function varargout = MINSURFO(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MINSURFO
%    *********
% 
%    Find the surface with minimal area, given boundary conditions, 
%    and above an obstacle.
% 
%    This is problem 17 in the COPS (Version 2) collection of 
%    E. Dolan and J. More'
%    see "Benchmarking Optimization Software with COPS"
%    Argonne National Labs Technical Report ANL/MCS-246 (2000)
% 
%    SIF input: Nick Gould, December 2000
% 
%    classification = 'C-COBR2-AN-V-V'
% 
%  grid points in x direction (fixed at 50 in COPS)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MINSURFO';

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
        v_('NX') = 5;
        v_('NY') = 10;
        v_('0') = 0;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('NX+1') = 1+v_('NX');
        v_('NY+1') = 1+v_('NY');
        v_('RNX+1') = v_('NX+1');
        v_('RNY+1') = v_('NY+1');
        v_('HX') = 1.0/v_('RNX+1');
        v_('HY') = 1.0/v_('RNY+1');
        v_('AREA') = v_('HX')*v_('HY');
        v_('AREA') = 0.5*v_('AREA');
        v_('1/AREA') = 1.0/v_('AREA');
        v_('1/HX') = 1.0/v_('HX');
        v_('1/HX2') = v_('1/HX')*v_('1/HX');
        v_('1/HY') = 1.0/v_('HY');
        v_('1/HY2') = v_('1/HY')*v_('1/HY');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('NX+1')
            for J=v_('0'):v_('NY+1')
                [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('0'):v_('NX')
            for J=v_('0'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('1/AREA');
            end
        end
        for I=v_('1'):v_('NX+1')
            for J=v_('1'):v_('NY+1')
                [ig,ig_] = s2mpjlib('ii',['B',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('1/AREA');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('0'):v_('NX')
            for J=v_('0'):v_('NY')
                pbm.gconst(ig_(['A',int2str(I),',',int2str(J)])) = -1.0;
            end
        end
        for I=v_('1'):v_('NX+1')
            for J=v_('1'):v_('NY+1')
                pbm.gconst(ig_(['B',int2str(I),',',int2str(J)])) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        v_('1/4HX') = 0.25/v_('HX');
        v_('3/4HX') = 0.75/v_('HX');
        v_('1/4HY') = 0.25/v_('HY');
        v_('3/4HY') = 0.75/v_('HY');
        v_('3/4HX') = 0.9999999999+v_('3/4HX');
        v_('3/4HY') = 0.9999999999+v_('3/4HY');
        v_('1/4HX') = fix(v_('1/4HX'));
        v_('1/4HY') = fix(v_('1/4HY'));
        v_('3/4HX') = fix(v_('3/4HX'));
        v_('3/4HY') = fix(v_('3/4HY'));
        for I=v_('1/4HX'):v_('3/4HX')
            for J=v_('1/4HY'):v_('3/4HY')
                pb.xlower(ix_(['V',int2str(I),',',int2str(J)]),1) = 1.0;
            end
        end
        for J=v_('0'):v_('NY+1')
            pb.xlower(ix_(['V',int2str(round(v_('0'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['V',int2str(round(v_('0'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['V',int2str(round(v_('NX+1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['V',int2str(round(v_('NX+1'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('0'):v_('NX+1')
            v_('I') = I;
            v_('VIJ') = 2.0*I;
            v_('VIJ') = v_('VIJ')*v_('HX');
            v_('VIJ') = -1.0+v_('VIJ');
            v_('VIJ') = v_('VIJ')*v_('VIJ');
            v_('VIJ') = v_('ONE')-v_('VIJ');
            pb.xlower(ix_(['V',int2str(I),',',int2str(round(v_('0')))]),1) = v_('VIJ');
            pb.xupper(ix_(['V',int2str(I),',',int2str(round(v_('0')))]),1) = v_('VIJ');
            pb.xlower(ix_(['V',int2str(I),',',int2str(round(v_('NY+1')))]),1) =...
                  v_('VIJ');
            pb.xupper(ix_(['V',int2str(I),',',int2str(round(v_('NY+1')))]),1) =...
                  v_('VIJ');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('0'):v_('NX+1')
            v_('I') = I;
            v_('VIJ') = 2.0*I;
            v_('VIJ') = v_('VIJ')*v_('HX');
            v_('VIJ') = -1.0+v_('VIJ');
            v_('VIJ') = v_('VIJ')*v_('VIJ');
            v_('VIJ') = v_('ONE')-v_('VIJ');
            for J=v_('0'):v_('NY+1')
                pb.x0(ix_(['V',int2str(I),',',int2str(J)]),1) = v_('VIJ');
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
        for I=v_('0'):v_('NX')
            v_('I+1') = 1+I;
            for J=v_('0'):v_('NY')
                v_('J+1') = 1+J;
                ename = ['I',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['V',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['J',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['V',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for J=v_('0'):v_('NY+1')
            v_('J1') = 1+J;
            ename = ['J',int2str(round(v_('NX+1'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eISQ';
            ielftype(ie) = iet_('eISQ');
            ename = ['J',int2str(round(v_('NX+1'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(round(v_('NX+1'))),',',int2str(round(v_('J1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['J',int2str(round(v_('NX+1'))),',',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(round(v_('NX+1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('0'):v_('NX+1')
            v_('I1') = 1+I;
            ename = ['I',int2str(I),',',int2str(round(v_('NY+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eISQ';
            ielftype(ie) = iet_('eISQ');
            ename = ['I',int2str(I),',',int2str(round(v_('NY+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(round(v_('I1'))),',',int2str(round(v_('NY+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['I',int2str(I),',',int2str(round(v_('NY+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['V',int2str(I),',',int2str(round(v_('NY+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gROOT',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('0'):v_('NX')
            for J=v_('0'):v_('NY')
                ig = ig_(['A',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gROOT';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['I',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('1/HX2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['J',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('1/HY2');
            end
        end
        for I=v_('1'):v_('NX+1')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('NY+1')
                v_('J-1') = -1+J;
                ig = ig_(['B',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gROOT';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['I',int2str(round(v_('I-1'))),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('1/HX2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['J',int2str(I),',',int2str(round(v_('J-1')))]);
                pbm.grelw{ig}(posel) = v_('1/HY2');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLUTION            2.51948D+00    $ (NX=50,NY=25)
% LO SOLUTION            2.51488D+00    $ (NX=50,NY=50)
% LO SOLUTION            2.50568D+00    $ (NX=50,NY=75)
% LO SOLUTION            2.50694D+00    $ (NX=50,NY=100)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COBR2-AN-V-V';
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

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gROOT'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        ROOTAL = sqrt(GVAR_);
        varargout{1} = ROOTAL;
        if(nargout>1)
            g_ = 0.5/ROOTAL;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -0.25/ROOTAL^3;
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

