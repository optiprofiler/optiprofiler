function varargout = ELEC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ELEC
%    *********
% 
%    Given np electrons, find the equilibrium state distribution of minimal
%    Columb potential of the electrons positioned on a conducting sphere
% 
%    This is problem 2 in the COPS (Version 2) collection of 
%    E. Dolan and J. More'
%    see "Benchmarking Optimization Software with COPS"
%    Argonne National Labs Technical Report ANL/MCS-246 (2000)
% 
%    SIF input: Nick Gould, November 2000
% 
%    classification = 'C-COOR2-AN-V-V'
% 
%    The number of electrons
% 
%       Alternative values for the SIF file parameters:
% IE NP                  25             $-PARAMETER
% IE NP                  50             $-PARAMETER
% IE NP                  100            $-PARAMETER
% IE NP                  200            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ELEC';

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
            v_('NP') = 25;  %  SIF file default value
        else
            v_('NP') = varargin{1};
        end
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('2PI') = 2.0*v_('PI');
        v_('0') = 0;
        v_('1') = 1;
        v_('NP-1') = -1+v_('NP');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('NP')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','PE',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('NP')
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['B',int2str(I)];
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
        for I=v_('1'):v_('NP')
            pbm.gconst(ig_(['B',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('N') = v_('NP');
        v_('1/N') = 1/v_('N');
        for I=v_('1'):v_('NP')
            v_('I') = I;
            v_('U') = I/v_('N');
            v_('THETA') = v_('2PI')*v_('U');
            v_(['THETA',int2str(I)]) = v_('THETA');
            v_('U') = v_('U')-v_('1/N');
            v_('PHI') = v_('PI')*v_('U');
            v_(['PHI',int2str(I)]) = v_('PHI');
        end
        for I=v_('1'):v_('NP')
            v_('COSTHETA') = cos(v_(['THETA',int2str(I)]));
            v_('SINTHETA') = sin(v_(['THETA',int2str(I)]));
            v_('SINPHI') = sin(v_(['PHI',int2str(I)]));
            v_('COSPHI') = cos(v_(['PHI',int2str(I)]));
            v_('X') = v_('COSTHETA')*v_('SINPHI');
            v_('Y') = v_('SINTHETA')*v_('SINPHI');
            v_('Z') = v_('COSPHI');
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_('X');
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_('X');
            end
            if(isKey(ix_,['Y',int2str(I)]))
                pb.x0(ix_(['Y',int2str(I)]),1) = v_('Y');
            else
                pb.y0(find(pbm.congrps==ig_(['Y',int2str(I)])),1) = v_('Y');
            end
            if(isKey(ix_,['Z',int2str(I)]))
                pb.x0(ix_(['Z',int2str(I)]),1) = v_('Z');
            else
                pb.y0(find(pbm.congrps==ig_(['Z',int2str(I)])),1) = v_('Z');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePE',iet_);
        elftv{it}{1} = 'XI';
        elftv{it}{2} = 'XJ';
        elftv{it}{3} = 'YI';
        elftv{it}{4} = 'YJ';
        elftv{it}{5} = 'ZI';
        elftv{it}{6} = 'ZJ';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NP-1')
            v_('I+1') = 1+I;
            for J=v_('I+1'):v_('NP')
                ename = ['P',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePE';
                ielftype(ie) = iet_('ePE');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('YJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ZI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Z',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ZJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('1'):v_('NP')
            ename = ['X',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Z',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NP-1')
            v_('I+1') = 1+I;
            for J=v_('I+1'):v_('NP')
                ig = ig_('PE');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for I=v_('1'):v_('NP')
            ig = ig_(['B',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['X',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Y',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Z',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLUTION             2.43812D+02   $ (NH=25)
% LO SOLUTION             1.05518D+03   $ (NH=50)
% LO SOLUTION             4.44836D+03   $ (NH=100)
% LO SOLUTION             1.84389D+04   $ (NH=200)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-V-V';
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

    case 'ePE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,6);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)-1;
        U_(3,5) = U_(3,5)+1;
        U_(3,6) = U_(3,6)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        SS = IV_(1)^2+IV_(2)^2+IV_(3)^2;
        ROOTSS = sqrt(SS);
        varargout{1} = 1.0/ROOTSS;
        if(nargout>1)
            g_(1,1) = -IV_(1)/ROOTSS^3;
            g_(2,1) = -IV_(2)/ROOTSS^3;
            g_(3,1) = -IV_(3)/ROOTSS^3;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 3.0*IV_(1)^2/ROOTSS^5-1.0/ROOTSS^3;
                H_(1,2) = 3.0*IV_(1)*IV_(2)/ROOTSS^5;
                H_(2,1) = H_(1,2);
                H_(1,3) = 3.0*IV_(1)*IV_(3)/ROOTSS^5;
                H_(3,1) = H_(1,3);
                H_(2,2) = 3.0*IV_(2)^2/ROOTSS^5-1.0/ROOTSS^3;
                H_(2,3) = 3.0*IV_(2)*IV_(3)/ROOTSS^5;
                H_(3,2) = H_(2,3);
                H_(3,3) = 3.0*IV_(3)^2/ROOTSS^5-1.0/ROOTSS^3;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eSQR'

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

