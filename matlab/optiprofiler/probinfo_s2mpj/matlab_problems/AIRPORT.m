function varargout = AIRPORT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    This problem is concerned with the localisation of airports in Brazil.
%    We consider  m  balls in the real plane, whose centers are the coordinates
%    of some Brazilian  cities and whose  radius were chosen such that the balls are
%    disjoint. The problem is to find one point  (xi, yi) on  each ball, i=1,..,m,
%    such that  SUM(||(xi,yi) - (xj,yj)||)  is  minimum, where the sum involves all
%    the pairs (i,j) such that 1 <= i <= m, 1 <= j <= m and i <> j.
% 
%    For this problem instance, we have m =  42 cities and n = 84 points, 
%    i.e, 42 nonlinear inequalities constraints and 84 variables.
% 
%    Source:
%    Contribution from a LANCELOT user.
% 
%    SIF input : Rodrigo de Barros Nabholz & Maria Aparecida Diniz Ehrhardt
%                November 1994, DMA - IMECC- UNICAMP
%    Adaptation for CUTE: Ph. Toint, November 1994.
% 
%    classification = 'SQR2-MN-84-42'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'AIRPORT';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 42;
        v_('N-1') = 41;
        v_('1') = 1;
        v_('R1') = 0.09;
        v_('R2') = 0.3;
        v_('R3') = 0.09;
        v_('R4') = 0.45;
        v_('R5') = 0.5;
        v_('R6') = 0.04;
        v_('R7') = 0.1;
        v_('R8') = 0.02;
        v_('R9') = 0.02;
        v_('R10') = 0.07;
        v_('R11') = 0.4;
        v_('R12') = 0.045;
        v_('R13') = 0.05;
        v_('R14') = 0.056;
        v_('R15') = 0.36;
        v_('R16') = 0.08;
        v_('R17') = 0.07;
        v_('R18') = 0.36;
        v_('R19') = 0.67;
        v_('R20') = 0.38;
        v_('R21') = 0.37;
        v_('R22') = 0.05;
        v_('R23') = 0.4;
        v_('R24') = 0.66;
        v_('R25') = 0.05;
        v_('R26') = 0.07;
        v_('R27') = 0.08;
        v_('R28') = 0.3;
        v_('R29') = 0.31;
        v_('R30') = 0.49;
        v_('R31') = 0.09;
        v_('R32') = 0.46;
        v_('R33') = 0.12;
        v_('R34') = 0.07;
        v_('R35') = 0.07;
        v_('R36') = 0.09;
        v_('R37') = 0.05;
        v_('R38') = 0.13;
        v_('R39') = 0.16;
        v_('R40') = 0.46;
        v_('R41') = 0.25;
        v_('R42') = 0.1;
        v_('CX1') = -6.3;
        v_('CX2') = -7.8;
        v_('CX3') = -9.0;
        v_('CX4') = -7.2;
        v_('CX5') = -5.7;
        v_('CX6') = -1.9;
        v_('CX7') = -3.5;
        v_('CX8') = -0.5;
        v_('CX9') = 1.4;
        v_('CX10') = 4.0;
        v_('CX11') = 2.1;
        v_('CX12') = 5.5;
        v_('CX13') = 5.7;
        v_('CX14') = 5.7;
        v_('CX15') = 3.8;
        v_('CX16') = 5.3;
        v_('CX17') = 4.7;
        v_('CX18') = 3.3;
        v_('CX19') = 0.0;
        v_('CX20') = -1.0;
        v_('CX21') = -0.4;
        v_('CX22') = 4.2;
        v_('CX23') = 3.2;
        v_('CX24') = 1.7;
        v_('CX25') = 3.3;
        v_('CX26') = 2.0;
        v_('CX27') = 0.7;
        v_('CX28') = 0.1;
        v_('CX29') = -0.1;
        v_('CX30') = -3.5;
        v_('CX31') = -4.0;
        v_('CX32') = -2.7;
        v_('CX33') = -0.5;
        v_('CX34') = -2.9;
        v_('CX35') = -1.2;
        v_('CX36') = -0.4;
        v_('CX37') = -0.1;
        v_('CX38') = -1.0;
        v_('CX39') = -1.7;
        v_('CX40') = -2.1;
        v_('CX41') = -1.8;
        v_('CX42') = 0.0;
        v_('CY1') = 8.0;
        v_('CY2') = 5.1;
        v_('CY3') = 2.0;
        v_('CY4') = 2.6;
        v_('CY5') = 5.5;
        v_('CY6') = 7.1;
        v_('CY7') = 5.9;
        v_('CY8') = 6.6;
        v_('CY9') = 6.1;
        v_('CY10') = 5.6;
        v_('CY11') = 4.9;
        v_('CY12') = 4.7;
        v_('CY13') = 4.3;
        v_('CY14') = 3.6;
        v_('CY15') = 4.1;
        v_('CY16') = 3.0;
        v_('CY17') = 2.4;
        v_('CY18') = 3.0;
        v_('CY19') = 4.7;
        v_('CY20') = 3.4;
        v_('CY21') = 2.3;
        v_('CY22') = 1.5;
        v_('CY23') = 0.5;
        v_('CY24') = -1.7;
        v_('CY25') = -2.0;
        v_('CY26') = -3.1;
        v_('CY27') = -3.5;
        v_('CY28') = -2.4;
        v_('CY29') = -1.3;
        v_('CY30') = 0.0;
        v_('CY31') = -1.7;
        v_('CY32') = -2.1;
        v_('CY33') = -0.4;
        v_('CY34') = -2.9;
        v_('CY35') = -3.4;
        v_('CY36') = -4.3;
        v_('CY37') = -5.2;
        v_('CY38') = -6.5;
        v_('CY39') = -7.5;
        v_('CY40') = -6.4;
        v_('CY41') = -5.1;
        v_('CY42') = 0.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-1')
            v_('I+1') = I+v_('1');
            for J=v_('I+1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['OBJ1',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2mpjlib('ii',['OBJ2',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['Y',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['CONS',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['CONS',int2str(I)];
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['CONS',int2str(I)])) = v_(['R',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = -10;
            pb.xupper(ix_(['X',int2str(I)])) = 10;
            pb.xlower(ix_(['Y',int2str(I)]),1) = -10;
            pb.xupper(ix_(['Y',int2str(I)])) = 10;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eDIFSQR',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDIFSQR';
            ielftype(ie) = iet_('eDIFSQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['CX',int2str(I)]);
        end
        for I=v_('1'):v_('N')
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDIFSQR';
            ielftype(ie) = iet_('eDIFSQR');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('W',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['CY',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gSQUARE',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N-1')
            v_('I+1') = I+v_('1');
            for J=v_('I+1'):v_('N')
                ig = ig_(['OBJ1',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gSQUARE';
                ig = ig_(['OBJ2',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gSQUARE';
            end
        end
        for I=v_('1'):v_('N')
            ig = ig_(['CONS',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = .0;
%    Solution
% LO SOLTN              47952.695811
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'SQR2-MN-84-42';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eDIFSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DIF = EV_(1)-pbm.elpar{iel_}(1);
        varargout{1} = DIF*DIF;
        if(nargout>1)
            g_(1,1) = 2.0*DIF;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQUARE'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

