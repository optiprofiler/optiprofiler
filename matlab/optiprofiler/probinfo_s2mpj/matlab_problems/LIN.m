function varargout = LIN(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LIN
%    *********
% 
%    A non-convex global optimization chemical equilibrium problem from the 
%    thesis of W.J. Lin.
%    It has a nonlinear objective and linear constraints.
% 
%    Source: illustrative example (section 4.6) in
%    C.M. McDonald and C.A. Floudas, "Global optimization for the phase 
%    and chemical equilibrium problem: application to the NRTL equation",
%    Computers & Chemical Engineering, (submitted), 1994.
% 
%    SIF input: Marcel Mongeau, 9 February 1994.
% 
%    classification = 'OLR2-AY-4-2'
% 
%    PARAMETERS likely to be changed for different problems:
% 
%    Number of variable sets (# of phases)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LIN';

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
        v_('P') = 2;
        v_('C') = 2;
        v_('1') = 1;
        v_('2') = 2;
        v_(['TAU',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 3.00498;
        v_(['TAU',int2str(round(v_('2'))),',',int2str(round(v_('1')))]) = 4.69071;
        v_(['ALF',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 0.391965;
        v_(['ALF',int2str(round(v_('2'))),',',int2str(round(v_('1')))]) = 0.391965;
        v_(['INIT',int2str(round(v_('1')))]) = 0.5;
        v_(['INIT',int2str(round(v_('2')))]) = 0.5;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('C')
            for K=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(K)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(K)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('C')
            for K=v_('1'):v_('P')
                [ig,ig_] = s2mpjlib('ii',['MB',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['MB',int2str(I)];
                iv = ix_(['X',int2str(I),',',int2str(K)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
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
        for I=v_('1'):v_('C')
            pbm.gconst(ig_(['MB',int2str(I)])) = v_(['INIT',int2str(I)]);
        end
        v_('ZERO') = 0.0;
        v_('ONE') = 1.0;
        for I=v_('1'):v_('C')
            v_(['ALF',int2str(I),',',int2str(I)]) = v_('ZERO');
            v_(['TAU',int2str(I),',',int2str(I)]) = v_('ZERO');
        end
        for I=v_('1'):v_('C')
            for J=v_('1'):v_('C')
                v_('MALF') = -1.0*v_(['ALF',int2str(I),',',int2str(J)]);
                v_('PROD') = v_('MALF')*v_(['TAU',int2str(I),',',int2str(J)]);
                v_(['G',int2str(I),',',int2str(J)]) = exp(v_('PROD'));
            end
        end
        for I=v_('1'):v_('C')
            v_(['G',int2str(I),',',int2str(I)]) = v_('ONE');
        end
        for I=v_('1'):v_('C')
            for J=v_('1'):v_('C')
                v_(['M',int2str(I),',',int2str(J)]) = v_(['G',int2str(I),',',int2str(J)])*...
                     v_(['TAU',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 1.e-12*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))])) =...
              v_('INIT1');
        pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))])) =...
              v_('INIT1');
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))])) =...
              v_('INIT2');
        pb.xupper(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))])) =...
              v_('INIT2');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) =...
              0.5;
        pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) =...
              0.0;
        pb.x0(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('1')))]),1) =...
              0.0;
        pb.x0(ix_(['X',int2str(round(v_('2'))),',',int2str(round(v_('2')))]),1) =...
              0.5;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXTAUG1',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftp{it}{1} = 'G11';
        elftp{it}{2} = 'G12';
        elftp{it}{3} = 'G21';
        elftp{it}{4} = 'G22';
        elftp{it}{5} = 'M11';
        elftp{it}{6} = 'M12';
        elftp{it}{7} = 'M21';
        elftp{it}{8} = 'M22';
        [it,iet_] = s2mpjlib( 'ii', 'eXTAUG2',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftp{it}{1} = 'G11';
        elftp{it}{2} = 'G12';
        elftp{it}{3} = 'G21';
        elftp{it}{4} = 'G22';
        elftp{it}{5} = 'M11';
        elftp{it}{6} = 'M12';
        elftp{it}{7} = 'M21';
        elftp{it}{8} = 'M22';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGXC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y1';
        elftv{it}{3} = 'Y2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for K=v_('1'):v_('P')
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXTAUG1';
            ielftype(ie) = iet_('eXTAUG1');
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('1'))),',',int2str(K)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('2'))),',',int2str(K)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G11',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G12',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G21',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G22',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M11',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M12',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M21',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('1'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M22',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXTAUG2';
            ielftype(ie) = iet_('eXTAUG2');
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('1'))),',',int2str(K)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('2'))),',',int2str(K)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G11',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G12',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G21',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('G22',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['G',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M11',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('1'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M12',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('1'))),',',int2str(round(v_('2')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M21',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('2'))),',',int2str(round(v_('1')))]);
            ename = ['A',int2str(round(v_('2'))),',',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('M22',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) =...
                  v_(['M',int2str(round(v_('2'))),',',int2str(round(v_('2')))]);
        end
        for K=v_('1'):v_('P')
            for I=v_('1'):v_('C')
                ename = ['B',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXLOGX';
                ielftype(ie) = iet_('eXLOGX');
                vname = ['X',int2str(I),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for K=v_('1'):v_('P')
            for I=v_('1'):v_('C')
                ename = ['C',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eXLOGXC';
                ielftype(ie) = iet_('eXLOGXC');
                vname = ['X',int2str(I),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('1'))),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
                posev = find(strcmp('Y1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('2'))),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.e-12,[],[]);
                posev = find(strcmp('Y2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for K=v_('1'):v_('P')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('1'))),',',int2str(K)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(round(v_('2'))),',',int2str(K)]);
            pbm.grelw{ig}(posel) = 1.;
            for I=v_('1'):v_('C')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = -1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% Global minimum:        -0.02020 : 
%  XV LIN       X(1,1)    0.00071
%  XV LIN       X(1,2)    0.49929
%  XV LIN       X(2,1)    0.15588
%  XV LIN       X(2,2)    0.34412
% local minimum:         -0.01961 :
%  XV LIN       X(1,1)    0.00213
%  XV LIN       X(1,2)    0.49787
%  XV LIN       X(2,1)    0.46547
%  XV LIN       X(2,2)    0.03453
% local maximum:         -0.01730 :
%  XV LIN       X(1,1)    0.00173
%  XV LIN       X(1,2)    0.49827
%  XV LIN       X(2,1)    0.37544
%  XV LIN       X(2,2)    0.12456
% LO SOLTN               -0.02020
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OLR2-AY-4-2';
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

    case 'eXTAUG1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        INSUM1 = EV_(1)*pbm.elpar{iel_}(1)+EV_(2)*pbm.elpar{iel_}(3);
        INSUM2 = EV_(1)*pbm.elpar{iel_}(2)+EV_(2)*pbm.elpar{iel_}(4);
        RATIO1 = pbm.elpar{iel_}(5)/INSUM1;
        RATIO2 = pbm.elpar{iel_}(6)/INSUM2;
        TERM1 = EV_(1)*RATIO1;
        TERM2 = EV_(2)*RATIO2;
        SUM = TERM1+TERM2;
        SQ1 = TERM1/INSUM1;
        SQ2 = TERM2/INSUM2;
        SQ11 = SQ1*pbm.elpar{iel_}(1);
        SQ12 = SQ2*pbm.elpar{iel_}(2);
        SQ21 = SQ1*pbm.elpar{iel_}(3);
        SQ22 = SQ2*pbm.elpar{iel_}(4);
        TRI1 = RATIO1-SQ11-SQ12;
        TRI2 = RATIO2-SQ21-SQ22;
        CUB1 = SQ11/INSUM1;
        CUB2 = SQ12/INSUM2;
        CUB11 = CUB1*pbm.elpar{iel_}(1);
        CUB12 = CUB2*pbm.elpar{iel_}(2);
        CUBM21 = CUB1*pbm.elpar{iel_}(3);
        CUBM22 = CUB2*pbm.elpar{iel_}(4);
        CUB21 = SQ21*pbm.elpar{iel_}(3)/INSUM1;
        CUB22 = SQ22*pbm.elpar{iel_}(4)/INSUM2;
        H1 = RATIO2-SQ22-2*SQ21;
        H2 = pbm.elpar{iel_}(6)*pbm.elpar{iel_}(2)/INSUM2^2;
        H3 = SQ1*pbm.elpar{iel_}(3)^2/INSUM1;
        H4 = pbm.elpar{iel_}(6)*pbm.elpar{iel_}(4)/INSUM2^2;
        H5 = SQ2*pbm.elpar{iel_}(4)^2/INSUM2;
        varargout{1} = EV_(1)*SUM;
        if(nargout>1)
            g_(1,1) = SUM+EV_(1)*TRI1;
            g_(2,1) = EV_(1)*TRI2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2*(TRI1+EV_(1)*(-SQ11+CUB11+CUB12));
                H_(1,2) = H1+EV_(1)*(-H2+2*(CUBM21+CUBM22));
                H_(2,1) = H_(1,2);
                H_(2,2) = 2*EV_(1)*(H3-H4+H5);
                varargout{3} = H_;
            end
        end

    case 'eXTAUG2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        INSUM1 = EV_(1)*pbm.elpar{iel_}(1)+EV_(2)*pbm.elpar{iel_}(3);
        INSUM2 = EV_(1)*pbm.elpar{iel_}(2)+EV_(2)*pbm.elpar{iel_}(4);
        RATIO1 = pbm.elpar{iel_}(7)/INSUM1;
        RATIO2 = pbm.elpar{iel_}(8)/INSUM2;
        TERM1 = EV_(1)*RATIO1;
        TERM2 = EV_(2)*RATIO2;
        SUM = TERM1+TERM2;
        SQ1 = TERM1/INSUM1;
        SQ2 = TERM2/INSUM2;
        SQ11 = SQ1*pbm.elpar{iel_}(1);
        SQ12 = SQ2*pbm.elpar{iel_}(2);
        SQ21 = SQ1*pbm.elpar{iel_}(3);
        SQ22 = SQ2*pbm.elpar{iel_}(4);
        TRI1 = RATIO1-SQ11-SQ12;
        TRI2 = RATIO2-SQ21-SQ22;
        CUB1 = SQ11/INSUM1;
        CUB2 = SQ12/INSUM2;
        CUB11 = CUB1*pbm.elpar{iel_}(1);
        CUB12 = CUB2*pbm.elpar{iel_}(2);
        CUBM21 = CUB1*pbm.elpar{iel_}(3);
        CUBM22 = CUB2*pbm.elpar{iel_}(4);
        CUB21 = SQ21*pbm.elpar{iel_}(3)/INSUM1;
        CUB22 = SQ22*pbm.elpar{iel_}(4)/INSUM2;
        H1 = RATIO1-SQ11-2*SQ12;
        H2 = pbm.elpar{iel_}(7)*pbm.elpar{iel_}(3)/INSUM1^2;
        H3 = SQ1*pbm.elpar{iel_}(1)^2/INSUM1;
        H4 = pbm.elpar{iel_}(7)*pbm.elpar{iel_}(1)/INSUM1^2;
        H5 = SQ2*pbm.elpar{iel_}(2)^2/INSUM2;
        varargout{1} = EV_(2)*SUM;
        if(nargout>1)
            g_(1,1) = EV_(2)*TRI1;
            g_(2,1) = SUM+EV_(2)*TRI2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2*EV_(2)*(H3-H4+H5);
                H_(1,2) = H1+EV_(2)*(-H2+2*(CUBM21+CUBM22));
                H_(2,1) = H_(1,2);
                H_(2,2) = 2*(TRI2+EV_(2)*(-SQ22+CUB21+CUB22));
                varargout{3} = H_;
            end
        end

    case 'eXLOGX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LOGX = log(EV_(1));
        varargout{1} = EV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX+1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0/EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eXLOGXC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(2,2) = U_(2,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(1,1) = U_(1,1)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        LOGX = log(IV_(2));
        varargout{1} = IV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = LOGX;
            g_(2,1) = IV_(1)/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0/IV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = -IV_(1)/IV_(2)^2;
                varargout{3} = U_.'*H_*U_;
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

