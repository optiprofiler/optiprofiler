function varargout = DTOC1NA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC1NA
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, NX control variables and NY state variables.
%    The nonlinearity parameter mu is set to 0.005.
% 
%    The problem is not convex.
% 
%    Sources: problem 1 in
%    T.F. Coleman and A. Liao,
%    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
%    Control Problems",
%    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    L.Z. Liao and C.A. Shoemaker,
%    "Advantages of differential dynamic programming over Newton's method for
%    discrete-time optimal control problems",
%    Tech. Report ctc92tr97, Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    SIF input: Ph. Toint, August 1993
% 
%    classification = 'C-COQR2-AN-V-V'
% 
%    Problem variants: they are identified by the values of
%    the parameter vector ( N, NX, NY )
% 
%    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
%    and (N-1)*NY constraints
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER # periods  } original value
% IE NX                  2              $-PARAMETER # controls } n=   58, m=  36
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   50             $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n=  298, m= 196
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   100            $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n=  598, m= 396
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   500            $-PARAMETER # periods  }
% IE NX                  2              $-PARAMETER # controls } n= 2998, m=1996
% IE NY                  4              $-PARAMETER # states   }
% 
% IE N                   1000           $-PARAMETER # periods  }
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DTOC1NA';

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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE NX                  2              $-PARAMETER # controls } n= 5998, m=3996
        if(nargs<2)
            v_('NX') = 2;  %  SIF file default value
        else
            v_('NX') = varargin{2};
        end
% IE NY                  4              $-PARAMETER # states   }
        if(nargs<3)
            v_('NY') = 4;  %  SIF file default value
        else
            v_('NY') = varargin{3};
        end
% IE N                   10             $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=  145, m=  90
% IE NY                  10             $-PARAMETER # states   }
% IE N                   50             $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=  745, m= 490
% IE NY                  10             $-PARAMETER # states   }
% IE N                   100            $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n= 1495, m= 990
% IE NY                  10             $-PARAMETER # states   }
% IE N                   500            $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n= 7495, m=4990
% IE NY                  10             $-PARAMETER # states   }
% IE N                   1000           $-PARAMETER # periods  }
% IE NX                  5              $-PARAMETER # controls } n=14995, m=9990
% IE NY                  10             $-PARAMETER # states   }
        if(nargs<4)
            v_('MU') = 0.005;  %  SIF file default value
        else
            v_('MU') = varargin{4};
        end
        v_('N-1') = -1+v_('N');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('NY-1') = -1+v_('NY');
        v_('NX+NY') = v_('NX')+v_('NY');
        v_('RXY') = v_('NX+NY');
        v_('1/RXY') = 1.0/v_('RXY');
        v_('MU/RXY') = v_('MU')*v_('1/RXY');
        v_('NYNX') = v_('NX')*v_('NY');
        v_('NYNX-1') = -1+v_('NYNX');
        for J=v_('1'):v_('NX')
            for I=v_('1'):v_('NY')
                v_('I-J') = I-J;
                v_('RI-J') = v_('I-J');
                v_(['B',int2str(I),',',int2str(J)]) = v_('RI-J')*v_('1/RXY');
                v_('I+J') = I+J;
                v_('RI+J') = v_('I+J');
                v_(['C',int2str(I),',',int2str(J)]) = v_('RI+J')*v_('MU/RXY');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(T),',',int2str(I)],ix_);
                pb.xnames{iv} = ['X',int2str(T),',',int2str(I)];
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(T),',',int2str(I)],ix_);
                pb.xnames{iv} = ['Y',int2str(T),',',int2str(I)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                [ig,ig_] = s2mpjlib('ii',['OX',int2str(T),',',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['OY',int2str(T),',',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['Y',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.5;
            end
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.25+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.25;
            end
            for I=v_('1'):v_('NX')
                [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) =...
                          v_(['B',int2str(round(v_('1'))),',',int2str(I)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('1'))),',',int2str(I)]);
                end
            end
            for J=v_('2'):v_('NY-1')
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(J)];
                iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['Y',int2str(T),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 0.5;
                end
                iv = ix_(['Y',int2str(T),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -0.25+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -0.25;
                end
                iv = ix_(['Y',int2str(T),',',int2str(round(v_('J+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 0.25+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 0.25;
                end
                for I=v_('1'):v_('NX')
                    [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(J)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['TT',int2str(T),',',int2str(J)];
                    iv = ix_(['X',int2str(T),',',int2str(I)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_(['B',int2str(J),',',int2str(I)])+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_(['B',int2str(J),',',int2str(I)]);
                    end
                end
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('NY')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] =...
                  s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('NY')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 0.5+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 0.5;
            end
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('NY-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -0.25+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -0.25;
            end
            for I=v_('1'):v_('NX')
                [ig,ig_] =...
                      s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('NY')))],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('NY')))];
                iv = ix_(['X',int2str(T),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('NY'))),',',int2str(I)])+...
                         pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['B',int2str(round(v_('NY'))),',',int2str(I)]);
                end
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
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                pbm.gconst(ig_(['OX',int2str(T),',',int2str(I)])) = -0.5;
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                pbm.gconst(ig_(['OY',int2str(T),',',int2str(I)])) = -0.25;
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('1'):v_('NY')
            pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'MUC';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for T=v_('1'):v_('N-1')
            for K=v_('0'):v_('NYNX-1')
                v_('I') = fix(K/v_('NX'));
                v_('INX') = v_('I')*v_('NX');
                v_('J') = K-v_('INX');
                v_('I') = 1+v_('I');
                v_('J') = 1+v_('J');
                ename = ['E',int2str(T),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePR';
                ielftype(ie) = iet_('ePR');
                vname = ['Y',int2str(T),',',int2str(round(v_('I')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(T),',',int2str(round(v_('J')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('MUC',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) =...
                      v_(['C',int2str(round(v_('I'))),',',int2str(round(v_('J')))]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for T=v_('1'):v_('N-1')
            for I=v_('1'):v_('NX')
                ig = ig_(['OX',int2str(T),',',int2str(I)]);
                pbm.grftype{ig} = 'gL4';
            end
        end
        for T=v_('1'):v_('N')
            for I=v_('1'):v_('NY')
                ig = ig_(['OY',int2str(T),',',int2str(I)]);
                pbm.grftype{ig} = 'gL4';
            end
        end
        for T=v_('1'):v_('N-1')
            for J=v_('1'):v_('NY')
                for K=v_('0'):v_('NYNX-1')
                    ig = ig_(['TT',int2str(T),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(T),',',int2str(K)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO S(  10,2, 4)        0.07531263517
% LO S(  50,2, 4)        0.23949933830
% LO S( 100,2, 4)        0.44472956731
% LO S( 500,2, 4)        2.08660731072
% LO S(1000,2, 4)        4.13895890755
% LO S(  10,5,10)        1.16782170926
% LO S(  50,5,10)        6.29429519092
% LO S( 100,5,10)        12.7020120922
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COQR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
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

    case 'ePR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2);
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
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

