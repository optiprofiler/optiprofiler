function varargout = DTOC2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC2
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, 2 control variables and 4 state variables.
% 
%    The problem is not convex.
% 
%    Sources: problem 2 in
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
%    classification = 'OOR2-AN-V-V'
% 
%    Problem variants: they are identified by the value of the parameter N.
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

persistent pbm;

name = 'DTOC2';

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
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('NY-1') = -1+v_('NY');
        v_('2NY') = v_('NY')+v_('NY');
        v_('R2NY') = v_('2NY');
        v_('1/2NY') = 1.0/v_('R2NY');
        for J=v_('1'):v_('NX')
            for I=v_('1'):v_('NY')
                v_('I+J') = I+J;
                v_('RI+J') = v_('I+J');
                v_(['C',int2str(I),',',int2str(J)]) = v_('RI+J')*v_('1/2NY');
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
        for T=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(T)],ig_);
            gtype{ig} = '<>';
        end
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            for J=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['TT',int2str(T),',',int2str(J)];
                iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for I=v_('1'):v_('NY')
            v_('RI') = I;
            v_('TMP') = v_('RI')*v_('1/2NY');
            pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = v_('TMP');
            pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = v_('TMP');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('NY')
            v_('RI') = I;
            v_('TMP') = v_('RI')*v_('1/2NY');
            pb.x0(ix_(['Y',int2str(round(v_('1'))),',',int2str(I)]),1) = v_('TMP');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOEL',iet_);
        elftv{it}{1} = 'YY1';
        elftv{it}{2} = 'YY2';
        elftv{it}{3} = 'YY3';
        elftv{it}{4} = 'YY4';
        elftv{it}{5} = 'XX1';
        elftv{it}{6} = 'XX2';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'YY';
        [it,iet_] = s2mpjlib( 'ii', 'eSINE',iet_);
        elftv{it}{1} = 'ZZ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for T=v_('1'):v_('N-1')
            ename = ['EO',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eOEL';
            ielftype(ie) = iet_('eOEL');
            vname = ['Y',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(T),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(T),',',int2str(round(v_('3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(T),',',int2str(round(v_('4')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XX1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(T),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XX2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            for J=v_('1'):v_('NY')
                ename = ['SY',int2str(T),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSINE';
                ielftype(ie) = iet_('eSINE');
                vname = ['Y',int2str(T),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ZZ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for I=v_('1'):v_('NX')
                ename = ['SX',int2str(T),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSINE';
                ielftype(ie) = iet_('eSINE');
                vname = ['X',int2str(T),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ZZ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for J=v_('1'):v_('NY')
            ename = ['YNSQ',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['Y',int2str(round(v_('N'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('YY',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for T=v_('1'):v_('N-1')
            ig = ig_(['O',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EO',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for J=v_('1'):v_('NY')
            ig = ig_(['O',int2str(round(v_('N')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['YNSQ',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for T=v_('1'):v_('N-1')
            for J=v_('1'):v_('NY')
                ig = ig_(['TT',int2str(T),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['SY',int2str(T),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                for I=v_('1'):v_('NX')
                    ig = ig_(['TT',int2str(T),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['SX',int2str(T),',',int2str(I)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_(['C',int2str(J),',',int2str(I)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION(  10)      0.485983918948
% LO SOLUTION(  20)      0.486212213803
% LO SOLUTION(  30)      0.486383392574
% LO SOLUTION(  40)      0.486572686778
% LO SOLUTION(  50)      0.486884900389
% LO SOLUTION( 100)      0.487532342563
% LO SOLUTION( 500)      0.490996540460
% LO SOLUTION(1000)      0.490200910983
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-V-V';
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

    case 'eSINE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SZ = sin(EV_(1));
        varargout{1} = SZ;
        if(nargout>1)
            g_(1,1) = cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -SZ;
                varargout{3} = H_;
            end
        end

    case 'eOEL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XN2 = EV_(5)*EV_(5)+EV_(6)*EV_(6);
        YN2 = EV_(1)*EV_(1)+EV_(2)*EV_(2)+EV_(3)*EV_(3)+EV_(4)*EV_(4);
        SZ = sin(0.5*XN2);
        CZ = cos(0.5*XN2);
        SZ2 = SZ*SZ+1.0;
        SC = SZ*CZ;
        CCSS = CZ*CZ-SZ*SZ;
        varargout{1} = YN2*SZ2;
        if(nargout>1)
            g_(5,1) = 2.0*YN2*SC*EV_(5);
            g_(6,1) = 2.0*YN2*SC*EV_(6);
            g_(1,1) = 2.0*EV_(1)*SZ2;
            g_(2,1) = 2.0*EV_(2)*SZ2;
            g_(3,1) = 2.0*EV_(3)*SZ2;
            g_(4,1) = 2.0*EV_(4)*SZ2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(6,6);
                H_(5,5) = 2.0*YN2*(SC+EV_(5)*EV_(5)*CCSS);
                H_(5,6) = 2.0*YN2*EV_(5)*EV_(6)*CCSS;
                H_(6,5) = H_(5,6);
                H_(5,1) = 4.0*EV_(1)*SC*EV_(5);
                H_(1,5) = H_(5,1);
                H_(5,2) = 4.0*EV_(2)*SC*EV_(5);
                H_(2,5) = H_(5,2);
                H_(5,3) = 4.0*EV_(3)*SC*EV_(5);
                H_(3,5) = H_(5,3);
                H_(5,4) = 4.0*EV_(4)*SC*EV_(5);
                H_(4,5) = H_(5,4);
                H_(6,6) = 2.0*YN2*(SC+EV_(6)*EV_(6)*CCSS);
                H_(6,1) = 4.0*EV_(1)*SC*EV_(6);
                H_(1,6) = H_(6,1);
                H_(6,2) = 4.0*EV_(2)*SC*EV_(6);
                H_(2,6) = H_(6,2);
                H_(6,3) = 4.0*EV_(3)*SC*EV_(6);
                H_(3,6) = H_(6,3);
                H_(6,4) = 4.0*EV_(4)*SC*EV_(6);
                H_(4,6) = H_(6,4);
                H_(1,1) = 2.0*SZ2;
                H_(2,2) = 2.0*SZ2;
                H_(3,3) = 2.0*SZ2;
                H_(4,4) = 2.0*SZ2;
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

