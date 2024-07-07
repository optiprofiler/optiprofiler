function varargout = DIXMAANM1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DIXMAANM1
%    *********
%    A variant on the Dixon-Maany test problem (version I) but removing 
%    elements/groups of type 2 since the parameter beta=0
% 
%    Source:
%    L. Luksan, C. Matonoha and J. Vlcek  
%    Modified CUTE problems for sparse unconstraoined optimization
%    Technical Report 1081
%    Institute of Computer Science
%    Academy of Science of the Czech Republic
% 
%    (problem 19) based on
% 
%    L.C.W. Dixon and Z. Maany,
%    "A family of test problems with sparse Hessians for unconstrained
%    optimization",
%    TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
% 
%    SIF input: Ph. Toint, Dec 1989.
%               correction by Ph. Shott, January 1995.
%               this version Nick Gould, June, 2013
%               update Nick Gould, August 2022, to remove beta=0 terms.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    M is equal to the third of the number of variables
% 
%       Alternative values for the SIF file parameters:
% IE M                   5              $-PARAMETER n = 15  original value 
% IE M                   30             $-PARAMETER n = 90
% IE M                   100            $-PARAMETER n = 300
% IE M                   500            $-PARAMETER n = 1500
% IE M                   1000           $-PARAMETER n = 3000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DIXMAANM1';

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
            v_('M') = 5;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
% IE M                   3000           $-PARAMETER n = 9000
        v_('N') = 3*v_('M');
        v_('ALPHA') = 1.0;
        v_('BETA') = 0.0;
        v_('GAMMA') = 0.125;
        v_('DELTA') = 0.125;
        v_('K1') = 2;
        v_('K3') = 1;
        v_('K4') = 2;
        v_('RN') = v_('N');
        v_('N-1') = -1+v_('N');
        v_('2M') = v_('M')+v_('M');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','GA',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','GC',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','GD',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('GA')) = -1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 2.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('2M')
            v_('I+M') = I+v_('M');
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQC';
            ielftype(ie) = iet_('eSQC');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+M')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('M')
            v_('I+2M') = I+v_('2M');
            ename = ['D',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+2M')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('I/N') = v_('RI')/v_('RN');
            v_('TMP') = 1.0;
            for J=v_('1'):v_('K1')
                v_('TMP') = v_('TMP')*v_('I/N');
            end
            v_('AI') = v_('TMP')*v_('ALPHA');
            ig = ig_('GA');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('AI');
        end
        for I=v_('1'):v_('2M')
            v_('RI') = I;
            v_('I/N') = v_('RI')/v_('RN');
            v_('TMP') = 1.0;
            for J=v_('1'):v_('K3')
                v_('TMP') = v_('TMP')*v_('I/N');
            end
            v_('CI') = v_('TMP')*v_('GAMMA');
            ig = ig_('GC');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('CI');
        end
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('I/N') = v_('RI')/v_('RN');
            v_('TMP') = 1.0;
            for J=v_('1'):v_('K4')
                v_('TMP') = v_('TMP')*v_('I/N');
            end
            v_('DI') = v_('TMP')*v_('DELTA');
            ig = ig_('GD');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('DI');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               1.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-V-0';
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

    case 'eSQC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        F1 = EV_(1)*EV_(1);
        F2 = EV_(2)^4;
        varargout{1} = F1*F2;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*F2;
            g_(2,1) = 4.0*F1*EV_(2)^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*F2;
                H_(1,2) = 8.0*EV_(1)*EV_(2)^3;
                H_(2,1) = H_(1,2);
                H_(2,2) = 12.0*F1*EV_(2)^2;
                varargout{3} = H_;
            end
        end

    case 'en2PR'

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

