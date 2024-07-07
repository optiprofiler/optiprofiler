function varargout = AIRCRFTA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : AIRCRFTA
%    *********
% 
%    The aircraft stability problem by Rheinboldt, as a function
%    of the elevator, aileron and rudder deflection controls.
% 
%    Source: Problem 9 in
%    J.J. More',"A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer Seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'NOR2-RN-8-5'
% 
%    Values for the controls
%    1) Elevator
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'AIRCRFTA';

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
        v_('ELVVAL') = 0.1;
        v_('AILVAL') = 0.0;
        v_('RUDVAL') = 0.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','ROLLRATE',ix_);
        pb.xnames{iv} = 'ROLLRATE';
        [iv,ix_] = s2mpjlib('ii','PITCHRAT',ix_);
        pb.xnames{iv} = 'PITCHRAT';
        [iv,ix_] = s2mpjlib('ii','YAWRATE',ix_);
        pb.xnames{iv} = 'YAWRATE';
        [iv,ix_] = s2mpjlib('ii','ATTCKANG',ix_);
        pb.xnames{iv} = 'ATTCKANG';
        [iv,ix_] = s2mpjlib('ii','SSLIPANG',ix_);
        pb.xnames{iv} = 'SSLIPANG';
        [iv,ix_] = s2mpjlib('ii','ELEVATOR',ix_);
        pb.xnames{iv} = 'ELEVATOR';
        [iv,ix_] = s2mpjlib('ii','AILERON',ix_);
        pb.xnames{iv} = 'AILERON';
        [iv,ix_] = s2mpjlib('ii','RUDDERDF',ix_);
        pb.xnames{iv} = 'RUDDERDF';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','G1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G1';
        iv = ix_('ROLLRATE');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.933+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.933;
        end
        iv = ix_('PITCHRAT');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.107+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.107;
        end
        iv = ix_('YAWRATE');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.126+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.126;
        end
        iv = ix_('SSLIPANG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.99+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.99;
        end
        iv = ix_('AILERON');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -45.83+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -45.83;
        end
        iv = ix_('RUDDERDF');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.64+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.64;
        end
        [ig,ig_] = s2mpjlib('ii','G2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G2';
        iv = ix_('PITCHRAT');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.987+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.987;
        end
        iv = ix_('ATTCKANG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -22.95+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -22.95;
        end
        iv = ix_('ELEVATOR');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -28.37+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -28.37;
        end
        [ig,ig_] = s2mpjlib('ii','G3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G3';
        iv = ix_('ROLLRATE');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.002;
        end
        iv = ix_('YAWRATE');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.235+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.235;
        end
        iv = ix_('SSLIPANG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.67;
        end
        iv = ix_('AILERON');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.921+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.921;
        end
        iv = ix_('RUDDERDF');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.51+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.51;
        end
        [ig,ig_] = s2mpjlib('ii','G4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G4';
        iv = ix_('PITCHRAT');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('ATTCKANG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('ELEVATOR');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.168+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.168;
        end
        [ig,ig_] = s2mpjlib('ii','G5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G5';
        iv = ix_('YAWRATE');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('SSLIPANG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.196+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.196;
        end
        iv = ix_('AILERON');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.0071+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.0071;
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
        pb.xlower(ix_('ELEVATOR'),1) = v_('ELVVAL');
        pb.xupper(ix_('ELEVATOR'),1) = v_('ELVVAL');
        pb.xlower(ix_('AILERON'),1) = v_('AILVAL');
        pb.xupper(ix_('AILERON'),1) = v_('AILVAL');
        pb.xlower(ix_('RUDDERDF'),1) = v_('RUDVAL');
        pb.xupper(ix_('RUDDERDF'),1) = v_('RUDVAL');
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        if(isKey(ix_,'ELEVATOR'))
            pb.x0(ix_('ELEVATOR'),1) = v_('ELVVAL');
        else
            pb.y0(find(pbm.congrps==ig_('ELEVATOR')),1) = v_('ELVVAL');
        end
        if(isKey(ix_,'AILERON'))
            pb.x0(ix_('AILERON'),1) = v_('AILVAL');
        else
            pb.y0(find(pbm.congrps==ig_('AILERON')),1) = v_('AILVAL');
        end
        if(isKey(ix_,'RUDDERDF'))
            pb.x0(ix_('RUDDERDF'),1) = v_('RUDVAL');
        else
            pb.y0(find(pbm.congrps==ig_('RUDDERDF')),1) = v_('RUDVAL');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1A';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'PITCHRAT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'YAWRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E1B';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'YAWRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'ATTCKANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E1C';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'ATTCKANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SSLIPANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E1D';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'PITCHRAT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'ATTCKANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2A';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'ROLLRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'YAWRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2B';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'ROLLRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'SSLIPANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3A';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'ROLLRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PITCHRAT';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3B';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = 'ROLLRATE';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'ATTCKANG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('G1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1A');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.727;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E1B');
        pbm.grelw{ig}(posel) = 8.39;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -684.4;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E1D');
        pbm.grelw{ig}(posel) = 63.5;
        ig = ig_('G2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2A');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.949;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2B');
        pbm.grelw{ig}(posel) = 0.173;
        ig = ig_('G3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3A');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.716;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E3B');
        pbm.grelw{ig}(posel) = -1.578;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1D');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.132;
        ig = ig_('G4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2B');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('G5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3B');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-RN-8-5';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

