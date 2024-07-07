function varargout = GENHUMPS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : GENHUMPS
%    *********
% 
%    A multi-dimensional variant of HUMPS, a two dimensional function
%    with a lot of humps. The density of humps increases with the
%    parameter ZETA, making the problem more difficult.
% 
%    The problem is nonconvex.
% 
%    Source:
%    Ph. Toint, private communication, 1997.
% 
%    SDIF input: N. Gould and Ph. Toint, November 1997.
% 
%    classification = 'OUR2-AN-V-0'
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   5              $-PARAMETER
% IE N                   10             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER    original value
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GENHUMPS';

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
% RE ZETA                2.0            $-PARAMETER
        if(nargs<2)
            v_('ZETA') = 20.0;  %  SIF file default value
        else
            v_('ZETA') = varargin{2};
        end
        v_('1') = 1;
        v_('N-1') = -1+v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = -506.2*ones(pb.n,1);
        pb.x0(ix_(['X',int2str(round(v_('1')))]),1) = -506.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eHMP',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'A';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            ename = ['Q',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHMP';
            ielftype(ie) = iet_('eHMP');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-506.2);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-506.2);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('ZETA');
            ename = ['SX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-506.2);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['SY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-506.2);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SX',int2str(I)]);
            pbm.grelw{ig}(posel) = 0.05;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['SY',int2str(I)]);
            pbm.grelw{ig}(posel) = 0.05;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLTN               0.0
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

    case 'eHMP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SAX = sin(pbm.elpar{iel_}(1)*EV_(1));
        SAY = sin(pbm.elpar{iel_}(1)*EV_(2));
        CAX = cos(pbm.elpar{iel_}(1)*EV_(1));
        CAY = cos(pbm.elpar{iel_}(1)*EV_(2));
        AA = pbm.elpar{iel_}(1)+pbm.elpar{iel_}(1);
        AAA = pbm.elpar{iel_}(1)*AA;
        varargout{1} = (SAX*SAY)^2;
        if(nargout>1)
            g_(1,1) = AA*SAX*CAX*SAY^2;
            g_(2,1) = AA*SAX^2*CAY*SAY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = AAA*SAY^2*(CAX^2-SAX^2);
                H_(1,2) = AA*AA*SAX*CAX*SAY*CAY;
                H_(2,1) = H_(1,2);
                H_(2,2) = AAA*SAX^2*(CAY^2-SAY^2);
                varargout{3} = H_;
            end
        end

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^2;
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

