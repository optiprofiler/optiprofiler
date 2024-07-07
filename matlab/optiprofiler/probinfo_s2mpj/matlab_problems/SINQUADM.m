function varargout = SINQUAD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SINQUADM
%    *********
% 
%    Another function with nontrivial groups and
%    repetitious elements.
% 
%    PhT thinks that the objective function is given by
%    f(x1,...,x_n) = (x_1-1)^4 + sum_{i=2}^{n-1}(x_i^2-x_1^2+sin(x_i-x_n))^2
%                    + (x_n^2 - x_1^2)^2
% 
%    Source:
%    N. Gould, private communication.
% 
%    SIF input: N. Gould, Dec 1989.
% 
%    classification = 'OUR2-AY-V-0'
% 
%    number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SINQUADM';

switch(action)

    case 'setup'

    pb.name      = 'SINQUADM';
    pb.sifpbname = 'SINQUADM';
    pbm.name     = 'SINQUADM';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('NM1') = -1+v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','G1',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('2'):v_('N')
            [ig,ig_] = s2xlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('G1')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.1*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'SQ',iet_);
        elftv{it}{1} = 'V1';
        [it,iet_] = s2xlib( 'ii', 'SINE',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQ';
        ielftype(ie) = iet_('SQ');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('2'):v_('N')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQ';
            ielftype(ie) = iet_('SQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.1);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('2'):v_('NM1')
            ename = ['S',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SINE';
            ielftype(ie) = iet_('SINE');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.1);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(v_('N'))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.1);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','L2',igt_);
        [it,igt_] = s2xlib('ii','L4',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('G1');
        pbm.grftype{ig} = 'L4';
        for I=v_('2'):v_('N')
            ig = ig_(['G',int2str(I)]);
%            pbm.grftype{ig} = 'L2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(v_('1'))]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        for I=v_('2'):v_('NM1')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            pbm.grelw{ig}(posel) = 1;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = 'OUR2-AY-V-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'SQ'

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

    case 'SINE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = sin(IV_(1));
        if(nargout>1)
            g_(1,1) = cos(IV_(1));
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = -sin(IV_(1));
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'L2'

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

    case 'L4'

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
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',pb.name,' with action = setup'])
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

