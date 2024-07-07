function varargout = BLOCKQP2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BLOCKQP2
%    *********
% 
%    A non-convex quadratic program with some structure.
% 
%    The objective function is of the form
%       sum (i=1,n) x_i y_i + sum(j=1,b) z_j^2
%    There are n equality constraints of the form
%         x_i - y_i + sum (j=1,b) z_j = b   
%    There is an inequality constraint of the form
%       sum(i=1,n) x_i + y_i + sum(j=1,b) z_j >= b + 1
%    Finally, there are simple bounds
%          1 <= x_i, y_i <= 1    (i=1,n)
%          0 <= z_j <= 2         (j=1,b)
% 
%    SIF input: Nick Gould, June 1994
% 
%    classification = 'QLR2-AN-V-V'
% 
%    The number of equality constraints
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BLOCKQP2';

switch(action)

    case 'setup'

    pb.name      = 'BLOCKQP2';
    pb.sifpbname = 'BLOCKQP2';
    pbm.name     = 'BLOCKQP2';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        if(nargin<3)
            v_('B') = 5;  %  SIF file default value
        else
            v_('B') = varargin{2};
        end
        v_('1') = 1;
        v_('RB') = v_('B');
        v_('RB+1') = 1+v_('RB');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        for J=v_('1'):v_('B')
            [iv,ix_] = s2xlib('ii',['Z',int2str(J)],ix_);
            pb.xnames{iv} = ['Z',int2str(J)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N')
            [ig,ig_] = s2xlib('ii',I,ig_);
            gtype{ig}  = '>=';
            cnames{ig} = I;
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2xlib('ii',['E',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for J=v_('1'):v_('B')
            [ig,ig_] = s2xlib('ii','I',ig_);
            gtype{ig}  = '>=';
            cnames{ig} = 'I';
            iv = ix_(['Z',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            for I=v_('1'):v_('N')
                [ig,ig_] = s2xlib('ii',['E',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['E',int2str(I)];
                iv = ix_(['Z',int2str(J)]);
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['E',int2str(I)])) = v_('RB');
        end
        pbm.gconst(ig_('I')) = v_('RB+1');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['X',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['X',int2str(I)])) = 1.0;
            pb.xlower(ix_(['Y',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['Y',int2str(I)])) = 1.0;
        end
        for J=v_('1'):v_('B')
            pb.xlower(ix_(['Z',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['Z',int2str(J)])) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'SQ',iet_);
        elftv{it}{1} = 'Z';
        [it,iet_] = s2xlib( 'ii', 'PROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['P',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'PROD';
            ielftype(ie) = iet_('PROD');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('B')
            ename = ['S',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQ';
            ielftype(ie) = iet_('SQ');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
        end
        for J=v_('1'):v_('B')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QLR2-AN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'SQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.5*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0;
                varargout{3} = H_;
            end
        end

    case 'PROD'

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
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
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

