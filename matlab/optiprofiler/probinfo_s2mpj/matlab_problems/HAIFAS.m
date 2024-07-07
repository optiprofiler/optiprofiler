function varargout = HAIFAS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HAIFAS
%    *********
% 
%   Truss Topology Design (t6-9)
% 
%   Source: M. Tsibulevsky, Optimization Laboratory,
%           Faculty of Industrial Engineering, Technion,
%           Haifa, 32000, Israel.
% 
%   SIF input: Conn, Gould and Toint, May, 1992
%              minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'LQR2-AN-13-9'
% 
%   2 * Number of nodes
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HAIFAS';

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
        v_('N') = 12;
        v_('M') = 9;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','Z',ix_);
        pb.xnames{iv} = 'Z';
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_('Z');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            v_('J') = 10;
            iv = ix_(['X',int2str(round(v_('J')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.00000+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.00000;
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
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
        v_('I1') = 4;
        v_('I2') = 4;
        v_('L') = 1;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 5;
        v_('I2') = 5;
        v_('L') = 2;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 5;
        v_('I2') = 11;
        v_('L') = 3;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 11;
        v_('I2') = 11;
        v_('L') = 4;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 10;
        v_('I2') = 10;
        v_('L') = 5;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 10;
        v_('I2') = 11;
        v_('L') = 6;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 11;
        v_('I2') = 11;
        v_('L') = 7;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 4;
        v_('I2') = 4;
        v_('L') = 8;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 4;
        v_('I2') = 10;
        v_('L') = 9;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 10;
        v_('I2') = 10;
        v_('L') = 10;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 5;
        v_('I2') = 5;
        v_('L') = 11;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 6;
        v_('I2') = 6;
        v_('L') = 12;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 6;
        v_('I2') = 12;
        v_('L') = 13;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 12;
        v_('I2') = 12;
        v_('L') = 14;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 11;
        v_('I2') = 11;
        v_('L') = 15;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 11;
        v_('I2') = 12;
        v_('L') = 16;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 12;
        v_('I2') = 12;
        v_('L') = 17;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 5;
        v_('I2') = 5;
        v_('L') = 18;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 5;
        v_('I2') = 11;
        v_('L') = 19;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 11;
        v_('I2') = 11;
        v_('L') = 20;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I1') = 6;
        v_('I2') = 6;
        v_('L') = 21;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('L')))];
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
        end
        vname = ['X',int2str(round(v_('I2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        v_('J') = 1;
        v_('L') = 1;
        ig = ig_(['G',int2str(round(v_('J')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.00000e+01;
        v_('J') = 2;
        v_('L') = 2;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 2;
        v_('L') = 3;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 2;
        v_('L') = 4;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.60000;
        v_('J') = 3;
        v_('L') = 5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.00000e+01;
        v_('J') = 3;
        v_('L') = 6;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -8.00000e+01;
        v_('J') = 3;
        v_('L') = 7;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.00000e+01;
        v_('J') = 4;
        v_('L') = 8;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 4;
        v_('L') = 9;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -6.40000;
        v_('J') = 4;
        v_('L') = 10;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.60000;
        v_('J') = 5;
        v_('L') = 11;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.00000e+01;
        v_('J') = 6;
        v_('L') = 12;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 6;
        v_('L') = 13;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 6;
        v_('L') = 14;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.60000;
        v_('J') = 7;
        v_('L') = 15;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.00000e+01;
        v_('J') = 7;
        v_('L') = 16;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -8.00000e+01;
        v_('J') = 7;
        v_('L') = 17;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.00000e+01;
        v_('J') = 8;
        v_('L') = 18;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.40000;
        v_('J') = 8;
        v_('L') = 19;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -6.40000;
        v_('J') = 8;
        v_('L') = 20;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.60000;
        v_('J') = 9;
        v_('L') = 21;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('L')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.00000e+01;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-AN-13-9';
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

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.5*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 0.5*EV_(2);
            g_(2,1) = 0.5*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(2,1) = 0.5;
                H_(1,2) = H_(2,1);
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

