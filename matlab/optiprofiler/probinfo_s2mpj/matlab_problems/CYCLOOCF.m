function varargout = CYCLOOCF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    The cyclooctane molecule is comprised of eight carbon atoms aligned
%    in an equally spaced ring. When they take a position of minimum
%    potential energy so that next-neighbours are equally spaced.
% 
%    Given positions v_1, ..., v_p in R^3 (with p = 8 for cyclooctane),
%    and given a spacing c^2 we have that
% 
%       ||v_i - v_i+1,mod p||^2 = c^2 for i = 1,..,p, and
%       ||v_i - v_i+2,mod p||^2 = 2p/(p-2) c^3
% 
%    where (arbitrarily) we have v_1 = 0 and component 1 of v_2 = 0
% 
%    Source:
%    an extension of the cyclooctane molecule configuration space as
%    described in (for example)
% 
%     E. Coutsias, S. Martin, A. Thompson & J. Watson
%     "Topology of cyclooctane energy landscape"
%     J. Chem. Phys. 132-234115 (2010)
% 
%    SIF input: Nick Gould, Feb 2020.
% 
%    This is a version of CYCLOOPT.SIF without the fixed variables
% 
%    classification = 'NQR2-MN-V-V'
% 
%    The number of molecules
% 
%       Alternative values for the SIF file parameters:
% IE P                   8              $-PARAMETER     original value
% IE P                   100            $-PARAMETER
% IE P                   1000           $-PARAMETER
% IE P                   10000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CYCLOOCF';

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
            v_('P') = 8;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
% IE P                   100000         $-PARAMETER
        v_('C') = 1.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('P-1') = -1+v_('P');
        v_('P-2') = -2+v_('P');
        v_('THREE') = 3.0;
        v_('RP') = v_('P');
        v_('2RP') = 2.0*v_('RP');
        v_('RP-2') = -2.0+v_('RP');
        v_('2RP/RP-2') = v_('2RP')/v_('RP-2');
        v_('C2') = v_('C')*v_('C');
        v_('SC2') = v_('2RP/RP-2')*v_('C2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','Y2',ix_);
        pb.xnames{iv} = 'Y2';
        [iv,ix_] = s2mpjlib('ii','Z2',ix_);
        pb.xnames{iv} = 'Z2';
        for I=v_('3'):v_('P')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('P')
            [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(I)];
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['B',int2str(I)];
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
        for I=v_('1'):v_('P')
            pbm.gconst(ig_(['A',int2str(I)])) = v_('C2');
            pbm.gconst(ig_(['B',int2str(I)])) = v_('SC2');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('RI') = 0.0;
        v_('START') = v_('RI')/v_('RP');
        pb.x0(ix_('Y2'),1) = v_('START');
        pb.x0(ix_('Z3'),1) = v_('START');
        for I=v_('3'):v_('P')
            v_('RI') = I;
            v_('START') = v_('RI')/v_('RP');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('START');
            pb.x0(ix_(['Y',int2str(I)]),1) = v_('START');
            pb.x0(ix_(['Z',int2str(I)]),1) = v_('START');
        end
        v_('RI') = 2.0;
        v_('START') = v_('RI')/v_('RP');
        for I=v_('3'):v_('P')
            v_('RI') = I;
            v_('START') = v_('RI')/v_('RP');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eSQRDIF',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        v_('I') = 1;
        v_('I+1') = 1+v_('I');
        ename = ['AY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I') = 2;
        v_('I+1') = 1+v_('I');
        ename = ['AX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('I+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['AY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['AZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I+1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('3'):v_('P-1')
            v_('I+1') = 1+I;
            ename = ['AX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['AY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['AZ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['AX',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AX',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AY',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AY',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['AZ',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['AZ',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I') = 1;
        v_('I+2') = 2+v_('I');
        ename = ['BX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        v_('I') = 2;
        v_('I+2') = 2+v_('I');
        ename = ['BX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BX',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['BY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['BZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('I')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('I+2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('3'):v_('P-2')
            v_('I+2') = 2+I;
            ename = ['BX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['BY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(round(v_('I+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['BZ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQRDIF';
            ielftype(ie) = iet_('eSQRDIF');
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(round(v_('I+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['BX',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BX',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('P-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BY',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('P-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BZ',int2str(round(v_('P-1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('P-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BX',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['BX',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['BY',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BY',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQRDIF';
        ielftype(ie) = iet_('eSQRDIF');
        ename = ['BZ',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('P')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['BZ',int2str(round(v_('P')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Z',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        v_('I') = 1;
        ig = ig_(['A',int2str(round(v_('I')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['AY',int2str(round(v_('I')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['AZ',int2str(round(v_('I')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['B',int2str(round(v_('I')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['BX',int2str(round(v_('I')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['BY',int2str(round(v_('I')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['BZ',int2str(round(v_('I')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('2'):v_('P')
            ig = ig_(['A',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['AX',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['AY',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['AZ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['B',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['BX',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['BY',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['BZ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION             0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NQR2-MN-V-V';
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

    case 'eSQR'

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

    case 'eSQRDIF'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

