function varargout = DRUGDISE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DRUGDISE
%    *********
% 
%    This is a variant of the drug displacement problem DRUGDIS where the
%    state equations have been Expanded in term of more intermediate
%    functions, each one of them being less nonlinear.
% 
%    The problem is based on the kinetic model of Aarons and Rowland which
%    simulates the interaction of the two drugs (warfarin and phenylnutazone)
%    in a patient bloodstream.  The state variable are the concentrations of
%    unbound warfarin (w) and phenylbutazone (p).  The problem is to control
%    the rate of injection (u) of the pain-killing phenylbutazone so that both
%    drugs reach a specified steady-state in minimum time and the concentration
%    of warfarin does not rise above a toxicity level.
% 
%    The problem is discretized using the trapeziodal rule.  It is non-convex.
% 
%    Source:
%    H. Maurer and M. Wiegand,
%    "Numerical solution of a drug displacement problem with bounded state
%    variables",
%    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'LOR2-MY-V-V'
% 
%    Discretization: specify the number of interior points + 1
% 
%       Alternative values for the SIF file parameters:
% IE NI                  10             $-PARAMETER n=63, m=50 
% IE NI                  100            $-PARAMETER n=603, m=500   original value
% IE NI                  100            $-PARAMETER n=6003, m=5000 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DRUGDISE';

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
            v_('NI') = 10;  %  SIF file default value
        else
            v_('NI') = varargin{1};
        end
        if(nargs<2)
            v_('TOXIC') = 0.026;  %  SIF file default value
        else
            v_('TOXIC') = varargin{2};
        end
        if(nargs<3)
            v_('WSS') = 0.02;  %  SIF file default value
        else
            v_('WSS') = varargin{3};
        end
        if(nargs<4)
            v_('UMAX') = 8.0;  %  SIF file default value
        else
            v_('UMAX') = varargin{4};
        end
        if(nargs<5)
            v_('PSTART') = 0.0;  %  SIF file default value
        else
            v_('PSTART') = varargin{5};
        end
        if(nargs<6)
            v_('PFINAL') = 2.0;  %  SIF file default value
        else
            v_('PFINAL') = varargin{6};
        end
        if(nargs<7)
            v_('Z') = 46.4;  %  SIF file default value
        else
            v_('Z') = varargin{7};
        end
        v_('AVP') = v_('PSTART')+v_('PFINAL');
        v_('AVP') = 0.5*v_('AVP');
        v_('-Z') = -1.0*v_('Z');
        v_('-ZZ') = v_('Z')*v_('-Z');
        v_('NI-1') = -1+v_('NI');
        v_('RNI') = v_('NI');
        v_('-1/NI') = -1.0/v_('RNI');
        v_('-Z/NI') = v_('Z')*v_('-1/NI');
        v_('0') = 0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','TF',ix_);
        pb.xnames{iv} = 'TF';
        pb.xscale(iv,1) = 200.0;
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I)],ix_);
            pb.xnames{iv} = ['W',int2str(I)];
            pb.xscale(iv,1) = 0.02;
        end
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for I=v_('0'):v_('NI-1')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        for I=v_('0'):v_('NI-1')
            [iv,ix_] = s2mpjlib('ii',['A',int2str(I)],ix_);
            pb.xnames{iv} = ['A',int2str(I)];
            pb.xscale(iv,1) = 200.0;
        end
        for I=v_('0'):v_('NI-1')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
            pb.xscale(iv,1) = 200.0;
        end
        for I=v_('0'):v_('NI-1')
            [iv,ix_] = s2mpjlib('ii',['C',int2str(I)],ix_);
            pb.xnames{iv} = ['C',int2str(I)];
            pb.xscale(iv,1) = 0.0000001;
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','TFINAL',ig_);
        gtype{ig} = '<>';
        iv = ix_('TF');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 100.0;
        for I=v_('0'):v_('NI-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['EW',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EW',int2str(I)];
            iv = ix_(['W',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['W',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            pbm.gscale(ig,1) = 0.02;
            [ig,ig_] = s2mpjlib('ii',['EP',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EP',int2str(I)];
            iv = ix_(['P',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['P',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['EA',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EA',int2str(I)];
            iv = ix_(['A',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['P',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-Z')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-Z');
            end
            pbm.gscale(ig,1) = 200.0;
            [ig,ig_] = s2mpjlib('ii',['EB',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EB',int2str(I)];
            iv = ix_(['B',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['W',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-Z')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-Z');
            end
            pbm.gscale(ig,1) = 200.0;
            [ig,ig_] = s2mpjlib('ii',['EC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EC',int2str(I)];
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
        for I=v_('0'):v_('NI-1')
            pbm.gconst(ig_(['EA',int2str(I)])) = 232.0;
            pbm.gconst(ig_(['EB',int2str(I)])) = 232.0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('0'):v_('NI-1')
            pb.xlower(ix_(['C',int2str(I)])) = -Inf;
            pb.xupper(ix_(['C',int2str(I)]),1) = +Inf;
        end
        pb.xlower(ix_('TF'),1) = 200.0;
        for I=v_('0'):v_('NI')
            pb.xupper(ix_(['W',int2str(I)])) = v_('TOXIC');
        end
        for I=v_('0'):v_('NI-1')
            pb.xupper(ix_(['U',int2str(I)])) = v_('UMAX');
        end
        pb.xlower(ix_(['W',int2str(round(v_('0')))]),1) = v_('WSS');
        pb.xupper(ix_(['W',int2str(round(v_('0')))]),1) = v_('WSS');
        pb.xlower(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.xupper(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.xlower(ix_(['P',int2str(round(v_('0')))]),1) = v_('PSTART');
        pb.xupper(ix_(['P',int2str(round(v_('0')))]),1) = v_('PSTART');
        pb.xlower(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        pb.xupper(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('2W/10') = 0.2*v_('WSS');
        v_('2P/10') = 0.2*v_('AVP');
        v_('2(W+P)/10') = v_('2W/10')+v_('2P/10');
        v_('D') = 1.0+v_('2(W+P)/10');
        v_('DD') = v_('D')*v_('D');
        v_('ZP') = v_('AVP')*v_('Z');
        v_('ZW') = v_('WSS')*v_('Z');
        v_('AA') = v_('DD')+v_('ZP');
        v_('AA') = 232.0+v_('AA');
        v_('BB') = v_('DD')+v_('ZW');
        v_('BB') = 232.0+v_('BB');
        v_('AB') = v_('AA')*v_('BB');
        v_('WP') = v_('WSS')*v_('AVP');
        v_('-ZZWP') = v_('WP')*v_('-ZZ');
        v_('CD') = v_('AB')+v_('-ZZWP');
        v_('CC') = v_('DD')/v_('CD');
        for I=v_('0'):v_('NI-1')
            pb.x0(ix_(['W',int2str(I)]),1) = v_('WSS');
            pb.x0(ix_(['P',int2str(I)]),1) = v_('AVP');
            pb.x0(ix_(['U',int2str(I)]),1) = v_('UMAX');
            pb.x0(ix_(['A',int2str(I)]),1) = v_('AA');
            pb.x0(ix_(['B',int2str(I)]),1) = v_('BB');
            pb.x0(ix_(['C',int2str(I)]),1) = v_('CC');
        end
        pb.x0(ix_('TF'),1) = 240.0;
        pb.x0(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.x0(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en3S',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        [it,iet_] = s2mpjlib( 'ii', 'en3D2',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        [it,iet_] = s2mpjlib( 'ii', 'eDSQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'en3PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('NI-1')
            ename = ['WA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3S';
            ielftype(ie) = iet_('en3S');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['A',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['WB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3D2';
            ielftype(ie) = iet_('en3D2');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['PA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3D2';
            ielftype(ie) = iet_('en3D2');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['B',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['PB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3S';
            ielftype(ie) = iet_('en3S');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['DD',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDSQ';
            ielftype(ie) = iet_('eDSQ');
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['CA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3PR';
            ielftype(ie) = iet_('en3PR');
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['A',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['B',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['CB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en3PR';
            ielftype(ie) = iet_('en3PR');
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('NI-1')
            ig = ig_(['EW',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/NI');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-Z/NI');
            ig = ig_(['EP',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/NI');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-Z/NI');
            ig = ig_(['EA',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['DD',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['EB',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['DD',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['EC',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['DD',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-ZZ');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 200.0;
%    Solution
% LO SOLTN               ????
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MY-V-V';
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 0.02;
        varargout{1} = pbm;

    case 'en3S'

        EV_  = varargin{1};
        iel_ = varargin{2};
        WSSMV4 = pbm.efpar(1)-EV_(4);
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*WSSMV4;
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*WSSMV4;
            g_(2,1) = EV_(1)*EV_(3)*WSSMV4;
            g_(3,1) = EV_(1)*EV_(2)*WSSMV4;
            g_(4,1) = -EV_(1)*EV_(2)*EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = EV_(3)*WSSMV4;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*WSSMV4;
                H_(3,1) = H_(1,3);
                H_(1,4) = -EV_(2)*EV_(3);
                H_(4,1) = H_(1,4);
                H_(2,3) = EV_(1)*WSSMV4;
                H_(3,2) = H_(2,3);
                H_(2,4) = -EV_(1)*EV_(3);
                H_(4,2) = H_(2,4);
                H_(3,4) = -EV_(1)*EV_(2);
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    case 'en3D2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(4,5);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(3,3) = U_(3,3)+1;
        U_(4,4) = U_(4,4)+1;
        U_(4,5) = U_(4,5)-2;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        IV_(4) = U_(4,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*IV_(3)*IV_(4);
        if(nargout>1)
            g_(1,1) = IV_(2)*IV_(3)*IV_(4);
            g_(2,1) = IV_(1)*IV_(3)*IV_(4);
            g_(3,1) = IV_(1)*IV_(2)*IV_(4);
            g_(4,1) = IV_(1)*IV_(2)*IV_(3);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = IV_(3)*IV_(4);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2)*IV_(4);
                H_(3,1) = H_(1,3);
                H_(1,4) = IV_(2)*IV_(3);
                H_(4,1) = H_(1,4);
                H_(2,3) = IV_(1)*IV_(4);
                H_(3,2) = H_(2,3);
                H_(2,4) = IV_(1)*IV_(3);
                H_(4,2) = H_(2,4);
                H_(3,4) = IV_(1)*IV_(2);
                H_(4,3) = H_(3,4);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eDSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+2.000000e-01;
        U_(1,2) = U_(1,2)+2.000000e-01;
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

    case 'en3PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
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

