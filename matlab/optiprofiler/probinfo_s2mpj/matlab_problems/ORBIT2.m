function varargout = ORBIT2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%        A reformulation of the discretized optimal control problem
%        ORBIT
%        Consider minimising the time needed for a spacecraft to
%        move from one circular orbit around the earth to another.
% 
%        This problem can be put in the following form:
%        let (x1,x2,x3) be the position and (x4,x5,x6) the velocity
%        let (u1,u2,u3) be the driving force vector
%        let q be the time required, then our problem becomes:
%        MINIMISE        q
%        with the equations of motion:
%                        dx1/dt = hv*q*x4
%                        dx2/dt = hv*q*x5
%                        dx3/dt = hv*q*x6
%                        dx4/dt = hv*q*(hg*x1/r^3-hf*u(1)/m)     (1)
%                        dx5/dt = hv*q*(hg*x1/r^3-hf*u(2)/m)             
%                        dx6/dt = hv*q*(hg*x1/r^3-hf*u(3)/m)
%        
%        with    m = m0-hm*q*t   ( - the mass variation)
%                        
%                r = sqrt( x1^2 + x2^2 + x3^2 )  -  dist. from the center
%                                                   of the earth
%        't' is a rescaled time varying between 0 and 1, and
%        'hv,hf,hg,hm' are scaling constants.
%        the driving force is bounded by 
%                        u1^2 + u2^2 + u3^2 <= 1                 (2)
%        (the rather arbitrary no. '1' representing the max. power
%        of the spacecraft).
%        We choose the initial conditions:
%                x1 = x2 = 0   ,   x3 = 1   -  initial position
%                x5 = x6 = 0   ,   x4 = Vorb    -   corresponding orbital
%                                                   speed
%        and the final conditions are:
%        x1^2 + x2^2 + x3^2 = Rf^2  -   final orbit radius
%        x4^2 + x5^2 + x6^2 = Vforb^2 -  corresponding orbital
%                                                 speed
%        x1*x4  + x2*x5 + x3*x6 = 0  -  direction must be parallel
%                                        to the velocity
%        we have chosen the constants hg,hf,hv,hm so that the x,u vectors
%        are of order one. These correspond to an initial orbit at 150 km
%        above the earth's surface, and a final orbit at 250 km above the 
%        earth's surface.
%        The reduction to an NLP is done in that same way as for CAR.SIF
%        The time taken should be:
%                        q = 315 secs
% 
%    SIF input: Thomas Felici, Universite de Nancy, France, 
%               October 1993, with modifications to compute derivatives
%               by Nick Gould
%      SAVEs removed December 3rd 2014
% 
%    classification = 'LOR1-RN-V-V'
% 
% *******************************************
% 
%        M = Number of time nodes.
% 
%        Change this for different resolution
% 
% IE M                              3   $-PARAMETER n= 25, m= 18
% IE M                             10   $-PARAMETER n= 88, m= 67
% IE M                             30   $-PARAMETER n=268, m=207 original value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORBIT2';

switch(action)

    case 'setup'

    pb.name      = 'ORBIT2';
    pb.sifpbname = 'ORBIT2';
    pbm.name     = 'ORBIT2';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('M') = 300;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
        v_('NEQ') = 6;
        v_('NCON') = 1;
        v_('NLIM') = 3;
        v_('NLINEQ') = 0;
        v_('NLINCON') = 0;
        v_('NLINLIM') = 0;
        v_('NTEQ') = v_('NEQ')+v_('NLINEQ');
        v_('NTCON') = v_('NCON')+v_('NLINCON');
        v_('NTLIM') = v_('NLIM')+v_('NLINLIM');
        v_('NX') = 6;
        v_('NU') = 3;
        v_('NQ') = 1;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('M-1') = v_('M')-v_('1');
        v_('RM') = v_('M-1');
        for I=v_('1'):v_('M')
            v_('I-1') = I-v_('1');
            v_('RI1') = v_('I-1');
            v_(['T',int2str(I)]) = v_('RI1')/v_('RM');
        end
        v_('RT') = 6371.0e0;
        v_('RT*RT') = v_('RT')*v_('RT');
        v_('MU') = 9.81e-3*v_('RT*RT');
        v_('ONE') = 1.0e+0;
        v_('PI/4') = atan(v_('ONE'));
        v_('PI') = 4.0e+0*v_('PI/4');
        v_('R0') = 1.5e+2+v_('RT');
        v_('R0*R0') = v_('R0')*v_('R0');
        v_('MU/R0') = v_('MU')/v_('R0');
        v_('VORB') = sqrt(v_('MU/R0'));
        v_('RF') = 2.5e+2+v_('RT');
        v_('M0') = 3.0e+3;
        v_('QT') = 3.333e+0;
        v_('VG') = 3.0e+0;
        v_('TS') = 1.0e+0;
        v_('HG') = v_('VORB')*v_('R0*R0');
        v_('HG') = v_('MU')/v_('HG');
        v_('HG') = v_('TS')*v_('HG');
        v_('HV') = v_('VORB')/v_('R0');
        v_('HV') = v_('TS')*v_('HV');
        v_('HF') = v_('VORB')*v_('M0');
        v_('HF') = v_('VG')/v_('HF');
        v_('HF') = v_('QT')*v_('HF');
        v_('HF') = v_('TS')*v_('HF');
        v_('HM') = v_('QT')/v_('M0');
        v_('HM') = v_('TS')*v_('HM');
        v_('VF') = v_('MU')/v_('RF');
        v_('VF') = sqrt(v_('VF'));
        v_('VF') = v_('VF')/v_('VORB');
        v_('VFVF') = v_('VF')*v_('VF');
        v_('RF') = v_('RF')/v_('R0');
        v_('RFRF') = v_('RF')*v_('RF');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('NX')
                [iv,ix_] = s2xlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('M-1')
            for J=v_('1'):v_('NU')
                [iv,ix_] = s2xlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
            end
        end
        for J=v_('1'):v_('NQ')
            [iv,ix_] = s2xlib('ii',['Q',int2str(J)],ix_);
            pb.xnames{iv} = ['Q',int2str(J)];
            xscale(iv,1) = 1.0e+2;
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_(['Q',int2str(v_('1'))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('M-1')
            for J=v_('1'):v_('NTEQ')
                [ig,ig_] = s2xlib('ii',['K',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['K',int2str(I),',',int2str(J)];
            end
        end
        for J=v_('1'):v_('NTLIM')
            [ig,ig_] = s2xlib('ii',['L',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['L',int2str(J)];
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('NTCON')
                [ig,ig_] = s2xlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['G',int2str(I),',',int2str(J)];
                pbm.gscale(ig,1) = 1.0e+2;
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
        pbm.gconst(ig_(['L',int2str(v_('1'))])) = v_('RFRF');
        pbm.gconst(ig_(['L',int2str(v_('3'))])) = v_('VFVF');
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['G',int2str(I),',',int2str(v_('1'))])) = 1.0e+0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_(['Q',int2str(v_('1'))]),1) = .10000E+03;
        pb.xupper(ix_(['Q',int2str(v_('1'))])) = .40000E+26;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('1'))]),1) = .00000E+00;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('1'))]),1) = .00000E+00;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('2'))]),1) = .00000E+00;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('2'))]),1) = .00000E+00;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('3'))]),1) = .10000E+01;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('3'))]),1) = .10000E+01;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('4'))]),1) = .10000E+01;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('4'))]),1) = .10000E+01;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('5'))]),1) = .00000E+00;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('5'))]),1) = .00000E+00;
        pb.xlower(ix_(['X',int2str(v_('1')),',',int2str(v_('6'))]),1) = .00000E+00;
        pb.xupper(ix_(['X',int2str(v_('1')),',',int2str(v_('6'))]),1) = .00000E+00;
        for I=v_('2'):v_('M')
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('1'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('1'))])) = .10000E+26;
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('2'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('2'))])) = .10000E+26;
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('3'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('3'))])) = .10000E+26;
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('4'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('4'))])) = .10000E+26;
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('5'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('5'))])) = .10000E+26;
            pb.xlower(ix_(['X',int2str(I),',',int2str(v_('6'))]),1) = -.10000E+26;
            pb.xupper(ix_(['X',int2str(I),',',int2str(v_('6'))])) = .10000E+26;
        end
        for I=v_('1'):v_('M-1')
            pb.xlower(ix_(['U',int2str(I),',',int2str(v_('1'))]),1) = -.10000E+01;
            pb.xupper(ix_(['U',int2str(I),',',int2str(v_('1'))])) = .10000E+01;
            pb.xlower(ix_(['U',int2str(I),',',int2str(v_('2'))]),1) = -.10000E+01;
            pb.xupper(ix_(['U',int2str(I),',',int2str(v_('2'))])) = .10000E+01;
            pb.xlower(ix_(['U',int2str(I),',',int2str(v_('3'))]),1) = -.10000E+01;
            pb.xupper(ix_(['U',int2str(I),',',int2str(v_('3'))])) = .10000E+01;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,['Q',int2str(v_('1'))]))
            pb.x0(ix_(['Q',int2str(v_('1'))]),1) = .10000E+03;
        else
            pb.y0(find(pbm.congrps==ig_(['Q',int2str(v_('1'))])),1) = .10000E+03;
        end
        for I=v_('1'):v_('M')
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('1'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('1'))]),1) = .00000E+00;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('1'))])),1) =...
                      .00000E+00;
            end
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('2'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('2'))]),1) = .00000E+00;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('2'))])),1) =...
                      .00000E+00;
            end
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('3'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('3'))]),1) = .10000E+01;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('3'))])),1) =...
                      .10000E+01;
            end
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('4'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('4'))]),1) = .10000E+01;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('4'))])),1) =...
                      .10000E+01;
            end
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('5'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('5'))]),1) = .00000E+00;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('5'))])),1) =...
                      .00000E+00;
            end
            if(isKey(ix_,['X',int2str(I),',',int2str(v_('6'))]))
                pb.x0(ix_(['X',int2str(I),',',int2str(v_('6'))]),1) = .00000E+00;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I),',',int2str(v_('6'))])),1) =...
                      .00000E+00;
            end
        end
        for I=v_('1'):v_('M-1')
            if(isKey(ix_,['U',int2str(I),',',int2str(v_('1'))]))
                pb.x0(ix_(['U',int2str(I),',',int2str(v_('1'))]),1) = .10000E+01;
            else
                pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(v_('1'))])),1) =...
                      .10000E+01;
            end
            if(isKey(ix_,['U',int2str(I),',',int2str(v_('2'))]))
                pb.x0(ix_(['U',int2str(I),',',int2str(v_('2'))]),1) = .10000E+01;
            else
                pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(v_('2'))])),1) =...
                      .10000E+01;
            end
            if(isKey(ix_,['U',int2str(I),',',int2str(v_('3'))]))
                pb.x0(ix_(['U',int2str(I),',',int2str(v_('3'))]),1) = .10000E+01;
            else
                pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(v_('3'))])),1) =...
                      .10000E+01;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'SQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2xlib( 'ii', 'PROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2xlib( 'ii', 'Ktyp',iet_);
        elftv{it}{1} = 'XI1';
        elftv{it}{2} = 'XF1';
        elftv{it}{3} = 'XI2';
        elftv{it}{4} = 'XF2';
        elftv{it}{5} = 'XI3';
        elftv{it}{6} = 'XF3';
        elftv{it}{7} = 'XI4';
        elftv{it}{8} = 'XF4';
        elftv{it}{9} = 'XI5';
        elftv{it}{10} = 'XF5';
        elftv{it}{11} = 'XI6';
        elftv{it}{12} = 'XF6';
        elftv{it}{13} = 'UV1';
        elftv{it}{14} = 'UV2';
        elftv{it}{15} = 'UV3';
        elftv{it}{16} = 'QV1';
        elftp{it}{1} = 'TN';
        elftp{it}{2} = 'TN1';
        elftp{it}{3} = 'ID';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M-1')
            v_('S') = 1+I;
            for J=v_('1'):v_('NTEQ')
                v_('ReJ') = J;
                ename = ['Ke',int2str(I),',',int2str(J)];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'Ktyp';
                ielftype(ie) = iet_('Ktyp');
                [~,posep] = ismember('TN',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
                [~,posep] = ismember('TN1',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['T',int2str(v_('S'))]);
                [~,posep] = ismember('ID',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('ReJ');
                vname = ['X',int2str(I),',',int2str(v_('1'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('1'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(v_('2'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('2'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(v_('3'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('3'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(v_('4'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('4'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(v_('5'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI5',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('5'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF5',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(v_('6'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XI6',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(v_('S')),',',int2str(v_('6'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('XF6',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(v_('1'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('UV1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(v_('2'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('UV2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['U',int2str(I),',',int2str(v_('3'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('UV3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Q',int2str(v_('1'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('QV1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for J=v_('1'):v_('NTCON')
            for I=v_('1'):v_('M-1')
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('1'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'SQR';
                ielftype(ie) = iet_('SQR');
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('1'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                vname = ['U',int2str(I),',',int2str(v_('1'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('2'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'SQR';
                ielftype(ie) = iet_('SQR');
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('2'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                vname = ['U',int2str(I),',',int2str(v_('2'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('3'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                pbm.elftype{ie} = 'SQR';
                ielftype(ie) = iet_('SQR');
                ename = ['Ge',int2str(I),',',int2str(J),',',int2str(v_('3'))];
                [ie,ie_] = s2xlib('ii',ename,ie_);
                vname = ['U',int2str(I),',',int2str(v_('3'))];
                [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for J=v_('1'):v_('NTCON')
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('1'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQR';
            ielftype(ie) = iet_('SQR');
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('1'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['U',int2str(v_('M-1')),',',int2str(v_('1'))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('2'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQR';
            ielftype(ie) = iet_('SQR');
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('2'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['U',int2str(v_('M-1')),',',int2str(v_('2'))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('3'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'SQR';
            ielftype(ie) = iet_('SQR');
            ename = ['Ge',int2str(v_('M')),',',int2str(J),',',int2str(v_('3'))];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            vname = ['U',int2str(v_('M-1')),',',int2str(v_('3'))];
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['Le',int2str(v_('1')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('1')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('1'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('1')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('1')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('2'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('1')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('1')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('3'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PROD';
        ielftype(ie) = iet_('PROD');
        ename = ['Le',int2str(v_('2')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('1'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('4'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PROD';
        ielftype(ie) = iet_('PROD');
        ename = ['Le',int2str(v_('2')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('2'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('5'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PROD';
        ielftype(ie) = iet_('PROD');
        ename = ['Le',int2str(v_('2')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('3'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('2')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('6'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('3')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('3')),',',int2str(v_('1'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('4'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('3')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('3')),',',int2str(v_('2'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('5'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Le',int2str(v_('3')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'SQR';
        ielftype(ie) = iet_('SQR');
        ename = ['Le',int2str(v_('3')),',',int2str(v_('3'))];
        [ie,ie_] = s2xlib('ii',ename,ie_);
        vname = ['X',int2str(v_('M')),',',int2str(v_('6'))];
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M-1')
            for J=v_('1'):v_('NTEQ')
                ig = ig_(['K',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Ke',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1;
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('NTCON')
                ig = ig_(['G',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['Ge',int2str(I),',',int2str(J),',',int2str(v_('1'))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1;
                posel = posel+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['Ge',int2str(I),',',int2str(J),',',int2str(v_('2'))]);
                pbm.grelw{ig}(posel) = 1;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['Ge',int2str(I),',',int2str(J),',',int2str(v_('3'))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1;
            end
        end
        for I=v_('1'):v_('NTLIM')
            ig = ig_(['L',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Le',int2str(I),',',int2str(v_('1'))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Le',int2str(I),',',int2str(v_('2'))]);
            pbm.grelw{ig}(posel) = 1;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Le',int2str(I),',',int2str(v_('3'))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        sA2 = size(pbm.A,2);
        for j = 1:min([sA2,pb.n,length(xscale)])
            if ( xscale(j) ~= 0.0 && xscale(j) ~= 1.0 )
                for i = find(pbm.A(:,j))
                      pbm.A(i,j) = pbm.A(i,j)/xscale(j);
                end
            end
        end
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR1-RN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
