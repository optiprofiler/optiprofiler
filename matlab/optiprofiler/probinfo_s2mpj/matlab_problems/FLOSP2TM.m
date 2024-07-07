function varargout = FLOSP2TM(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : FLOSP2TM
%    *********
% 
%    A  two-dimensional base  flow  problem in an inclined enclosure.
% 
%    Temperature constant y = +/- 1 boundary conditions
%    Middling Reynold's number
% 
%    The flow is considered in a square of length 2,  centered on the
%    origin and aligned with the x-y axes. The square is divided into
%    4 n ** 2  sub-squares,  each of  length 1 / n.  The differential
%    equation is replaced by  discrete nonlinear equations at each of 
%    the grid points. 
% 
%    The differential equation relates the vorticity, temperature and
%    a stream function.
%    
%    Source: 
%    J. N. Shadid
%    "Experimental and computational study of the stability
%    of Natural convection flow in an inclined enclosure",
%    Ph. D. Thesis, University of Minnesota, 1989,
%    problem SP2 (pp.128-130), 
% 
%    SIF input: Nick Gould, August 1993.
% 
%    classification = 'NQR2-MY-V-V'
% 
%    Half the number of discretization intervals
%    Number of variables = 3(2M+1)**2 
% 
%       Alternative values for the SIF file parameters:
% IE M                   1              $-PARAMETER n=27
% IE M                   2              $-PARAMETER n=75
% IE M                   5              $-PARAMETER n=363     original value
% IE M                   8              $-PARAMETER n=867
% IE M                   10             $-PARAMETER n=1323
% IE M                   15             $-PARAMETER n=2883
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'FLOSP2TM';

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
            v_('M') = 1;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
        if(nargs<2)
            v_('RA') = 1.0e+5;  %  SIF file default value
        else
            v_('RA') = varargin{2};
        end
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('AX') = 1.0;
        v_('THETA') = 0.5*v_('PI');
        v_('A1') = 0.0;
        v_('A2') = 1.0;
        v_('A3') = 0.0;
        v_('B1') = 0.0;
        v_('B2') = 1.0;
        v_('B3') = 1.0;
        v_('F1') = 1.0;
        v_('F2') = 0.0;
        v_('F3') = 0.0;
        v_('G1') = 1.0;
        v_('G2') = 0.0;
        v_('G3') = 0.0;
        v_('M-1') = -1+v_('M');
        v_('-M') = -1*v_('M');
        v_('-M+1') = -1*v_('M-1');
        v_('1/H') = v_('M');
        v_('-1/H') = -1.0*v_('1/H');
        v_('2/H') = 2.0*v_('1/H');
        v_('-2/H') = -2.0*v_('1/H');
        v_('H') = 1.0/v_('1/H');
        v_('H2') = v_('H')*v_('H');
        v_('1/H2') = v_('1/H')*v_('1/H');
        v_('-2/H2') = -2.0*v_('1/H2');
        v_('1/2H') = 0.5*v_('1/H');
        v_('-1/2H') = -0.5*v_('1/H');
        v_('AXX') = v_('AX')*v_('AX');
        v_('SINTHETA') = sin(v_('THETA'));
        v_('COSTHETA') = cos(v_('THETA'));
        v_('PI1') = v_('AX')*v_('RA');
        v_('PI1') = v_('PI1')*v_('COSTHETA');
        v_('PI1') = -0.5*v_('PI1');
        v_('-PI1') = -1.0*v_('PI1');
        v_('PI2') = v_('AXX')*v_('RA');
        v_('PI2') = v_('PI2')*v_('SINTHETA');
        v_('PI2') = 0.5*v_('PI2');
        v_('-PI2') = -1.0*v_('PI2');
        v_('2A1') = 2.0*v_('A1');
        v_('2B1') = 2.0*v_('B1');
        v_('2F1') = 2.0*v_('F1');
        v_('2G1') = 2.0*v_('G1');
        v_('2F1/AX') = v_('2F1')/v_('AX');
        v_('2G1/AX') = v_('2G1')/v_('AX');
        v_('AX/2') = 0.5*v_('AX');
        v_('AXX/2') = 0.5*v_('AXX');
        v_('AXX/4') = 0.25*v_('AXX');
        v_('2AX') = 2.0*v_('AX');
        v_('2AXX') = 2.0*v_('AXX');
        v_('2/AX') = 2.0/v_('AX');
        v_('2/AXH') = v_('2/H')/v_('AX');
        v_('-2/AXH') = -1.0*v_('2/AXH');
        v_('PI1/2H') = v_('PI1')*v_('1/2H');
        v_('-PI1/2H') = v_('PI1')*v_('-1/2H');
        v_('PI2/2H') = v_('PI2')*v_('1/2H');
        v_('-PI2/2H') = v_('PI2')*v_('-1/2H');
        v_('2A1/H') = v_('2A1')*v_('1/H');
        v_('-2A1/H') = v_('2A1')*v_('-1/H');
        v_('2B1/H') = v_('2B1')*v_('1/H');
        v_('-2B1/H') = v_('2B1')*v_('-1/H');
        v_('2F1/AXH') = v_('2F1/AX')*v_('1/H');
        v_('-2F1/AXH') = v_('2F1/AX')*v_('-1/H');
        v_('2G1/AXH') = v_('2G1/AX')*v_('1/H');
        v_('-2G1/AXH') = v_('2G1/AX')*v_('-1/H');
        v_('AX/H2') = v_('AX')*v_('1/H2');
        v_('-AX/H2') = -1.0*v_('AX/H2');
        v_('AX/4H2') = 0.25*v_('AX/H2');
        v_('-AX/4H2') = -0.25*v_('AX/H2');
        v_('AXX/H2') = v_('AXX')*v_('1/H2');
        v_('-2AXX/H2') = -2.0*v_('AXX/H2');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('-M'):v_('M')
            for I=v_('-M'):v_('M')
                [iv,ix_] = s2mpjlib('ii',['OM',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['OM',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['PH',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['PH',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['PS',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['PS',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('-M+1'):v_('M-1')
            v_('J+') = 1+J;
            v_('J-') = -1+J;
            for I=v_('-M+1'):v_('M-1')
                v_('I+') = 1+I;
                v_('I-') = -1+I;
                [ig,ig_] = s2mpjlib('ii',['S',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['S',int2str(I),',',int2str(J)];
                iv = ix_(['OM',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/H2');
                end
                iv = ix_(['OM',int2str(round(v_('I+'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['OM',int2str(round(v_('I-'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['OM',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2AXX/H2');
                end
                iv = ix_(['OM',int2str(I),',',int2str(round(v_('J+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
                iv = ix_(['OM',int2str(I),',',int2str(round(v_('J-')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
                iv = ix_(['PH',int2str(round(v_('I+'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-PI1/2H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-PI1/2H');
                end
                iv = ix_(['PH',int2str(round(v_('I-'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('PI1/2H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('PI1/2H');
                end
                iv = ix_(['PH',int2str(I),',',int2str(round(v_('J+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-PI2/2H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-PI2/2H');
                end
                iv = ix_(['PH',int2str(I),',',int2str(round(v_('J-')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('PI2/2H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('PI2/2H');
                end
                [ig,ig_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['V',int2str(I),',',int2str(J)];
                iv = ix_(['PS',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/H2');
                end
                iv = ix_(['PS',int2str(round(v_('I+'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['PS',int2str(round(v_('I-'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['PS',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2AXX/H2');
                end
                iv = ix_(['PS',int2str(I),',',int2str(round(v_('J+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
                iv = ix_(['PS',int2str(I),',',int2str(round(v_('J-')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
                iv = ix_(['OM',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/4')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/4');
                end
                [ig,ig_] = s2mpjlib('ii',['E',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['E',int2str(I),',',int2str(J)];
                iv = ix_(['PH',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/H2');
                end
                iv = ix_(['PH',int2str(round(v_('I+'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['PH',int2str(round(v_('I-'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/H2');
                end
                iv = ix_(['PH',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2AXX/H2');
                end
                iv = ix_(['PH',int2str(I),',',int2str(round(v_('J+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
                iv = ix_(['PH',int2str(I),',',int2str(round(v_('J-')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('AXX/H2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('AXX/H2');
                end
            end
        end
        for K=v_('-M'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2A1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2A1/H');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('M-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2A1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2A1/H');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('A2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('A2');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('-M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('-M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('-M+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2B1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2B1/H');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('-M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('-M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('-M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2B1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2B1/H');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(K),',',int2str(round(v_('-M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(K),',',int2str(round(v_('-M')))];
            iv = ix_(['PH',int2str(K),',',int2str(round(v_('-M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('B2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('B2');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2F1/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2F1/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('M-1'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2F1/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2F1/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('F2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('F2');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('-M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('-M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('-M+1'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2G1/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2G1/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('-M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('-M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('-M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2G1/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2G1/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['T',int2str(round(v_('-M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['T',int2str(round(v_('-M'))),',',int2str(K)];
            iv = ix_(['PH',int2str(round(v_('-M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('G2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('G2');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(K),',',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(K),',',int2str(round(v_('M')))];
            iv = ix_(['PS',int2str(K),',',int2str(round(v_('M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2/H');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(K),',',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(K),',',int2str(round(v_('M')))];
            iv = ix_(['PS',int2str(K),',',int2str(round(v_('M-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2/H');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(K),',',int2str(round(v_('-M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(K),',',int2str(round(v_('-M')))];
            iv = ix_(['PS',int2str(K),',',int2str(round(v_('-M+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2/H');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(K),',',int2str(round(v_('-M')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(K),',',int2str(round(v_('-M')))];
            iv = ix_(['PS',int2str(K),',',int2str(round(v_('-M')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2/H');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(round(v_('M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(round(v_('M'))),',',int2str(K)];
            iv = ix_(['PS',int2str(round(v_('M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(round(v_('M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(round(v_('M'))),',',int2str(K)];
            iv = ix_(['PS',int2str(round(v_('M-1'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(round(v_('-M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(round(v_('-M'))),',',int2str(K)];
            iv = ix_(['PS',int2str(round(v_('-M+1'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('2/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('2/AXH');
            end
            [ig,ig_] = s2mpjlib('ii',['V',int2str(round(v_('-M'))),',',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['V',int2str(round(v_('-M'))),',',int2str(K)];
            iv = ix_(['PS',int2str(round(v_('-M'))),',',int2str(K)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-2/AXH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-2/AXH');
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
        for K=v_('-M'):v_('M')
            pbm.gconst(ig_(['T',int2str(K),',',int2str(round(v_('M')))])) = v_('A3');
            pbm.gconst(ig_(['T',int2str(K),',',int2str(round(v_('-M')))])) = v_('B3');
            pbm.gconst(ig_(['T',int2str(round(v_('M'))),',',int2str(K)])) = v_('F3');
            pbm.gconst(ig_(['T',int2str(round(v_('-M'))),',',int2str(K)])) = v_('G3');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for K=v_('-M'):v_('M')
            pb.xlower(ix_(['PS',int2str(K),',',int2str(round(v_('-M')))]),1) = 1.0;
            pb.xupper(ix_(['PS',int2str(K),',',int2str(round(v_('-M')))]),1) = 1.0;
            pb.xlower(ix_(['PS',int2str(round(v_('-M'))),',',int2str(K)]),1) = 1.0;
            pb.xupper(ix_(['PS',int2str(round(v_('-M'))),',',int2str(K)]),1) = 1.0;
            pb.xlower(ix_(['PS',int2str(K),',',int2str(round(v_('M')))]),1) = 1.0;
            pb.xupper(ix_(['PS',int2str(K),',',int2str(round(v_('M')))]),1) = 1.0;
            pb.xlower(ix_(['PS',int2str(round(v_('M'))),',',int2str(K)]),1) = 1.0;
            pb.xupper(ix_(['PS',int2str(round(v_('M'))),',',int2str(K)]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'PSIM';
        elftv{it}{2} = 'PSIP';
        elftv{it}{3} = 'PHIM';
        elftv{it}{4} = 'PHIP';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('-M+1'):v_('M-1')
            v_('J+') = 1+J;
            v_('J-') = -1+J;
            for I=v_('-M+1'):v_('M-1')
                v_('I+') = 1+I;
                v_('I-') = -1+I;
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['PS',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PSIP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PS',int2str(I),',',int2str(round(v_('J-')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PSIM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PH',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PHIP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PH',int2str(round(v_('I-'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PHIM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['F',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePROD';
                ielftype(ie) = iet_('ePROD');
                vname = ['PS',int2str(round(v_('I+'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PSIP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PS',int2str(round(v_('I-'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PSIM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PH',int2str(I),',',int2str(round(v_('J+')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PHIP',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['PH',int2str(I),',',int2str(round(v_('J-')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('PHIM',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for J=v_('-M+1'):v_('M-1')
            for I=v_('-M+1'):v_('M-1')
                ig = ig_(['E',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-AX/4H2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['F',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('AX/4H2');
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NQR2-MY-V-V';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,4);
        U_(1,2) = U_(1,2)+1;
        U_(1,1) = U_(1,1)-1;
        U_(2,4) = U_(2,4)+1;
        U_(2,3) = U_(2,3)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
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

