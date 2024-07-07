function varargout = NASH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NASH
%    *********
% 
%    A quadratic programming reformulation of a linear
%    complementarity problem arising from Nash equilibrium
%    provided by Michael Ferris
% 
%    classification = 'QLR2-AN-72-24'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NASH';

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
        v_('1') = 1;
        v_('N') = 72;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(J)],ix_);
            pb.xnames{iv} = ['X',int2str(J)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1000.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1000.0;
        end
        iv = ix_('X57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 500.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 500.0;
        end
        iv = ix_('X58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1000.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1000.0;
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        iv = ix_('X25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        iv = ix_('X26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
        iv = ix_('X27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C4';
        iv = ix_('X28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C5';
        iv = ix_('X29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C6';
        iv = ix_('X30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C7';
        iv = ix_('X31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02309+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02309;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.288626;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.263887;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.447486;
        end
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C8';
        iv = ix_('X32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C9';
        iv = ix_('X33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C10';
        iv = ix_('X34');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C11';
        iv = ix_('X35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C12';
        iv = ix_('X36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C13';
        iv = ix_('X37');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.33;
        end
        [ig,ig_] = s2mpjlib('ii','C14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C14';
        iv = ix_('X38');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.67;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33;
        end
        [ig,ig_] = s2mpjlib('ii','C15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C15';
        iv = ix_('X39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X63');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.67;
        end
        [ig,ig_] = s2mpjlib('ii','C16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C16';
        iv = ix_('X40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X64');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C17',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C17';
        iv = ix_('X41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C18';
        iv = ix_('X42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X66');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C19',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C19';
        iv = ix_('X43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X67');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.33;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.67;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.33;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C20';
        iv = ix_('X44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X68');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.33;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.33;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.67+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.67;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C21',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C21';
        iv = ix_('X45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X69');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C22',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C22';
        iv = ix_('X46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X70');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.288626+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.288626;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.892169+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.892169;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.298588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.298588;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.593581+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.593581;
        end
        [ig,ig_] = s2mpjlib('ii','C23',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C23';
        iv = ix_('X47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X71');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.263887+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.263887;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.298588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.298588;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.412719+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.412719;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.114131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.114131;
        end
        [ig,ig_] = s2mpjlib('ii','C24',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C24';
        iv = ix_('X48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X72');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.447486+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.447486;
        end
        iv = ix_('X21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.593581+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.593581;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.114131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.114131;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.707712+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.707712;
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
        pbm.gconst(ig_('C2')) = 35.100673;
        pbm.gconst(ig_('C3')) = 35.100673;
        pbm.gconst(ig_('C4')) = 35.100673;
        pbm.gconst(ig_('C5')) = 35.100673;
        pbm.gconst(ig_('C6')) = 35.100673;
        pbm.gconst(ig_('C7')) = 35.100673;
        pbm.gconst(ig_('C8')) = -15.0;
        pbm.gconst(ig_('C9')) = -15.0;
        pbm.gconst(ig_('C10')) = -20.0;
        pbm.gconst(ig_('C22')) = 61.241589;
        pbm.gconst(ig_('C23')) = -1.150548;
        pbm.gconst(ig_('C24')) = -60.091041;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.E+00*ones(pb.n,1);
        pb.xupper = 0.E+00*ones(pb.n,1);
        pb.xlower(ix_('X1')) = -Inf;
        pb.xupper(ix_('X1'),1) = +Inf;
        pb.xupper(ix_('X8')) = 1000.0;
        pb.xupper(ix_('X9')) = 500.0;
        pb.xupper(ix_('X10')) = 1000.0;
        pb.xlower(ix_('X11')) = -Inf;
        pb.xupper(ix_('X11'),1) = +Inf;
        pb.xlower(ix_('X12')) = -Inf;
        pb.xupper(ix_('X12'),1) = +Inf;
        pb.xlower(ix_('X13')) = -Inf;
        pb.xupper(ix_('X13'),1) = +Inf;
        pb.xlower(ix_('X14')) = -Inf;
        pb.xupper(ix_('X14'),1) = +Inf;
        pb.xlower(ix_('X15')) = -Inf;
        pb.xupper(ix_('X15'),1) = +Inf;
        pb.xlower(ix_('X16')) = -Inf;
        pb.xupper(ix_('X16'),1) = +Inf;
        pb.xlower(ix_('X17')) = -Inf;
        pb.xupper(ix_('X17'),1) = +Inf;
        pb.xlower(ix_('X18')) = -Inf;
        pb.xupper(ix_('X18'),1) = +Inf;
        pb.xlower(ix_('X19')) = -Inf;
        pb.xupper(ix_('X19'),1) = +Inf;
        pb.xlower(ix_('X20')) = -Inf;
        pb.xupper(ix_('X20'),1) = +Inf;
        pb.xlower(ix_('X21')) = -Inf;
        pb.xupper(ix_('X21'),1) = +Inf;
        pb.xlower(ix_('X22')) = -Inf;
        pb.xupper(ix_('X22'),1) = +Inf;
        pb.xlower(ix_('X23')) = -Inf;
        pb.xupper(ix_('X23'),1) = +Inf;
        pb.xlower(ix_('X24')) = -Inf;
        pb.xupper(ix_('X24'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = sparse( pb.n, pb.n );
        ix1 = ix_('X25');
        ix2 = ix_('X1');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X26');
        ix2 = ix_('X2');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X27');
        ix2 = ix_('X3');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X28');
        ix2 = ix_('X4');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X29');
        ix2 = ix_('X5');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X30');
        ix2 = ix_('X6');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X31');
        ix2 = ix_('X7');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X32');
        ix2 = ix_('X8');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X33');
        ix2 = ix_('X9');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X34');
        ix2 = ix_('X10');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X35');
        ix2 = ix_('X11');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X36');
        ix2 = ix_('X12');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X37');
        ix2 = ix_('X13');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X38');
        ix2 = ix_('X14');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X39');
        ix2 = ix_('X15');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X40');
        ix2 = ix_('X16');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X41');
        ix2 = ix_('X17');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X42');
        ix2 = ix_('X18');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X43');
        ix2 = ix_('X19');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X44');
        ix2 = ix_('X20');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X45');
        ix2 = ix_('X21');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X46');
        ix2 = ix_('X22');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X47');
        ix2 = ix_('X23');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X48');
        ix2 = ix_('X24');
        pbm.H(ix1,ix2) = 1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X49');
        ix2 = ix_('X1');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X50');
        ix2 = ix_('X2');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X51');
        ix2 = ix_('X3');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X52');
        ix2 = ix_('X4');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X53');
        ix2 = ix_('X5');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X54');
        ix2 = ix_('X6');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X55');
        ix2 = ix_('X7');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X56');
        ix2 = ix_('X8');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X57');
        ix2 = ix_('X9');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X58');
        ix2 = ix_('X10');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X59');
        ix2 = ix_('X11');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X60');
        ix2 = ix_('X12');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X61');
        ix2 = ix_('X13');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X62');
        ix2 = ix_('X14');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X63');
        ix2 = ix_('X15');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X64');
        ix2 = ix_('X16');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X65');
        ix2 = ix_('X17');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X66');
        ix2 = ix_('X18');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X67');
        ix2 = ix_('X19');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X68');
        ix2 = ix_('X20');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X69');
        ix2 = ix_('X21');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X70');
        ix2 = ix_('X22');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X71');
        ix2 = ix_('X23');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        ix1 = ix_('X72');
        ix2 = ix_('X24');
        pbm.H(ix1,ix2) = -1.0+pbm.H(ix1,ix2);
        pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'QLR2-AN-72-24';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

