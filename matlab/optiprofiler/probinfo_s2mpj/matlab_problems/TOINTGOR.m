function varargout = TOINTGOR(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TOINTGOR
%    *********
% 
%    Toint's  Operations Research problem
% 
%    Source:
%    Ph.L. Toint,
%    "Some numerical results using a sparse matrix updating formula in
%    unconstrained optimization",
%    Mathematics of Computation 32(1):839-852, 1978.
% 
%    See also Buckley#55 (p.94) (With a slightly lower optimal value?)
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-COUR2-MN-50-0'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TOINTGOR';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 50;
        v_('ALPH1') = 1.25;
        v_('ALPH2') = 1.40;
        v_('ALPH3') = 2.40;
        v_('ALPH4') = 1.40;
        v_('ALPH5') = 1.75;
        v_('ALPH6') = 1.20;
        v_('ALPH7') = 2.25;
        v_('ALPH8') = 1.20;
        v_('ALPH9') = 1.00;
        v_('ALPH10') = 1.10;
        v_('ALPH11') = 1.50;
        v_('ALPH12') = 1.60;
        v_('ALPH13') = 1.25;
        v_('ALPH14') = 1.25;
        v_('ALPH15') = 1.20;
        v_('ALPH16') = 1.20;
        v_('ALPH17') = 1.40;
        v_('ALPH18') = 0.50;
        v_('ALPH19') = 0.50;
        v_('ALPH20') = 1.25;
        v_('ALPH21') = 1.80;
        v_('ALPH22') = 0.75;
        v_('ALPH23') = 1.25;
        v_('ALPH24') = 1.40;
        v_('ALPH25') = 1.60;
        v_('ALPH26') = 2.00;
        v_('ALPH27') = 1.00;
        v_('ALPH28') = 1.60;
        v_('ALPH29') = 1.25;
        v_('ALPH30') = 2.75;
        v_('ALPH31') = 1.25;
        v_('ALPH32') = 1.25;
        v_('ALPH33') = 1.25;
        v_('ALPH34') = 3.00;
        v_('ALPH35') = 1.50;
        v_('ALPH36') = 2.00;
        v_('ALPH37') = 1.25;
        v_('ALPH38') = 1.40;
        v_('ALPH39') = 1.80;
        v_('ALPH40') = 1.50;
        v_('ALPH41') = 2.20;
        v_('ALPH42') = 1.40;
        v_('ALPH43') = 1.50;
        v_('ALPH44') = 1.25;
        v_('ALPH45') = 2.00;
        v_('ALPH46') = 1.50;
        v_('ALPH47') = 1.25;
        v_('ALPH48') = 1.40;
        v_('ALPH49') = 0.60;
        v_('ALPH50') = 1.50;
        v_('BETA1') = 1.0;
        v_('BETA2') = 1.5;
        v_('BETA3') = 1.0;
        v_('BETA4') = 0.1;
        v_('BETA5') = 1.5;
        v_('BETA6') = 2.0;
        v_('BETA7') = 1.0;
        v_('BETA8') = 1.5;
        v_('BETA9') = 3.0;
        v_('BETA10') = 2.0;
        v_('BETA11') = 1.0;
        v_('BETA12') = 3.0;
        v_('BETA13') = 0.1;
        v_('BETA14') = 1.5;
        v_('BETA15') = 0.15;
        v_('BETA16') = 2.0;
        v_('BETA17') = 1.0;
        v_('BETA18') = 0.1;
        v_('BETA19') = 3.0;
        v_('BETA20') = 0.1;
        v_('BETA21') = 1.2;
        v_('BETA22') = 1.0;
        v_('BETA23') = 0.1;
        v_('BETA24') = 2.0;
        v_('BETA25') = 1.2;
        v_('BETA26') = 3.0;
        v_('BETA27') = 1.5;
        v_('BETA28') = 3.0;
        v_('BETA29') = 2.0;
        v_('BETA30') = 1.0;
        v_('BETA31') = 1.2;
        v_('BETA32') = 2.0;
        v_('BETA33') = 1.0;
        v_('D1') = -5.0;
        v_('D2') = -5.0;
        v_('D3') = -5.0;
        v_('D4') = -2.5;
        v_('D5') = -6.0;
        v_('D6') = -6.0;
        v_('D7') = -5.0;
        v_('D8') = -6.0;
        v_('D9') = -10.0;
        v_('D10') = -6.0;
        v_('D11') = -5.0;
        v_('D12') = -9.0;
        v_('D13') = -2.0;
        v_('D14') = -7.0;
        v_('D15') = -2.5;
        v_('D16') = -6.0;
        v_('D17') = -5.0;
        v_('D18') = -2.0;
        v_('D19') = -9.0;
        v_('D20') = -2.0;
        v_('D21') = -5.0;
        v_('D22') = -5.0;
        v_('D23') = -2.5;
        v_('D24') = -5.0;
        v_('D25') = -6.0;
        v_('D26') = -10.0;
        v_('D27') = -7.0;
        v_('D28') = -10.0;
        v_('D29') = -6.0;
        v_('D30') = -5.0;
        v_('D31') = -4.0;
        v_('D32') = -4.0;
        v_('D33') = -4.0;
        v_('1') = 1;
        v_('33') = 33;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['GA',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
            v_('SCALE') = 1.0/v_(['ALPH',int2str(I)]);
            pbm.gscale(ig,1) = v_('SCALE');
        end
        [ig,ig_] = s2mpjlib('ii','GB1',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X31');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB2',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB3',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB4',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB5',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB6',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB7',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB8',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB9',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB10',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB11',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB12',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X21');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','GB13',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X22');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X23');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X24');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB14',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X23');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X25');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X26');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB15',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X25');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X27');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X28');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB16',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X28');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X29');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X30');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB17',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X29');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X31');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X32');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB18',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X32');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X33');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X34');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB19',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X33');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X35');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB20',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X35');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X21');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X36');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB21',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X36');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X37');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X38');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB22',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X30');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X37');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X39');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB23',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X38');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X39');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X40');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB24',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X40');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X41');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X42');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB25',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X41');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X43');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X44');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X50');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB26',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X44');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X45');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X46');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X47');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB27',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X46');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X48');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB28',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X42');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X45');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X48');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X50');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X49');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','GB29',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X26');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X34');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X43');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','GB30',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X24');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X47');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','GB31',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X49');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','GB32',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X22');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','GB33',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X27');
        valA(end+1) = -1.0;
        for I=v_('1'):v_('33')
            v_('SCALE') = 1.0/v_(['BETA',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['GB',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('SCALE');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('33')
            pbm.gconst(ig_(['GB',int2str(I)])) = v_(['D',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gACT',igt_);
        [it,igt_] = s2mpjlib('ii','gBBT',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['GA',int2str(I)]);
            pbm.grftype{ig} = 'gACT';
        end
        for I=v_('1'):v_('33')
            ig = ig_(['GB',int2str(I)]);
            pbm.grftype{ig} = 'gBBT';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN             1373.90546067
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COUR2-MN-50-0';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'g_globs'

        pbm = varargin{1};
        pbm.gfpar(1) = 1.0e0;    % this is ONE
        pbm.gfpar(2) = 0.0e0;    % this is ZERO
        varargout{1} = pbm;

    case 'gACT'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        AT = abs(GVAR_);
        AT1 = AT+pbm.gfpar(1);
        LAT = log(AT1);
        AA = AT/AT1;
        AG = AA+LAT;
        varargout{1} = AT*LAT;
        if(nargout>1)
            g_ = sign(GVAR_)*(AG);
            varargout{2} = g_;
            if(nargout>2)
                H_ = (2.0-AA)/AT1;
                varargout{3} = H_;
            end
        end

    case 'gBBT'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        AT = abs(GVAR_);
        AT1 = AT+pbm.gfpar(1);
        LAT = log(AT1);
        AA = AT/AT1;
        TPOS = max(sign(GVAR_)*(pbm.gfpar(1)),pbm.gfpar(2));
        TNEG = pbm.gfpar(1)-TPOS;
        varargout{1} = GVAR_*GVAR_*(TNEG+TPOS*LAT);
        if(nargout>1)
            g_ = TNEG*2.0*GVAR_+TPOS*GVAR_*(AA+2.0*LAT);
            varargout{2} = g_;
            if(nargout>2)
                H_ = TNEG*2.0+TPOS*(AA*(4.0-AA)+2.0*LAT);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,2];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

