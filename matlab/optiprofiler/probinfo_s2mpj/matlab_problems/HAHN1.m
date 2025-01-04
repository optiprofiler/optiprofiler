function varargout = HAHN1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HAHN1
%    *********
% 
%    NIST Data fitting problem HAHN1 given as an inconsistent set of
%    nonlinear equations.
% 
%    Fit: y = (b1 + b2*x + b3*x**2 + b4*x**3) / 
%             (1 + b5*x + b6*x**2 + b7*x**3) + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Hahn, T., NIST (197?). 
%      Copper Thermal Expansion Study.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'C-CNOR2-MN-7-37'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HAHN1';

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
        v_('M') = 236;
        v_('N') = 7;
        v_('1') = 1;
        v_('X1') = 24.41;
        v_('X2') = 34.82;
        v_('X3') = 44.09;
        v_('X4') = 45.07;
        v_('X5') = 54.98;
        v_('X6') = 65.51;
        v_('X7') = 70.53;
        v_('X8') = 75.70;
        v_('X9') = 89.57;
        v_('X10') = 91.14;
        v_('X11') = 96.40;
        v_('X12') = 97.19;
        v_('X13') = 114.26;
        v_('X14') = 120.25;
        v_('X15') = 127.08;
        v_('X16') = 133.55;
        v_('X17') = 133.61;
        v_('X18') = 158.67;
        v_('X19') = 172.74;
        v_('X20') = 171.31;
        v_('X21') = 202.14;
        v_('X22') = 220.55;
        v_('X23') = 221.05;
        v_('X24') = 221.39;
        v_('X25') = 250.99;
        v_('X26') = 268.99;
        v_('X27') = 271.80;
        v_('X28') = 271.97;
        v_('X29') = 321.31;
        v_('X30') = 321.69;
        v_('X31') = 330.14;
        v_('X32') = 333.03;
        v_('X33') = 333.47;
        v_('X34') = 340.77;
        v_('X35') = 345.65;
        v_('X36') = 373.11;
        v_('X37') = 373.79;
        v_('X38') = 411.82;
        v_('X39') = 419.51;
        v_('X40') = 421.59;
        v_('X41') = 422.02;
        v_('X42') = 422.47;
        v_('X43') = 422.61;
        v_('X44') = 441.75;
        v_('X45') = 447.41;
        v_('X46') = 448.7;
        v_('X47') = 472.89;
        v_('X48') = 476.69;
        v_('X49') = 522.47;
        v_('X50') = 522.62;
        v_('X51') = 524.43;
        v_('X52') = 546.75;
        v_('X53') = 549.53;
        v_('X54') = 575.29;
        v_('X55') = 576.00;
        v_('X56') = 625.55;
        v_('X57') = 20.15;
        v_('X58') = 28.78;
        v_('X59') = 29.57;
        v_('X60') = 37.41;
        v_('X61') = 39.12;
        v_('X62') = 50.24;
        v_('X63') = 61.38;
        v_('X64') = 66.25;
        v_('X65') = 73.42;
        v_('X66') = 95.52;
        v_('X67') = 107.32;
        v_('X68') = 122.04;
        v_('X69') = 134.03;
        v_('X70') = 163.19;
        v_('X71') = 163.48;
        v_('X72') = 175.70;
        v_('X73') = 179.86;
        v_('X74') = 211.27;
        v_('X75') = 217.78;
        v_('X76') = 219.14;
        v_('X77') = 262.52;
        v_('X78') = 268.01;
        v_('X79') = 268.62;
        v_('X80') = 336.25;
        v_('X81') = 337.23;
        v_('X82') = 339.33;
        v_('X83') = 427.38;
        v_('X84') = 428.58;
        v_('X85') = 432.68;
        v_('X86') = 528.99;
        v_('X87') = 531.08;
        v_('X88') = 628.34;
        v_('X89') = 253.24;
        v_('X90') = 273.13;
        v_('X91') = 273.66;
        v_('X92') = 282.10;
        v_('X93') = 346.62;
        v_('X94') = 347.19;
        v_('X95') = 348.78;
        v_('X96') = 351.18;
        v_('X97') = 450.10;
        v_('X98') = 450.35;
        v_('X99') = 451.92;
        v_('X100') = 455.56;
        v_('X101') = 552.22;
        v_('X102') = 553.56;
        v_('X103') = 555.74;
        v_('X104') = 652.59;
        v_('X105') = 656.20;
        v_('X106') = 14.13;
        v_('X107') = 20.41;
        v_('X108') = 31.30;
        v_('X109') = 33.84;
        v_('X110') = 39.70;
        v_('X111') = 48.83;
        v_('X112') = 54.50;
        v_('X113') = 60.41;
        v_('X114') = 72.77;
        v_('X115') = 75.25;
        v_('X116') = 86.84;
        v_('X117') = 94.88;
        v_('X118') = 96.40;
        v_('X119') = 117.37;
        v_('X120') = 139.08;
        v_('X121') = 147.73;
        v_('X122') = 158.63;
        v_('X123') = 161.84;
        v_('X124') = 192.11;
        v_('X125') = 206.76;
        v_('X126') = 209.07;
        v_('X127') = 213.32;
        v_('X128') = 226.44;
        v_('X129') = 237.12;
        v_('X130') = 330.90;
        v_('X131') = 358.72;
        v_('X132') = 370.77;
        v_('X133') = 372.72;
        v_('X134') = 396.24;
        v_('X135') = 416.59;
        v_('X136') = 484.02;
        v_('X137') = 495.47;
        v_('X138') = 514.78;
        v_('X139') = 515.65;
        v_('X140') = 519.47;
        v_('X141') = 544.47;
        v_('X142') = 560.11;
        v_('X143') = 620.77;
        v_('X144') = 18.97;
        v_('X145') = 28.93;
        v_('X146') = 33.91;
        v_('X147') = 40.03;
        v_('X148') = 44.66;
        v_('X149') = 49.87;
        v_('X150') = 55.16;
        v_('X151') = 60.90;
        v_('X152') = 72.08;
        v_('X153') = 85.15;
        v_('X154') = 97.06;
        v_('X155') = 119.63;
        v_('X156') = 133.27;
        v_('X157') = 143.84;
        v_('X158') = 161.91;
        v_('X159') = 180.67;
        v_('X160') = 198.44;
        v_('X161') = 226.86;
        v_('X162') = 229.65;
        v_('X163') = 258.27;
        v_('X164') = 273.77;
        v_('X165') = 339.15;
        v_('X166') = 350.13;
        v_('X167') = 362.75;
        v_('X168') = 371.03;
        v_('X169') = 393.32;
        v_('X170') = 448.53;
        v_('X171') = 473.78;
        v_('X172') = 511.12;
        v_('X173') = 524.70;
        v_('X174') = 548.75;
        v_('X175') = 551.64;
        v_('X176') = 574.02;
        v_('X177') = 623.86;
        v_('X178') = 21.46;
        v_('X179') = 24.33;
        v_('X180') = 33.43;
        v_('X181') = 39.22;
        v_('X182') = 44.18;
        v_('X183') = 55.02;
        v_('X184') = 94.33;
        v_('X185') = 96.44;
        v_('X186') = 118.82;
        v_('X187') = 128.48;
        v_('X188') = 141.94;
        v_('X189') = 156.92;
        v_('X190') = 171.65;
        v_('X191') = 190.00;
        v_('X192') = 223.26;
        v_('X193') = 223.88;
        v_('X194') = 231.50;
        v_('X195') = 265.05;
        v_('X196') = 269.44;
        v_('X197') = 271.78;
        v_('X198') = 273.46;
        v_('X199') = 334.61;
        v_('X200') = 339.79;
        v_('X201') = 349.52;
        v_('X202') = 358.18;
        v_('X203') = 377.98;
        v_('X204') = 394.77;
        v_('X205') = 429.66;
        v_('X206') = 468.22;
        v_('X207') = 487.27;
        v_('X208') = 519.54;
        v_('X209') = 523.03;
        v_('X210') = 612.99;
        v_('X211') = 638.59;
        v_('X212') = 641.36;
        v_('X213') = 622.05;
        v_('X214') = 631.50;
        v_('X215') = 663.97;
        v_('X216') = 646.9;
        v_('X217') = 748.29;
        v_('X218') = 749.21;
        v_('X219') = 750.14;
        v_('X220') = 647.04;
        v_('X221') = 646.89;
        v_('X222') = 746.9;
        v_('X223') = 748.43;
        v_('X224') = 747.35;
        v_('X225') = 749.27;
        v_('X226') = 647.61;
        v_('X227') = 747.78;
        v_('X228') = 750.51;
        v_('X229') = 851.37;
        v_('X230') = 845.97;
        v_('X231') = 847.54;
        v_('X232') = 849.93;
        v_('X233') = 851.61;
        v_('X234') = 849.75;
        v_('X235') = 850.98;
        v_('X236') = 848.23;
        v_('Y1') = 0.591;
        v_('Y2') = 1.547;
        v_('Y3') = 2.902;
        v_('Y4') = 2.894;
        v_('Y5') = 4.703;
        v_('Y6') = 6.307;
        v_('Y7') = 7.03;
        v_('Y8') = 7.898;
        v_('Y9') = 9.470;
        v_('Y10') = 9.484;
        v_('Y11') = 10.072;
        v_('Y12') = 10.163;
        v_('Y13') = 11.615;
        v_('Y14') = 12.005;
        v_('Y15') = 12.478;
        v_('Y16') = 12.982;
        v_('Y17') = 12.970;
        v_('Y18') = 13.926;
        v_('Y19') = 14.452;
        v_('Y20') = 14.404;
        v_('Y21') = 15.190;
        v_('Y22') = 15.550;
        v_('Y23') = 15.528;
        v_('Y24') = 15.499;
        v_('Y25') = 16.131;
        v_('Y26') = 16.438;
        v_('Y27') = 16.387;
        v_('Y28') = 16.549;
        v_('Y29') = 16.872;
        v_('Y30') = 16.830;
        v_('Y31') = 16.926;
        v_('Y32') = 16.907;
        v_('Y33') = 16.966;
        v_('Y34') = 17.060;
        v_('Y35') = 17.122;
        v_('Y36') = 17.311;
        v_('Y37') = 17.355;
        v_('Y38') = 17.668;
        v_('Y39') = 17.767;
        v_('Y40') = 17.803;
        v_('Y41') = 17.765;
        v_('Y42') = 17.768;
        v_('Y43') = 17.736;
        v_('Y44') = 17.858;
        v_('Y45') = 17.877;
        v_('Y46') = 17.912;
        v_('Y47') = 18.046;
        v_('Y48') = 18.085;
        v_('Y49') = 18.291;
        v_('Y50') = 18.357;
        v_('Y51') = 18.426;
        v_('Y52') = 18.584;
        v_('Y53') = 18.610;
        v_('Y54') = 18.870;
        v_('Y55') = 18.795;
        v_('Y56') = 19.111;
        v_('Y57') = 0.367;
        v_('Y58') = 0.796;
        v_('Y59') = 0.892;
        v_('Y60') = 1.903;
        v_('Y61') = 2.150;
        v_('Y62') = 3.697;
        v_('Y63') = 5.870;
        v_('Y64') = 6.421;
        v_('Y65') = 7.422;
        v_('Y66') = 9.944;
        v_('Y67') = 11.023;
        v_('Y68') = 11.87;
        v_('Y69') = 12.786;
        v_('Y70') = 14.067;
        v_('Y71') = 13.974;
        v_('Y72') = 14.462;
        v_('Y73') = 14.464;
        v_('Y74') = 15.381;
        v_('Y75') = 15.483;
        v_('Y76') = 15.59;
        v_('Y77') = 16.075;
        v_('Y78') = 16.347;
        v_('Y79') = 16.181;
        v_('Y80') = 16.915;
        v_('Y81') = 17.003;
        v_('Y82') = 16.978;
        v_('Y83') = 17.756;
        v_('Y84') = 17.808;
        v_('Y85') = 17.868;
        v_('Y86') = 18.481;
        v_('Y87') = 18.486;
        v_('Y88') = 19.090;
        v_('Y89') = 16.062;
        v_('Y90') = 16.337;
        v_('Y91') = 16.345;
        v_('Y92') = 16.388;
        v_('Y93') = 17.159;
        v_('Y94') = 17.116;
        v_('Y95') = 17.164;
        v_('Y96') = 17.123;
        v_('Y97') = 17.979;
        v_('Y98') = 17.974;
        v_('Y99') = 18.007;
        v_('Y100') = 17.993;
        v_('Y101') = 18.523;
        v_('Y102') = 18.669;
        v_('Y103') = 18.617;
        v_('Y104') = 19.371;
        v_('Y105') = 19.330;
        v_('Y106') = 0.080;
        v_('Y107') = 0.248;
        v_('Y108') = 1.089;
        v_('Y109') = 1.418;
        v_('Y110') = 2.278;
        v_('Y111') = 3.624;
        v_('Y112') = 4.574;
        v_('Y113') = 5.556;
        v_('Y114') = 7.267;
        v_('Y115') = 7.695;
        v_('Y116') = 9.136;
        v_('Y117') = 9.959;
        v_('Y118') = 9.957;
        v_('Y119') = 11.600;
        v_('Y120') = 13.138;
        v_('Y121') = 13.564;
        v_('Y122') = 13.871;
        v_('Y123') = 13.994;
        v_('Y124') = 14.947;
        v_('Y125') = 15.473;
        v_('Y126') = 15.379;
        v_('Y127') = 15.455;
        v_('Y128') = 15.908;
        v_('Y129') = 16.114;
        v_('Y130') = 17.071;
        v_('Y131') = 17.135;
        v_('Y132') = 17.282;
        v_('Y133') = 17.368;
        v_('Y134') = 17.483;
        v_('Y135') = 17.764;
        v_('Y136') = 18.185;
        v_('Y137') = 18.271;
        v_('Y138') = 18.236;
        v_('Y139') = 18.237;
        v_('Y140') = 18.523;
        v_('Y141') = 18.627;
        v_('Y142') = 18.665;
        v_('Y143') = 19.086;
        v_('Y144') = 0.214;
        v_('Y145') = 0.943;
        v_('Y146') = 1.429;
        v_('Y147') = 2.241;
        v_('Y148') = 2.951;
        v_('Y149') = 3.782;
        v_('Y150') = 4.757;
        v_('Y151') = 5.602;
        v_('Y152') = 7.169;
        v_('Y153') = 8.920;
        v_('Y154') = 10.055;
        v_('Y155') = 12.035;
        v_('Y156') = 12.861;
        v_('Y157') = 13.436;
        v_('Y158') = 14.167;
        v_('Y159') = 14.755;
        v_('Y160') = 15.168;
        v_('Y161') = 15.651;
        v_('Y162') = 15.746;
        v_('Y163') = 16.216;
        v_('Y164') = 16.445;
        v_('Y165') = 16.965;
        v_('Y166') = 17.121;
        v_('Y167') = 17.206;
        v_('Y168') = 17.250;
        v_('Y169') = 17.339;
        v_('Y170') = 17.793;
        v_('Y171') = 18.123;
        v_('Y172') = 18.49;
        v_('Y173') = 18.566;
        v_('Y174') = 18.645;
        v_('Y175') = 18.706;
        v_('Y176') = 18.924;
        v_('Y177') = 19.1;
        v_('Y178') = 0.375;
        v_('Y179') = 0.471;
        v_('Y180') = 1.504;
        v_('Y181') = 2.204;
        v_('Y182') = 2.813;
        v_('Y183') = 4.765;
        v_('Y184') = 9.835;
        v_('Y185') = 10.040;
        v_('Y186') = 11.946;
        v_('Y187') = 12.596;
        v_('Y188') = 13.303;
        v_('Y189') = 13.922;
        v_('Y190') = 14.440;
        v_('Y191') = 14.951;
        v_('Y192') = 15.627;
        v_('Y193') = 15.639;
        v_('Y194') = 15.814;
        v_('Y195') = 16.315;
        v_('Y196') = 16.334;
        v_('Y197') = 16.430;
        v_('Y198') = 16.423;
        v_('Y199') = 17.024;
        v_('Y200') = 17.009;
        v_('Y201') = 17.165;
        v_('Y202') = 17.134;
        v_('Y203') = 17.349;
        v_('Y204') = 17.576;
        v_('Y205') = 17.848;
        v_('Y206') = 18.090;
        v_('Y207') = 18.276;
        v_('Y208') = 18.404;
        v_('Y209') = 18.519;
        v_('Y210') = 19.133;
        v_('Y211') = 19.074;
        v_('Y212') = 19.239;
        v_('Y213') = 19.280;
        v_('Y214') = 19.101;
        v_('Y215') = 19.398;
        v_('Y216') = 19.252;
        v_('Y217') = 19.89;
        v_('Y218') = 20.007;
        v_('Y219') = 19.929;
        v_('Y220') = 19.268;
        v_('Y221') = 19.324;
        v_('Y222') = 20.049;
        v_('Y223') = 20.107;
        v_('Y224') = 20.062;
        v_('Y225') = 20.065;
        v_('Y226') = 19.286;
        v_('Y227') = 19.972;
        v_('Y228') = 20.088;
        v_('Y229') = 20.743;
        v_('Y230') = 20.83;
        v_('Y231') = 20.935;
        v_('Y232') = 21.035;
        v_('Y233') = 20.93;
        v_('Y234') = 21.074;
        v_('Y235') = 21.085;
        v_('Y236') = 20.935;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['F',int2str(I)];
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'B1'))
            pb.x0(ix_('B1'),1) = 10.0;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 10.0;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = -1.0;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = -1.0;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = 0.05;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = 0.05;
        end
        if(isKey(ix_,'B4'))
            pb.x0(ix_('B4'),1) = -0.00001;
        else
            pb.y0(find(pbm.congrps==ig_('B4')),1) = -0.00001;
        end
        if(isKey(ix_,'B5'))
            pb.x0(ix_('B5'),1) = -0.05;
        else
            pb.y0(find(pbm.congrps==ig_('B5')),1) = -0.05;
        end
        if(isKey(ix_,'B6'))
            pb.x0(ix_('B6'),1) = 0.001;
        else
            pb.y0(find(pbm.congrps==ig_('B6')),1) = 0.001;
        end
        if(isKey(ix_,'B7'))
            pb.x0(ix_('B7'),1) = -0.000001;
        else
            pb.y0(find(pbm.congrps==ig_('B7')),1) = -0.000001;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE19',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE19';
            ielftype(ie) = iet_('eE19');
            vname = 'B1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-7-37';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eE19'

        EV_  = varargin{1};
        iel_ = varargin{2};
        X2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        X3 = X2*pbm.elpar{iel_}(1);
        X4 = X3*pbm.elpar{iel_}(1);
        X5 = X4*pbm.elpar{iel_}(1);
        X6 = X5*pbm.elpar{iel_}(1);
        T = EV_(1)+EV_(2)*pbm.elpar{iel_}(1)+EV_(3)*X2+EV_(4)*X3;
        D = 1.0e0+EV_(5)*pbm.elpar{iel_}(1)+EV_(6)*X2+EV_(7)*X3;
        D2 = D*D;
        TD3 = 0.5e0*D2*D;
        varargout{1} = T/D;
        if(nargout>1)
            g_(1,1) = 1.0e0/D;
            g_(2,1) = pbm.elpar{iel_}(1)/D;
            g_(3,1) = X2/D;
            g_(4,1) = X3/D;
            g_(5,1) = -pbm.elpar{iel_}(1)*T/D2;
            g_(6,1) = -X2*T/D2;
            g_(7,1) = -X3*T/D2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(7,7);
                H_(1,5) = -pbm.elpar{iel_}(1)/D2;
                H_(5,1) = H_(1,5);
                H_(1,6) = -X2/D2;
                H_(6,1) = H_(1,6);
                H_(1,7) = -X3/D2;
                H_(7,1) = H_(1,7);
                H_(2,5) = -X2/D2;
                H_(5,2) = H_(2,5);
                H_(2,6) = -X3/D2;
                H_(6,2) = H_(2,6);
                H_(2,7) = -X4/D2;
                H_(7,2) = H_(2,7);
                H_(3,5) = -X3/D2;
                H_(5,3) = H_(3,5);
                H_(3,6) = -X4/D2;
                H_(6,3) = H_(3,6);
                H_(3,7) = -X5/D2;
                H_(7,3) = H_(3,7);
                H_(4,5) = -X4/D2;
                H_(5,4) = H_(4,5);
                H_(4,6) = -X5/D2;
                H_(6,4) = H_(4,6);
                H_(4,7) = -X6/D2;
                H_(7,4) = H_(4,7);
                H_(5,5) = X2*T/TD3;
                H_(5,6) = X3*T/TD3;
                H_(6,5) = H_(5,6);
                H_(5,7) = X4*T/TD3;
                H_(7,5) = H_(5,7);
                H_(6,6) = X4*T/TD3;
                H_(6,7) = X5*T/TD3;
                H_(7,6) = H_(6,7);
                H_(7,7) = X6*T/TD3;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
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

