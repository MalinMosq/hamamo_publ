% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

%% EXPERIMENT 1 & 2: Enhanced localisation. Convergence w.r.t. H, keeping l fixed


% Test case EXPERIMENT 1: 4 interfaces, Omega = [0,1]^2
% A0 = 1, A1 = 1, B1 = 1
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = 5; % (1/32)
% % titlestring = {'Error (FEM) vs H for different $\ell$ with enhanced basis','4 interfaces, h = 5 (1/32)'};
% energy_error_fem_enh_loc = [... % ell = 1       ell = 2         ell = 3         ell = 4         ell = 5
%    0.060254061058756   0.032326829415923   0.032377284884709  0.032377119022522 % H = 0
%    0.057246664635360   0.013465623282883   0.008328132953719  0.007987110844053 % H = 1
%    0.065944305132824   0.018299225253942   0.004963821130320  0.002415002852364 % H = 2
%    0.056708261649998   0.018675500362774   0.006255905090371  0.001973001423631 % H = 3
%         ]
h = 6;
titlestring = {'Error (FEM) vs H for different $\ell$','4 interfaces, h = 6 (1/64)','$A_i = 1$', '$A_j = 1$','$f_i = \sin(\pi x)\sin(\pi y)$, $f_j = x + 2y$'};
energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
   0.059842420607469   0.032383192571294   0.032433647860598   0.032433481153816 % H = 0
   0.056901293431663   0.013470390061628   0.008380031515889   0.008044422399445 % H = 1
   0.065162833314657   0.018070212548787   0.004896393971657   0.002434228247662 % H = 2
   0.054244565367396   0.017575441946848   0.005668379081384   0.001743181365697 % H = 3
   0.042774906725813   0.014424467297749   0.004885013943305   0.001583813608686 % H = 4
];
% % Non-stabilised method
% h = 6;
% titlestring = {'Error (FEM) vs H for different $\ell$','4 interfaces, h = 6 (1/64)','$A_i = 1$', '$A_j = 1$','$f_i = \sin(\pi x)\sin(\pi y)$, $f_j = x + 2y$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.272445199164623   0.038124179905599   0.032435357422754   0.032433481241416 % H = 0
%    0.662053577304605   0.183375046724699   0.034238559329295   0.009443470962587 % H = 1
%    0.808795809166761   0.455433272968067   0.130900957904597   0.029716480574561 % H = 2
%    0.871598634681145   0.693708032822455   0.288003240358335   0.080075897869762 % H = 3
%    0.893892890587483   0.831146876552423   0.545650424964277   0.193201314091882 % H = 4
% ]
% % Non-stabilised method
% % h = 5;
% energy_error_fem_enh_loc = [
%    0.272901613397024   0.038093427558349   0.032378981654272   0.032377119112561 % H = 0
%    0.663100623851102   0.184434538320044   0.034495228846799   0.009425110996793 % H = 1
%    0.810345334051655   0.462613158815768   0.135105827967941   0.031012472662278 % H = 2
%    0.874011741437058   0.714320822307931   0.321998745969311   0.094752894507525 % H = 3
%    ];

% Test case: 4 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(x)*sin(y)+2
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 5; % 1/32
% titlestring = {'Error (FEM) vs H for different $\ell$ with enhanced basis','4 interfaces, h = 5 (1/32)','$A_i$ pw const, random in $[0.1,1]$', '$A_j = \sin(x)*\sin(y)+2$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.042721497870477   0.023759024666190   0.023784545296760   0.023784268802542 % H = 0
%    0.038315674838965   0.009594406778065   0.005999526708050   0.005753386025663 % H = 1
%    0.046155078035831   0.013264547117893   0.003583297003502   0.001736465821372 % H = 2
%    0.039943488426333   0.013547713060959   0.004549805261015   0.001434264041882 % H = 3
%    ]
% h = 6; % 1/64
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.042376602522081   0.023767746758773   0.023792405239541   0.023792150272122 % H = 0
%    0.038304693087549   0.009596082489252   0.006029473240831   0.005786850427452 % H = 1
%    0.045635622761912   0.013087545470626   0.003532834953418   0.001748079753632 % H = 2
%    0.038066695625866   0.012693910048722   0.004114456837290   0.001265180777384 % H = 3
%    0.030278321580690   0.010385971200318   0.003523358039801   0.001144030521781 % H = 4
%    ]

% Test case: 4 interfaces, Omega = [0,1]^2
% A0 = 1, A1 = pw const random
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 5;
% titlestring = {'Error (FEM) vs H for different $\ell$ with enhanced basis','4 interfaces, h = 5 (1/32)','$A_i = 1$', '$A_j$ pw const, random in $[0.1,1]$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.111570876466790   0.063311485258688   0.063442917236906   0.063442939023135 % H = 0
%    0.118384079733076   0.020247892908863   0.016158890693326   0.015888556353343 % H = 1
%    0.141918195616736   0.037845252872512   0.008605098498574   0.004508148236132 % H = 2
%    0.234352849924444   0.075381609051871   0.019112638774879   0.006331418695641 % H = 3
%    ]


% Test case: 4 interfaces, Omega = [0,1]^2
% A0 = 1, A1 = sin(30pi x)*sin(30pi y)+2
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6; % VECTOR VERSION (same as function handle)
% titlestring = {'Error (FEM) vs H for different $\ell$ with enhanced basis','4 interfaces, h = 6 (1/64)','$A_i = 1$', '$A_j = \sin(30\pi x)*\sin(30\pi y)+2$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.044175060857907   0.024500311648570   0.024525656027651   0.024525518924115 % H = 0
%    0.042949929743188   0.010464915558149   0.006282684522438   0.005999248555966 % H = 1
%    0.049659838413790   0.014048287427466   0.003684494856480   0.001814370752233 % H = 2
%    0.041870562444875   0.014588905566341   0.004510592726523   0.001323867942310 % H = 3
%    0.039255352270228   0.013404513572805   0.004225900129345   0.001316868194793 % H = 4
%    ]

% Test case: 4 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(30pi x)*sin(30pi y)+2, B1 = 1
% f0 = 1, f1 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6;
% titlestring = {'Error (FEM) vs H for different $\ell$','4 interfaces, h = 6 (1/64)','$A_i =$ pw random', '$A_j = \sin(30\pi x) \sin(30\pi y) + 2$','$f_i = 1$, $f_j = 1$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.029728125150632   0.001253291821954   0.000004078442647   0.000000012904784 % H = 0
%    0.034876680154556   0.005643874749911   0.001020213188508   0.000206238098591 % H = 1
%    0.037144638303999   0.009620155418142   0.002238407585459   0.000724476534101 % H = 2
%    0.031428754126410   0.010104452278005   0.003060085156600   0.000856771491470 % H = 3
%    0.029726692632799   0.009394508680421   0.002893434161373   0.000892268966405 % H = 4
%    ]
% h = 5;
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.030181224218483   0.001256355547568   0.000004216758756   0.000000032476328 % H = 0
%    0.034518073804039   0.005663423583267   0.001025189935198   0.000207579931120 % H = 1
%    0.037626404023336   0.009754140294389   0.002292597259763   0.000751485456395 % H = 2
%    0.033775610114949   0.010849843467338   0.003412756699328   0.000996880746591 % H = 3
% ]

% % Test case EXPERIMENT 2: 4 interfaces, Omega = [0,1]^2
% % A0 = pw const random, A1 = sin(30pi x)*sin(30pi y)+2, B1 = 1
% % f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = 5;
% % energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
% %    0.045526568678546   0.025657703913366   0.025687747338252   0.025687461636259 % H = 0
% %    0.042598430202008   0.010786471900048   0.006499356524997   0.006203780066151 % H = 1
% %    0.051017276620897   0.014657818567069   0.003862813921012   0.001868806864947 % H = 2
% %    0.045097611263711   0.016051568830745   0.005168102407379   0.001560960042530 % H = 3
% %    ]
% h = 6;
% titlestring = {'Error (FEM) vs H for different $\ell$ with enhanced basis','4 interfaces, h = 6 (1/64)','$A_i =$ p.w. random in $[.01,1]$ on fine grid', '$A_j = \sin(30\pi x)*\sin(30\pi y)+2$'};
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.045335547395428   0.025683444007209   0.025712709054594   0.025712443151336 % H = 0
%    0.042767106292551   0.010794335600730   0.006531002214135   0.006242434650653 % H = 1
%    0.050279626830323   0.014441572198531   0.003802192779816   0.001879541291834 % H = 2
%    0.042627873241417   0.015005383615950   0.004642232883428   0.001363050747599 % H = 3
%    0.040583937600623   0.013828634550154   0.004350696801523   0.001355202781147 % H = 4
%    ];



% Test case: 8 interfaces, Omega = [0,1]^2
% A0 = 1, A1 = 1, B1 = 1
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6;
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.091448897488995   0.017365314125374   0.014354090310365   0.014466748645432 % H = 0
%    0.080488297754164   0.021742095921485   0.006337435653437   0.003709355470378 % H = 1
%    0.092589738273061   0.024939412260124   0.006819532483270   0.002070736751460 % H = 2
%    0.072276389244285   0.024265461767220   0.007559892565467   0.002071161322879 % H = 3
%    0.055752840088301   0.019297539919030   0.006394396795526   0.002038401784879 % H = 4
% ]

% Test case: 8 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(30pi x)*sin(30pi y) + 1.1, B1 = 1
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6;
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.119675713064805   0.019685695048591   0.016731585940802   0.016859975793654 % H = 0
%    0.118827860813359   0.034403581087306   0.007428848976001   0.004362967385655 % H = 1
%    0.162726192536235   0.048937681439517   0.011037953781432   0.003667643299811 % H = 2
%    0.154131990721052   0.046139456494564   0.012740118625685   0.003958762517010 % H = 3
%    0.114032935548316   0.041563327128211   0.011696365386961   0.003425259033119 % H = 4
% ]

% Test case: 8 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(30pi x)*sin(30pi y) + 1.1, B1 = 1
% f = pw const rand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6;
% energy_error_fem_enh_loc = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.037840636936484   0.002918470725004   0.000489522673571   0.000037523543471 % H = 0
%    0.048487984591454   0.011634021648081   0.002248630657336   0.000556926108602 % H = 1
%    0.058314384917203   0.016938645582902   0.003707821700482   0.001155981776026 % H = 2
%    0.054370731660914   0.016428805196591   0.004885908703894   0.001436080546087 % H = 3
%    0.041396654679542   0.015181178518958   0.004158736908687   0.001201355931801 % H = 4
% ]



figure('Color','white')
set(gcf,'Position',[50,0,700,700])
H = 1./2.^(0:(size(energy_error_fem_enh_loc,1)-1));
h_fig = loglog(H,energy_error_fem_enh_loc,'-o');%;ell6;ell7]','-o') %ell1;ell2;ell3; % 1;ell2;ell3;ell4;ell5
col = get(h_fig,'color');
hold on
loglog(H,.03*H.^2,'--','Color',col{4})
% loglog(H,.03*H,'--','Color',col{3})
ylim([1e-3 1])
legendstring = {[repmat('$\ell = $',size(energy_error_fem_enh_loc,2),1), num2str((1:size(energy_error_fem_enh_loc,2))'); '     $H^2$']};
legend(legendstring,'interpreter','latex','Location','southeast');
xlabel({'$H$'},'Interpreter','latex','fontsize',16)
ylabel('$\|u_h - \tilde{u}^\ell_{H,h}\|_a$','Interpreter','latex','fontsize',16)
xticks(flip(H))
xticklabels(flip({'1','1/2','1/4','1/8','1/16'}))

% title(titlestring,'Interpreter','latex')
% matlab2tikz('figures/nonstab_Ai1_Aj1_Bj1_fisin_fjx2y.tex')

% Image = getframe(gcf);
% exportgraphics(gcf, 'figures/Ai1_Aj1_Bj1_fisin_fjx2y.png');
% Aipwrand_Aj30sin_Bj1_fi1_fj1
% Ai1_Aj1_Bj1_fisin_fjx2y
% Aipwrand_Aj30sin_Bj1_fisin_fjx2y.png




%% EXPERIMENT 3:  Non-convex coarse triangles (fine grid interface)
% A0 = 1
e_fem_lod = [0.643639937180619   0.152282335964427   0.024014935419862   0.003119456671311   1.0e-03*0.391544298090456   1.0e-03*0.044931126430044];

% A0 pw random
e_fem_lod = [   0.524951357613358   0.126788150645751   0.019701416608526   0.002568974026859   0.000307441986633   0.000037131984985];

ell = 1:length(e_fem_lod);

% subplot(1,2,1)
figure('Color','white')
set(gcf,'Position',[50,50,700,700])
h = semilogy(ell, e_fem_lod,'-o'); hold on
col = get(h,'color');
%semilogy(ell, 50*1./10.^ell,'--','color',col)
semilogy(ell, exp(-2*ell),'--','color',col)
% semilogy(ell, 50*1./10.^(ell/2),'--')
% semilogy(ell, 50*1./10.^(2*ell),'--')
xticks(ell)
xlabel('$\ell$','Interpreter','latex','fontsize',16)
ylabel('$\|u_h - \tilde{u}_{H,h}^\ell\|_a$','fontsize',16,'Interpreter','Latex')

legend('$\|u_h - \tilde{u}_{H,h}^\ell\|_a$','$\exp(-2\ell)$','Interpreter','Latex')

% matlab2tikz('figures/nonconvex_domain.tex')
% subplot(1,2,2)
% plot(ell, e_fem_lod,'-o')
% xticks(ell)
% xlabel('$\ell$','Interpreter','latex')
% ylabel('Energy error FEM - LOD')




%% EXPERIMENT 4: Intricate example with different smoothness on RHS. ell = log2(1/H)



% Test case: 8 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(30pi x)*sin(30pi y) + 1.1, B1 = 1
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 6;
b_smooth = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
   0.119675713064805   0.019685695048591   0.016731585940802   0.016859975793654
   0.118827860813359   0.034403581087306   0.007428848976001   0.004362967385655
   0.162726192536235   0.048937681439517   0.011037953781432   0.003667643299811
   0.154131990721052   0.046139456494564   0.012740118625685   0.003958762517010
   0.114032935548316   0.041563327128211   0.011696365386961   0.003425259033119
]

% Test case: 8 interfaces, Omega = [0,1]^2
% A0 = pw const random, A1 = sin(30pi x)*sin(30pi y) + 1.1, B1 = 1
% f0 = f1 = sin(30pi x)sin(30pi y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = 6;
b_rapid = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
   0.029818135066014   0.028999108537120   0.028969307101192  0.028971329359575 % H = 0
   0.020761954519147   0.019097082258440   0.019022273641885  0.019056957146597 % H = 1
   0.013946674058934   0.012063336804951   0.012103307556698  0.012086404439087 % H = 2
                   0                   0   0.007258597699250  0                 % H = 3
                   0                   0                   0   0.002652453302586 % H = 4
]

% % Test case: 8 interfaces, Omega = [0,1]^2
% % A0 = pw const random, A1 = sin(30pi x)*sin(30pi y) + 1.1, B1 = 1
% % f0 = f1 = pw const random
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = 6;
% b_rapid = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
%    0.037840636936484   0.002918470725004   0.000489522673571   0.000037523543471 % H = 0
%    0.048487984591454   0.011634021648081   0.002248630657336   0.000556926108602 % H = 1
%    0.058314384917203   0.016938645582902   0.003707821700482   0.001155981776026 % H = 2
%    0.054370731660914   0.016428805196591   0.004885908703894   0.001436080546087 % H = 3
%    0.041396654679542   0.015181178518958   0.004158736908687   0.001201355931801] % H = 4
% Gives too good result. Unclear why, did something go wrong?

% Smooth RHS
temp = diag(b_smooth(2:5,1:4));
H = 1./2.^(0:(size(temp,1)-1));
loglog(H,temp,'-*'), hold on

% Rough RHS
b_pw_const = diag(b_rapid(2:5,1:4));
loglog(H(1:4),b_pw_const,'-*')

% Reference lines
loglog(H,.3*H.^2,'--')
loglog(H,.015*H,'--')

ylim([1e-3 1])
legendstring = {'$\ell = \log_2 \frac{1}{H}$',...
    '$f = \sin(30\pi x)\sin(\pi y)$','$H^2$','$H$'};
legend(legendstring,'interpreter','latex','Location','southeast');
xlabel({'$H$'},'Interpreter','latex','fontsize',16)
ylabel('$\|u_h - \tilde{u}^\ell_{H,h}\|_a$','Interpreter','latex','fontsize',16)
xticks(flip(H))
xticklabels(flip({'1','1/2','1/4','1/8','1/16'}))
% title(titlestring,'Interpreter','latex')
% matlab2tikz('figures/convergence_new_complicated_H_H2_1.1_h6.tex')






















%% Localisation experiment

% f_i = 0, f_j = 1
energy_error_lod10_lod = [0.400454599074059   0.123310577478648   0.023980381522457   0.003937358531811 0.000604332615483   0.000078270151687   0.000000102157214];
energy_error_lod10_enh_loc = [0.361366187653110   0.100956716758286   0.018141369501705   0.003150457717176 0.000460003518814   0.000060693702238   0.000000060068501];


ell = 1:7;
semilogy(ell,energy_error_lod10_lod,'-o'), hold on %;ell6;ell7]','-o') %ell1;ell2;ell3; % 1;ell2;ell3;ell4;ell5
semilogy(ell,energy_error_lod10_enh_loc,'-o')
legend('LOD', 'enh_loc basis','interpreter','latex') %'$\ell = 1$', '$\ell = 2$', '$\ell = 3$',  '$\ell = 6$', '$\ell = 7$',
xlabel({'$\ell$'},'Interpreter','latex','fontsize',16)
ylabel('$\|u_h - \tilde{u}^\ell_{H,h}\|_a$','Interpreter','latex','fontsize',16)
% xticks(flip(H))
% xticklabels(flip({'1','2','3','4','5'}))
title('Error (FEM) vs H for different $\ell$ with enh_loc basis (maybe faulty)','Interpreter','latex')









%% Error FEM - LOD
clc

% 4 interfaces, h = 6, Omega = [0,1]^2
% A0 = rand, A1 = sin(30pi x)*sin(30pi y)+2, B1 = 1
% f0 = sin(pi*x).*sin(pi*y), f1 = x + 2*y
titlestring = {'Error (FEM) vs H for different $\ell$ with enh_loc basis','4 interfaces, h = 6 (1/64)','$A_i =$ p.w. random in $[.01,1]$ on fine grid', '$A_j = \sin(30\pi x)*\sin(30\pi y)+2$'};
energy_error_lod = [... % ell = 1         ell = 2             ell = 3             ell = 4         ell = 5
   0.045335547395428   0.025683444007209   0.025712709054594   0.025712443151336 % H = 0
   0.042767106292551   0.010794335600730   0.006531002214135   0.006242434650653 % H = 1
   0.050279626830323   0.014441572198531   0.003802192779816   0.001879541291834 % H = 2
   0.042627873241417   0.015005383615950   0.004642232883428   0.001363050747599 % H = 3
   0.040583937600623   0.013828634550154   0.004350696801523   0.001355202781147 % H = 4
   ];


energy_error_fem = [0.1 %false
    0.1 %false
    0.091387725080591  
    0.162018695715850   
    0.080067462667971
    %0.004972757991356
    ];

figure('Color','white')
set(gcf,'Position',[50,0,700,700])
H = 1./2.^(0:(size(energy_error_lod,1)-1));
h_fig = loglog(H,energy_error_lod,'-o');%;ell6;ell7]','-o') %ell1;ell2;ell3; % 1;ell2;ell3;ell4;ell5
col = get(h_fig,'color');
hold on
h_fig = loglog(H,energy_error_fem,'-o');

loglog(H,.03*H.^2,'--','Color',col{4})
loglog(H,.03*H,'--','Color',col{3})
ylim([1e-4 .5])
legendstring = {[repmat('$\ell = $',size(energy_error_lod,2),1), num2str((1:size(energy_error_lod,2))'); '     $H^2$']};
legend(legendstring,'interpreter','latex','Location','southeast');
xlabel({'H','(as proportion of subdomain edges)'},'Interpreter','latex','fontsize',16)
ylabel('Energy norm error','Interpreter','latex','fontsize',16)
xticks(flip(H))
xticklabels(flip({'1','1/2','1/4','1/8','1/16'}))