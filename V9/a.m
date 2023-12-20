%function a(cal);

for cal=1:2;
% cal 0:help  1: 1 D   2:2 D  3:tests
switch(cal);

case 0;
help exemple_general_1D;help exemple1_1D;help exemple_1D_pertes;help exemple2_1D;help exemple3_1D;help exemple5_1D;help exemple6_1D;help exemple9_1D;help exemple10_1D;help exemple11_1D;help exemple_1D_metal;
help exemple_general_conique;help exemple_conique_pertes;help exemple1_conique;
help exemple_general_2D;help exemple_2D_pertes;help exemple1_2D;help exemple2_2D;help exemple3_2D;help exemple5_2D;help exemple6_2D;help exemple7_2D;help exemple9_2D;help exemple10_2D;help exemple11_2D;help exemple8;
%help test1D2D;help test1D;

case 1;
exemple_general_1D;retfig;
exemple_1D_pertes;retfig;
exemple1_1D;retfig;
exemple2_1D;retfig;
exemple3_1D;retfig;
exemple5_1D;retfig;
exemple6_1D;retfig;
exemple9_1D;retfig;
exemple10_1D;retfig;
exemple11_1D;retfig;
exemple_1D_metal;retfig;

exemple_general_conique;retfig;
exemple1_conique;retfig;

exemple_V9_0D_anisotrope
exemple_V9_1D_anisotrope
exemple_V9_conique_anisotrope

case 2;
exemple_general_2D;retfig;
exemple1_2D;retfig;
exemple2_2D;retfig;
exemple3_2D;retfig;
exemple5_2D;retfig;
exemple6_2D;retfig;
exemple7_2D;retfig;
exemple9_2D;retfig;
exemple10_2D;retfig;
exemple11_2D;retfig;
exemple_2D_pertes;retfig;

exemple8;retfig;
exemple_V9_2D_anisotrope

end

end;% cal