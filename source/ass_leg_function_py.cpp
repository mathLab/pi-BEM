
#include <math.h>

#include "../include/ass_leg_function.h"


AssLegFunction::AssLegFunction()
{
  leg_pointers.resize(21);
  for (unsigned int i = 0; i < 21; ++i)
    {
      leg_pointers[i].resize(i + 1);
    }

  leg_pointers[0][0]   = &AssLegFunction::P_0_0;
  leg_pointers[1][0]   = &AssLegFunction::P_1_0;
  leg_pointers[1][1]   = &AssLegFunction::P_1_1;
  leg_pointers[2][0]   = &AssLegFunction::P_2_0;
  leg_pointers[2][1]   = &AssLegFunction::P_2_1;
  leg_pointers[2][2]   = &AssLegFunction::P_2_2;
  leg_pointers[3][0]   = &AssLegFunction::P_3_0;
  leg_pointers[3][1]   = &AssLegFunction::P_3_1;
  leg_pointers[3][2]   = &AssLegFunction::P_3_2;
  leg_pointers[3][3]   = &AssLegFunction::P_3_3;
  leg_pointers[4][0]   = &AssLegFunction::P_4_0;
  leg_pointers[4][1]   = &AssLegFunction::P_4_1;
  leg_pointers[4][2]   = &AssLegFunction::P_4_2;
  leg_pointers[4][3]   = &AssLegFunction::P_4_3;
  leg_pointers[4][4]   = &AssLegFunction::P_4_4;
  leg_pointers[5][0]   = &AssLegFunction::P_5_0;
  leg_pointers[5][1]   = &AssLegFunction::P_5_1;
  leg_pointers[5][2]   = &AssLegFunction::P_5_2;
  leg_pointers[5][3]   = &AssLegFunction::P_5_3;
  leg_pointers[5][4]   = &AssLegFunction::P_5_4;
  leg_pointers[5][5]   = &AssLegFunction::P_5_5;
  leg_pointers[6][0]   = &AssLegFunction::P_6_0;
  leg_pointers[6][1]   = &AssLegFunction::P_6_1;
  leg_pointers[6][2]   = &AssLegFunction::P_6_2;
  leg_pointers[6][3]   = &AssLegFunction::P_6_3;
  leg_pointers[6][4]   = &AssLegFunction::P_6_4;
  leg_pointers[6][5]   = &AssLegFunction::P_6_5;
  leg_pointers[6][6]   = &AssLegFunction::P_6_6;
  leg_pointers[7][0]   = &AssLegFunction::P_7_0;
  leg_pointers[7][1]   = &AssLegFunction::P_7_1;
  leg_pointers[7][2]   = &AssLegFunction::P_7_2;
  leg_pointers[7][3]   = &AssLegFunction::P_7_3;
  leg_pointers[7][4]   = &AssLegFunction::P_7_4;
  leg_pointers[7][5]   = &AssLegFunction::P_7_5;
  leg_pointers[7][6]   = &AssLegFunction::P_7_6;
  leg_pointers[7][7]   = &AssLegFunction::P_7_7;
  leg_pointers[8][0]   = &AssLegFunction::P_8_0;
  leg_pointers[8][1]   = &AssLegFunction::P_8_1;
  leg_pointers[8][2]   = &AssLegFunction::P_8_2;
  leg_pointers[8][3]   = &AssLegFunction::P_8_3;
  leg_pointers[8][4]   = &AssLegFunction::P_8_4;
  leg_pointers[8][5]   = &AssLegFunction::P_8_5;
  leg_pointers[8][6]   = &AssLegFunction::P_8_6;
  leg_pointers[8][7]   = &AssLegFunction::P_8_7;
  leg_pointers[8][8]   = &AssLegFunction::P_8_8;
  leg_pointers[9][0]   = &AssLegFunction::P_9_0;
  leg_pointers[9][1]   = &AssLegFunction::P_9_1;
  leg_pointers[9][2]   = &AssLegFunction::P_9_2;
  leg_pointers[9][3]   = &AssLegFunction::P_9_3;
  leg_pointers[9][4]   = &AssLegFunction::P_9_4;
  leg_pointers[9][5]   = &AssLegFunction::P_9_5;
  leg_pointers[9][6]   = &AssLegFunction::P_9_6;
  leg_pointers[9][7]   = &AssLegFunction::P_9_7;
  leg_pointers[9][8]   = &AssLegFunction::P_9_8;
  leg_pointers[9][9]   = &AssLegFunction::P_9_9;
  leg_pointers[10][0]  = &AssLegFunction::P_10_0;
  leg_pointers[10][1]  = &AssLegFunction::P_10_1;
  leg_pointers[10][2]  = &AssLegFunction::P_10_2;
  leg_pointers[10][3]  = &AssLegFunction::P_10_3;
  leg_pointers[10][4]  = &AssLegFunction::P_10_4;
  leg_pointers[10][5]  = &AssLegFunction::P_10_5;
  leg_pointers[10][6]  = &AssLegFunction::P_10_6;
  leg_pointers[10][7]  = &AssLegFunction::P_10_7;
  leg_pointers[10][8]  = &AssLegFunction::P_10_8;
  leg_pointers[10][9]  = &AssLegFunction::P_10_9;
  leg_pointers[10][10] = &AssLegFunction::P_10_10;
  leg_pointers[11][0]  = &AssLegFunction::P_11_0;
  leg_pointers[11][1]  = &AssLegFunction::P_11_1;
  leg_pointers[11][2]  = &AssLegFunction::P_11_2;
  leg_pointers[11][3]  = &AssLegFunction::P_11_3;
  leg_pointers[11][4]  = &AssLegFunction::P_11_4;
  leg_pointers[11][5]  = &AssLegFunction::P_11_5;
  leg_pointers[11][6]  = &AssLegFunction::P_11_6;
  leg_pointers[11][7]  = &AssLegFunction::P_11_7;
  leg_pointers[11][8]  = &AssLegFunction::P_11_8;
  leg_pointers[11][9]  = &AssLegFunction::P_11_9;
  leg_pointers[11][10] = &AssLegFunction::P_11_10;
  leg_pointers[11][11] = &AssLegFunction::P_11_11;
  leg_pointers[12][0]  = &AssLegFunction::P_12_0;
  leg_pointers[12][1]  = &AssLegFunction::P_12_1;
  leg_pointers[12][2]  = &AssLegFunction::P_12_2;
  leg_pointers[12][3]  = &AssLegFunction::P_12_3;
  leg_pointers[12][4]  = &AssLegFunction::P_12_4;
  leg_pointers[12][5]  = &AssLegFunction::P_12_5;
  leg_pointers[12][6]  = &AssLegFunction::P_12_6;
  leg_pointers[12][7]  = &AssLegFunction::P_12_7;
  leg_pointers[12][8]  = &AssLegFunction::P_12_8;
  leg_pointers[12][9]  = &AssLegFunction::P_12_9;
  leg_pointers[12][10] = &AssLegFunction::P_12_10;
  leg_pointers[12][11] = &AssLegFunction::P_12_11;
  leg_pointers[12][12] = &AssLegFunction::P_12_12;
  leg_pointers[13][0]  = &AssLegFunction::P_13_0;
  leg_pointers[13][1]  = &AssLegFunction::P_13_1;
  leg_pointers[13][2]  = &AssLegFunction::P_13_2;
  leg_pointers[13][3]  = &AssLegFunction::P_13_3;
  leg_pointers[13][4]  = &AssLegFunction::P_13_4;
  leg_pointers[13][5]  = &AssLegFunction::P_13_5;
  leg_pointers[13][6]  = &AssLegFunction::P_13_6;
  leg_pointers[13][7]  = &AssLegFunction::P_13_7;
  leg_pointers[13][8]  = &AssLegFunction::P_13_8;
  leg_pointers[13][9]  = &AssLegFunction::P_13_9;
  leg_pointers[13][10] = &AssLegFunction::P_13_10;
  leg_pointers[13][11] = &AssLegFunction::P_13_11;
  leg_pointers[13][12] = &AssLegFunction::P_13_12;
  leg_pointers[13][13] = &AssLegFunction::P_13_13;
  leg_pointers[14][0]  = &AssLegFunction::P_14_0;
  leg_pointers[14][1]  = &AssLegFunction::P_14_1;
  leg_pointers[14][2]  = &AssLegFunction::P_14_2;
  leg_pointers[14][3]  = &AssLegFunction::P_14_3;
  leg_pointers[14][4]  = &AssLegFunction::P_14_4;
  leg_pointers[14][5]  = &AssLegFunction::P_14_5;
  leg_pointers[14][6]  = &AssLegFunction::P_14_6;
  leg_pointers[14][7]  = &AssLegFunction::P_14_7;
  leg_pointers[14][8]  = &AssLegFunction::P_14_8;
  leg_pointers[14][9]  = &AssLegFunction::P_14_9;
  leg_pointers[14][10] = &AssLegFunction::P_14_10;
  leg_pointers[14][11] = &AssLegFunction::P_14_11;
  leg_pointers[14][12] = &AssLegFunction::P_14_12;
  leg_pointers[14][13] = &AssLegFunction::P_14_13;
  leg_pointers[14][14] = &AssLegFunction::P_14_14;
  leg_pointers[15][0]  = &AssLegFunction::P_15_0;
  leg_pointers[15][1]  = &AssLegFunction::P_15_1;
  leg_pointers[15][2]  = &AssLegFunction::P_15_2;
  leg_pointers[15][3]  = &AssLegFunction::P_15_3;
  leg_pointers[15][4]  = &AssLegFunction::P_15_4;
  leg_pointers[15][5]  = &AssLegFunction::P_15_5;
  leg_pointers[15][6]  = &AssLegFunction::P_15_6;
  leg_pointers[15][7]  = &AssLegFunction::P_15_7;
  leg_pointers[15][8]  = &AssLegFunction::P_15_8;
  leg_pointers[15][9]  = &AssLegFunction::P_15_9;
  leg_pointers[15][10] = &AssLegFunction::P_15_10;
  leg_pointers[15][11] = &AssLegFunction::P_15_11;
  leg_pointers[15][12] = &AssLegFunction::P_15_12;
  leg_pointers[15][13] = &AssLegFunction::P_15_13;
  leg_pointers[15][14] = &AssLegFunction::P_15_14;
  leg_pointers[15][15] = &AssLegFunction::P_15_15;
  leg_pointers[16][0]  = &AssLegFunction::P_16_0;
  leg_pointers[16][1]  = &AssLegFunction::P_16_1;
  leg_pointers[16][2]  = &AssLegFunction::P_16_2;
  leg_pointers[16][3]  = &AssLegFunction::P_16_3;
  leg_pointers[16][4]  = &AssLegFunction::P_16_4;
  leg_pointers[16][5]  = &AssLegFunction::P_16_5;
  leg_pointers[16][6]  = &AssLegFunction::P_16_6;
  leg_pointers[16][7]  = &AssLegFunction::P_16_7;
  leg_pointers[16][8]  = &AssLegFunction::P_16_8;
  leg_pointers[16][9]  = &AssLegFunction::P_16_9;
  leg_pointers[16][10] = &AssLegFunction::P_16_10;
  leg_pointers[16][11] = &AssLegFunction::P_16_11;
  leg_pointers[16][12] = &AssLegFunction::P_16_12;
  leg_pointers[16][13] = &AssLegFunction::P_16_13;
  leg_pointers[16][14] = &AssLegFunction::P_16_14;
  leg_pointers[16][15] = &AssLegFunction::P_16_15;
  leg_pointers[16][16] = &AssLegFunction::P_16_16;
  leg_pointers[17][0]  = &AssLegFunction::P_17_0;
  leg_pointers[17][1]  = &AssLegFunction::P_17_1;
  leg_pointers[17][2]  = &AssLegFunction::P_17_2;
  leg_pointers[17][3]  = &AssLegFunction::P_17_3;
  leg_pointers[17][4]  = &AssLegFunction::P_17_4;
  leg_pointers[17][5]  = &AssLegFunction::P_17_5;
  leg_pointers[17][6]  = &AssLegFunction::P_17_6;
  leg_pointers[17][7]  = &AssLegFunction::P_17_7;
  leg_pointers[17][8]  = &AssLegFunction::P_17_8;
  leg_pointers[17][9]  = &AssLegFunction::P_17_9;
  leg_pointers[17][10] = &AssLegFunction::P_17_10;
  leg_pointers[17][11] = &AssLegFunction::P_17_11;
  leg_pointers[17][12] = &AssLegFunction::P_17_12;
  leg_pointers[17][13] = &AssLegFunction::P_17_13;
  leg_pointers[17][14] = &AssLegFunction::P_17_14;
  leg_pointers[17][15] = &AssLegFunction::P_17_15;
  leg_pointers[17][16] = &AssLegFunction::P_17_16;
  leg_pointers[17][17] = &AssLegFunction::P_17_17;
  leg_pointers[18][0]  = &AssLegFunction::P_18_0;
  leg_pointers[18][1]  = &AssLegFunction::P_18_1;
  leg_pointers[18][2]  = &AssLegFunction::P_18_2;
  leg_pointers[18][3]  = &AssLegFunction::P_18_3;
  leg_pointers[18][4]  = &AssLegFunction::P_18_4;
  leg_pointers[18][5]  = &AssLegFunction::P_18_5;
  leg_pointers[18][6]  = &AssLegFunction::P_18_6;
  leg_pointers[18][7]  = &AssLegFunction::P_18_7;
  leg_pointers[18][8]  = &AssLegFunction::P_18_8;
  leg_pointers[18][9]  = &AssLegFunction::P_18_9;
  leg_pointers[18][10] = &AssLegFunction::P_18_10;
  leg_pointers[18][11] = &AssLegFunction::P_18_11;
  leg_pointers[18][12] = &AssLegFunction::P_18_12;
  leg_pointers[18][13] = &AssLegFunction::P_18_13;
  leg_pointers[18][14] = &AssLegFunction::P_18_14;
  leg_pointers[18][15] = &AssLegFunction::P_18_15;
  leg_pointers[18][16] = &AssLegFunction::P_18_16;
  leg_pointers[18][17] = &AssLegFunction::P_18_17;
  leg_pointers[18][18] = &AssLegFunction::P_18_18;
  leg_pointers[19][0]  = &AssLegFunction::P_19_0;
  leg_pointers[19][1]  = &AssLegFunction::P_19_1;
  leg_pointers[19][2]  = &AssLegFunction::P_19_2;
  leg_pointers[19][3]  = &AssLegFunction::P_19_3;
  leg_pointers[19][4]  = &AssLegFunction::P_19_4;
  leg_pointers[19][5]  = &AssLegFunction::P_19_5;
  leg_pointers[19][6]  = &AssLegFunction::P_19_6;
  leg_pointers[19][7]  = &AssLegFunction::P_19_7;
  leg_pointers[19][8]  = &AssLegFunction::P_19_8;
  leg_pointers[19][9]  = &AssLegFunction::P_19_9;
  leg_pointers[19][10] = &AssLegFunction::P_19_10;
  leg_pointers[19][11] = &AssLegFunction::P_19_11;
  leg_pointers[19][12] = &AssLegFunction::P_19_12;
  leg_pointers[19][13] = &AssLegFunction::P_19_13;
  leg_pointers[19][14] = &AssLegFunction::P_19_14;
  leg_pointers[19][15] = &AssLegFunction::P_19_15;
  leg_pointers[19][16] = &AssLegFunction::P_19_16;
  leg_pointers[19][17] = &AssLegFunction::P_19_17;
  leg_pointers[19][18] = &AssLegFunction::P_19_18;
  leg_pointers[19][19] = &AssLegFunction::P_19_19;
  leg_pointers[20][0]  = &AssLegFunction::P_20_0;
  leg_pointers[20][1]  = &AssLegFunction::P_20_1;
  leg_pointers[20][2]  = &AssLegFunction::P_20_2;
  leg_pointers[20][3]  = &AssLegFunction::P_20_3;
  leg_pointers[20][4]  = &AssLegFunction::P_20_4;
  leg_pointers[20][5]  = &AssLegFunction::P_20_5;
  leg_pointers[20][6]  = &AssLegFunction::P_20_6;
  leg_pointers[20][7]  = &AssLegFunction::P_20_7;
  leg_pointers[20][8]  = &AssLegFunction::P_20_8;
  leg_pointers[20][9]  = &AssLegFunction::P_20_9;
  leg_pointers[20][10] = &AssLegFunction::P_20_10;
  leg_pointers[20][11] = &AssLegFunction::P_20_11;
  leg_pointers[20][12] = &AssLegFunction::P_20_12;
  leg_pointers[20][13] = &AssLegFunction::P_20_13;
  leg_pointers[20][14] = &AssLegFunction::P_20_14;
  leg_pointers[20][15] = &AssLegFunction::P_20_15;
  leg_pointers[20][16] = &AssLegFunction::P_20_16;
  leg_pointers[20][17] = &AssLegFunction::P_20_17;
  leg_pointers[20][18] = &AssLegFunction::P_20_18;
  leg_pointers[20][19] = &AssLegFunction::P_20_19;
  leg_pointers[20][20] = &AssLegFunction::P_20_20;

  leg_der_pointers.resize(21);
  for (unsigned int i = 0; i < 21; ++i)
    {
      leg_der_pointers[i].resize(i + 1);
    }

  leg_der_pointers[0][0]   = &AssLegFunction::P_0_0_Deriv;
  leg_der_pointers[1][0]   = &AssLegFunction::P_1_0_Deriv;
  leg_der_pointers[1][1]   = &AssLegFunction::P_1_1_Deriv;
  leg_der_pointers[2][0]   = &AssLegFunction::P_2_0_Deriv;
  leg_der_pointers[2][1]   = &AssLegFunction::P_2_1_Deriv;
  leg_der_pointers[2][2]   = &AssLegFunction::P_2_2_Deriv;
  leg_der_pointers[3][0]   = &AssLegFunction::P_3_0_Deriv;
  leg_der_pointers[3][1]   = &AssLegFunction::P_3_1_Deriv;
  leg_der_pointers[3][2]   = &AssLegFunction::P_3_2_Deriv;
  leg_der_pointers[3][3]   = &AssLegFunction::P_3_3_Deriv;
  leg_der_pointers[4][0]   = &AssLegFunction::P_4_0_Deriv;
  leg_der_pointers[4][1]   = &AssLegFunction::P_4_1_Deriv;
  leg_der_pointers[4][2]   = &AssLegFunction::P_4_2_Deriv;
  leg_der_pointers[4][3]   = &AssLegFunction::P_4_3_Deriv;
  leg_der_pointers[4][4]   = &AssLegFunction::P_4_4_Deriv;
  leg_der_pointers[5][0]   = &AssLegFunction::P_5_0_Deriv;
  leg_der_pointers[5][1]   = &AssLegFunction::P_5_1_Deriv;
  leg_der_pointers[5][2]   = &AssLegFunction::P_5_2_Deriv;
  leg_der_pointers[5][3]   = &AssLegFunction::P_5_3_Deriv;
  leg_der_pointers[5][4]   = &AssLegFunction::P_5_4_Deriv;
  leg_der_pointers[5][5]   = &AssLegFunction::P_5_5_Deriv;
  leg_der_pointers[6][0]   = &AssLegFunction::P_6_0_Deriv;
  leg_der_pointers[6][1]   = &AssLegFunction::P_6_1_Deriv;
  leg_der_pointers[6][2]   = &AssLegFunction::P_6_2_Deriv;
  leg_der_pointers[6][3]   = &AssLegFunction::P_6_3_Deriv;
  leg_der_pointers[6][4]   = &AssLegFunction::P_6_4_Deriv;
  leg_der_pointers[6][5]   = &AssLegFunction::P_6_5_Deriv;
  leg_der_pointers[6][6]   = &AssLegFunction::P_6_6_Deriv;
  leg_der_pointers[7][0]   = &AssLegFunction::P_7_0_Deriv;
  leg_der_pointers[7][1]   = &AssLegFunction::P_7_1_Deriv;
  leg_der_pointers[7][2]   = &AssLegFunction::P_7_2_Deriv;
  leg_der_pointers[7][3]   = &AssLegFunction::P_7_3_Deriv;
  leg_der_pointers[7][4]   = &AssLegFunction::P_7_4_Deriv;
  leg_der_pointers[7][5]   = &AssLegFunction::P_7_5_Deriv;
  leg_der_pointers[7][6]   = &AssLegFunction::P_7_6_Deriv;
  leg_der_pointers[7][7]   = &AssLegFunction::P_7_7_Deriv;
  leg_der_pointers[8][0]   = &AssLegFunction::P_8_0_Deriv;
  leg_der_pointers[8][1]   = &AssLegFunction::P_8_1_Deriv;
  leg_der_pointers[8][2]   = &AssLegFunction::P_8_2_Deriv;
  leg_der_pointers[8][3]   = &AssLegFunction::P_8_3_Deriv;
  leg_der_pointers[8][4]   = &AssLegFunction::P_8_4_Deriv;
  leg_der_pointers[8][5]   = &AssLegFunction::P_8_5_Deriv;
  leg_der_pointers[8][6]   = &AssLegFunction::P_8_6_Deriv;
  leg_der_pointers[8][7]   = &AssLegFunction::P_8_7_Deriv;
  leg_der_pointers[8][8]   = &AssLegFunction::P_8_8_Deriv;
  leg_der_pointers[9][0]   = &AssLegFunction::P_9_0_Deriv;
  leg_der_pointers[9][1]   = &AssLegFunction::P_9_1_Deriv;
  leg_der_pointers[9][2]   = &AssLegFunction::P_9_2_Deriv;
  leg_der_pointers[9][3]   = &AssLegFunction::P_9_3_Deriv;
  leg_der_pointers[9][4]   = &AssLegFunction::P_9_4_Deriv;
  leg_der_pointers[9][5]   = &AssLegFunction::P_9_5_Deriv;
  leg_der_pointers[9][6]   = &AssLegFunction::P_9_6_Deriv;
  leg_der_pointers[9][7]   = &AssLegFunction::P_9_7_Deriv;
  leg_der_pointers[9][8]   = &AssLegFunction::P_9_8_Deriv;
  leg_der_pointers[9][9]   = &AssLegFunction::P_9_9_Deriv;
  leg_der_pointers[10][0]  = &AssLegFunction::P_10_0_Deriv;
  leg_der_pointers[10][1]  = &AssLegFunction::P_10_1_Deriv;
  leg_der_pointers[10][2]  = &AssLegFunction::P_10_2_Deriv;
  leg_der_pointers[10][3]  = &AssLegFunction::P_10_3_Deriv;
  leg_der_pointers[10][4]  = &AssLegFunction::P_10_4_Deriv;
  leg_der_pointers[10][5]  = &AssLegFunction::P_10_5_Deriv;
  leg_der_pointers[10][6]  = &AssLegFunction::P_10_6_Deriv;
  leg_der_pointers[10][7]  = &AssLegFunction::P_10_7_Deriv;
  leg_der_pointers[10][8]  = &AssLegFunction::P_10_8_Deriv;
  leg_der_pointers[10][9]  = &AssLegFunction::P_10_9_Deriv;
  leg_der_pointers[10][10] = &AssLegFunction::P_10_10_Deriv;
  leg_der_pointers[11][0]  = &AssLegFunction::P_11_0_Deriv;
  leg_der_pointers[11][1]  = &AssLegFunction::P_11_1_Deriv;
  leg_der_pointers[11][2]  = &AssLegFunction::P_11_2_Deriv;
  leg_der_pointers[11][3]  = &AssLegFunction::P_11_3_Deriv;
  leg_der_pointers[11][4]  = &AssLegFunction::P_11_4_Deriv;
  leg_der_pointers[11][5]  = &AssLegFunction::P_11_5_Deriv;
  leg_der_pointers[11][6]  = &AssLegFunction::P_11_6_Deriv;
  leg_der_pointers[11][7]  = &AssLegFunction::P_11_7_Deriv;
  leg_der_pointers[11][8]  = &AssLegFunction::P_11_8_Deriv;
  leg_der_pointers[11][9]  = &AssLegFunction::P_11_9_Deriv;
  leg_der_pointers[11][10] = &AssLegFunction::P_11_10_Deriv;
  leg_der_pointers[11][11] = &AssLegFunction::P_11_11_Deriv;
  leg_der_pointers[12][0]  = &AssLegFunction::P_12_0_Deriv;
  leg_der_pointers[12][1]  = &AssLegFunction::P_12_1_Deriv;
  leg_der_pointers[12][2]  = &AssLegFunction::P_12_2_Deriv;
  leg_der_pointers[12][3]  = &AssLegFunction::P_12_3_Deriv;
  leg_der_pointers[12][4]  = &AssLegFunction::P_12_4_Deriv;
  leg_der_pointers[12][5]  = &AssLegFunction::P_12_5_Deriv;
  leg_der_pointers[12][6]  = &AssLegFunction::P_12_6_Deriv;
  leg_der_pointers[12][7]  = &AssLegFunction::P_12_7_Deriv;
  leg_der_pointers[12][8]  = &AssLegFunction::P_12_8_Deriv;
  leg_der_pointers[12][9]  = &AssLegFunction::P_12_9_Deriv;
  leg_der_pointers[12][10] = &AssLegFunction::P_12_10_Deriv;
  leg_der_pointers[12][11] = &AssLegFunction::P_12_11_Deriv;
  leg_der_pointers[12][12] = &AssLegFunction::P_12_12_Deriv;
  leg_der_pointers[13][0]  = &AssLegFunction::P_13_0_Deriv;
  leg_der_pointers[13][1]  = &AssLegFunction::P_13_1_Deriv;
  leg_der_pointers[13][2]  = &AssLegFunction::P_13_2_Deriv;
  leg_der_pointers[13][3]  = &AssLegFunction::P_13_3_Deriv;
  leg_der_pointers[13][4]  = &AssLegFunction::P_13_4_Deriv;
  leg_der_pointers[13][5]  = &AssLegFunction::P_13_5_Deriv;
  leg_der_pointers[13][6]  = &AssLegFunction::P_13_6_Deriv;
  leg_der_pointers[13][7]  = &AssLegFunction::P_13_7_Deriv;
  leg_der_pointers[13][8]  = &AssLegFunction::P_13_8_Deriv;
  leg_der_pointers[13][9]  = &AssLegFunction::P_13_9_Deriv;
  leg_der_pointers[13][10] = &AssLegFunction::P_13_10_Deriv;
  leg_der_pointers[13][11] = &AssLegFunction::P_13_11_Deriv;
  leg_der_pointers[13][12] = &AssLegFunction::P_13_12_Deriv;
  leg_der_pointers[13][13] = &AssLegFunction::P_13_13_Deriv;
  leg_der_pointers[14][0]  = &AssLegFunction::P_14_0_Deriv;
  leg_der_pointers[14][1]  = &AssLegFunction::P_14_1_Deriv;
  leg_der_pointers[14][2]  = &AssLegFunction::P_14_2_Deriv;
  leg_der_pointers[14][3]  = &AssLegFunction::P_14_3_Deriv;
  leg_der_pointers[14][4]  = &AssLegFunction::P_14_4_Deriv;
  leg_der_pointers[14][5]  = &AssLegFunction::P_14_5_Deriv;
  leg_der_pointers[14][6]  = &AssLegFunction::P_14_6_Deriv;
  leg_der_pointers[14][7]  = &AssLegFunction::P_14_7_Deriv;
  leg_der_pointers[14][8]  = &AssLegFunction::P_14_8_Deriv;
  leg_der_pointers[14][9]  = &AssLegFunction::P_14_9_Deriv;
  leg_der_pointers[14][10] = &AssLegFunction::P_14_10_Deriv;
  leg_der_pointers[14][11] = &AssLegFunction::P_14_11_Deriv;
  leg_der_pointers[14][12] = &AssLegFunction::P_14_12_Deriv;
  leg_der_pointers[14][13] = &AssLegFunction::P_14_13_Deriv;
  leg_der_pointers[14][14] = &AssLegFunction::P_14_14_Deriv;
  leg_der_pointers[15][0]  = &AssLegFunction::P_15_0_Deriv;
  leg_der_pointers[15][1]  = &AssLegFunction::P_15_1_Deriv;
  leg_der_pointers[15][2]  = &AssLegFunction::P_15_2_Deriv;
  leg_der_pointers[15][3]  = &AssLegFunction::P_15_3_Deriv;
  leg_der_pointers[15][4]  = &AssLegFunction::P_15_4_Deriv;
  leg_der_pointers[15][5]  = &AssLegFunction::P_15_5_Deriv;
  leg_der_pointers[15][6]  = &AssLegFunction::P_15_6_Deriv;
  leg_der_pointers[15][7]  = &AssLegFunction::P_15_7_Deriv;
  leg_der_pointers[15][8]  = &AssLegFunction::P_15_8_Deriv;
  leg_der_pointers[15][9]  = &AssLegFunction::P_15_9_Deriv;
  leg_der_pointers[15][10] = &AssLegFunction::P_15_10_Deriv;
  leg_der_pointers[15][11] = &AssLegFunction::P_15_11_Deriv;
  leg_der_pointers[15][12] = &AssLegFunction::P_15_12_Deriv;
  leg_der_pointers[15][13] = &AssLegFunction::P_15_13_Deriv;
  leg_der_pointers[15][14] = &AssLegFunction::P_15_14_Deriv;
  leg_der_pointers[15][15] = &AssLegFunction::P_15_15_Deriv;
  leg_der_pointers[16][0]  = &AssLegFunction::P_16_0_Deriv;
  leg_der_pointers[16][1]  = &AssLegFunction::P_16_1_Deriv;
  leg_der_pointers[16][2]  = &AssLegFunction::P_16_2_Deriv;
  leg_der_pointers[16][3]  = &AssLegFunction::P_16_3_Deriv;
  leg_der_pointers[16][4]  = &AssLegFunction::P_16_4_Deriv;
  leg_der_pointers[16][5]  = &AssLegFunction::P_16_5_Deriv;
  leg_der_pointers[16][6]  = &AssLegFunction::P_16_6_Deriv;
  leg_der_pointers[16][7]  = &AssLegFunction::P_16_7_Deriv;
  leg_der_pointers[16][8]  = &AssLegFunction::P_16_8_Deriv;
  leg_der_pointers[16][9]  = &AssLegFunction::P_16_9_Deriv;
  leg_der_pointers[16][10] = &AssLegFunction::P_16_10_Deriv;
  leg_der_pointers[16][11] = &AssLegFunction::P_16_11_Deriv;
  leg_der_pointers[16][12] = &AssLegFunction::P_16_12_Deriv;
  leg_der_pointers[16][13] = &AssLegFunction::P_16_13_Deriv;
  leg_der_pointers[16][14] = &AssLegFunction::P_16_14_Deriv;
  leg_der_pointers[16][15] = &AssLegFunction::P_16_15_Deriv;
  leg_der_pointers[16][16] = &AssLegFunction::P_16_16_Deriv;
  leg_der_pointers[17][0]  = &AssLegFunction::P_17_0_Deriv;
  leg_der_pointers[17][1]  = &AssLegFunction::P_17_1_Deriv;
  leg_der_pointers[17][2]  = &AssLegFunction::P_17_2_Deriv;
  leg_der_pointers[17][3]  = &AssLegFunction::P_17_3_Deriv;
  leg_der_pointers[17][4]  = &AssLegFunction::P_17_4_Deriv;
  leg_der_pointers[17][5]  = &AssLegFunction::P_17_5_Deriv;
  leg_der_pointers[17][6]  = &AssLegFunction::P_17_6_Deriv;
  leg_der_pointers[17][7]  = &AssLegFunction::P_17_7_Deriv;
  leg_der_pointers[17][8]  = &AssLegFunction::P_17_8_Deriv;
  leg_der_pointers[17][9]  = &AssLegFunction::P_17_9_Deriv;
  leg_der_pointers[17][10] = &AssLegFunction::P_17_10_Deriv;
  leg_der_pointers[17][11] = &AssLegFunction::P_17_11_Deriv;
  leg_der_pointers[17][12] = &AssLegFunction::P_17_12_Deriv;
  leg_der_pointers[17][13] = &AssLegFunction::P_17_13_Deriv;
  leg_der_pointers[17][14] = &AssLegFunction::P_17_14_Deriv;
  leg_der_pointers[17][15] = &AssLegFunction::P_17_15_Deriv;
  leg_der_pointers[17][16] = &AssLegFunction::P_17_16_Deriv;
  leg_der_pointers[17][17] = &AssLegFunction::P_17_17_Deriv;
  leg_der_pointers[18][0]  = &AssLegFunction::P_18_0_Deriv;
  leg_der_pointers[18][1]  = &AssLegFunction::P_18_1_Deriv;
  leg_der_pointers[18][2]  = &AssLegFunction::P_18_2_Deriv;
  leg_der_pointers[18][3]  = &AssLegFunction::P_18_3_Deriv;
  leg_der_pointers[18][4]  = &AssLegFunction::P_18_4_Deriv;
  leg_der_pointers[18][5]  = &AssLegFunction::P_18_5_Deriv;
  leg_der_pointers[18][6]  = &AssLegFunction::P_18_6_Deriv;
  leg_der_pointers[18][7]  = &AssLegFunction::P_18_7_Deriv;
  leg_der_pointers[18][8]  = &AssLegFunction::P_18_8_Deriv;
  leg_der_pointers[18][9]  = &AssLegFunction::P_18_9_Deriv;
  leg_der_pointers[18][10] = &AssLegFunction::P_18_10_Deriv;
  leg_der_pointers[18][11] = &AssLegFunction::P_18_11_Deriv;
  leg_der_pointers[18][12] = &AssLegFunction::P_18_12_Deriv;
  leg_der_pointers[18][13] = &AssLegFunction::P_18_13_Deriv;
  leg_der_pointers[18][14] = &AssLegFunction::P_18_14_Deriv;
  leg_der_pointers[18][15] = &AssLegFunction::P_18_15_Deriv;
  leg_der_pointers[18][16] = &AssLegFunction::P_18_16_Deriv;
  leg_der_pointers[18][17] = &AssLegFunction::P_18_17_Deriv;
  leg_der_pointers[18][18] = &AssLegFunction::P_18_18_Deriv;
  leg_der_pointers[19][0]  = &AssLegFunction::P_19_0_Deriv;
  leg_der_pointers[19][1]  = &AssLegFunction::P_19_1_Deriv;
  leg_der_pointers[19][2]  = &AssLegFunction::P_19_2_Deriv;
  leg_der_pointers[19][3]  = &AssLegFunction::P_19_3_Deriv;
  leg_der_pointers[19][4]  = &AssLegFunction::P_19_4_Deriv;
  leg_der_pointers[19][5]  = &AssLegFunction::P_19_5_Deriv;
  leg_der_pointers[19][6]  = &AssLegFunction::P_19_6_Deriv;
  leg_der_pointers[19][7]  = &AssLegFunction::P_19_7_Deriv;
  leg_der_pointers[19][8]  = &AssLegFunction::P_19_8_Deriv;
  leg_der_pointers[19][9]  = &AssLegFunction::P_19_9_Deriv;
  leg_der_pointers[19][10] = &AssLegFunction::P_19_10_Deriv;
  leg_der_pointers[19][11] = &AssLegFunction::P_19_11_Deriv;
  leg_der_pointers[19][12] = &AssLegFunction::P_19_12_Deriv;
  leg_der_pointers[19][13] = &AssLegFunction::P_19_13_Deriv;
  leg_der_pointers[19][14] = &AssLegFunction::P_19_14_Deriv;
  leg_der_pointers[19][15] = &AssLegFunction::P_19_15_Deriv;
  leg_der_pointers[19][16] = &AssLegFunction::P_19_16_Deriv;
  leg_der_pointers[19][17] = &AssLegFunction::P_19_17_Deriv;
  leg_der_pointers[19][18] = &AssLegFunction::P_19_18_Deriv;
  leg_der_pointers[19][19] = &AssLegFunction::P_19_19_Deriv;
  leg_der_pointers[20][0]  = &AssLegFunction::P_20_0_Deriv;
  leg_der_pointers[20][1]  = &AssLegFunction::P_20_1_Deriv;
  leg_der_pointers[20][2]  = &AssLegFunction::P_20_2_Deriv;
  leg_der_pointers[20][3]  = &AssLegFunction::P_20_3_Deriv;
  leg_der_pointers[20][4]  = &AssLegFunction::P_20_4_Deriv;
  leg_der_pointers[20][5]  = &AssLegFunction::P_20_5_Deriv;
  leg_der_pointers[20][6]  = &AssLegFunction::P_20_6_Deriv;
  leg_der_pointers[20][7]  = &AssLegFunction::P_20_7_Deriv;
  leg_der_pointers[20][8]  = &AssLegFunction::P_20_8_Deriv;
  leg_der_pointers[20][9]  = &AssLegFunction::P_20_9_Deriv;
  leg_der_pointers[20][10] = &AssLegFunction::P_20_10_Deriv;
  leg_der_pointers[20][11] = &AssLegFunction::P_20_11_Deriv;
  leg_der_pointers[20][12] = &AssLegFunction::P_20_12_Deriv;
  leg_der_pointers[20][13] = &AssLegFunction::P_20_13_Deriv;
  leg_der_pointers[20][14] = &AssLegFunction::P_20_14_Deriv;
  leg_der_pointers[20][15] = &AssLegFunction::P_20_15_Deriv;
  leg_der_pointers[20][16] = &AssLegFunction::P_20_16_Deriv;
  leg_der_pointers[20][17] = &AssLegFunction::P_20_17_Deriv;
  leg_der_pointers[20][18] = &AssLegFunction::P_20_18_Deriv;
  leg_der_pointers[20][19] = &AssLegFunction::P_20_19_Deriv;
  leg_der_pointers[20][20] = &AssLegFunction::P_20_20_Deriv;
}

AssLegFunction::~AssLegFunction()
{}

void
AssLegFunction::AssLegFunSph(const unsigned int p,
                             const unsigned int m,
                             double             x,
                             double             y[])
{
  for (unsigned int n = m; n < p + 1; n++)
    {
      y[n - m] = (this->*leg_pointers[n][m])(x);
    }
}

void
AssLegFunction::AssLegFunSphDeriv(const unsigned int p,
                                  const unsigned int m,
                                  double             x,
                                  double             y[],
                                  double             dy[])
{
  for (unsigned int n = m; n < p + 1; n++)
    {
      y[n - m]  = (this->*leg_pointers[n][m])(x);
      dy[n - m] = (this->*leg_der_pointers[n][m])(x);
    }
}

double
AssLegFunction::GetAssLegFunSph(const unsigned int n,
                                const unsigned int m,
                                double             x) const
{
  return (this->*leg_pointers[n][m])(x);
}

double
AssLegFunction::GetAssLegFunSphDeriv(const unsigned int n,
                                     const unsigned int m,
                                     double             x) const
{
  return (this->*leg_der_pointers[n][m])(x);
}

double
AssLegFunction::P_0_0(const double x) const
{
  return 1;
}

double
AssLegFunction::P_0_0_Deriv(const double x) const
{
  return 0;
}

double
AssLegFunction::P_1_0(const double x) const
{
  return x;
}

double
AssLegFunction::P_1_0_Deriv(const double x) const
{
  return 1;
}

double
AssLegFunction::P_1_1(const double x) const
{
  return -sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_1_1_Deriv(const double x) const
{
  return x / sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_2_0(const double x) const
{
  return (3.0 / 2.0) * pow(x, 2) - 1.0 / 2.0;
}

double
AssLegFunction::P_2_0_Deriv(const double x) const
{
  return 3 * x;
}

double
AssLegFunction::P_2_1(const double x) const
{
  return -3 * x * sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_2_1_Deriv(const double x) const
{
  return 3 * (2 * pow(x, 2) - 1) / sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_2_2(const double x) const
{
  return 3 - 3 * pow(x, 2);
}

double
AssLegFunction::P_2_2_Deriv(const double x) const
{
  return -6 * x;
}

double
AssLegFunction::P_3_0(const double x) const
{
  return (1.0 / 2.0) * x * (5 * pow(x, 2) - 3);
}

double
AssLegFunction::P_3_0_Deriv(const double x) const
{
  return (15.0 / 2.0) * pow(x, 2) - 3.0 / 2.0;
}

double
AssLegFunction::P_3_1(const double x) const
{
  return (3.0 / 2.0) * (1 - 5 * pow(x, 2)) * sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_3_1_Deriv(const double x) const
{
  return (1.0 / 2.0) * (45 * pow(x, 3) - 33 * x) / sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_3_2(const double x) const
{
  return 15 * x * (1 - pow(x, 2));
}

double
AssLegFunction::P_3_2_Deriv(const double x) const
{
  return 15 - 45 * pow(x, 2);
}

double
AssLegFunction::P_3_3(const double x) const
{
  return -15 * pow(1 - pow(x, 2), 3.0 / 2.0);
}

double
AssLegFunction::P_3_3_Deriv(const double x) const
{
  return 45 * x * sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_4_0(const double x) const
{
  return (35.0 / 8.0) * pow(x, 4) - 15.0 / 4.0 * pow(x, 2) + 3.0 / 8.0;
}

double
AssLegFunction::P_4_0_Deriv(const double x) const
{
  return (5.0 / 2.0) * x * (7 * pow(x, 2) - 3);
}

double
AssLegFunction::P_4_1(const double x) const
{
  return (5.0 / 2.0) * x * sqrt(1 - pow(x, 2)) * (3 - 7 * pow(x, 2));
}

double
AssLegFunction::P_4_1_Deriv(const double x) const
{
  return (5.0 / 2.0) * (28 * pow(x, 4) - 27 * pow(x, 2) + 3) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_4_2(const double x) const
{
  return -105.0 / 2.0 * pow(x, 4) + 60 * pow(x, 2) - 15.0 / 2.0;
}

double
AssLegFunction::P_4_2_Deriv(const double x) const
{
  return -210 * pow(x, 3) + 120 * x;
}

double
AssLegFunction::P_4_3(const double x) const
{
  return -105 * x * pow(1 - pow(x, 2), 3.0 / 2.0);
}

double
AssLegFunction::P_4_3_Deriv(const double x) const
{
  return sqrt(1 - pow(x, 2)) * (420 * pow(x, 2) - 105);
}

double
AssLegFunction::P_4_4(const double x) const
{
  return 105 * pow(pow(x, 2) - 1, 2);
}

double
AssLegFunction::P_4_4_Deriv(const double x) const
{
  return 420 * x * (pow(x, 2) - 1);
}

double
AssLegFunction::P_5_0(const double x) const
{
  return (1.0 / 8.0) * x * (63 * pow(x, 4) - 70 * pow(x, 2) + 15);
}

double
AssLegFunction::P_5_0_Deriv(const double x) const
{
  return (315.0 / 8.0) * pow(x, 4) - 105.0 / 4.0 * pow(x, 2) + 15.0 / 8.0;
}

double
AssLegFunction::P_5_1(const double x) const
{
  return (15.0 / 8.0) * sqrt(1 - pow(x, 2)) *
         (-21 * pow(x, 4) + 14 * pow(x, 2) - 1);
}

double
AssLegFunction::P_5_1_Deriv(const double x) const
{
  return (1.0 / 8.0) * (1575 * pow(x, 5) - 1890 * pow(x, 3) + 435 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_5_2(const double x) const
{
  return (105.0 / 2.0) * x * (-3 * pow(x, 4) + 4 * pow(x, 2) - 1);
}

double
AssLegFunction::P_5_2_Deriv(const double x) const
{
  return -1575.0 / 2.0 * pow(x, 4) + 630 * pow(x, 2) - 105.0 / 2.0;
}

double
AssLegFunction::P_5_3(const double x) const
{
  return (105.0 / 2.0) * (1 - 9 * pow(x, 2)) * pow(1 - pow(x, 2), 3.0 / 2.0);
}

double
AssLegFunction::P_5_3_Deriv(const double x) const
{
  return (315.0 / 2.0) * x * sqrt(1 - pow(x, 2)) * (15 * pow(x, 2) - 7);
}

double
AssLegFunction::P_5_4(const double x) const
{
  return 945 * x * pow(pow(x, 2) - 1, 2);
}

double
AssLegFunction::P_5_4_Deriv(const double x) const
{
  return 4725 * pow(x, 4) - 5670 * pow(x, 2) + 945;
}

double
AssLegFunction::P_5_5(const double x) const
{
  return -945 * pow(1 - pow(x, 2), 5.0 / 2.0);
}

double
AssLegFunction::P_5_5_Deriv(const double x) const
{
  return 4725 * x * pow(1 - pow(x, 2), 3.0 / 2.0);
}

double
AssLegFunction::P_6_0(const double x) const
{
  return (231.0 / 16.0) * pow(x, 6) - 315.0 / 16.0 * pow(x, 4) +
         (105.0 / 16.0) * pow(x, 2) - 5.0 / 16.0;
}

double
AssLegFunction::P_6_0_Deriv(const double x) const
{
  return (21.0 / 8.0) * x * (33 * pow(x, 4) - 30 * pow(x, 2) + 5);
}

double
AssLegFunction::P_6_1(const double x) const
{
  return (21.0 / 8.0) * x * sqrt(1 - pow(x, 2)) *
         (-33 * pow(x, 4) + 30 * pow(x, 2) - 5);
}

double
AssLegFunction::P_6_1_Deriv(const double x) const
{
  return (21.0 / 8.0) *
         (198 * pow(x, 6) - 285 * pow(x, 4) + 100 * pow(x, 2) - 5) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_6_2(const double x) const
{
  return -3465.0 / 8.0 * pow(x, 6) + (5355.0 / 8.0) * pow(x, 4) -
         1995.0 / 8.0 * pow(x, 2) + 105.0 / 8.0;
}

double
AssLegFunction::P_6_2_Deriv(const double x) const
{
  return (105.0 / 4.0) * x * (-99 * pow(x, 4) + 102 * pow(x, 2) - 19);
}

double
AssLegFunction::P_6_3(const double x) const
{
  return (315.0 / 2.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (3 - 11 * pow(x, 2));
}

double
AssLegFunction::P_6_3_Deriv(const double x) const
{
  return (945.0 / 2.0) * sqrt(1 - pow(x, 2)) *
         (22 * pow(x, 4) - 15 * pow(x, 2) + 1);
}

double
AssLegFunction::P_6_4(const double x) const
{
  return (945.0 / 2.0) * pow(pow(x, 2) - 1, 2) * (11 * pow(x, 2) - 1);
}

double
AssLegFunction::P_6_4_Deriv(const double x) const
{
  return 31185 * pow(x, 5) - 43470 * pow(x, 3) + 12285 * x;
}

double
AssLegFunction::P_6_5(const double x) const
{
  return -10395 * x * pow(1 - pow(x, 2), 5.0 / 2.0);
}

double
AssLegFunction::P_6_5_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 3.0 / 2.0) * (62370 * pow(x, 2) - 10395);
}

double
AssLegFunction::P_6_6(const double x) const
{
  return -10395 * pow(pow(x, 2) - 1, 3);
}

double
AssLegFunction::P_6_6_Deriv(const double x) const
{
  return -62370 * x * pow(pow(x, 2) - 1, 2);
}

double
AssLegFunction::P_7_0(const double x) const
{
  return (1.0 / 16.0) * x *
         (429 * pow(x, 6) - 693 * pow(x, 4) + 315 * pow(x, 2) - 35);
}

double
AssLegFunction::P_7_0_Deriv(const double x) const
{
  return (3003.0 / 16.0) * pow(x, 6) - 3465.0 / 16.0 * pow(x, 4) +
         (945.0 / 16.0) * pow(x, 2) - 35.0 / 16.0;
}

double
AssLegFunction::P_7_1(const double x) const
{
  return (7.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (-429 * pow(x, 6) + 495 * pow(x, 4) - 135 * pow(x, 2) + 5);
}

double
AssLegFunction::P_7_1_Deriv(const double x) const
{
  return (1.0 / 16.0) *
         (21021 * pow(x, 7) - 35343 * pow(x, 5) + 16695 * pow(x, 3) -
          1925 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_7_2(const double x) const
{
  return (63.0 / 8.0) * x *
         (-143 * pow(x, 6) + 253 * pow(x, 4) - 125 * pow(x, 2) + 15);
}

double
AssLegFunction::P_7_2_Deriv(const double x) const
{
  return -63063.0 / 8.0 * pow(x, 6) + (79695.0 / 8.0) * pow(x, 4) -
         23625.0 / 8.0 * pow(x, 2) + 945.0 / 8.0;
}

double
AssLegFunction::P_7_3(const double x) const
{
  return (315.0 / 8.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-143 * pow(x, 4) + 66 * pow(x, 2) - 3);
}

double
AssLegFunction::P_7_3_Deriv(const double x) const
{
  return (315.0 / 8.0) * x * sqrt(1 - pow(x, 2)) *
         (1001 * pow(x, 4) - 902 * pow(x, 2) + 141);
}

double
AssLegFunction::P_7_4(const double x) const
{
  return (3465.0 / 2.0) * x * pow(pow(x, 2) - 1, 2) * (13 * pow(x, 2) - 3);
}

double
AssLegFunction::P_7_4_Deriv(const double x) const
{
  return (315315.0 / 2.0) * pow(x, 6) - 502425.0 / 2.0 * pow(x, 4) +
         (197505.0 / 2.0) * pow(x, 2) - 10395.0 / 2.0;
}

double
AssLegFunction::P_7_5(const double x) const
{
  return (10395.0 / 2.0) * (1 - 13 * pow(x, 2)) * pow(1 - pow(x, 2), 5.0 / 2.0);
}

double
AssLegFunction::P_7_5_Deriv(const double x) const
{
  return (10395.0 / 2.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (91 * pow(x, 2) - 31);
}

double
AssLegFunction::P_7_6(const double x) const
{
  return -135135 * x * pow(pow(x, 2) - 1, 3);
}

double
AssLegFunction::P_7_6_Deriv(const double x) const
{
  return (135135 - 945945 * pow(x, 2)) * pow(pow(x, 2) - 1, 2);
}

double
AssLegFunction::P_7_7(const double x) const
{
  return -135135 * pow(1 - pow(x, 2), 7.0 / 2.0);
}

double
AssLegFunction::P_7_7_Deriv(const double x) const
{
  return 945945 * x * pow(1 - pow(x, 2), 5.0 / 2.0);
}

double
AssLegFunction::P_8_0(const double x) const
{
  return (6435.0 / 128.0) * pow(x, 8) - 3003.0 / 32.0 * pow(x, 6) +
         (3465.0 / 64.0) * pow(x, 4) - 315.0 / 32.0 * pow(x, 2) + 35.0 / 128.0;
}

double
AssLegFunction::P_8_0_Deriv(const double x) const
{
  return (9.0 / 16.0) * x *
         (715 * pow(x, 6) - 1001 * pow(x, 4) + 385 * pow(x, 2) - 35);
}

double
AssLegFunction::P_8_1(const double x) const
{
  return (9.0 / 16.0) * x * sqrt(1 - pow(x, 2)) *
         (-715 * pow(x, 6) + 1001 * pow(x, 4) - 385 * pow(x, 2) + 35);
}

double
AssLegFunction::P_8_1_Deriv(const double x) const
{
  return (9.0 / 16.0) *
         (5720 * pow(x, 8) - 11011 * pow(x, 6) + 6545 * pow(x, 4) -
          1225 * pow(x, 2) + 35) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_8_2(const double x) const
{
  return -45045.0 / 16.0 * pow(x, 8) + (45045.0 / 8.0) * pow(x, 6) -
         3465 * pow(x, 4) + (5355.0 / 8.0) * pow(x, 2) - 315.0 / 16.0;
}

double
AssLegFunction::P_8_2_Deriv(const double x) const
{
  return (315.0 / 4.0) * x *
         (-286 * pow(x, 6) + 429 * pow(x, 4) - 176 * pow(x, 2) + 17);
}

double
AssLegFunction::P_8_3(const double x) const
{
  return (3465.0 / 8.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-39 * pow(x, 4) + 26 * pow(x, 2) - 3);
}

double
AssLegFunction::P_8_3_Deriv(const double x) const
{
  return (10395.0 / 8.0) * sqrt(1 - pow(x, 2)) *
         (104 * pow(x, 6) - 117 * pow(x, 4) + 30 * pow(x, 2) - 1);
}

double
AssLegFunction::P_8_4(const double x) const
{
  return (10395.0 / 8.0) * pow(pow(x, 2) - 1, 2) *
         (65 * pow(x, 4) - 26 * pow(x, 2) + 1);
}

double
AssLegFunction::P_8_4_Deriv(const double x) const
{
  return 675675 * pow(x, 7) - 1216215 * pow(x, 5) + 613305 * pow(x, 3) -
         72765 * x;
}

double
AssLegFunction::P_8_5(const double x) const
{
  return (135135.0 / 2.0) * x * (1 - 5 * pow(x, 2)) *
         pow(1 - pow(x, 2), 5.0 / 2.0);
}

double
AssLegFunction::P_8_5_Deriv(const double x) const
{
  return (135135.0 / 2.0) * sqrt(1 - pow(x, 2)) *
         (-40 * pow(x, 6) + 61 * pow(x, 4) - 22 * pow(x, 2) + 1);
}

double
AssLegFunction::P_8_6(const double x) const
{
  return (135135.0 / 2.0) * (1 - 15 * pow(x, 2)) * pow(pow(x, 2) - 1, 3);
}

double
AssLegFunction::P_8_6_Deriv(const double x) const
{
  return 810810 * x * (3 - 10 * pow(x, 2)) * pow(pow(x, 2) - 1, 2);
}

double
AssLegFunction::P_8_7(const double x) const
{
  return -2027025 * x * pow(1 - pow(x, 2), 7.0 / 2.0);
}

double
AssLegFunction::P_8_7_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 5.0 / 2.0) * (16216200 * pow(x, 2) - 2027025);
}

double
AssLegFunction::P_8_8(const double x) const
{
  return 2027025 * pow(pow(x, 2) - 1, 4);
}

double
AssLegFunction::P_8_8_Deriv(const double x) const
{
  return 16216200 * x * pow(pow(x, 2) - 1, 3);
}

double
AssLegFunction::P_9_0(const double x) const
{
  return (1.0 / 128.0) * x *
         (12155 * pow(x, 8) - 25740 * pow(x, 6) + 18018 * pow(x, 4) -
          4620 * pow(x, 2) + 315);
}

double
AssLegFunction::P_9_0_Deriv(const double x) const
{
  return (109395.0 / 128.0) * pow(x, 8) - 45045.0 / 32.0 * pow(x, 6) +
         (45045.0 / 64.0) * pow(x, 4) - 3465.0 / 32.0 * pow(x, 2) +
         315.0 / 128.0;
}

double
AssLegFunction::P_9_1(const double x) const
{
  return (45.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (-2431 * pow(x, 8) + 4004 * pow(x, 6) - 2002 * pow(x, 4) +
          308 * pow(x, 2) - 7);
}

double
AssLegFunction::P_9_1_Deriv(const double x) const
{
  return (1.0 / 128.0) *
         (984555 * pow(x, 9) - 2136420 * pow(x, 7) + 1531530 * pow(x, 5) -
          401940 * pow(x, 3) + 28035 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_9_2(const double x) const
{
  return (495.0 / 16.0) * x *
         (-221 * pow(x, 8) + 494 * pow(x, 6) - 364 * pow(x, 4) +
          98 * pow(x, 2) - 7);
}

double
AssLegFunction::P_9_2_Deriv(const double x) const
{
  return -984555.0 / 16.0 * pow(x, 8) + (855855.0 / 8.0) * pow(x, 6) -
         225225.0 / 4.0 * pow(x, 4) + (72765.0 / 8.0) * pow(x, 2) -
         3465.0 / 16.0;
}

double
AssLegFunction::P_9_3(const double x) const
{
  return (3465.0 / 16.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-221 * pow(x, 6) + 195 * pow(x, 4) - 39 * pow(x, 2) + 1);
}

double
AssLegFunction::P_9_3_Deriv(const double x) const
{
  return (10395.0 / 16.0) * x * sqrt(1 - pow(x, 2)) *
         (663 * pow(x, 6) - 897 * pow(x, 4) + 325 * pow(x, 2) - 27);
}

double
AssLegFunction::P_9_4(const double x) const
{
  return (135135.0 / 8.0) * x * pow(pow(x, 2) - 1, 2) *
         (17 * pow(x, 4) - 10 * pow(x, 2) + 1);
}

double
AssLegFunction::P_9_4_Deriv(const double x) const
{
  return (20675655.0 / 8.0) * pow(x, 8) - 10405395.0 / 2.0 * pow(x, 6) +
         (12837825.0 / 4.0) * pow(x, 4) - 1216215.0 / 2.0 * pow(x, 2) +
         135135.0 / 8.0;
}

double
AssLegFunction::P_9_5(const double x) const
{
  return (135135.0 / 8.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-85 * pow(x, 4) + 30 * pow(x, 2) - 1);
}

double
AssLegFunction::P_9_5_Deriv(const double x) const
{
  return (675675.0 / 8.0) * x * sqrt(1 - pow(x, 2)) *
         (-153 * pow(x, 6) + 263 * pow(x, 4) - 123 * pow(x, 2) + 13);
}

double
AssLegFunction::P_9_6(const double x) const
{
  return (675675.0 / 2.0) * x * (3 - 17 * pow(x, 2)) * pow(pow(x, 2) - 1, 3);
}

double
AssLegFunction::P_9_6_Deriv(const double x) const
{
  return -103378275.0 / 2.0 * pow(x, 8) + 127702575.0 * pow(x, 6) -
         101351250.0 * pow(x, 4) + 26351325 * pow(x, 2) - 2027025.0 / 2.0;
}

double
AssLegFunction::P_9_7(const double x) const
{
  return (2027025.0 / 2.0) * (1 - 17 * pow(x, 2)) *
         pow(1 - pow(x, 2), 7.0 / 2.0);
}

double
AssLegFunction::P_9_7_Deriv(const double x) const
{
  return (2027025.0 / 2.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (153 * pow(x, 2) - 41);
}

double
AssLegFunction::P_9_8(const double x) const
{
  return 34459425 * x * pow(pow(x, 2) - 1, 4);
}

double
AssLegFunction::P_9_8_Deriv(const double x) const
{
  return pow(pow(x, 2) - 1, 3) * (310134825.0 * pow(x, 2) - 34459425);
}

double
AssLegFunction::P_9_9(const double x) const
{
  return -34459425 * pow(1 - pow(x, 2), 9.0 / 2.0);
}

double
AssLegFunction::P_9_9_Deriv(const double x) const
{
  return 310134825.0 * x * pow(1 - pow(x, 2), 7.0 / 2.0);
}

double
AssLegFunction::P_10_0(const double x) const
{
  return (46189.0 / 256.0) * pow(x, 10) - 109395.0 / 256.0 * pow(x, 8) +
         (45045.0 / 128.0) * pow(x, 6) - 15015.0 / 128.0 * pow(x, 4) +
         (3465.0 / 256.0) * pow(x, 2) - 63.0 / 256.0;
}

double
AssLegFunction::P_10_0_Deriv(const double x) const
{
  return (55.0 / 128.0) * x *
         (4199 * pow(x, 8) - 7956 * pow(x, 6) + 4914 * pow(x, 4) -
          1092 * pow(x, 2) + 63);
}

double
AssLegFunction::P_10_1(const double x) const
{
  return (55.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (-4199 * pow(x, 8) + 7956 * pow(x, 6) - 4914 * pow(x, 4) +
          1092 * pow(x, 2) - 63);
}

double
AssLegFunction::P_10_1_Deriv(const double x) const
{
  return (55.0 / 128.0) *
         (41990 * pow(x, 10) - 101439 * pow(x, 8) + 85176 * pow(x, 6) -
          28938 * pow(x, 4) + 3402 * pow(x, 2) - 63) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_10_2(const double x) const
{
  return -2078505.0 / 128.0 * pow(x, 10) + (5141565.0 / 128.0) * pow(x, 8) -
         2207205.0 / 64.0 * pow(x, 6) + (765765.0 / 64.0) * pow(x, 4) -
         183645.0 / 128.0 * pow(x, 2) + 3465.0 / 128.0;
}

double
AssLegFunction::P_10_2_Deriv(const double x) const
{
  return (495.0 / 64.0) * x *
         (-20995 * pow(x, 8) + 41548 * pow(x, 6) - 26754 * pow(x, 4) +
          6188 * pow(x, 2) - 371);
}

double
AssLegFunction::P_10_3(const double x) const
{
  return (6435.0 / 16.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-323 * pow(x, 6) + 357 * pow(x, 4) - 105 * pow(x, 2) + 7);
}

double
AssLegFunction::P_10_3_Deriv(const double x) const
{
  return (6435.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (3230 * pow(x, 8) - 5117 * pow(x, 6) + 2415 * pow(x, 4) -
          343 * pow(x, 2) + 7);
}

double
AssLegFunction::P_10_4(const double x) const
{
  return (45045.0 / 16.0) * pow(pow(x, 2) - 1, 2) *
         (323 * pow(x, 6) - 255 * pow(x, 4) + 45 * pow(x, 2) - 1);
}

double
AssLegFunction::P_10_4_Deriv(const double x) const
{
  return (45045.0 / 8.0) * x *
         (1615 * pow(x, 8) - 3604 * pow(x, 6) + 2634 * pow(x, 4) -
          692 * pow(x, 2) + 47);
}

double
AssLegFunction::P_10_5(const double x) const
{
  return (135135.0 / 8.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-323 * pow(x, 4) + 170 * pow(x, 2) - 15);
}

double
AssLegFunction::P_10_5_Deriv(const double x) const
{
  return (675675.0 / 8.0) * sqrt(1 - pow(x, 2)) *
         (-646 * pow(x, 8) + 1241 * pow(x, 6) - 715 * pow(x, 4) +
          123 * pow(x, 2) - 3);
}

double
AssLegFunction::P_10_6(const double x) const
{
  return (675675.0 / 8.0) * pow(pow(x, 2) - 1, 3) *
         (-323 * pow(x, 4) + 102 * pow(x, 2) - 3);
}

double
AssLegFunction::P_10_6_Deriv(const double x) const
{
  return (675675.0 / 4.0) * x *
         (-1615 * pow(x, 8) + 4284 * pow(x, 6) - 3834 * pow(x, 4) +
          1276 * pow(x, 2) - 111);
}

double
AssLegFunction::P_10_7(const double x) const
{
  return (11486475.0 / 2.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (3 - 19 * pow(x, 2));
}

double
AssLegFunction::P_10_7_Deriv(const double x) const
{
  return (11486475.0 / 2.0) * sqrt(1 - pow(x, 2)) *
         (190 * pow(x, 8) - 461 * pow(x, 6) + 355 * pow(x, 4) - 87 * pow(x, 2) +
          3);
}

double
AssLegFunction::P_10_8(const double x) const
{
  return (34459425.0 / 2.0) * pow(pow(x, 2) - 1, 4) * (19 * pow(x, 2) - 1);
}

double
AssLegFunction::P_10_8_Deriv(const double x) const
{
  return 34459425 * x * pow(pow(x, 2) - 1, 3) * (95 * pow(x, 2) - 23);
}

double
AssLegFunction::P_10_9(const double x) const
{
  return -654729075.0 * x * pow(1 - pow(x, 2), 9.0 / 2.0);
}

double
AssLegFunction::P_10_9_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 7.0 / 2.0) *
         (6547290750.0 * pow(x, 2) - 654729075.0);
}

double
AssLegFunction::P_10_10(const double x) const
{
  return -654729075.0 * pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_10_10_Deriv(const double x) const
{
  return -6547290750.0 * x * pow(pow(x, 2) - 1, 4);
}

double
AssLegFunction::P_11_0(const double x) const
{
  return (1.0 / 256.0) * x *
         (88179 * pow(x, 10) - 230945 * pow(x, 8) + 218790 * pow(x, 6) -
          90090 * pow(x, 4) + 15015 * pow(x, 2) - 693);
}

double
AssLegFunction::P_11_0_Deriv(const double x) const
{
  return (969969.0 / 256.0) * pow(x, 10) - 2078505.0 / 256.0 * pow(x, 8) +
         (765765.0 / 128.0) * pow(x, 6) - 225225.0 / 128.0 * pow(x, 4) +
         (45045.0 / 256.0) * pow(x, 2) - 693.0 / 256.0;
}

double
AssLegFunction::P_11_1(const double x) const
{
  return (33.0 / 256.0) * sqrt(1 - pow(x, 2)) *
         (-29393 * pow(x, 10) + 62985 * pow(x, 8) - 46410 * pow(x, 6) +
          13650 * pow(x, 4) - 1365 * pow(x, 2) + 21);
}

double
AssLegFunction::P_11_1_Deriv(const double x) const
{
  return (1.0 / 256.0) *
         (10669659 * pow(x, 11) - 28406235 * pow(x, 9) + 27348750 * pow(x, 7) -
          11441430 * pow(x, 5) + 1936935 * pow(x, 3) - 90783 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_11_2(const double x) const
{
  return (2145.0 / 128.0) * x *
         (-2261 * pow(x, 10) + 6137 * pow(x, 8) - 6018 * pow(x, 6) +
          2562 * pow(x, 4) - 441 * pow(x, 2) + 21);
}

double
AssLegFunction::P_11_2_Deriv(const double x) const
{
  return -53348295.0 / 128.0 * pow(x, 10) + (118474785.0 / 128.0) * pow(x, 8) -
         45180135.0 / 64.0 * pow(x, 6) + (13738725.0 / 64.0) * pow(x, 4) -
         2837835.0 / 128.0 * pow(x, 2) + 45045.0 / 128.0;
}

double
AssLegFunction::P_11_3(const double x) const
{
  return (45045.0 / 128.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-969 * pow(x, 8) + 1292 * pow(x, 6) - 510 * pow(x, 4) +
          60 * pow(x, 2) - 1);
}

double
AssLegFunction::P_11_3_Deriv(const double x) const
{
  return (135135.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (3553 * pow(x, 8) - 6460 * pow(x, 6) + 3774 * pow(x, 4) -
          780 * pow(x, 2) + 41);
}

double
AssLegFunction::P_11_4(const double x) const
{
  return (135135.0 / 16.0) * x * pow(pow(x, 2) - 1, 2) *
         (323 * pow(x, 6) - 323 * pow(x, 4) + 85 * pow(x, 2) - 5);
}

double
AssLegFunction::P_11_4_Deriv(const double x) const
{
  return (480134655.0 / 16.0) * pow(x, 10) - 1178512335.0 / 16.0 * pow(x, 8) +
         (498513015.0 / 8.0) * pow(x, 6) - 168243075.0 / 8.0 * pow(x, 4) +
         (38513475.0 / 16.0) * pow(x, 2) - 675675.0 / 16.0;
}

double
AssLegFunction::P_11_5(const double x) const
{
  return (135135.0 / 16.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-2261 * pow(x, 6) + 1615 * pow(x, 4) - 255 * pow(x, 2) + 5);
}

double
AssLegFunction::P_11_5_Deriv(const double x) const
{
  return (135135.0 / 16.0) * x * sqrt(1 - pow(x, 2)) *
         (-24871 * pow(x, 8) + 52972 * pow(x, 6) - 36346 * pow(x, 4) +
          8780 * pow(x, 2) - 535);
}

double
AssLegFunction::P_11_6(const double x) const
{
  return (2297295.0 / 8.0) * x * pow(pow(x, 2) - 1, 3) *
         (-399 * pow(x, 4) + 190 * pow(x, 2) - 15);
}

double
AssLegFunction::P_11_6_Deriv(const double x) const
{
  return -10082827755.0 / 8.0 * pow(x, 10) + (28677133485.0 / 8.0) * pow(x, 8) -
         14328228915.0 / 4.0 * pow(x, 6) + (5823642825.0 / 4.0) * pow(x, 4) -
         1619592975.0 / 8.0 * pow(x, 2) + 34459425.0 / 8.0;
}

double
AssLegFunction::P_11_7(const double x) const
{
  return (34459425.0 / 8.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-133 * pow(x, 4) + 38 * pow(x, 2) - 1);
}

double
AssLegFunction::P_11_7_Deriv(const double x) const
{
  return (34459425.0 / 8.0) * x * sqrt(1 - pow(x, 2)) *
         (1463 * pow(x, 8) - 3800 * pow(x, 6) + 3294 * pow(x, 4) -
          1040 * pow(x, 2) + 83);
}

double
AssLegFunction::P_11_8(const double x) const
{
  return (654729075.0 / 2.0) * x * pow(pow(x, 2) - 1, 4) * (7 * pow(x, 2) - 1);
}

double
AssLegFunction::P_11_8_Deriv(const double x) const
{
  return (50414138775.0 / 2.0) * pow(x, 10) - 170884288575.0 / 2.0 * pow(x, 8) +
         105411381075.0 * pow(x, 6) - 55651971375.0 * pow(x, 4) +
         (21606059475.0 / 2.0) * pow(x, 2) - 654729075.0 / 2.0;
}

double
AssLegFunction::P_11_9(const double x) const
{
  return (654729075.0 / 2.0) * (1 - 21 * pow(x, 2)) *
         pow(1 - pow(x, 2), 9.0 / 2.0);
}

double
AssLegFunction::P_11_9_Deriv(const double x) const
{
  return (1964187225.0 / 2.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (77 * pow(x, 2) - 17);
}

double
AssLegFunction::P_11_10(const double x) const
{
  return -13749310575.0 * x * pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_11_10_Deriv(const double x) const
{
  return (13749310575.0 - 151242416325.0 * pow(x, 2)) * pow(pow(x, 2) - 1, 4);
}

double
AssLegFunction::P_11_11(const double x) const
{
  return -13749310575.0 * pow(1 - pow(x, 2), 11.0 / 2.0);
}

double
AssLegFunction::P_11_11_Deriv(const double x) const
{
  return 151242416325.0 * x * pow(1 - pow(x, 2), 9.0 / 2.0);
}

double
AssLegFunction::P_12_0(const double x) const
{
  return (676039.0 / 1024.0) * pow(x, 12) - 969969.0 / 512.0 * pow(x, 10) +
         (2078505.0 / 1024.0) * pow(x, 8) - 255255.0 / 256.0 * pow(x, 6) +
         (225225.0 / 1024.0) * pow(x, 4) - 9009.0 / 512.0 * pow(x, 2) +
         231.0 / 1024.0;
}

double
AssLegFunction::P_12_0_Deriv(const double x) const
{
  return (39.0 / 256.0) * x *
         (52003 * pow(x, 10) - 124355 * pow(x, 8) + 106590 * pow(x, 6) -
          39270 * pow(x, 4) + 5775 * pow(x, 2) - 231);
}

double
AssLegFunction::P_12_1(const double x) const
{
  return (39.0 / 256.0) * x * sqrt(1 - pow(x, 2)) *
         (-52003 * pow(x, 10) + 124355 * pow(x, 8) - 106590 * pow(x, 6) +
          39270 * pow(x, 4) - 5775 * pow(x, 2) + 231);
}

double
AssLegFunction::P_12_1_Deriv(const double x) const
{
  return (39.0 / 256.0) *
         (624036 * pow(x, 12) - 1815583 * pow(x, 10) + 1971915 * pow(x, 8) -
          981750 * pow(x, 6) + 219450 * pow(x, 4) - 17787 * pow(x, 2) + 231) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_12_2(const double x) const
{
  return -22309287.0 / 256.0 * pow(x, 12) + (16489473.0 / 64.0) * pow(x, 10) -
         72747675.0 / 256.0 * pow(x, 8) + (2297295.0 / 16.0) * pow(x, 6) -
         8333325.0 / 256.0 * pow(x, 4) + (171171.0 / 64.0) * pow(x, 2) -
         9009.0 / 256.0;
}

double
AssLegFunction::P_12_2_Deriv(const double x) const
{
  return (3003.0 / 64.0) * x *
         (-22287 * pow(x, 10) + 54910 * pow(x, 8) - 48450 * pow(x, 6) +
          18360 * pow(x, 4) - 2775 * pow(x, 2) + 114);
}

double
AssLegFunction::P_12_3(const double x) const
{
  return (15015.0 / 128.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-7429 * pow(x, 8) + 11628 * pow(x, 6) - 5814 * pow(x, 4) +
          1020 * pow(x, 2) - 45);
}

double
AssLegFunction::P_12_3_Deriv(const double x) const
{
  return (45045.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (29716 * pow(x, 10) - 61047 * pow(x, 8) + 42636 * pow(x, 6) -
          11730 * pow(x, 4) + 1080 * pow(x, 2) - 15);
}

double
AssLegFunction::P_12_4(const double x) const
{
  return (135135.0 / 128.0) * pow(pow(x, 2) - 1, 2) *
         (7429 * pow(x, 8) - 9044 * pow(x, 6) + 3230 * pow(x, 4) -
          340 * pow(x, 2) + 5);
}

double
AssLegFunction::P_12_4_Deriv(const double x) const
{
  return (135135.0 / 32.0) * x *
         (22287 * pow(x, 10) - 59755 * pow(x, 8) + 57494 * pow(x, 6) -
          23766 * pow(x, 4) + 3915 * pow(x, 2) - 175);
}

double
AssLegFunction::P_12_5(const double x) const
{
  return (2297295.0 / 16.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-437 * pow(x, 6) + 399 * pow(x, 4) - 95 * pow(x, 2) + 5);
}

double
AssLegFunction::P_12_5_Deriv(const double x) const
{
  return (2297295.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (-5244 * pow(x, 10) + 12293 * pow(x, 8) - 9804 * pow(x, 6) +
          3070 * pow(x, 4) - 320 * pow(x, 2) + 5);
}

double
AssLegFunction::P_12_6(const double x) const
{
  return (2297295.0 / 16.0) * pow(pow(x, 2) - 1, 3) *
         (-3059 * pow(x, 6) + 1995 * pow(x, 4) - 285 * pow(x, 2) + 5);
}

double
AssLegFunction::P_12_6_Deriv(const double x) const
{
  return (6891885.0 / 4.0) * x *
         (-3059 * pow(x, 10) + 9310 * pow(x, 8) - 10298 * pow(x, 6) +
          4952 * pow(x, 4) - 955 * pow(x, 2) + 50);
}

double
AssLegFunction::P_12_7(const double x) const
{
  return (130945815.0 / 8.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-161 * pow(x, 4) + 70 * pow(x, 2) - 5);
}

double
AssLegFunction::P_12_7_Deriv(const double x) const
{
  return (130945815.0 / 8.0) * sqrt(1 - pow(x, 2)) *
         (1932 * pow(x, 10) - 5369 * pow(x, 8) + 5192 * pow(x, 6) -
          2010 * pow(x, 4) + 260 * pow(x, 2) - 5);
}

double
AssLegFunction::P_12_8(const double x) const
{
  return (654729075.0 / 8.0) * pow(pow(x, 2) - 1, 4) *
         (161 * pow(x, 4) - 42 * pow(x, 2) + 1);
}

double
AssLegFunction::P_12_8_Deriv(const double x) const
{
  return (654729075.0 / 2.0) * x *
         (483 * pow(x, 10) - 1715 * pow(x, 8) + 2270 * pow(x, 6) -
          1350 * pow(x, 4) + 335 * pow(x, 2) - 23);
}

double
AssLegFunction::P_12_9(const double x) const
{
  return (4583103525.0 / 2.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (3 - 23 * pow(x, 2));
}

double
AssLegFunction::P_12_9_Deriv(const double x) const
{
  return (13749310575.0 / 2.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (3 * pow(x, 2) * (23 * pow(x, 2) - 3) +
          (1 - 23 * pow(x, 2)) * (1 - pow(x, 2)));
}

double
AssLegFunction::P_12_10(const double x) const
{
  return (13749310575.0 / 2.0) * (1 - 23 * pow(x, 2)) * pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_12_10_Deriv(const double x) const
{
  return 27498621150.0 * x * (14 - 69 * pow(x, 2)) * pow(pow(x, 2) - 1, 4);
}

double
AssLegFunction::P_12_11(const double x) const
{
  return -316234143225.0 * x * pow(1 - pow(x, 2), 11.0 / 2.0);
}

double
AssLegFunction::P_12_11_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 9.0 / 2.0) *
         (3794809718700.0 * pow(x, 2) - 316234143225.0);
}

double
AssLegFunction::P_12_12(const double x) const
{
  return 316234143225.0 * pow(pow(x, 2) - 1, 6);
}

double
AssLegFunction::P_12_12_Deriv(const double x) const
{
  return 3794809718700.0 * x * pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_13_0(const double x) const
{
  return (1.0 / 1024.0) * x *
         (1300075 * pow(x, 12) - 4056234 * pow(x, 10) + 4849845 * pow(x, 8) -
          2771340 * pow(x, 6) + 765765 * pow(x, 4) - 90090 * pow(x, 2) + 3003);
}

double
AssLegFunction::P_13_0_Deriv(const double x) const
{
  return (16900975.0 / 1024.0) * pow(x, 12) - 22309287.0 / 512.0 * pow(x, 10) +
         (43648605.0 / 1024.0) * pow(x, 8) - 4849845.0 / 256.0 * pow(x, 6) +
         (3828825.0 / 1024.0) * pow(x, 4) - 135135.0 / 512.0 * pow(x, 2) +
         3003.0 / 1024.0;
}

double
AssLegFunction::P_13_1(const double x) const
{
  return (91.0 / 1024.0) * sqrt(1 - pow(x, 2)) *
         (-185725 * pow(x, 12) + 490314 * pow(x, 10) - 479655 * pow(x, 8) +
          213180 * pow(x, 6) - 42075 * pow(x, 4) + 2970 * pow(x, 2) - 33);
}

double
AssLegFunction::P_13_1_Deriv(const double x) const
{
  return (1.0 / 1024.0) *
         (219712675.0 * pow(x, 13) - 693616014.0 * pow(x, 11) +
          839023185.0 * pow(x, 9) - 484984500.0 * pow(x, 7) +
          135540405.0 * pow(x, 5) - 16126110 * pow(x, 3) + 543543 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_13_2(const double x) const
{
  return (1365.0 / 256.0) * x *
         (-37145 * pow(x, 12) + 118864 * pow(x, 10) - 145673 * pow(x, 8) +
          85272 * pow(x, 6) - 24123 * pow(x, 4) + 2904 * pow(x, 2) - 99);
}

double
AssLegFunction::P_13_2_Deriv(const double x) const
{
  return -659138025.0 / 256.0 * pow(x, 12) + (111546435.0 / 16.0) * pow(x, 10) -
         1789592805.0 / 256.0 * pow(x, 8) + (101846745.0 / 32.0) * pow(x, 6) -
         164639475.0 / 256.0 * pow(x, 4) + (1486485.0 / 32.0) * pow(x, 2) -
         135135.0 / 256.0;
}

double
AssLegFunction::P_13_3(const double x) const
{
  return (15015.0 / 256.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-37145 * pow(x, 10) + 66861 * pow(x, 8) - 40698 * pow(x, 6) +
          9690 * pow(x, 4) - 765 * pow(x, 2) + 9);
}

double
AssLegFunction::P_13_3_Deriv(const double x) const
{
  return (15015.0 / 256.0) * x * sqrt(1 - pow(x, 2)) *
         (482885 * pow(x, 10) - 1106921 * pow(x, 8) + 901170 * pow(x, 6) -
          312018 * pow(x, 4) + 42585 * pow(x, 2) - 1557);
}

double
AssLegFunction::P_13_4(const double x) const
{
  return (255255.0 / 128.0) * x * pow(pow(x, 2) - 1, 2) *
         (10925 * pow(x, 8) - 15732 * pow(x, 6) + 7182 * pow(x, 4) -
          1140 * pow(x, 2) + 45);
}

double
AssLegFunction::P_13_4_Deriv(const double x) const
{
  return (36252591375.0 / 128.0) * pow(x, 12) -
         52761463755.0 / 64.0 * pow(x, 10) +
         (113879210445.0 / 128.0) * pow(x, 8) -
         13953004065.0 / 32.0 * pow(x, 6) +
         (12133546425.0 / 128.0) * pow(x, 4) - 470945475.0 / 64.0 * pow(x, 2) +
         11486475.0 / 128.0;
}

double
AssLegFunction::P_13_5(const double x) const
{
  return (2297295.0 / 128.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-10925 * pow(x, 8) + 12236 * pow(x, 6) - 3990 * pow(x, 4) +
          380 * pow(x, 2) - 5);
}

double
AssLegFunction::P_13_5_Deriv(const double x) const
{
  return (2297295.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (-142025 * pow(x, 10) + 364021 * pow(x, 8) - 331322 * pow(x, 6) +
          127946 * pow(x, 4) - 19405 * pow(x, 2) + 785);
}

double
AssLegFunction::P_13_6(const double x) const
{
  return (43648605.0 / 16.0) * x * pow(pow(x, 2) - 1, 3) *
         (-575 * pow(x, 6) + 483 * pow(x, 4) - 105 * pow(x, 2) + 5);
}

double
AssLegFunction::P_13_6_Deriv(const double x) const
{
  return -326273322375.0 / 16.0 * pow(x, 12) + 66258582390.0 * pow(x, 10) -
         1288113982155.0 / 16.0 * pow(x, 8) +
         (89523288855.0 / 2.0) * pow(x, 6) - 177431579325.0 / 16.0 * pow(x, 4) +
         (1964187225.0 / 2.0) * pow(x, 2) - 218243025.0 / 16.0;
}

double
AssLegFunction::P_13_7(const double x) const
{
  return (218243025.0 / 16.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-805 * pow(x, 6) + 483 * pow(x, 4) - 63 * pow(x, 2) + 1);
}

double
AssLegFunction::P_13_7_Deriv(const double x) const
{
  return (1527701175.0 / 16.0) * x * sqrt(1 - pow(x, 2)) *
         (1495 * pow(x, 10) - 4439 * pow(x, 8) + 4750 * pow(x, 6) -
          2182 * pow(x, 4) + 395 * pow(x, 2) - 19);
}

double
AssLegFunction::P_13_8(const double x) const
{
  return (4583103525.0 / 8.0) * x * pow(pow(x, 2) - 1, 4) *
         (115 * pow(x, 4) - 46 * pow(x, 2) + 3);
}

double
AssLegFunction::P_13_8_Deriv(const double x) const
{
  return (6851739769875.0 / 8.0) * pow(x, 12) -
         12754777110075.0 / 4.0 * pow(x, 10) +
         (36174436122825.0 / 8.0) * pow(x, 8) -
         5999282514225.0 / 2.0 * pow(x, 6) +
         (7264219087125.0 / 8.0) * pow(x, 4) -
         398730006675.0 / 4.0 * pow(x, 2) + 13749310575.0 / 8.0;
}

double
AssLegFunction::P_13_9(const double x) const
{
  return (4583103525.0 / 8.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-575 * pow(x, 4) + 138 * pow(x, 2) - 3);
}

double
AssLegFunction::P_13_9_Deriv(const double x) const
{
  return (4583103525.0 / 8.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (5175 * pow(x, 4) - 1242 * pow(x, 2) +
          92 * (1 - pow(x, 2)) * (3 - 25 * pow(x, 2)) + 27);
}

double
AssLegFunction::P_13_10(const double x) const
{
  return (105411381075.0 / 2.0) * x * (3 - 25 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_13_10_Deriv(const double x) const
{
  return -34258698849375.0 / 2.0 * pow(x, 12) + 74209612276800.0 * pow(x, 10) -
         251406143863875.0 / 2.0 * pow(x, 8) + 103303153453500.0 * pow(x, 6) -
         81693820333125.0 / 2.0 * pow(x, 4) + 6324682864500.0 * pow(x, 2) -
         316234143225.0 / 2.0;
}

double
AssLegFunction::P_13_11(const double x) const
{
  return (316234143225.0 / 2.0) * (1 - 25 * pow(x, 2)) *
         pow(1 - pow(x, 2), 11.0 / 2.0);
}

double
AssLegFunction::P_13_11_Deriv(const double x) const
{
  return (316234143225.0 / 2.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (325 * pow(x, 2) - 61);
}

double
AssLegFunction::P_13_12(const double x) const
{
  return 7905853580625.0 * x * pow(pow(x, 2) - 1, 6);
}

double
AssLegFunction::P_13_12_Deriv(const double x) const
{
  return pow(pow(x, 2) - 1, 5) *
         (102776096548125.0 * pow(x, 2) - 7905853580625.0);
}

double
AssLegFunction::P_13_13(const double x) const
{
  return -7905853580625.0 * pow(1 - pow(x, 2), 13.0 / 2.0);
}

double
AssLegFunction::P_13_13_Deriv(const double x) const
{
  return 102776096548125.0 * x * pow(1 - pow(x, 2), 11.0 / 2.0);
}

double
AssLegFunction::P_14_0(const double x) const
{
  return (5014575.0 / 2048.0) * pow(x, 14) - 16900975.0 / 2048.0 * pow(x, 12) +
         (22309287.0 / 2048.0) * pow(x, 10) - 14549535.0 / 2048.0 * pow(x, 8) +
         (4849845.0 / 2048.0) * pow(x, 6) - 765765.0 / 2048.0 * pow(x, 4) +
         (45045.0 / 2048.0) * pow(x, 2) - 429.0 / 2048.0;
}

double
AssLegFunction::P_14_0_Deriv(const double x) const
{
  return (105.0 / 1024.0) * x *
         (334305 * pow(x, 12) - 965770 * pow(x, 10) + 1062347 * pow(x, 8) -
          554268 * pow(x, 6) + 138567 * pow(x, 4) - 14586 * pow(x, 2) + 429);
}

double
AssLegFunction::P_14_1(const double x) const
{
  return (105.0 / 1024.0) * x * sqrt(1 - pow(x, 2)) *
         (-334305 * pow(x, 12) + 965770 * pow(x, 10) - 1062347 * pow(x, 8) +
          554268 * pow(x, 6) - 138567 * pow(x, 4) + 14586 * pow(x, 2) - 429);
}

double
AssLegFunction::P_14_1_Deriv(const double x) const
{
  return (105.0 / 1024.0) *
         (4680270 * pow(x, 14) - 15935205 * pow(x, 12) + 21246940 * pow(x, 10) -
          13995267 * pow(x, 8) + 4711278 * pow(x, 6) - 751179 * pow(x, 4) +
          44616 * pow(x, 2) - 429) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_14_2(const double x) const
{
  return -456326325.0 / 1024.0 * pow(x, 14) +
         (1571790675.0 / 1024.0) * pow(x, 12) -
         2119382265.0 / 1024.0 * pow(x, 10) +
         (1411304895.0 / 1024.0) * pow(x, 8) -
         480134655.0 / 1024.0 * pow(x, 6) + (77342265.0 / 1024.0) * pow(x, 4) -
         4639635.0 / 1024.0 * pow(x, 2) + 45045.0 / 1024.0;
}

double
AssLegFunction::P_14_2_Deriv(const double x) const
{
  return (1365.0 / 512.0) * x *
         (-2340135 * pow(x, 12) + 6908970 * pow(x, 10) - 7763305 * pow(x, 8) +
          4135692 * pow(x, 6) - 1055241 * pow(x, 4) + 113322 * pow(x, 2) -
          3399);
}

double
AssLegFunction::P_14_3(const double x) const
{
  return (23205.0 / 256.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-58995 * pow(x, 10) + 120175 * pow(x, 8) - 86526 * pow(x, 6) +
          26334 * pow(x, 4) - 3135 * pow(x, 2) + 99);
}

double
AssLegFunction::P_14_3_Deriv(const double x) const
{
  return (69615.0 / 256.0) * sqrt(1 - pow(x, 2)) *
         (275310 * pow(x, 12) - 697015 * pow(x, 10) + 648945 * pow(x, 8) -
          272118 * pow(x, 6) + 50160 * pow(x, 4) - 3267 * pow(x, 2) + 33);
}

double
AssLegFunction::P_14_4(const double x) const
{
  return (2297295.0 / 256.0) * pow(pow(x, 2) - 1, 2) *
         (6555 * pow(x, 10) - 10925 * pow(x, 8) + 6118 * pow(x, 6) -
          1330 * pow(x, 4) + 95 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_4_Deriv(const double x) const
{
  return (2297295.0 / 128.0) * x *
         (45885 * pow(x, 12) - 144210 * pow(x, 10) + 172615 * pow(x, 8) -
          97964 * pow(x, 6) + 26619 * pow(x, 4) - 3042 * pow(x, 2) + 97);
}

double
AssLegFunction::P_14_5(const double x) const
{
  return (43648605.0 / 128.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-1725 * pow(x, 8) + 2300 * pow(x, 6) - 966 * pow(x, 4) +
          140 * pow(x, 2) - 5);
}

double
AssLegFunction::P_14_5_Deriv(const double x) const
{
  return (218243025.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (-4830 * pow(x, 12) + 13455 * pow(x, 10) - 13777 * pow(x, 8) +
          6342 * pow(x, 6) - 1280 * pow(x, 4) + 91 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_6(const double x) const
{
  return (218243025.0 / 128.0) * pow(pow(x, 2) - 1, 3) *
         (-3105 * pow(x, 8) + 3220 * pow(x, 6) - 966 * pow(x, 4) +
          84 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_6_Deriv(const double x) const
{
  return (654729075.0 / 64.0) * x *
         (-7245 * pow(x, 12) + 25070 * pow(x, 10) - 33235 * pow(x, 8) +
          20996 * pow(x, 6) - 6371 * pow(x, 4) + 814 * pow(x, 2) - 29);
}

double
AssLegFunction::P_14_7(const double x) const
{
  return (654729075.0 / 16.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-1035 * pow(x, 6) + 805 * pow(x, 4) - 161 * pow(x, 2) + 7);
}

double
AssLegFunction::P_14_7_Deriv(const double x) const
{
  return (4583103525.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (2070 * pow(x, 12) - 6555 * pow(x, 10) + 7705 * pow(x, 8) -
          4102 * pow(x, 6) + 960 * pow(x, 4) - 79 * pow(x, 2) + 1);
}

double
AssLegFunction::P_14_8(const double x) const
{
  return (4583103525.0 / 16.0) * pow(pow(x, 2) - 1, 4) *
         (1035 * pow(x, 6) - 575 * pow(x, 4) + 69 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_8_Deriv(const double x) const
{
  return (4583103525.0 / 8.0) * x *
         (7245 * pow(x, 12) - 28290 * pow(x, 10) + 42895 * pow(x, 8) -
          31468 * pow(x, 6) + 11259 * pow(x, 4) - 1714 * pow(x, 2) + 73);
}

double
AssLegFunction::P_14_9(const double x) const
{
  return (105411381075.0 / 8.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-135 * pow(x, 4) + 50 * pow(x, 2) - 3);
}

double
AssLegFunction::P_14_9_Deriv(const double x) const
{
  return (316234143225.0 / 8.0) * sqrt(1 - pow(x, 2)) *
         (-630 * pow(x, 12) + 2315 * pow(x, 10) - 3225 * pow(x, 8) +
          2086 * pow(x, 6) - 608 * pow(x, 4) + 63 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_10(const double x) const
{
  return (316234143225.0 / 8.0) * pow(pow(x, 2) - 1, 5) *
         (-225 * pow(x, 4) + 50 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_10_Deriv(const double x) const
{
  return (1581170716125.0 / 4.0) * x * pow(pow(x, 2) - 1, 4) *
         (-225 * pow(x, 4) + 50 * pow(x, 2) +
          10 * (1 - 9 * pow(x, 2)) * (pow(x, 2) - 1) - 1);
}

double
AssLegFunction::P_14_11(const double x) const
{
  return (7905853580625.0 / 2.0) * x * (1 - 9 * pow(x, 2)) *
         pow(1 - pow(x, 2), 11.0 / 2.0);
}

double
AssLegFunction::P_14_11_Deriv(const double x) const
{
  return (7905853580625.0 / 2.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (11 * pow(x, 2) * (9 * pow(x, 2) - 1) +
          (1 - 27 * pow(x, 2)) * (1 - pow(x, 2)));
}

double
AssLegFunction::P_14_12(const double x) const
{
  return (7905853580625.0 / 2.0) * pow(pow(x, 2) - 1, 6) * (27 * pow(x, 2) - 1);
}

double
AssLegFunction::P_14_12_Deriv(const double x) const
{
  return 23717560741875.0 * x * pow(pow(x, 2) - 1, 5) * (63 * pow(x, 2) - 11);
}

double
AssLegFunction::P_14_13(const double x) const
{
  return -213458046676875.0 * x * pow(1 - pow(x, 2), 13.0 / 2.0);
}

double
AssLegFunction::P_14_13_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 11.0 / 2.0) *
         (2988412653476250.0 * pow(x, 2) - 213458046676875.0);
}

double
AssLegFunction::P_14_14(const double x) const
{
  return -213458046676875.0 * pow(pow(x, 2) - 1, 7);
}

double
AssLegFunction::P_14_14_Deriv(const double x) const
{
  return -2988412653476250.0 * x * pow(pow(x, 2) - 1, 6);
}

double
AssLegFunction::P_15_0(const double x) const
{
  return (1.0 / 2048.0) * x *
         (9694845 * pow(x, 14) - 35102025 * pow(x, 12) + 50702925 * pow(x, 10) -
          37182145 * pow(x, 8) + 14549535 * pow(x, 6) - 2909907 * pow(x, 4) +
          255255 * pow(x, 2) - 6435);
}

double
AssLegFunction::P_15_0_Deriv(const double x) const
{
  return (145422675.0 / 2048.0) * pow(x, 14) -
         456326325.0 / 2048.0 * pow(x, 12) +
         (557732175.0 / 2048.0) * pow(x, 10) -
         334639305.0 / 2048.0 * pow(x, 8) + (101846745.0 / 2048.0) * pow(x, 6) -
         14549535.0 / 2048.0 * pow(x, 4) + (765765.0 / 2048.0) * pow(x, 2) -
         6435.0 / 2048.0;
}

double
AssLegFunction::P_15_1(const double x) const
{
  return (15.0 / 2048.0) * sqrt(1 - pow(x, 2)) *
         (-9694845 * pow(x, 14) + 30421755 * pow(x, 12) -
          37182145 * pow(x, 10) + 22309287 * pow(x, 8) - 6789783 * pow(x, 6) +
          969969 * pow(x, 4) - 51051 * pow(x, 2) + 429);
}

double
AssLegFunction::P_15_1_Deriv(const double x) const
{
  return (1.0 / 2048.0) *
         (2181340125.0 * pow(x, 15) - 7968159675.0 * pow(x, 13) +
          11610969825.0 * pow(x, 11) - 8589075495.0 * pow(x, 9) +
          3390041655.0 * pow(x, 7) - 683828145.0 * pow(x, 5) +
          60495435 * pow(x, 3) - 1537965 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_15_2(const double x) const
{
  return (1785.0 / 1024.0) * x *
         (-570285 * pow(x, 14) + 2104155 * pow(x, 12) - 3096145 * pow(x, 10) +
          2312167 * pow(x, 8) - 921063 * pow(x, 6) + 187473 * pow(x, 4) -
          16731 * pow(x, 2) + 429);
}

double
AssLegFunction::P_15_2_Deriv(const double x) const
{
  return -15269380875.0 / 1024.0 * pow(x, 14) +
         (48826916775.0 / 1024.0) * pow(x, 12) -
         60792807075.0 / 1024.0 * pow(x, 10) +
         (37144962855.0 / 1024.0) * pow(x, 8) -
         11508682185.0 / 1024.0 * pow(x, 6) +
         (1673196525.0 / 1024.0) * pow(x, 4) - 89594505.0 / 1024.0 * pow(x, 2) +
         765765.0 / 1024.0;
}

double
AssLegFunction::P_15_3(const double x) const
{
  return (69615.0 / 1024.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-190095 * pow(x, 12) + 432630 * pow(x, 10) - 360525 * pow(x, 8) +
          134596 * pow(x, 6) - 21945 * pow(x, 4) + 1254 * pow(x, 2) - 11);
}

double
AssLegFunction::P_15_3_Deriv(const double x) const
{
  return (208845.0 / 1024.0) * x * sqrt(1 - pow(x, 2)) *
         (950475 * pow(x, 12) - 2635110 * pow(x, 10) + 2764025 * pow(x, 8) -
          1365188 * pow(x, 6) + 320397 * pow(x, 4) - 31350 * pow(x, 2) + 847);
}

double
AssLegFunction::P_15_4(const double x) const
{
  return (3968055.0 / 256.0) * x * pow(pow(x, 2) - 1, 2) *
         (10005 * pow(x, 10) - 18975 * pow(x, 8) + 12650 * pow(x, 6) -
          3542 * pow(x, 4) + 385 * pow(x, 2) - 11);
}

double
AssLegFunction::P_15_4_Deriv(const double x) const
{
  return (595505854125.0 / 256.0) * pow(x, 14) -
         2011030114275.0 / 256.0 * pow(x, 12) +
         (2645323706025.0 / 256.0) * pow(x, 10) -
         1707664373415.0 / 256.0 * pow(x, 8) +
         (558833089815.0 / 256.0) * pow(x, 6) -
         85769508825.0 / 256.0 * pow(x, 4) +
         (4844995155.0 / 256.0) * pow(x, 2) - 43648605.0 / 256.0;
}

double
AssLegFunction::P_15_5(const double x) const
{
  return (43648605.0 / 256.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-10005 * pow(x, 10) + 15525 * pow(x, 8) - 8050 * pow(x, 6) +
          1610 * pow(x, 4) - 105 * pow(x, 2) + 1);
}

double
AssLegFunction::P_15_5_Deriv(const double x) const
{
  return (218243025.0 / 256.0) * x * sqrt(1 - pow(x, 2)) *
         (-30015 * pow(x, 12) + 90390 * pow(x, 10) - 102925 * pow(x, 8) +
          55108 * pow(x, 6) - 13993 * pow(x, 4) + 1478 * pow(x, 2) - 43);
}

double
AssLegFunction::P_15_6(const double x) const
{
  return (218243025.0 / 128.0) * x * pow(pow(x, 2) - 1, 3) *
         (-10005 * pow(x, 8) + 12420 * pow(x, 6) - 4830 * pow(x, 4) +
          644 * pow(x, 2) - 21);
}

double
AssLegFunction::P_15_6_Deriv(const double x) const
{
  return -32752821976875.0 / 128.0 * pow(x, 14) +
         (120394855956375.0 / 128.0) * pow(x, 12) -
         173100546493875.0 / 128.0 * pow(x, 10) +
         (122563318652775.0 / 128.0) * pow(x, 8) -
         44094039014025.0 / 128.0 * pow(x, 6) +
         (7447543228125.0 / 128.0) * pow(x, 4) -
         462893456025.0 / 128.0 * pow(x, 2) + 4583103525.0 / 128.0;
}

double
AssLegFunction::P_15_7(const double x) const
{
  return (654729075.0 / 128.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-30015 * pow(x, 8) + 28980 * pow(x, 6) - 8050 * pow(x, 4) +
          644 * pow(x, 2) - 7);
}

double
AssLegFunction::P_15_7_Deriv(const double x) const
{
  return (654729075.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (450225 * pow(x, 12) - 1517310 * pow(x, 10) + 1946375 * pow(x, 8) -
          1179716 * pow(x, 6) + 339759 * pow(x, 4) - 40670 * pow(x, 2) + 1337);
}

double
AssLegFunction::P_15_8(const double x) const
{
  return (15058768725.0 / 16.0) * x * pow(pow(x, 2) - 1, 4) *
         (1305 * pow(x, 6) - 945 * pow(x, 4) + 175 * pow(x, 2) - 7);
}

double
AssLegFunction::P_15_8_Deriv(const double x) const
{
  return (294775397791875.0 / 16.0) * pow(x, 14) -
         1206885019465125.0 / 16.0 * pow(x, 12) +
         (1952143483665375.0 / 16.0) * pow(x, 10) -
         1571728868134425.0 / 16.0 * pow(x, 8) +
         (649650341565225.0 / 16.0) * pow(x, 6) -
         127020714195375.0 / 16.0 * pow(x, 4) +
         (9170790153525.0 / 16.0) * pow(x, 2) - 105411381075.0 / 16.0;
}

double
AssLegFunction::P_15_9(const double x) const
{
  return (105411381075.0 / 16.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-1305 * pow(x, 6) + 675 * pow(x, 4) - 75 * pow(x, 2) + 1);
}

double
AssLegFunction::P_15_9_Deriv(const double x) const
{
  return (316234143225.0 / 16.0) * x * sqrt(1 - pow(x, 2)) *
         (-6525 * pow(x, 12) + 25110 * pow(x, 10) - 37355 * pow(x, 8) +
          26708 * pow(x, 6) - 9219 * pow(x, 4) + 1334 * pow(x, 2) - 53);
}

double
AssLegFunction::P_15_10(const double x) const
{
  return (1581170716125.0 / 8.0) * x * pow(pow(x, 2) - 1, 5) *
         (-261 * pow(x, 4) + 90 * pow(x, 2) - 5);
}

double
AssLegFunction::P_15_10_Deriv(const double x) const
{
  return -6190283353629375.0 / 8.0 * pow(x, 14) +
         (28674530936926875.0 / 8.0) * pow(x, 12) -
         53309170694154375.0 / 8.0 * pow(x, 10) +
         (50304946333516875.0 / 8.0) * pow(x, 8) -
         24958779754033125.0 / 8.0 * pow(x, 6) +
         (6016354574855625.0 / 8.0) * pow(x, 4) -
         545503897063125.0 / 8.0 * pow(x, 2) + 7905853580625.0 / 8.0;
}

double
AssLegFunction::P_15_11(const double x) const
{
  return (7905853580625.0 / 8.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-261 * pow(x, 4) + 54 * pow(x, 2) - 1);
}

double
AssLegFunction::P_15_11_Deriv(const double x) const
{
  return (7905853580625.0 / 8.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (2871 * pow(x, 4) - 594 * pow(x, 2) +
          36 * (1 - pow(x, 2)) * (3 - 29 * pow(x, 2)) + 11);
}

double
AssLegFunction::P_15_12(const double x) const
{
  return (71152682225625.0 / 2.0) * x * pow(pow(x, 2) - 1, 6) *
         (29 * pow(x, 2) - 3);
}

double
AssLegFunction::P_15_12_Deriv(const double x) const
{
  return (213458046676875.0 / 2.0) * pow(pow(x, 2) - 1, 5) *
         (4 * pow(x, 2) * (29 * pow(x, 2) - 3) +
          (pow(x, 2) - 1) * (29 * pow(x, 2) - 1));
}

double
AssLegFunction::P_15_13(const double x) const
{
  return (213458046676875.0 / 2.0) * (1 - 29 * pow(x, 2)) *
         pow(1 - pow(x, 2), 13.0 / 2.0);
}

double
AssLegFunction::P_15_13_Deriv(const double x) const
{
  return (213458046676875.0 / 2.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (435 * pow(x, 2) - 71);
}

double
AssLegFunction::P_15_14(const double x) const
{
  return -6190283353629375.0 * x * pow(pow(x, 2) - 1, 7);
}

double
AssLegFunction::P_15_14_Deriv(const double x) const
{
  return (6190283353629375.0 - 92854250304440625.0 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 6);
}

double
AssLegFunction::P_15_15(const double x) const
{
  return -6190283353629375.0 * pow(1 - pow(x, 2), 15.0 / 2.0);
}

double
AssLegFunction::P_15_15_Deriv(const double x) const
{
  return 92854250304440625.0 * x * pow(1 - pow(x, 2), 13.0 / 2.0);
}

double
AssLegFunction::P_16_0(const double x) const
{
  return (300540195.0 / 32768.0) * pow(x, 16) -
         145422675.0 / 4096.0 * pow(x, 14) +
         (456326325.0 / 8192.0) * pow(x, 12) -
         185910725.0 / 4096.0 * pow(x, 10) +
         (334639305.0 / 16384.0) * pow(x, 8) - 20369349.0 / 4096.0 * pow(x, 6) +
         (4849845.0 / 8192.0) * pow(x, 4) - 109395.0 / 4096.0 * pow(x, 2) +
         6435.0 / 32768.0;
}

double
AssLegFunction::P_16_0_Deriv(const double x) const
{
  return (17.0 / 2048.0) * x *
         (17678835 * pow(x, 14) - 59879925 * pow(x, 12) +
          80528175 * pow(x, 10) - 54679625 * pow(x, 8) + 19684665 * pow(x, 6) -
          3594591 * pow(x, 4) + 285285 * pow(x, 2) - 6435);
}

double
AssLegFunction::P_16_1(const double x) const
{
  return (17.0 / 2048.0) * x * sqrt(1 - pow(x, 2)) *
         (-17678835 * pow(x, 14) + 59879925 * pow(x, 12) -
          80528175 * pow(x, 10) + 54679625 * pow(x, 8) - 19684665 * pow(x, 6) +
          3594591 * pow(x, 4) - 285285 * pow(x, 2) + 6435);
}

double
AssLegFunction::P_16_1_Deriv(const double x) const
{
  return (17.0 / 2048.0) *
         (282861360.0 * pow(x, 16) - 1103501475.0 * pow(x, 14) +
          1744777125.0 * pow(x, 12) - 1432606175.0 * pow(x, 10) +
          649593945.0 * pow(x, 8) - 159360201.0 * pow(x, 6) +
          19114095 * pow(x, 4) - 868725 * pow(x, 2) + 6435) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_16_2(const double x) const
{
  return -4508102925.0 / 2048.0 * pow(x, 16) +
         (8870783175.0 / 1024.0) * pow(x, 14) -
         14146116075.0 / 1024.0 * pow(x, 12) +
         (11712375675.0 / 1024.0) * pow(x, 10) -
         334639305.0 / 64.0 * pow(x, 8) + (1324007685.0 / 1024.0) * pow(x, 6) -
         160044885.0 / 1024.0 * pow(x, 4) + (7329465.0 / 1024.0) * pow(x, 2) -
         109395.0 / 2048.0;
}

double
AssLegFunction::P_16_2_Deriv(const double x) const
{
  return (765.0 / 512.0) * x *
         (-23571780 * pow(x, 14) + 81170565 * pow(x, 12) -
          110949930.0 * pow(x, 10) + 76551475 * pow(x, 8) -
          27995968 * pow(x, 6) + 5192187 * pow(x, 4) - 418418 * pow(x, 2) +
          9581);
}

double
AssLegFunction::P_16_3(const double x) const
{
  return (101745.0 / 1024.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-310155 * pow(x, 12) + 780390 * pow(x, 10) - 740025 * pow(x, 8) +
          328900 * pow(x, 6) - 69069 * pow(x, 4) + 6006 * pow(x, 2) - 143);
}

double
AssLegFunction::P_16_3_Deriv(const double x) const
{
  return (101745.0 / 1024.0) * sqrt(1 - pow(x, 2)) *
         (4962480 * pow(x, 14) - 14957475 * pow(x, 12) + 17464590 * pow(x, 10) -
          9949225 * pow(x, 8) + 2854852 * pow(x, 6) - 381381 * pow(x, 4) +
          18590 * pow(x, 2) - 143);
}

double
AssLegFunction::P_16_4(const double x) const
{
  return (1322685.0 / 1024.0) * pow(pow(x, 2) - 1, 2) *
         (310155 * pow(x, 12) - 660330 * pow(x, 10) + 512325 * pow(x, 8) -
          177100 * pow(x, 6) + 26565 * pow(x, 4) - 1386 * pow(x, 2) + 11);
}

double
AssLegFunction::P_16_4_Deriv(const double x) const
{
  return (1322685.0 / 64.0) * x *
         (310155 * pow(x, 14) - 1120560 * pow(x, 12) + 1607355 * pow(x, 10) -
          1163800 * pow(x, 8) + 446545 * pow(x, 6) - 86856 * pow(x, 4) +
          7337 * pow(x, 2) - 176);
}

double
AssLegFunction::P_16_5(const double x) const
{
  return (3968055.0 / 256.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-310155 * pow(x, 10) + 550275 * pow(x, 8) - 341550 * pow(x, 6) +
          88550 * pow(x, 4) - 8855 * pow(x, 2) + 231);
}

double
AssLegFunction::P_16_5_Deriv(const double x) const
{
  return (3968055.0 / 256.0) * sqrt(1 - pow(x, 2)) *
         (-4962480 * pow(x, 14) + 16078035 * pow(x, 12) -
          20166630 * pow(x, 10) + 12327425 * pow(x, 8) - 3789940 * pow(x, 6) +
          541541 * pow(x, 4) - 28182 * pow(x, 2) + 231);
}

double
AssLegFunction::P_16_6(const double x) const
{
  return (43648605.0 / 256.0) * pow(pow(x, 2) - 1, 3) *
         (-310155 * pow(x, 10) + 450225 * pow(x, 8) - 217350 * pow(x, 6) +
          40250 * pow(x, 4) - 2415 * pow(x, 2) + 21);
}

double
AssLegFunction::P_16_6_Deriv(const double x) const
{
  return (43648605.0 / 64.0) * x *
         (-1240620 * pow(x, 14) + 4832415 * pow(x, 12) - 7495470 * pow(x, 10) +
          5882825 * pow(x, 8) - 2450880 * pow(x, 6) + 518049 * pow(x, 4) -
          47558 * pow(x, 2) + 1239);
}

double
AssLegFunction::P_16_7(const double x) const
{
  return (5019589575.0 / 128.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-13485 * pow(x, 8) + 15660 * pow(x, 6) - 5670 * pow(x, 4) +
          700 * pow(x, 2) - 21);
}

double
AssLegFunction::P_16_7_Deriv(const double x) const
{
  return (5019589575.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (215760 * pow(x, 14) - 772125 * pow(x, 12) + 1074630 * pow(x, 10) -
          731275 * pow(x, 8) + 250628 * pow(x, 6) - 39907 * pow(x, 4) +
          2310 * pow(x, 2) - 21);
}

double
AssLegFunction::P_16_8(const double x) const
{
  return (15058768725.0 / 128.0) * pow(pow(x, 2) - 1, 4) *
         (40455 * pow(x, 8) - 36540 * pow(x, 6) + 9450 * pow(x, 4) -
          700 * pow(x, 2) + 7);
}

double
AssLegFunction::P_16_8_Deriv(const double x) const
{
  return (15058768725.0 / 8.0) * x *
         (40455 * pow(x, 14) - 173565 * pow(x, 12) + 298755 * pow(x, 10) -
          262225 * pow(x, 8) + 123061 * pow(x, 6) - 29463 * pow(x, 4) +
          3073 * pow(x, 2) - 91);
}

double
AssLegFunction::P_16_9(const double x) const
{
  return (75293843625.0 / 16.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-8091 * pow(x, 6) + 5481 * pow(x, 4) - 945 * pow(x, 2) + 35);
}

double
AssLegFunction::P_16_9_Deriv(const double x) const
{
  return (75293843625.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (-129456 * pow(x, 14) + 521739 * pow(x, 12) - 827226 * pow(x, 10) +
          648989 * pow(x, 8) - 259196 * pow(x, 6) + 48405 * pow(x, 4) -
          3290 * pow(x, 2) + 35);
}

double
AssLegFunction::P_16_10(const double x) const
{
  return (527056905375.0 / 16.0) * pow(pow(x, 2) - 1, 5) *
         (-8091 * pow(x, 6) + 3915 * pow(x, 4) - 405 * pow(x, 2) + 5);
}

double
AssLegFunction::P_16_10_Deriv(const double x) const
{
  return (527056905375.0 / 4.0) * x *
         (-32364 * pow(x, 14) + 155295 * pow(x, 12) - 302670 * pow(x, 10) +
          305225 * pow(x, 8) - 167360 * pow(x, 6) + 47649 * pow(x, 4) -
          5990 * pow(x, 2) + 215);
}

double
AssLegFunction::P_16_11(const double x) const
{
  return (14230536445125.0 / 8.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-899 * pow(x, 4) + 290 * pow(x, 2) - 15);
}

double
AssLegFunction::P_16_11_Deriv(const double x) const
{
  return (14230536445125.0 / 8.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (11 * pow(x, 2) * (899 * pow(x, 4) - 290 * pow(x, 2) + 15) +
          5 * (1 - pow(x, 2)) * (-899 * pow(x, 4) + 174 * pow(x, 2) - 3));
}

double
AssLegFunction::P_16_12(const double x) const
{
  return (71152682225625.0 / 8.0) * pow(pow(x, 2) - 1, 6) *
         (899 * pow(x, 4) - 174 * pow(x, 2) + 3);
}

double
AssLegFunction::P_16_12_Deriv(const double x) const
{
  return 127932522641673750.0 * pow(x, 5) * pow(pow(x, 2) - 1, 5) -
         53649122398121250.0 * pow(x, 3) * pow(pow(x, 2) - 1, 5) +
         3415328746830000.0 * x * pow(pow(x, 2) - 1, 5);
}

double
AssLegFunction::P_16_13(const double x) const
{
  return (2063427784543125.0 / 2.0) * x * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (3 - 31 * pow(x, 2));
}

double
AssLegFunction::P_16_13_Deriv(const double x) const
{
  return (2063427784543125.0 / 2.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (13 * pow(x, 2) * (31 * pow(x, 2) - 3) +
          3 * (1 - 31 * pow(x, 2)) * (1 - pow(x, 2)));
}

double
AssLegFunction::P_16_14(const double x) const
{
  return (6190283353629375.0 / 2.0) * (1 - 31 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 7);
}

double
AssLegFunction::P_16_14_Deriv(const double x) const
{
  return 12380566707258750.0 * x * (19 - 124 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 6);
}

double
AssLegFunction::P_16_15(const double x) const
{
  return -191898783962510625.0 * x * pow(1 - pow(x, 2), 15.0 / 2.0);
}

double
AssLegFunction::P_16_15_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 13.0 / 2.0) *
         (3070380543400170000.0 * pow(x, 2) - 191898783962510625.0);
}

double
AssLegFunction::P_16_16(const double x) const
{
  return 191898783962510625.0 * pow(pow(x, 2) - 1, 8);
}

double
AssLegFunction::P_16_16_Deriv(const double x) const
{
  return 3070380543400170000.0 * x * pow(pow(x, 2) - 1, 7);
}

double
AssLegFunction::P_17_0(const double x) const
{
  return (1.0 / 32768.0) * x *
         (583401555.0 * pow(x, 16) - 2404321560.0 * pow(x, 14) +
          4071834900.0 * pow(x, 12) - 3650610600.0 * pow(x, 10) +
          1859107250.0 * pow(x, 8) - 535422888.0 * pow(x, 6) +
          81477396 * pow(x, 4) - 5542680 * pow(x, 2) + 109395);
}

double
AssLegFunction::P_17_0_Deriv(const double x) const
{
  return (9917826435.0 / 32768.0) * pow(x, 16) -
         4508102925.0 / 4096.0 * pow(x, 14) +
         (13233463425.0 / 8192.0) * pow(x, 12) -
         5019589575.0 / 4096.0 * pow(x, 10) +
         (8365982625.0 / 16384.0) * pow(x, 8) -
         468495027.0 / 4096.0 * pow(x, 6) + (101846745.0 / 8192.0) * pow(x, 4) -
         2078505.0 / 4096.0 * pow(x, 2) + 109395.0 / 32768.0;
}

double
AssLegFunction::P_17_1(const double x) const
{
  return (153.0 / 32768.0) * sqrt(1 - pow(x, 2)) *
         (-64822395 * pow(x, 16) + 235717800.0 * pow(x, 14) -
          345972900.0 * pow(x, 12) + 262462200.0 * pow(x, 10) -
          109359250.0 * pow(x, 8) + 24496472 * pow(x, 6) - 2662660 * pow(x, 4) +
          108680 * pow(x, 2) - 715);
}

double
AssLegFunction::P_17_1_Deriv(const double x) const
{
  return (1.0 / 32768.0) *
         (168603049395.0 * pow(x, 17) - 699657573960.0 * pow(x, 15) +
          1193047625700.0 * pow(x, 13) - 1076930127000.0 * pow(x, 11) +
          552154853250.0 * pow(x, 9) - 160091443512.0 * pow(x, 7) +
          24524696196.0 * pow(x, 5) - 1679432040.0 * pow(x, 3) + 33365475 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_17_2(const double x) const
{
  return (2907.0 / 2048.0) * x *
         (-3411705 * pow(x, 16) + 14267130 * pow(x, 14) -
          24512250 * pow(x, 12) + 22290450 * pow(x, 10) - 11511500 * pow(x, 8) +
          3361358 * pow(x, 6) - 518518 * pow(x, 4) + 35750 * pow(x, 2) - 715);
}

double
AssLegFunction::P_17_2_Deriv(const double x) const
{
  return -168603049395.0 / 2048.0 * pow(x, 16) +
         (311059101825.0 / 1024.0) * pow(x, 14) -
         463171219875.0 / 1024.0 * pow(x, 12) +
         (356390859825.0 / 1024.0) * pow(x, 10) -
         75293843625.0 / 512.0 * pow(x, 8) +
         (34200136971.0 / 1024.0) * pow(x, 6) -
         3768329565.0 / 1024.0 * pow(x, 4) +
         (155887875.0 / 1024.0) * pow(x, 2) - 2078505.0 / 2048.0;
}

double
AssLegFunction::P_17_3(const double x) const
{
  return (14535.0 / 2048.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-10235115 * pow(x, 14) + 28224105 * pow(x, 12) -
          30045015 * pow(x, 10) + 15540525 * pow(x, 8) - 4029025 * pow(x, 6) +
          483483 * pow(x, 4) - 21021 * pow(x, 2) + 143);
}

double
AssLegFunction::P_17_3_Deriv(const double x) const
{
  return (43605.0 / 2048.0) * x * sqrt(1 - pow(x, 2)) *
         (57998985 * pow(x, 14) - 188884395.0 * pow(x, 12) +
          243091485.0 * pow(x, 10) - 157131975.0 * pow(x, 8) +
          53528475 * pow(x, 6) - 9186177 * pow(x, 4) + 679679 * pow(x, 2) -
          14157);
}

double
AssLegFunction::P_17_4(const double x) const
{
  return (305235.0 / 1024.0) * x * pow(pow(x, 2) - 1, 2) *
         (3411705 * pow(x, 12) - 8064030 * pow(x, 10) + 7153575 * pow(x, 8) -
          2960100 * pow(x, 6) + 575575 * pow(x, 4) - 46046 * pow(x, 2) + 1001);
}

double
AssLegFunction::P_17_4_Deriv(const double x) const
{
  return (17703320186475.0 / 1024.0) * pow(x, 16) -
         4260157264125.0 / 64.0 * pow(x, 14) +
         (26480160313425.0 / 256.0) * pow(x, 12) -
         5315745359925.0 / 64.0 * pow(x, 10) +
         (18748167062625.0 / 512.0) * pow(x, 8) -
         555166606995.0 / 64.0 * pow(x, 6) +
         (255126096225.0 / 256.0) * pow(x, 4) -
         2749862115.0 / 64.0 * pow(x, 2) + 305540235.0 / 1024.0;
}

double
AssLegFunction::P_17_5(const double x) const
{
  return (43648605.0 / 1024.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-310155 * pow(x, 12) + 620310 * pow(x, 10) - 450225 * pow(x, 8) +
          144900 * pow(x, 6) - 20125 * pow(x, 4) + 966 * pow(x, 2) - 7);
}

double
AssLegFunction::P_17_5_Deriv(const double x) const
{
  return (43648605.0 / 1024.0) * x * sqrt(1 - pow(x, 2)) *
         (-5272635 * pow(x, 14) + 18299145 * pow(x, 12) -
          25082535 * pow(x, 10) + 17251725 * pow(x, 8) - 6246225 * pow(x, 6) +
          1137787 * pow(x, 4) - 89229 * pow(x, 2) + 1967);
}

double
AssLegFunction::P_17_6(const double x) const
{
  return (1003917915.0 / 256.0) * x * pow(pow(x, 2) - 1, 3) *
         (-40455 * pow(x, 10) + 67425 * pow(x, 8) - 39150 * pow(x, 6) +
          9450 * pow(x, 4) - 875 * pow(x, 2) + 21);
}

double
AssLegFunction::P_17_6_Deriv(const double x) const
{
  return -690429487272525.0 / 256.0 * pow(x, 16) +
         (1421472473796375.0 / 128.0) * pow(x, 14) -
         2367373972488525.0 / 128.0 * pow(x, 12) +
         (2040929984067975.0 / 128.0) * pow(x, 10) -
         483612357603375.0 / 64.0 * pow(x, 8) +
         (246472891229565.0 / 128.0) * pow(x, 6) -
         30463889130675.0 / 128.0 * pow(x, 4) +
         (1412512506405.0 / 128.0) * pow(x, 2) - 21082276215.0 / 256.0;
}

double
AssLegFunction::P_17_7(const double x) const
{
  return (3011753745.0 / 256.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-148335 * pow(x, 10) + 202275 * pow(x, 8) - 91350 * pow(x, 6) +
          15750 * pow(x, 4) - 875 * pow(x, 2) + 7);
}

double
AssLegFunction::P_17_7_Deriv(const double x) const
{
  return (3011753745.0 / 256.0) * x * sqrt(1 - pow(x, 2)) *
         (2521695 * pow(x, 14) - 9560865 * pow(x, 12) + 14362395 * pow(x, 10) -
          10850325 * pow(x, 8) + 4319325 * pow(x, 6) - 864899 * pow(x, 4) +
          74473 * pow(x, 2) - 1799);
}

double
AssLegFunction::P_17_8(const double x) const
{
  return (75293843625.0 / 128.0) * x * pow(pow(x, 2) - 1, 4) *
         (29667 * pow(x, 8) - 32364 * pow(x, 6) + 10962 * pow(x, 4) -
          1260 * pow(x, 2) + 35);
}

double
AssLegFunction::P_17_8_Deriv(const double x) const
{
  return (37973621799988875.0 / 128.0) * pow(x, 16) -
         21322087106945625.0 / 16.0 * pow(x, 14) +
         (77918963482985625.0 / 32.0) * pow(x, 12) -
         37059253363006875.0 / 16.0 * pow(x, 10) +
         (77918963482985625.0 / 64.0) * pow(x, 8) -
         5528299880478375.0 / 16.0 * pow(x, 6) +
         (1525829741060625.0 / 32.0) * pow(x, 4) -
         39529267903125.0 / 16.0 * pow(x, 2) + 2635284526875.0 / 128.0;
}

double
AssLegFunction::P_17_9(const double x) const
{
  return (75293843625.0 / 128.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-267003 * pow(x, 8) + 226548 * pow(x, 6) - 54810 * pow(x, 4) +
          3780 * pow(x, 2) - 35);
}

double
AssLegFunction::P_17_9_Deriv(const double x) const
{
  return (677644592625.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (-504339 * pow(x, 14) + 2127933 * pow(x, 12) - 3587967 * pow(x, 10) +
          3068673 * pow(x, 8) - 1393337 * pow(x, 6) + 319767 * pow(x, 4) -
          31605 * pow(x, 2) + 875);
}

double
AssLegFunction::P_17_10(const double x) const
{
  return (2032933777875.0 / 16.0) * x * pow(pow(x, 2) - 1, 5) *
         (-9889 * pow(x, 6) + 6293 * pow(x, 4) - 1015 * pow(x, 2) + 35);
}

double
AssLegFunction::P_17_10_Deriv(const double x) const
{
  return -341762596199899875.0 / 16.0 * pow(x, 16) +
         (849837471833975625.0 / 8.0) * pow(x, 14) -
         1735932317596351875.0 / 8.0 * pow(x, 12) +
         (1866466995473705625.0 / 8.0) * pow(x, 10) -
         561242192726840625.0 / 4.0 * pow(x, 8) +
         (368955118412755875.0 / 8.0) * pow(x, 6) -
         59554795022848125.0 / 8.0 * pow(x, 4) +
         (3628786793506875.0 / 8.0) * pow(x, 2) - 71152682225625.0 / 16.0;
}

double
AssLegFunction::P_17_11(const double x) const
{
  return (14230536445125.0 / 16.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-9889 * pow(x, 6) + 4495 * pow(x, 4) - 435 * pow(x, 2) + 5);
}

double
AssLegFunction::P_17_11_Deriv(const double x) const
{
  return (14230536445125.0 / 16.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (108779 * pow(x, 6) - 49445 * pow(x, 4) + 4785 * pow(x, 2) +
          58 * (1 - pow(x, 2)) * (-1023 * pow(x, 4) + 310 * pow(x, 2) - 15) -
          55);
}

double
AssLegFunction::P_17_12(const double x) const
{
  return (412685556908625.0 / 8.0) * x * pow(pow(x, 2) - 1, 6) *
         (1023 * pow(x, 4) - 310 * pow(x, 2) + 15);
}

double
AssLegFunction::P_17_12_Deriv(const double x) const
{
  return (7177014520197897375.0 / 8.0) * pow(x, 16) -
         4989368383025276250.0 * pow(x, 14) +
         (23095947192391198125.0 / 2.0) * pow(x, 12) -
         14299554546883856250.0 * pow(x, 10) +
         (40428740582553448125.0 / 4.0) * pow(x, 8) -
         4003875273127479750.0 * pow(x, 6) +
         (1603283388590008125.0 / 2.0) * pow(x, 4) -
         61902833536293750.0 * pow(x, 2) + 6190283353629375.0 / 8.0;
}

double
AssLegFunction::P_17_13(const double x) const
{
  return (6190283353629375.0 / 8.0) * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (-341 * pow(x, 4) + 62 * pow(x, 2) - 1);
}

double
AssLegFunction::P_17_13_Deriv(const double x) const
{
  return (6190283353629375.0 / 8.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (4433 * pow(x, 4) - 806 * pow(x, 2) +
          124 * (1 - 11 * pow(x, 2)) * (1 - pow(x, 2)) + 13);
}

double
AssLegFunction::P_17_14(const double x) const
{
  return (191898783962510625.0 / 2.0) * x * (1 - 11 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 7);
}

double
AssLegFunction::P_17_14_Deriv(const double x) const
{
  return (191898783962510625.0 / 2.0) * pow(pow(x, 2) - 1, 6) *
         (-14 * pow(x, 2) * (11 * pow(x, 2) - 1) +
          (1 - 33 * pow(x, 2)) * (pow(x, 2) - 1));
}

double
AssLegFunction::P_17_15(const double x) const
{
  return (191898783962510625.0 / 2.0) * (1 - 33 * pow(x, 2)) *
         pow(1 - pow(x, 2), 15.0 / 2.0);
}

double
AssLegFunction::P_17_15_Deriv(const double x) const
{
  return (575696351887531875.0 / 2.0) * x * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (187 * pow(x, 2) - 27);
}

double
AssLegFunction::P_17_16(const double x) const
{
  return 6332659870762850625.0 * x * pow(pow(x, 2) - 1, 8);
}

double
AssLegFunction::P_17_16_Deriv(const double x) const
{
  return pow(pow(x, 2) - 1, 7) *
         (107655217802968460625.0 * pow(x, 2) - 6332659870762850625.0);
}

double
AssLegFunction::P_17_17(const double x) const
{
  return -6332659870762850625.0 * pow(1 - pow(x, 2), 17.0 / 2.0);
}

double
AssLegFunction::P_17_17_Deriv(const double x) const
{
  return 107655217802968460625.0 * x * pow(1 - pow(x, 2), 15.0 / 2.0);
}

double
AssLegFunction::P_18_0(const double x) const
{
  return (2268783825.0 / 65536.0) * pow(x, 18) -
         9917826435.0 / 65536.0 * pow(x, 16) +
         (4508102925.0 / 16384.0) * pow(x, 14) -
         4411154475.0 / 16384.0 * pow(x, 12) +
         (5019589575.0 / 32768.0) * pow(x, 10) -
         1673196525.0 / 32768.0 * pow(x, 8) +
         (156165009.0 / 16384.0) * pow(x, 6) -
         14549535.0 / 16384.0 * pow(x, 4) + (2078505.0 / 65536.0) * pow(x, 2) -
         12155.0 / 65536.0;
}

double
AssLegFunction::P_18_0_Deriv(const double x) const
{
  return (171.0 / 32768.0) * x *
         (119409675.0 * pow(x, 16) - 463991880.0 * pow(x, 14) +
          738168900.0 * pow(x, 12) - 619109400.0 * pow(x, 10) +
          293543250.0 * pow(x, 8) - 78278200 * pow(x, 6) +
          10958948 * pow(x, 4) - 680680 * pow(x, 2) + 12155);
}

double
AssLegFunction::P_18_1(const double x) const
{
  return (171.0 / 32768.0) * x * sqrt(1 - pow(x, 2)) *
         (-119409675.0 * pow(x, 16) + 463991880.0 * pow(x, 14) -
          738168900.0 * pow(x, 12) + 619109400.0 * pow(x, 10) -
          293543250.0 * pow(x, 8) + 78278200 * pow(x, 6) -
          10958948 * pow(x, 4) + 680680 * pow(x, 2) - 12155);
}

double
AssLegFunction::P_18_1_Deriv(const double x) const
{
  return (171.0 / 32768.0) *
         (2149374150.0 * pow(x, 18) - 9453834555.0 * pow(x, 16) +
          17294242800.0 * pow(x, 14) - 17025508500.0 * pow(x, 12) +
          9745635900.0 * pow(x, 10) - 3268114850.0 * pow(x, 8) +
          613701088.0 * pow(x, 6) - 57517460 * pow(x, 4) + 2066350 * pow(x, 2) -
          12155) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_18_2(const double x) const
{
  return -347123925225.0 / 32768.0 * pow(x, 18) +
         (1537263097425.0 / 32768.0) * pow(x, 16) -
         707772159225.0 / 8192.0 * pow(x, 14) +
         (701373561525.0 / 8192.0) * pow(x, 12) -
         808153921575.0 / 16384.0 * pow(x, 10) +
         (272731033575.0 / 16384.0) * pow(x, 8) -
         25767226485.0 / 8192.0 * pow(x, 6) +
         (2429772345.0 / 8192.0) * pow(x, 4) -
         351267345.0 / 32768.0 * pow(x, 2) + 2078505.0 / 32768.0;
}

double
AssLegFunction::P_18_2_Deriv(const double x) const
{
  return (14535.0 / 16384.0) * x *
         (-214937415.0 * pow(x, 16) + 846102840.0 * pow(x, 14) -
          1363441380.0 * pow(x, 12) + 1158098760.0 * pow(x, 10) -
          556005450.0 * pow(x, 8) + 150109960.0 * pow(x, 6) -
          21273252 * pow(x, 4) + 1337336 * pow(x, 2) - 24167);
}

double
AssLegFunction::P_18_3(const double x) const
{
  return (101745.0 / 2048.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-3411705 * pow(x, 14) + 10235115 * pow(x, 12) -
          12096045 * pow(x, 10) + 7153575 * pow(x, 8) - 2220075 * pow(x, 6) +
          345345 * pow(x, 4) - 23023 * pow(x, 2) + 429);
}

double
AssLegFunction::P_18_3_Deriv(const double x) const
{
  return (305235.0 / 2048.0) * sqrt(1 - pow(x, 2)) *
         (20470230 * pow(x, 16) - 71645805 * pow(x, 14) +
          100800375.0 * pow(x, 12) - 72966465 * pow(x, 10) +
          28860975 * pow(x, 8) - 6101095 * pow(x, 6) + 621621 * pow(x, 4) -
          23595 * pow(x, 2) + 143);
}

double
AssLegFunction::P_18_4(const double x) const
{
  return (3357585.0 / 2048.0) * pow(pow(x, 2) - 1, 2) *
         (1550775 * pow(x, 14) - 4032015 * pow(x, 12) + 4032015 * pow(x, 10) -
          1950975 * pow(x, 8) + 470925 * pow(x, 6) - 52325 * pow(x, 4) +
          2093 * pow(x, 2) - 13);
}

double
AssLegFunction::P_18_4_Deriv(const double x) const
{
  return (3357585.0 / 1024.0) * x *
         (13956975 * pow(x, 16) - 57068520 * pow(x, 14) +
          95527740 * pow(x, 12) - 84282120 * pow(x, 10) + 42024450 * pow(x, 8) -
          11780600 * pow(x, 6) + 1733004 * pow(x, 4) - 113048 * pow(x, 2) +
          2119);
}

double
AssLegFunction::P_18_5(const double x) const
{
  return (77224455.0 / 1024.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-471975 * pow(x, 12) + 1051830 * pow(x, 10) - 876525 * pow(x, 8) +
          339300 * pow(x, 6) - 61425 * pow(x, 4) + 4550 * pow(x, 2) - 91);
}

double
AssLegFunction::P_18_5_Deriv(const double x) const
{
  return (77224455.0 / 1024.0) * sqrt(1 - pow(x, 2)) *
         (-8495550 * pow(x, 16) + 31460505 * pow(x, 14) -
          46806435 * pow(x, 12) + 35801805 * pow(x, 10) - 14949675 * pow(x, 8) +
          3332875 * pow(x, 6) - 357721 * pow(x, 4) + 14287 * pow(x, 2) - 91);
}

double
AssLegFunction::P_18_6(const double x) const
{
  return (1003917915.0 / 1024.0) * pow(pow(x, 2) - 1, 3) *
         (-471975 * pow(x, 12) + 890010 * pow(x, 10) - 606825 * pow(x, 8) +
          182700 * pow(x, 6) - 23625 * pow(x, 4) + 1050 * pow(x, 2) - 7);
}

double
AssLegFunction::P_18_6_Deriv(const double x) const
{
  return (3011753745.0 / 512.0) * x *
         (-1415925 * pow(x, 16) + 6149160 * pow(x, 14) - 10949820 * pow(x, 12) +
          10290360 * pow(x, 10) - 5470350 * pow(x, 8) + 1635800 * pow(x, 6) -
          256732 * pow(x, 4) + 17864 * pow(x, 2) - 357);
}

double
AssLegFunction::P_18_7(const double x) const
{
  return (75293843625.0 / 256.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-18879 * pow(x, 10) + 29667 * pow(x, 8) - 16182 * pow(x, 6) +
          3654 * pow(x, 4) - 315 * pow(x, 2) + 7);
}

double
AssLegFunction::P_18_7_Deriv(const double x) const
{
  return (75293843625.0 / 256.0) * sqrt(1 - pow(x, 2)) *
         (339822 * pow(x, 16) - 1361985 * pow(x, 14) + 2198055 * pow(x, 12) -
          1826565 * pow(x, 10) + 829215 * pow(x, 8) - 200963 * pow(x, 6) +
          23429 * pow(x, 4) - 1015 * pow(x, 2) + 7);
}

double
AssLegFunction::P_18_8(const double x) const
{
  return (75293843625.0 / 256.0) * pow(pow(x, 2) - 1, 4) *
         (207669 * pow(x, 10) - 267003 * pow(x, 8) + 113274 * pow(x, 6) -
          18270 * pow(x, 4) + 945 * pow(x, 2) - 7);
}

double
AssLegFunction::P_18_8_Deriv(const double x) const
{
  return (75293843625.0 / 128.0) * x *
         (1869021 * pow(x, 16) - 8781432 * pow(x, 14) + 16991100 * pow(x, 12) -
          17424360 * pow(x, 10) + 10146750 * pow(x, 8) - 3334024 * pow(x, 6) +
          576156 * pow(x, 4) - 44184 * pow(x, 2) + 973);
}

double
AssLegFunction::P_18_9(const double x) const
{
  return (225881530875.0 / 128.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-346115 * pow(x, 8) + 356004 * pow(x, 6) - 113274 * pow(x, 4) +
          12180 * pow(x, 2) - 315);
}

double
AssLegFunction::P_18_9_Deriv(const double x) const
{
  return (2032933777875.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (-692230 * pow(x, 16) + 3055701 * pow(x, 14) - 5466819 * pow(x, 12) +
          5067721 * pow(x, 10) - 2580219 * pow(x, 8) + 703871 * pow(x, 6) -
          92505 * pow(x, 4) + 4515 * pow(x, 2) - 35);
}

double
AssLegFunction::P_18_10(const double x) const
{
  return (14230536445125.0 / 128.0) * pow(pow(x, 2) - 1, 5) *
         (-49445 * pow(x, 8) + 39556 * pow(x, 6) - 8990 * pow(x, 4) +
          580 * pow(x, 2) - 5);
}

double
AssLegFunction::P_18_10_Deriv(const double x) const
{
  return (14230536445125.0 / 64.0) * x *
         (-445005 * pow(x, 16) + 2294248 * pow(x, 14) - 4908540 * pow(x, 12) +
          5613240 * pow(x, 10) - 3677950 * pow(x, 8) + 1371800 * pow(x, 6) -
          271068 * pow(x, 4) + 23880 * pow(x, 2) - 605);
}

double
AssLegFunction::P_18_11(const double x) const
{
  return (412685556908625.0 / 16.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-1705 * pow(x, 6) + 1023 * pow(x, 4) - 155 * pow(x, 2) + 5);
}

double
AssLegFunction::P_18_11_Deriv(const double x) const
{
  return (412685556908625.0 / 16.0) * sqrt(1 - pow(x, 2)) *
         (30690 * pow(x, 16) - 151063 * pow(x, 14) + 304637 * pow(x, 12) -
          322243 * pow(x, 10) + 189717 * pow(x, 8) - 60613 * pow(x, 6) +
          9415 * pow(x, 4) - 545 * pow(x, 2) + 5);
}

double
AssLegFunction::P_18_12(const double x) const
{
  return (2063427784543125.0 / 16.0) * pow(pow(x, 2) - 1, 6) *
         (2387 * pow(x, 6) - 1023 * pow(x, 4) + 93 * pow(x, 2) - 1);
}

double
AssLegFunction::P_18_12_Deriv(const double x) const
{
  return (6190283353629375.0 / 8.0) * x * pow(pow(x, 2) - 1, 5) *
         (4774 * pow(x, 6) - 2046 * pow(x, 4) + 186 * pow(x, 2) +
          31 * (pow(x, 2) - 1) * (77 * pow(x, 4) - 22 * pow(x, 2) + 1) - 2);
}

double
AssLegFunction::P_18_13(const double x) const
{
  return (191898783962510625.0 / 8.0) * x * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (-77 * pow(x, 4) + 22 * pow(x, 2) - 1);
}

double
AssLegFunction::P_18_13_Deriv(const double x) const
{
  return (191898783962510625.0 / 8.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (13 * pow(x, 2) * (77 * pow(x, 4) - 22 * pow(x, 2) + 1) +
          (1 - pow(x, 2)) * (-385 * pow(x, 4) + 66 * pow(x, 2) - 1));
}

double
AssLegFunction::P_18_14(const double x) const
{
  return (191898783962510625.0 / 8.0) * pow(pow(x, 2) - 1, 7) *
         (-385 * pow(x, 4) + 66 * pow(x, 2) - 1);
}

double
AssLegFunction::P_18_14_Deriv(const double x) const
{
  return (191898783962510625.0 / 4.0) * x * pow(pow(x, 2) - 1, 6) *
         (-2695 * pow(x, 4) + 462 * pow(x, 2) +
          22 * (3 - 35 * pow(x, 2)) * (pow(x, 2) - 1) - 7);
}

double
AssLegFunction::P_18_15(const double x) const
{
  return (2110886623587616875.0 / 2.0) * x * pow(1 - pow(x, 2), 15.0 / 2.0) *
         (3 - 35 * pow(x, 2));
}

double
AssLegFunction::P_18_15_Deriv(const double x) const
{
  return (6332659870762850625.0 / 2.0) * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (5 * pow(x, 2) * (35 * pow(x, 2) - 3) +
          (1 - 35 * pow(x, 2)) * (1 - pow(x, 2)));
}

double
AssLegFunction::P_18_16(const double x) const
{
  return (6332659870762850625.0 / 2.0) * pow(pow(x, 2) - 1, 8) *
         (35 * pow(x, 2) - 1);
}

double
AssLegFunction::P_18_16_Deriv(const double x) const
{
  return 6332659870762850625.0 * x * pow(pow(x, 2) - 1, 7) *
         (315 * pow(x, 2) - 43);
}

double
AssLegFunction::P_18_17(const double x) const
{
  return -221643095476699771875.0 * x * pow(1 - pow(x, 2), 17.0 / 2.0);
}

double
AssLegFunction::P_18_17_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 15.0 / 2.0) *
         (3989575718580595893750.0 * pow(x, 2) - 221643095476699771875.0);
}

double
AssLegFunction::P_18_18(const double x) const
{
  return -221643095476699771875.0 * pow(pow(x, 2) - 1, 9);
}

double
AssLegFunction::P_18_18_Deriv(const double x) const
{
  return -3989575718580595893750.0 * x * pow(pow(x, 2) - 1, 8);
}

double
AssLegFunction::P_19_0(const double x) const
{
  return (1.0 / 65536.0) * x *
         (4418157975.0 * pow(x, 18) - 20419054425.0 * pow(x, 16) +
          39671305740.0 * pow(x, 14) - 42075627300.0 * pow(x, 12) +
          26466926850.0 * pow(x, 10) - 10039179150.0 * pow(x, 8) +
          2230928700.0 * pow(x, 6) - 267711444.0 * pow(x, 4) +
          14549535 * pow(x, 2) - 230945);
}

double
AssLegFunction::P_19_0_Deriv(const double x) const
{
  return (83945001525.0 / 65536.0) * pow(x, 18) -
         347123925225.0 / 65536.0 * pow(x, 16) +
         (148767396525.0 / 16384.0) * pow(x, 14) -
         136745788725.0 / 16384.0 * pow(x, 12) +
         (145568097675.0 / 32768.0) * pow(x, 10) -
         45176306175.0 / 32768.0 * pow(x, 8) +
         (3904125225.0 / 16384.0) * pow(x, 6) -
         334639305.0 / 16384.0 * pow(x, 4) +
         (43648605.0 / 65536.0) * pow(x, 2) - 230945.0 / 65536.0;
}

double
AssLegFunction::P_19_1(const double x) const
{
  return (95.0 / 65536.0) * sqrt(1 - pow(x, 2)) *
         (-883631595.0 * pow(x, 18) + 3653936055.0 * pow(x, 16) -
          6263890380.0 * pow(x, 14) + 5757717420.0 * pow(x, 12) -
          3064591530.0 * pow(x, 10) + 951080130.0 * pow(x, 8) -
          164384220.0 * pow(x, 6) + 14090076 * pow(x, 4) - 459459 * pow(x, 2) +
          2431);
}

double
AssLegFunction::P_19_1_Deriv(const double x) const
{
  return (1.0 / 65536.0) *
         (1594955028975.0 * pow(x, 19) - 7412116756275.0 * pow(x, 17) +
          14480026595100.0 * pow(x, 15) - 15441755219100.0 * pow(x, 13) +
          9766296007650.0 * pow(x, 11) - 3724535464650.0 * pow(x, 9) +
          832136405100.0 * pow(x, 7) - 100391791500.0 * pow(x, 5) +
          5485174695.0 * pow(x, 3) - 87528155 * x) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_19_2(const double x) const
{
  return (5985.0 / 32768.0) * x *
         (-126233085.0 * pow(x, 18) + 590224965.0 * pow(x, 16) -
          1159979700.0 * pow(x, 14) + 1244341860.0 * pow(x, 12) -
          791575590.0 * pow(x, 10) + 303607590.0 * pow(x, 8) -
          68213860 * pow(x, 6) + 8275124 * pow(x, 4) - 454597 * pow(x, 2) +
          7293);
}

double
AssLegFunction::P_19_2_Deriv(const double x) const
{
  return -14354595260775.0 / 32768.0 * pow(x, 18) +
         (60052439063925.0 / 32768.0) * pow(x, 16) -
         26034294391875.0 / 8192.0 * pow(x, 14) +
         (24204004604325.0 / 8192.0) * pow(x, 12) -
         26056689483825.0 / 16384.0 * pow(x, 10) +
         (8176911417675.0 / 16384.0) * pow(x, 8) -
         714454916175.0 / 8192.0 * pow(x, 6) +
         (61908271425.0 / 8192.0) * pow(x, 4) -
         8162289135.0 / 32768.0 * pow(x, 2) + 43648605.0 / 32768.0;
}

double
AssLegFunction::P_19_3(const double x) const
{
  return (1119195.0 / 32768.0) * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-11475735 * pow(x, 16) + 37218600 * pow(x, 14) -
          48384180 * pow(x, 12) + 32256120 * pow(x, 10) - 11705850 * pow(x, 8) +
          2260440 * pow(x, 6) - 209300 * pow(x, 4) + 7176 * pow(x, 2) - 39);
}

double
AssLegFunction::P_19_3_Deriv(const double x) const
{
  return (1119195.0 / 32768.0) * x * sqrt(1 - pow(x, 2)) *
         (218038965.0 * pow(x, 16) - 816327960.0 * pow(x, 14) +
          1246823100.0 * pow(x, 12) - 999939720.0 * pow(x, 10) +
          451325550.0 * pow(x, 8) - 113990760.0 * pow(x, 6) +
          15027740 * pow(x, 4) - 873080 * pow(x, 2) + 14469);
}

double
AssLegFunction::P_19_4(const double x) const
{
  return (25741485.0 / 2048.0) * x * pow(pow(x, 2) - 1, 2) *
         (498945 * pow(x, 14) - 1415925 * pow(x, 12) + 1577745 * pow(x, 10) -
          876525 * pow(x, 8) + 254475 * pow(x, 6) - 36855 * pow(x, 4) +
          2275 * pow(x, 2) - 39);
}

double
AssLegFunction::P_19_4_Deriv(const double x) const
{
  return (244028119433175.0 / 2048.0) * pow(x, 18) -
         1056298104459675.0 / 2048.0 * pow(x, 16) +
         (473824157932125.0 / 512.0) * pow(x, 14) -
         455773713820425.0 / 512.0 * pow(x, 12) +
         (507595956592725.0 / 1024.0) * pow(x, 10) -
         164757988620225.0 / 1024.0 * pow(x, 8) +
         (14886429482925.0 / 512.0) * pow(x, 6) -
         1333537630425.0 / 512.0 * pow(x, 4) +
         (181709142615.0 / 2048.0) * pow(x, 2) - 1003917915.0 / 2048.0;
}

double
AssLegFunction::P_19_5(const double x) const
{
  return (77224455.0 / 2048.0) * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-2494725 * pow(x, 14) + 6135675 * pow(x, 12) - 5785065 * pow(x, 10) +
          2629575 * pow(x, 8) - 593775 * pow(x, 6) + 61425 * pow(x, 4) -
          2275 * pow(x, 2) + 13);
}

double
AssLegFunction::P_19_5_Deriv(const double x) const
{
  return (386122275.0 / 2048.0) * x * sqrt(1 - pow(x, 2)) *
         (-9479955 * pow(x, 16) + 37326480 * pow(x, 14) -
          59927340 * pow(x, 12) + 50487840 * pow(x, 10) - 23920650 * pow(x, 8) +
          6336720 * pow(x, 6) - 875420 * pow(x, 4) + 53248 * pow(x, 2) - 923);
}

double
AssLegFunction::P_19_6(const double x) const
{
  return (1930611375.0 / 1024.0) * x * pow(pow(x, 2) - 1, 3) *
         (-698523 * pow(x, 12) + 1472562 * pow(x, 10) - 1157013 * pow(x, 8) +
          420732 * pow(x, 6) - 71253 * pow(x, 4) + 4914 * pow(x, 2) - 91);
}

double
AssLegFunction::P_19_6_Deriv(const double x) const
{
  return -25622952540483375.0 / 1024.0 * pow(x, 18) +
         (117107463033532125.0 / 1024.0) * pow(x, 16) -
         55531149937869375.0 / 256.0 * pow(x, 14) +
         (56520453124760625.0 / 256.0) * pow(x, 12) -
         66651992722940625.0 / 512.0 * pow(x, 10) +
         (22917262477984875.0 / 512.0) * pow(x, 8) -
         2193786525805875.0 / 256.0 * pow(x, 6) +
         (208187477623125.0 / 256.0) * pow(x, 4) -
         30042243606375.0 / 1024.0 * pow(x, 2) + 175685635125.0 / 1024.0;
}

double
AssLegFunction::P_19_7(const double x) const
{
  return (25097947875.0 / 1024.0) * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-698523 * pow(x, 12) + 1246014 * pow(x, 10) - 801009 * pow(x, 8) +
          226548 * pow(x, 6) - 27405 * pow(x, 4) + 1134 * pow(x, 2) - 7);
}

double
AssLegFunction::P_19_7_Deriv(const double x) const
{
  return (25097947875.0 / 1024.0) * x * sqrt(1 - pow(x, 2)) *
         (13271937 * pow(x, 16) - 56108388 * pow(x, 14) +
          96876240 * pow(x, 12) - 87868260 * pow(x, 10) + 44842410 * pow(x, 8) -
          12794508 * pow(x, 6) + 1902712 * pow(x, 4) - 124460 * pow(x, 2) +
          2317);
}

double
AssLegFunction::P_19_8(const double x) const
{
  return (225881530875.0 / 256.0) * x * pow(pow(x, 2) - 1, 4) *
         (232841 * pow(x, 10) - 346115 * pow(x, 8) + 178002 * pow(x, 6) -
          37758 * pow(x, 4) + 3045 * pow(x, 2) - 63);
}

double
AssLegFunction::P_19_8_Deriv(const double x) const
{
  return (999295149078851625.0 / 256.0) * pow(x, 18) -
         4905501507071290125.0 / 256.0 * pow(x, 16) +
         (2506868241288035625.0 / 64.0) * pow(x, 14) -
         2758671936646250625.0 / 64.0 * pow(x, 12) +
         (3527511679731414375.0 / 128.0) * pow(x, 10) -
         1318270138799488875.0 / 128.0 * pow(x, 8) +
         (137386342353385125.0 / 64.0) * pow(x, 6) -
         14206818884383125.0 / 64.0 * pow(x, 4) +
         (2234194221884625.0 / 256.0) * pow(x, 2) - 14230536445125.0 / 256.0;
}

double
AssLegFunction::P_19_9(const double x) const
{
  return (1581170716125.0 / 256.0) * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-365893 * pow(x, 10) + 445005 * pow(x, 8) - 178002 * pow(x, 6) +
          26970 * pow(x, 4) - 1305 * pow(x, 2) + 9);
}

double
AssLegFunction::P_19_9_Deriv(const double x) const
{
  return (1581170716125.0 / 256.0) * x * sqrt(1 - pow(x, 2)) *
         (-6951967 * pow(x, 16) + 32079916 * pow(x, 14) -
          60758016 * pow(x, 12) + 60732844 * pow(x, 10) - 34292326 * pow(x, 8) +
          10855332 * pow(x, 6) - 1793400 * pow(x, 4) + 130308 * pow(x, 2) -
          2691);
}

double
AssLegFunction::P_19_10(const double x) const
{
  return (45853950767625.0 / 128.0) * x * pow(pow(x, 2) - 1, 5) *
         (-63085 * pow(x, 8) + 61380 * pow(x, 6) - 18414 * pow(x, 4) +
          1860 * pow(x, 2) - 45);
}

double
AssLegFunction::P_19_10_Deriv(const double x) const
{
  return -54961233199336839375.0 / 128.0 * pow(x, 18) +
         (293725964622913948125.0 / 128.0) * pow(x, 16) -
         164414613681657714375.0 / 32.0 * pow(x, 14) +
         (199482339610214285625.0 / 32.0) * pow(x, 12) -
         283143789864761450625.0 / 64.0 * pow(x, 10) +
         (118224094915398346875.0 / 64.0) * pow(x, 8) -
         13842161387976796875.0 / 32.0 * pow(x, 6) +
         (1614288336774238125.0 / 32.0) * pow(x, 4) -
         286816462051494375.0 / 128.0 * pow(x, 2) + 2063427784543125.0 / 128.0;
}

double
AssLegFunction::P_19_11(const double x) const
{
  return (2063427784543125.0 / 128.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-12617 * pow(x, 8) + 9548 * pow(x, 6) - 2046 * pow(x, 4) +
          124 * pow(x, 2) - 1);
}

double
AssLegFunction::P_19_11_Deriv(const double x) const
{
  return (2063427784543125.0 / 128.0) * x * sqrt(1 - pow(x, 2)) *
         (239723 * pow(x, 16) - 1222144 * pow(x, 14) + 2579324 * pow(x, 12) -
          2900112 * pow(x, 10) + 1860042 * pow(x, 8) - 674976 * pow(x, 6) +
          128716 * pow(x, 4) - 10832 * pow(x, 2) + 259);
}

double
AssLegFunction::P_19_12(const double x) const
{
  return (63966261320836875.0 / 16.0) * x * pow(pow(x, 2) - 1, 6) *
         (407 * pow(x, 6) - 231 * pow(x, 4) + 33 * pow(x, 2) - 1);
}

double
AssLegFunction::P_19_12_Deriv(const double x) const
{
  return (494651098794031554375.0 / 16.0) * pow(x, 18) -
         2906690880680148436875.0 / 16.0 * pow(x, 16) +
         (1804808063167412428125.0 / 4.0) * pow(x, 14) -
         2453937683051265035625.0 / 4.0 * pow(x, 12) +
         (3949468872732431173125.0 / 8.0) * pow(x, 10) -
         1894616694061867400625.0 / 8.0 * pow(x, 8) +
         (258359729474860138125.0 / 4.0) * pow(x, 6) -
         35501275033064465625.0 / 4.0 * pow(x, 4) +
         (7484052574537914375.0 / 16.0) * pow(x, 2) -
         63966261320836875.0 / 16.0;
}

double
AssLegFunction::P_19_13(const double x) const
{
  return (63966261320836875.0 / 16.0) * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (-2849 * pow(x, 6) + 1155 * pow(x, 4) - 99 * pow(x, 2) + 1);
}

double
AssLegFunction::P_19_13_Deriv(const double x) const
{
  return (63966261320836875.0 / 16.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (37037 * pow(x, 6) - 15015 * pow(x, 4) + 1287 * pow(x, 2) +
          66 * (1 - pow(x, 2)) * (-259 * pow(x, 4) + 70 * pow(x, 2) - 3) - 13);
}

double
AssLegFunction::P_19_14(const double x) const
{
  return (2110886623587616875.0 / 8.0) * x * pow(pow(x, 2) - 1, 7) *
         (-259 * pow(x, 4) + 70 * pow(x, 2) - 3);
}

double
AssLegFunction::P_19_14_Deriv(const double x) const
{
  return (2110886623587616875.0 / 8.0) * pow(pow(x, 2) - 1, 6) *
         (-14 * pow(x, 2) * (259 * pow(x, 4) - 70 * pow(x, 2) + 3) +
          (pow(x, 2) - 1) * (-1295 * pow(x, 4) + 210 * pow(x, 2) - 3));
}

double
AssLegFunction::P_19_15(const double x) const
{
  return (2110886623587616875.0 / 8.0) * pow(1 - pow(x, 2), 15.0 / 2.0) *
         (-1295 * pow(x, 4) + 210 * pow(x, 2) - 3);
}

double
AssLegFunction::P_19_15_Deriv(const double x) const
{
  return (10554433117938084375.0 / 8.0) * x * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (3885 * pow(x, 4) - 630 * pow(x, 2) +
          28 * (1 - pow(x, 2)) * (3 - 37 * pow(x, 2)) + 9);
}

double
AssLegFunction::P_19_16(const double x) const
{
  return (73881031825566590625.0 / 2.0) * x * pow(pow(x, 2) - 1, 8) *
         (37 * pow(x, 2) - 3);
}

double
AssLegFunction::P_19_16_Deriv(const double x) const
{
  return (73881031825566590625.0 / 2.0) * pow(pow(x, 2) - 1, 7) *
         (16 * pow(x, 2) * (37 * pow(x, 2) - 3) +
          3 * (pow(x, 2) - 1) * (37 * pow(x, 2) - 1));
}

double
AssLegFunction::P_19_17(const double x) const
{
  return (221643095476699771875.0 / 2.0) * (1 - 37 * pow(x, 2)) *
         pow(1 - pow(x, 2), 17.0 / 2.0);
}

double
AssLegFunction::P_19_17_Deriv(const double x) const
{
  return (221643095476699771875.0 / 2.0) * x * pow(1 - pow(x, 2), 15.0 / 2.0) *
         (703 * pow(x, 2) - 91);
}

double
AssLegFunction::P_19_18(const double x) const
{
  return -8200794532637891559375.0 * x * pow(pow(x, 2) - 1, 9);
}

double
AssLegFunction::P_19_18_Deriv(const double x) const
{
  return (8200794532637891559375.0 - 155815096120119939628125.0 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 8);
}

double
AssLegFunction::P_19_19(const double x) const
{
  return -8200794532637891559375.0 * pow(1 - pow(x, 2), 19.0 / 2.0);
}

double
AssLegFunction::P_19_19_Deriv(const double x) const
{
  return 155815096120119939628125.0 * x * pow(1 - pow(x, 2), 17.0 / 2.0);
}

double
AssLegFunction::P_20_0(const double x) const
{
  return (34461632205.0 / 262144.0) * pow(x, 20) -
         83945001525.0 / 131072.0 * pow(x, 18) +
         (347123925225.0 / 262144.0) * pow(x, 16) -
         49589132175.0 / 32768.0 * pow(x, 14) +
         (136745788725.0 / 131072.0) * pow(x, 12) -
         29113619535.0 / 65536.0 * pow(x, 10) +
         (15058768725.0 / 131072.0) * pow(x, 8) -
         557732175.0 / 32768.0 * pow(x, 6) +
         (334639305.0 / 262144.0) * pow(x, 4) -
         4849845.0 / 131072.0 * pow(x, 2) + 46189.0 / 262144.0;
}

double
AssLegFunction::P_20_0_Deriv(const double x) const
{
  return (105.0 / 65536.0) * x *
         (1641030105.0 * pow(x, 18) - 7195285845.0 * pow(x, 16) +
          13223768580.0 * pow(x, 14) - 13223768580.0 * pow(x, 12) +
          7814045070.0 * pow(x, 10) - 2772725670.0 * pow(x, 8) +
          573667380.0 * pow(x, 6) - 63740820 * pow(x, 4) + 3187041 * pow(x, 2) -
          46189);
}

double
AssLegFunction::P_20_1(const double x) const
{
  return (105.0 / 65536.0) * x * sqrt(1 - pow(x, 2)) *
         (-1641030105.0 * pow(x, 18) + 7195285845.0 * pow(x, 16) -
          13223768580.0 * pow(x, 14) + 13223768580.0 * pow(x, 12) -
          7814045070.0 * pow(x, 10) + 2772725670.0 * pow(x, 8) -
          573667380.0 * pow(x, 6) + 63740820 * pow(x, 4) - 3187041 * pow(x, 2) +
          46189);
}

double
AssLegFunction::P_20_1_Deriv(const double x) const
{
  return (105.0 / 65536.0) *
         (32820602100.0 * pow(x, 20) - 160694717205.0 * pow(x, 18) +
          333900156645.0 * pow(x, 16) - 383489288820.0 * pow(x, 14) +
          265677532380.0 * pow(x, 12) - 113681752470.0 * pow(x, 10) +
          29543870070.0 * pow(x, 8) - 4398116580.0 * pow(x, 6) +
          331452264.0 * pow(x, 4) - 9653501 * pow(x, 2) + 46189) /
         sqrt(1 - pow(x, 2));
}

double
AssLegFunction::P_20_2(const double x) const
{
  return -3273855059475.0 / 65536.0 * pow(x, 20) +
         (251835004575.0 / 1024.0) * pow(x, 18) -
         33671020746825.0 / 65536.0 * pow(x, 16) +
         (2429867476575.0 / 4096.0) * pow(x, 14) -
         13537833083775.0 / 32768.0 * pow(x, 12) +
         (727840488375.0 / 4096.0) * pow(x, 10) -
         1520935641225.0 / 32768.0 * pow(x, 8) +
         (28444340925.0 / 4096.0) * pow(x, 6) -
         34467848415.0 / 65536.0 * pow(x, 4) +
         (63047985.0 / 4096.0) * pow(x, 2) - 4849845.0 / 65536.0;
}

double
AssLegFunction::P_20_2_Deriv(const double x) const
{
  return (21945.0 / 16384.0) * x *
         (-745922775.0 * pow(x, 18) + 3305011680.0 * pow(x, 16) -
          6137347140.0 * pow(x, 14) + 6200618760.0 * pow(x, 12) -
          3701389770.0 * pow(x, 10) + 1326663000.0 * pow(x, 8) -
          277226820.0 * pow(x, 6) + 31107960 * pow(x, 4) - 1570647 * pow(x, 2) +
          22984);
}

double
AssLegFunction::P_20_3(const double x) const
{
  return (1514205.0 / 32768.0) * x * pow(1 - pow(x, 2), 3.0 / 2.0) *
         (-19458855 * pow(x, 16) + 67856520 * pow(x, 14) -
          96282900 * pow(x, 12) + 71524440 * pow(x, 10) - 29801850 * pow(x, 8) +
          6921720 * pow(x, 6) - 835380 * pow(x, 4) + 44200 * pow(x, 2) - 663);
}

double
AssLegFunction::P_20_3_Deriv(const double x) const
{
  return (4542615.0 / 32768.0) * sqrt(1 - pow(x, 2)) *
         (129725700.0 * pow(x, 18) - 517405965.0 * pow(x, 16) +
          852791400.0 * pow(x, 14) - 751006620.0 * pow(x, 12) +
          381463680.0 * pow(x, 10) - 112477950.0 * pow(x, 8) +
          18378360 * pow(x, 6) - 1480700 * pow(x, 4) + 45084 * pow(x, 2) - 221);
}

double
AssLegFunction::P_20_4(const double x) const
{
  return (77224455.0 / 32768.0) * pow(pow(x, 2) - 1, 2) *
         (6486285 * pow(x, 16) - 19957800 * pow(x, 14) + 24542700 * pow(x, 12) -
          15426840 * pow(x, 10) + 5259150 * pow(x, 8) - 950040 * pow(x, 6) +
          81900 * pow(x, 4) - 2600 * pow(x, 2) + 13);
}

double
AssLegFunction::P_20_4_Deriv(const double x) const
{
  return (77224455.0 / 8192.0) * x *
         (32431425 * pow(x, 18) - 148186665.0 * pow(x, 16) +
          283778340.0 * pow(x, 14) - 295645140.0 * pow(x, 12) +
          181966590.0 * pow(x, 10) - 67237950 * pow(x, 8) +
          14482260 * pow(x, 6) - 1674660 * pow(x, 4) + 87113 * pow(x, 2) -
          1313);
}

double
AssLegFunction::P_20_5(const double x) const
{
  return (386122275.0 / 2048.0) * x * pow(1 - pow(x, 2), 5.0 / 2.0) *
         (-1297257 * pow(x, 14) + 3492615 * pow(x, 12) - 3681405 * pow(x, 10) +
          1928355 * pow(x, 8) - 525915 * pow(x, 6) + 71253 * pow(x, 4) -
          4095 * pow(x, 2) + 65);
}

double
AssLegFunction::P_20_5_Deriv(const double x) const
{
  return (1930611375.0 / 2048.0) * sqrt(1 - pow(x, 2)) *
         (-5189028 * pow(x, 18) + 21654213 * pow(x, 16) -
          37326480 * pow(x, 14) + 34359780 * pow(x, 12) -
          18231720 * pow(x, 10) + 5612022 * pow(x, 8) - 956592 * pow(x, 6) +
          80340 * pow(x, 4) - 2548 * pow(x, 2) + 13);
}

double
AssLegFunction::P_20_6(const double x) const
{
  return (25097947875.0 / 2048.0) * pow(pow(x, 2) - 1, 3) *
         (-299367 * pow(x, 14) + 698523 * pow(x, 12) - 623007 * pow(x, 10) +
          267003 * pow(x, 8) - 56637 * pow(x, 6) + 5481 * pow(x, 4) -
          189 * pow(x, 2) + 1);
}

double
AssLegFunction::P_20_6_Deriv(const double x) const
{
  return (75293843625.0 / 512.0) * x *
         (-498945 * pow(x, 18) + 2394936 * pow(x, 16) - 4822236 * pow(x, 14) +
          5286120 * pow(x, 12) - 3425190 * pow(x, 10) + 1332840 * pow(x, 8) -
          302364 * pow(x, 6) + 36824 * pow(x, 4) - 2017 * pow(x, 2) + 32);
}

double
AssLegFunction::P_20_7(const double x) const
{
  return (225881530875.0 / 1024.0) * x * pow(1 - pow(x, 2), 7.0 / 2.0) *
         (-232841 * pow(x, 12) + 465682 * pow(x, 10) - 346115 * pow(x, 8) +
          118668 * pow(x, 6) - 18879 * pow(x, 4) + 1218 * pow(x, 2) - 21);
}

double
AssLegFunction::P_20_7_Deriv(const double x) const
{
  return (1581170716125.0 / 1024.0) * sqrt(1 - pow(x, 2)) *
         (665260 * pow(x, 18) - 2960407 * pow(x, 16) + 5447940 * pow(x, 14) -
          5358040 * pow(x, 12) + 3038620 * pow(x, 10) - 999630 * pow(x, 8) +
          182028 * pow(x, 6) - 16320 * pow(x, 4) + 552 * pow(x, 2) - 3);
}

double
AssLegFunction::P_20_8(const double x) const
{
  return (1581170716125.0 / 1024.0) * pow(pow(x, 2) - 1, 4) *
         (432419 * pow(x, 12) - 731786 * pow(x, 10) + 445005 * pow(x, 8) -
          118668 * pow(x, 6) + 13485 * pow(x, 4) - 522 * pow(x, 2) + 3);
}

double
AssLegFunction::P_20_8_Deriv(const double x) const
{
  return (1581170716125.0 / 256.0) * x *
         (2162095 * pow(x, 18) - 11076579 * pow(x, 16) + 23866652 * pow(x, 14) -
          28066780 * pow(x, 12) + 19553250 * pow(x, 10) - 8195690 * pow(x, 8) +
          2005356 * pow(x, 6) - 263628 * pow(x, 4) + 15591 * pow(x, 2) - 267);
}

double
AssLegFunction::P_20_9(const double x) const
{
  return (45853950767625.0 / 256.0) * x * pow(1 - pow(x, 2), 9.0 / 2.0) *
         (-44733 * pow(x, 10) + 63085 * pow(x, 8) - 30690 * pow(x, 6) +
          6138 * pow(x, 4) - 465 * pow(x, 2) + 9);
}

double
AssLegFunction::P_20_9_Deriv(const double x) const
{
  return (137561852302875.0 / 256.0) * sqrt(1 - pow(x, 2)) *
         (-298220 * pow(x, 18) + 1437191 * pow(x, 16) - 2875188 * pow(x, 14) +
          3084872 * pow(x, 12) - 1914188 * pow(x, 10) + 690462 * pow(x, 8) -
          138012 * pow(x, 6) + 13584 * pow(x, 4) - 504 * pow(x, 2) + 3);
}

double
AssLegFunction::P_20_10(const double x) const
{
  return (137561852302875.0 / 256.0) * pow(pow(x, 2) - 1, 5) *
         (-164021 * pow(x, 10) + 189255 * pow(x, 8) - 71610 * pow(x, 6) +
          10230 * pow(x, 4) - 465 * pow(x, 2) + 3);
}

double
AssLegFunction::P_20_10_Deriv(const double x) const
{
  return (687809261514375.0 / 64.0) * x *
         (-164021 * pow(x, 18) + 908424 * pow(x, 16) - 2126476 * pow(x, 14) +
          2730728 * pow(x, 12) - 2088222 * pow(x, 10) + 965512 * pow(x, 8) -
          261708 * pow(x, 6) + 38232 * pow(x, 4) - 2517 * pow(x, 2) + 48);
}

double
AssLegFunction::P_20_11(const double x) const
{
  return (21322087106945625.0 / 128.0) * x * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (-5291 * pow(x, 8) + 4884 * pow(x, 6) - 1386 * pow(x, 4) +
          132 * pow(x, 2) - 3);
}

double
AssLegFunction::P_20_11_Deriv(const double x) const
{
  return (21322087106945625.0 / 128.0) * sqrt(1 - pow(x, 2)) *
         (105820 * pow(x, 18) - 558811 * pow(x, 16) + 1233408 * pow(x, 14) -
          1470700 * pow(x, 12) + 1021672 * pow(x, 10) - 415386 * pow(x, 8) +
          94080 * pow(x, 6) - 10524 * pow(x, 4) + 444 * pow(x, 2) - 3);
}

double
AssLegFunction::P_20_12(const double x) const
{
  return (63966261320836875.0 / 128.0) * pow(pow(x, 2) - 1, 6) *
         (15873 * pow(x, 8) - 11396 * pow(x, 6) + 2310 * pow(x, 4) -
          132 * pow(x, 2) + 1);
}

double
AssLegFunction::P_20_12_Deriv(const double x) const
{
  return (191898783962510625.0 / 32.0) * x *
         (26455 * pow(x, 18) - 159951 * pow(x, 16) + 411708 * pow(x, 14) -
          586124 * pow(x, 12) + 501458 * pow(x, 10) - 261970 * pow(x, 8) +
          81036 * pow(x, 6) - 13628 * pow(x, 4) + 1039 * pow(x, 2) - 23);
}

double
AssLegFunction::P_20_13(const double x) const
{
  return (2110886623587616875.0 / 16.0) * x * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (-481 * pow(x, 6) + 259 * pow(x, 4) - 35 * pow(x, 2) + 1);
}

double
AssLegFunction::P_20_13_Deriv(const double x) const
{
  return (2110886623587616875.0 / 16.0) * pow(1 - pow(x, 2), 11.0 / 2.0) *
         (13 * pow(x, 2) *
            (481 * pow(x, 6) - 259 * pow(x, 4) + 35 * pow(x, 2) - 1) +
          (1 - pow(x, 2)) *
            (-3367 * pow(x, 6) + 1295 * pow(x, 4) - 105 * pow(x, 2) + 1));
}

double
AssLegFunction::P_20_14(const double x) const
{
  return (2110886623587616875.0 / 16.0) * pow(pow(x, 2) - 1, 7) *
         (-3367 * pow(x, 6) + 1295 * pow(x, 4) - 105 * pow(x, 2) + 1);
}

double
AssLegFunction::P_20_14_Deriv(const double x) const
{
  return (14776206365113318125.0 / 8.0) * x * pow(pow(x, 2) - 1, 6) *
         (-3367 * pow(x, 6) + 1295 * pow(x, 4) - 105 * pow(x, 2) +
          (pow(x, 2) - 1) * (-1443 * pow(x, 4) + 370 * pow(x, 2) - 15) + 1);
}

double
AssLegFunction::P_20_15(const double x) const
{
  return (14776206365113318125.0 / 8.0) * x * pow(1 - pow(x, 2), 15.0 / 2.0) *
         (-1443 * pow(x, 4) + 370 * pow(x, 2) - 15);
}

double
AssLegFunction::P_20_15_Deriv(const double x) const
{
  return (221643095476699771875.0 / 8.0) * pow(1 - pow(x, 2), 13.0 / 2.0) *
         (pow(x, 2) * (1443 * pow(x, 4) - 370 * pow(x, 2) + 15) +
          (1 - pow(x, 2)) * (-481 * pow(x, 4) + 74 * pow(x, 2) - 1));
}

double
AssLegFunction::P_20_16(const double x) const
{
  return (221643095476699771875.0 / 8.0) * pow(pow(x, 2) - 1, 8) *
         (481 * pow(x, 4) - 74 * pow(x, 2) + 1);
}

double
AssLegFunction::P_20_16_Deriv(const double x) const
{
  return (221643095476699771875.0 / 2.0) * x * pow(pow(x, 2) - 1, 7) *
         (1924 * pow(x, 4) - 296 * pow(x, 2) +
          37 * (pow(x, 2) - 1) * (13 * pow(x, 2) - 1) + 4);
}

double
AssLegFunction::P_20_17(const double x) const
{
  return (8200794532637891559375.0 / 2.0) * x * (1 - 13 * pow(x, 2)) *
         pow(1 - pow(x, 2), 17.0 / 2.0);
}

double
AssLegFunction::P_20_17_Deriv(const double x) const
{
  return (8200794532637891559375.0 / 2.0) * pow(1 - pow(x, 2), 15.0 / 2.0) *
         (17 * pow(x, 2) * (13 * pow(x, 2) - 1) +
          (1 - 39 * pow(x, 2)) * (1 - pow(x, 2)));
}

double
AssLegFunction::P_20_18(const double x) const
{
  return (8200794532637891559375.0 / 2.0) * (1 - 39 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 9);
}

double
AssLegFunction::P_20_18_Deriv(const double x) const
{
  return 49204767195827349356250.0 * x * (8 - 65 * pow(x, 2)) *
         pow(pow(x, 2) - 1, 8);
}

double
AssLegFunction::P_20_19(const double x) const
{
  return -319830986772877770815625.0 * x * pow(1 - pow(x, 2), 19.0 / 2.0);
}

double
AssLegFunction::P_20_19_Deriv(const double x) const
{
  return pow(1 - pow(x, 2), 17.0 / 2.0) *
         (6396619735457555416312500.0 * pow(x, 2) - 319830986772877770815625.0);
}

double
AssLegFunction::P_20_20(const double x) const
{
  return 319830986772877770815625.0 * pow(pow(x, 2) - 1, 10);
}

double
AssLegFunction::P_20_20_Deriv(const double x) const
{
  return 6396619735457555416312500.0 * x * pow(pow(x, 2) - 1, 9);
}
