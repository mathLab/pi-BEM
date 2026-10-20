#include "../include/ass_leg_function.h"

#include <math.h>

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
AssLegFunction::P_0_0(const double /*x*/) const
{
  return (1);
}

double
AssLegFunction::P_0_0_Deriv(const double /*x*/) const
{
  return (0);
}

double
AssLegFunction::P_1_0(const double x) const
{
  return (x);
}

double
AssLegFunction::P_1_0_Deriv(const double /*x*/) const
{
  return (1);
}

double
AssLegFunction::P_1_1(const double x) const
{
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.2e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (-t1 * t4 / 0.2e1);
}



double
AssLegFunction::P_1_1_Deriv(const double x) const
{
  double t2;
  double t4;
  double t1;
  t1 = sqrt(0.2e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (t1 / t4 * x / 0.2e1);
}

double
AssLegFunction::P_2_0(const double x) const
{
  double t1;
  t1 = x * x;
  return (0.3e1 / 0.2e1 * t1 - 0.1e1 / 0.2e1);
}

double
AssLegFunction::P_2_0_Deriv(const double x) const
{
  return (0.3e1 * x);
}



double
AssLegFunction::P_2_1(const double x) const
{
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.6e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (-t1 * t4 * x / 0.2e1);
}



double
AssLegFunction::P_2_1_Deriv(const double x) const
{
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.6e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (t1 / t4 * t2 / 0.2e1 - t1 * t4 / 0.2e1);
}



double
AssLegFunction::P_2_2(const double x) const
{
  double t2;
  double t1;
  t1 = sqrt(0.6e1);
  t2 = x * x;
  return (t1 * (0.1e1 - t2) / 0.4e1);
}



double
AssLegFunction::P_2_2_Deriv(const double x) const
{
  double t1;
  t1 = sqrt(0.6e1);
  return (-t1 * x / 0.2e1);
}

double
AssLegFunction::P_3_0(const double x) const
{
  double t1;
  t1 = x * x;
  return (t1 * x + 0.3e1 / 0.2e1 * (t1 - 0.1e1) * x);
}

double
AssLegFunction::P_3_0_Deriv(const double x) const
{
  double t1;
  t1 = x * x;
  return (0.15e2 / 0.2e1 * t1 - 0.3e1 / 0.2e1);
}



double
AssLegFunction::P_3_1(const double x) const
{
  double t4;
  double t2;
  double t1;
  t1 = sqrt(0.3e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (-t1 * t4 * (0.360e3 * t2 - 0.72e2) / 0.288e3);
}



double
AssLegFunction::P_3_1_Deriv(const double x) const
{
  double t4;
  double t1;
  double t2;
  t1 = sqrt(0.3e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (t1 / t4 * (0.360e3 * t2 - 0.72e2) * x / 0.288e3 -
          0.5e1 / 0.2e1 * t1 * t4 * x);
}



double
AssLegFunction::P_3_2(const double x) const
{
  double t1;
  double t2;
  t1 = sqrt(0.30e2);
  t2 = x * x;
  return (t1 * (0.1e1 - t2) * x / 0.4e1);
}



double
AssLegFunction::P_3_2_Deriv(const double x) const
{
  double t1;
  double t2;
  t1 = sqrt(0.30e2);
  t2 = x * x;
  return (-t1 * t2 / 0.2e1 + t1 * (0.1e1 - t2) / 0.4e1);
}



double
AssLegFunction::P_3_3(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t1;
  t1 = sqrt(0.5e1);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (-t1 * t4 * t3 / 0.4e1);
}



double
AssLegFunction::P_3_3_Deriv(const double x) const
{
  double t2;
  double t4;
  double t1;
  t1 = sqrt(0.5e1);
  t2 = x * x;
  t4 = sqrt(0.1e1 - t2);
  return (0.3e1 / 0.4e1 * t1 * t4 * x);
}

double
AssLegFunction::P_4_0(const double x) const
{
  double t6;
  double t2;
  double t3;
  double t1;
  t1 = x * x;
  t2 = t1 * t1;
  t3 = t1 - 0.1e1;
  t6 = t3 * t3;
  return (t2 + 0.3e1 * t3 * t1 + 0.3e1 / 0.8e1 * t6);
}

double
AssLegFunction::P_4_0_Deriv(const double x) const
{
  double t1;
  t1 = x * x;
  return (0.10e2 * t1 * x + 0.15e2 / 0.2e1 * (t1 - 0.1e1) * x);
}



double
AssLegFunction::P_4_1(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.5e1);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (-t1 * t4 * (0.3840e4 * t2 * x - 0.2880e4 * t3 * x) / 0.3840e4);
}



double
AssLegFunction::P_4_1_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.5e1);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (t1 / t4 * (0.3840e4 * t2 * x - 0.2880e4 * t3 * x) * x / 0.3840e4 -
          t1 * t4 * (0.20160e5 * t2 - 0.2880e4) / 0.3840e4);
}



double
AssLegFunction::P_4_2(const double x) const
{
  double t2;
  double t1;
  t1 = sqrt(0.10e2);
  t2 = x * x;
  return (t1 * (0.1e1 - t2) * (0.20160e5 * t2 - 0.2880e4) / 0.23040e5);
}



double
AssLegFunction::P_4_2_Deriv(const double x) const
{
  double t3;
  double t1;
  t1 = sqrt(0.10e2);
  t3 = x * x;
  return (-t1 * x * (0.20160e5 * t3 - 0.2880e4) / 0.11520e5 +
          0.7e1 / 0.4e1 * t1 * (0.1e1 - t3) * x);
}



double
AssLegFunction::P_4_3(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.35e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (-t1 * t4 * t3 * x / 0.4e1);
}



double
AssLegFunction::P_4_3_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.35e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (0.3e1 / 0.4e1 * t1 * t4 * t2 - t1 * t4 * t3 / 0.4e1);
}



double
AssLegFunction::P_4_4(const double x) const
{
  double t2;
  double t1;
  double t4;
  t1 = sqrt(0.70e2);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  return (t1 * t4 / 0.16e2);
}



double
AssLegFunction::P_4_4_Deriv(const double x) const
{
  double t2;
  double t1;
  t1 = sqrt(0.70e2);
  t2 = x * x;
  return (-t1 * (0.1e1 - t2) * x / 0.4e1);
}

double
AssLegFunction::P_5_0(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t8;
  t1 = x * x;
  t2 = t1 * t1;
  t4 = t1 - 0.1e1;
  t8 = t4 * t4;
  return (t2 * x + 0.5e1 * t4 * t1 * x + 0.15e2 / 0.8e1 * t8 * x);
}

double
AssLegFunction::P_5_0_Deriv(const double x) const
{
  double t2;
  double t1;
  double t4;
  double t7;
  t1 = x * x;
  t2 = t1 * t1;
  t4 = t1 - 0.1e1;
  t7 = t4 * t4;
  return (0.15e2 * t2 + 0.45e2 / 0.2e1 * t1 * t4 + 0.15e2 / 0.8e1 * t7);
}



double
AssLegFunction::P_5_1(const double x) const
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  t1  = sqrt(0.30e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  return (-t1 * t4 * (0.57600e5 * t6 + 0.86400e5 * t8 * t2 + 0.7200e4 * t11) /
          0.115200e6);
}



double
AssLegFunction::P_5_1_Deriv(const double x) const
{
  double t9;
  double t7;
  double t1;
  double t4;
  double t12;
  double t3;
  double t2;
  t1  = sqrt(0.30e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (t1 / t4 * (0.57600e5 * t7 + 0.86400e5 * t9 * t2 + 0.7200e4 * t12) *
            x / 0.115200e6 -
          t1 * t4 * (0.403200e6 * t2 * x + 0.201600e6 * t9 * x) / 0.115200e6);
}



double
AssLegFunction::P_5_2(const double x) const
{
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.210e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  return (t1 * t3 * (0.403200e6 * t2 * x - 0.201600e6 * t3 * x) / 0.1612800e7);
}



double
AssLegFunction::P_5_2_Deriv(const double x) const
{
  double t6;
  double t1;
  double t3;
  t1 = sqrt(0.210e3);
  t3 = x * x;
  t6 = t3 - 0.1e1;
  return (-t1 * x * (0.403200e6 * t3 * x + 0.201600e6 * t6 * x) / 0.806400e6 -
          t1 * t6 * (0.1814400e7 * t3 - 0.201600e6) / 0.1612800e7);
}



double
AssLegFunction::P_5_3(const double x) const
{
  double t4;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.35e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (-t1 * t4 * t3 * (0.1814400e7 * t2 - 0.201600e6) / 0.3225600e7);
}



double
AssLegFunction::P_5_3_Deriv(const double x) const
{
  double t3;
  double t4;
  double t1;
  double t2;
  t1 = sqrt(0.35e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (t1 * t4 * (0.1814400e7 * t2 - 0.201600e6) * x / 0.1075200e7 -
          0.9e1 / 0.8e1 * t1 * t4 * t3 * x);
}



double
AssLegFunction::P_5_4(const double x) const
{
  double t2;
  double t4;
  double t1;
  t1 = sqrt(0.70e2);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  return (0.3e1 / 0.16e2 * t1 * t4 * x);
}



double
AssLegFunction::P_5_4_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t7;
  t1 = sqrt(0.70e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t7 = t3 * t3;
  return (-0.3e1 / 0.4e1 * t1 * t3 * t2 + 0.3e1 / 0.16e2 * t1 * t7);
}



double
AssLegFunction::P_5_5(const double x) const
{
  double t5;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.7e1);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (-0.3e1 / 0.16e2 * t1 * t5 * t4);
}



double
AssLegFunction::P_5_5_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.7e1);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (0.15e2 / 0.16e2 * t1 * t4 * t3 * x);
}

double
AssLegFunction::P_6_0(const double x) const
{
  double t7;
  double t4;
  double t1;
  double t2;
  t1 = x * x;
  t2 = t1 * t1;
  t4 = t1 - 0.1e1;
  t7 = t4 * t4;
  return (t2 * t1 + 0.15e2 / 0.2e1 * t4 * t2 + 0.45e2 / 0.8e1 * t7 * t1 +
          0.5e1 / 0.16e2 * t7 * t4);
}

double
AssLegFunction::P_6_0_Deriv(const double x) const
{
  double t1;
  double t2;
  double t5;
  double t9;
  t1 = x * x;
  t2 = t1 * t1;
  t5 = t1 - 0.1e1;
  t9 = t5 * t5;
  return (0.21e2 * t2 * x + 0.105e3 / 0.2e1 * t5 * t1 * x +
          0.105e3 / 0.8e1 * t9 * x);
}



double
AssLegFunction::P_6_1(const double x) const
{
  double t13;
  double t6;
  double t2;
  double t3;
  double t4;
  double t9;
  double t1;
  t1  = sqrt(0.42e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  return (
    -t1 * t4 *
    (0.967680e6 * t6 * x + 0.2419200e7 * t9 * t2 * x + 0.604800e6 * t13 * x) /
    0.1935360e7);
}



double
AssLegFunction::P_6_1_Deriv(const double x) const
{
  double t2;
  double t3;
  double t14;
  double t4;
  double t7;
  double t1;
  double t10;
  t1  = sqrt(0.42e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (t1 / t4 *
            (0.967680e6 * t7 * x + 0.2419200e7 * t10 * t2 * x +
             0.604800e6 * t14 * x) *
            x / 0.1935360e7 -
          t1 * t4 *
            (0.9676800e7 * t7 + 0.9676800e7 * t10 * t2 + 0.604800e6 * t14) /
            0.1935360e7);
}



double
AssLegFunction::P_6_2(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t5;
  double t7;
  double t10;
  t1  = sqrt(0.105e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t7  = -t3;
  t10 = t7 * t7;
  return (t1 * t3 *
          (0.9676800e7 * t5 + 0.9676800e7 * t7 * t2 + 0.604800e6 * t10) /
          0.19353600e8);
}



double
AssLegFunction::P_6_2_Deriv(const double x) const
{
  double t1;
  double t6;
  double t9;
  double t3;
  double t4;
  t1 = sqrt(0.105e3);
  t3 = x * x;
  t4 = t3 * t3;
  t6 = t3 - 0.1e1;
  t9 = t6 * t6;
  return (
    -t1 * x * (0.9676800e7 * t4 + 0.9676800e7 * t6 * t3 + 0.604800e6 * t9) /
      0.9676800e7 -
    t1 * t6 * (0.58060800e8 * t3 * x + 0.21772800e8 * t6 * x) / 0.19353600e8);
}



double
AssLegFunction::P_6_3(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.105e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (-t1 * t4 * t3 * (0.58060800e8 * t2 * x - 0.21772800e8 * t3 * x) /
          0.116121600e9);
}



double
AssLegFunction::P_6_3_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.105e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  return (t1 * t4 * (0.58060800e8 * t2 * x - 0.21772800e8 * t3 * x) * x /
            0.38707200e8 -
          t1 * t4 * t3 * (0.239500800e9 * t2 - 0.21772800e8) / 0.116121600e9);
}



double
AssLegFunction::P_6_4(const double x) const
{
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.14e2);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  return (t1 * t4 * (0.239500800e9 * t2 - 0.21772800e8) / 0.232243200e9);
}



double
AssLegFunction::P_6_4_Deriv(const double x) const
{
  double t10;
  double t2;
  double t3;
  double t1;
  t1  = sqrt(0.14e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t10 = t3 * t3;
  return (-t1 * t3 * (0.239500800e9 * t2 - 0.21772800e8) * x / 0.58060800e8 +
          0.33e2 / 0.16e2 * t1 * t10 * x);
}



double
AssLegFunction::P_6_5(const double x) const
{
  double t3;
  double t4;
  double t5;
  double t2;
  double t1;
  t1 = sqrt(0.77e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (-0.3e1 / 0.16e2 * t1 * t5 * t4 * x);
}



double
AssLegFunction::P_6_5_Deriv(const double x) const
{
  double t9;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.77e2);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = sqrt(t3);
  t9 = t3 * t3;
  return (0.15e2 / 0.16e2 * t1 * t4 * t3 * t2 - 0.3e1 / 0.16e2 * t1 * t4 * t9);
}



double
AssLegFunction::P_6_6(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t1;
  t1 = sqrt(0.231e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (t1 * t4 * t3 / 0.32e2);
}



double
AssLegFunction::P_6_6_Deriv(const double x) const
{
  double t2;
  double t4;
  double t1;
  t1 = sqrt(0.231e3);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  return (-0.3e1 / 0.16e2 * t1 * t4 * x);
}

double
AssLegFunction::P_7_0(const double x) const
{
  double t5;
  double t9;
  double t1;
  double t2;
  double t3;
  t1 = x * x;
  t2 = t1 * x;
  t3 = t1 * t1;
  t5 = t1 - 0.1e1;
  t9 = t5 * t5;
  return (t3 * t2 + 0.21e2 / 0.2e1 * t5 * t3 * x + 0.105e3 / 0.8e1 * t9 * t2 +
          0.35e2 / 0.16e2 * t9 * t5 * x);
}

double
AssLegFunction::P_7_0_Deriv(const double x) const
{
  double t1;
  double t2;
  double t5;
  double t8;
  t1 = x * x;
  t2 = t1 * t1;
  t5 = t1 - 0.1e1;
  t8 = t5 * t5;
  return (0.28e2 * t2 * t1 + 0.105e3 * t2 * t5 + 0.105e3 / 0.2e1 * t8 * t1 +
          0.35e2 / 0.16e2 * t8 * t5);
}



double
AssLegFunction::P_7_1(const double x) const
{
  double t9;
  double t2;
  double t3;
  double t4;
  double t6;
  double t12;
  double t1;
  t1  = sqrt(0.14e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (-t1 * t4 *
          (0.18063360e8 * t6 * t2 + 0.67737600e8 * t9 * t6 +
           0.33868800e8 * t12 * t2 + 0.1411200e7 * t12 * t9) /
          0.18063360e8);
}



double
AssLegFunction::P_7_1_Deriv(const double x) const
{
  double t10;
  double t7;
  double t13;
  double t1;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.14e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 / t4 *
            (0.18063360e8 * t2 * t7 + 0.67737600e8 * t10 * t7 +
             0.33868800e8 * t13 * t2 + 0.1411200e7 * t13 * t10) *
            x / 0.18063360e8 -
          t1 * t4 *
            (0.243855360e9 * t7 * x + 0.406425600e9 * t10 * t2 * x +
             0.76204800e8 * t13 * x) /
            0.18063360e8);
}



double
AssLegFunction::P_7_2(const double x) const
{
  double t5;
  double t8;
  double t2;
  double t3;
  double t1;
  double t12;
  t1  = sqrt(0.21e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t8  = -t3;
  t12 = t8 * t8;
  return (t1 * t3 *
          (0.243855360e9 * t5 * x + 0.406425600e9 * t8 * t2 * x +
           0.76204800e8 * t12 * x) /
          0.162570240e9);
}



double
AssLegFunction::P_7_2_Deriv(const double x) const
{
  double t1;
  double t3;
  double t4;
  double t7;
  double t11;
  t1  = sqrt(0.21e2);
  t3  = x * x;
  t4  = t3 * t3;
  t7  = t3 - 0.1e1;
  t11 = t7 * t7;
  return (
    -t1 * x *
      (0.243855360e9 * t4 * x + 0.406425600e9 * t7 * t3 * x +
       0.76204800e8 * t11 * x) /
      0.81285120e8 -
    t1 * t7 *
      (0.2032128000e10 * t4 + 0.1524096000e10 * t7 * t3 + 0.76204800e8 * t11) /
      0.162570240e9);
}



double
AssLegFunction::P_7_3(const double x) const
{
  double t1;
  double t12;
  double t7;
  double t9;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.42e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (
    -t1 * t4 * t3 *
    (0.2032128000e10 * t7 + 0.1524096000e10 * t9 * t2 + 0.76204800e8 * t12) /
    0.1625702400e10);
}



double
AssLegFunction::P_7_3_Deriv(const double x) const
{
  double t6;
  double t1;
  double t8;
  double t4;
  double t2;
  double t3;
  double t11;
  t1  = sqrt(0.42e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  return (
    t1 * t4 *
      (0.2032128000e10 * t6 + 0.1524096000e10 * t8 * t2 + 0.76204800e8 * t11) *
      x / 0.541900800e9 -
    t1 * t4 * t3 * (0.11176704000e11 * t2 * x + 0.3353011200e10 * t8 * x) /
      0.1625702400e10);
}



double
AssLegFunction::P_7_4(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.462e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (t1 * t4 * (0.11176704000e11 * t2 * x - 0.3353011200e10 * t3 * x) /
          0.35765452800e11);
}



double
AssLegFunction::P_7_4_Deriv(const double x) const
{
  double t14;
  double t2;
  double t3;
  double t1;
  t1  = sqrt(0.462e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t14 = t3 * t3;
  return (-t1 * t3 * (0.11176704000e11 * t2 * x - 0.3353011200e10 * t3 * x) *
            x / 0.8941363200e10 +
          t1 * t14 * (0.43589145600e11 * t2 - 0.3353011200e10) /
            0.35765452800e11);
}



double
AssLegFunction::P_7_5(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.462e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (-t1 * t5 * t4 * (0.43589145600e11 * t2 - 0.3353011200e10) /
          0.214592716800e12);
}



double
AssLegFunction::P_7_5_Deriv(const double x) const
{
  double t2;
  double t3;
  double t12;
  double t4;
  double t1;
  t1  = sqrt(0.462e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t12 = t3 * t3;
  return (t1 * t4 * t3 * (0.43589145600e11 * t2 - 0.3353011200e10) * x /
            0.42918543360e11 -
          0.13e2 / 0.32e2 * t1 * t4 * t12 * x);
}



double
AssLegFunction::P_7_6(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  t1 = sqrt(0.3003e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (t1 * t4 * t3 * x / 0.32e2);
}



double
AssLegFunction::P_7_6_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.3003e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (-0.3e1 / 0.16e2 * t1 * t4 * t2 + t1 * t4 * t3 / 0.32e2);
}



double
AssLegFunction::P_7_7(const double x) const
{
  double t6;
  double t1;
  double t3;
  double t2;
  double t4;
  t1 = sqrt(0.858e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = sqrt(t3);
  return (-t1 * t6 * t4 * t3 / 0.64e2);
}



double
AssLegFunction::P_7_7_Deriv(const double x) const
{
  double t2;
  double t5;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.858e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (0.7e1 / 0.64e2 * t1 * t5 * t4 * x);
}

double
AssLegFunction::P_8_0(const double x) const
{
  double t4;
  double t14;
  double t8;
  double t1;
  double t3;
  double t2;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t1 - 0.1e1;
  t8  = t4 * t4;
  t14 = t8 * t8;
  return (t3 + 0.14e2 * t4 * t2 * t1 + 0.105e3 / 0.4e1 * t8 * t2 +
          0.35e2 / 0.4e1 * t8 * t4 * t1 + 0.35e2 / 0.128e3 * t14);
}

double
AssLegFunction::P_8_0_Deriv(const double x) const
{
  double t3;
  double t2;
  double t10;
  double t6;
  double t1;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t6  = t1 - 0.1e1;
  t10 = t6 * t6;
  return (0.36e2 * t3 * t2 + 0.189e3 * t6 * t3 * x +
          0.315e3 / 0.2e1 * t2 * t10 + 0.315e3 / 0.16e2 * t10 * t6 * x);
}



double
AssLegFunction::P_8_1(const double x) const
{
  double t7;
  double t10;
  double t14;
  double t1;
  double t2;
  double t4;
  double t3;
  double t6;
  t1  = sqrt(0.2e1);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (-t1 * t4 *
          (0.371589120e9 * t7 * t6 + 0.1950842880e10 * t10 * t7 * x +
           0.1625702400e10 * t14 * t6 + 0.203212800e9 * t14 * t10 * x) /
          0.123863040e9);
}



double
AssLegFunction::P_8_1_Deriv(const double x) const
{
  double t18;
  double t7;
  double t8;
  double t11;
  double t1;
  double t2;
  double t4;
  double t3;
  double t15;
  t1  = sqrt(0.2e1);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  return (t1 / t4 *
            (0.371589120e9 * t7 * t8 + 0.1950842880e10 * t11 * t8 * x +
             0.1625702400e10 * t15 * t7 + 0.203212800e9 * t18 * x) *
            x / 0.123863040e9 -
          t1 * t4 *
            (0.6502809600e10 * t8 * t2 + 0.16257024000e11 * t11 * t8 +
             0.6096384000e10 * t15 * t2 + 0.203212800e9 * t18) /
            0.123863040e9);
}



double
AssLegFunction::P_8_2(const double x) const
{
  double t1;
  double t5;
  double t8;
  double t2;
  double t3;
  double t11;
  t1  = sqrt(0.35e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  return (t1 * t3 *
          (0.6502809600e10 * t5 * t2 + 0.16257024000e11 * t5 * t8 +
           0.6096384000e10 * t11 * t2 + 0.203212800e9 * t11 * t8) /
          0.4335206400e10);
}



double
AssLegFunction::P_8_2_Deriv(const double x) const
{
  double t10;
  double t1;
  double t3;
  double t4;
  double t7;
  t1  = sqrt(0.35e2);
  t3  = x * x;
  t4  = t3 * t3;
  t7  = t3 - 0.1e1;
  t10 = t7 * t7;
  return (-t1 * x *
            (0.6502809600e10 * t4 * t3 + 0.16257024000e11 * t7 * t4 +
             0.6096384000e10 * t10 * t3 + 0.203212800e9 * t10 * t7) /
            0.2167603200e10 -
          t1 * t7 *
            (0.71530905600e11 * t4 * x + 0.89413632000e11 * t7 * t3 * x +
             0.13412044800e11 * t10 * x) /
            0.4335206400e10);
}



double
AssLegFunction::P_8_3(const double x) const
{
  double t1;
  double t7;
  double t10;
  double t14;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.2310e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (-t1 * t4 * t3 *
          (0.71530905600e11 * t7 * x + 0.89413632000e11 * t10 * t2 * x +
           0.13412044800e11 * t14 * x) /
          0.286123622400e12);
}



double
AssLegFunction::P_8_3_Deriv(const double x) const
{
  double t6;
  double t13;
  double t2;
  double t3;
  double t9;
  double t4;
  double t1;
  t1  = sqrt(0.2310e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  return (t1 * t4 *
            (0.71530905600e11 * t6 * x + 0.89413632000e11 * t2 * t9 * x +
             0.13412044800e11 * t13 * x) *
            x / 0.95374540800e11 -
          t1 * t4 * t3 *
            (0.536481792000e12 * t6 + 0.321889075200e12 * t9 * t2 +
             0.13412044800e11 * t13) /
            0.286123622400e12);
}



double
AssLegFunction::P_8_4(const double x) const
{
  double t11;
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  t1  = sqrt(0.154e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  return (t1 * t4 *
          (0.536481792000e12 * t6 + 0.321889075200e12 * t8 * t2 +
           0.13412044800e11 * t11) /
          0.572247244800e12);
}



double
AssLegFunction::P_8_4_Deriv(const double x) const
{
  double t5;
  double t16;
  double t7;
  double t1;
  double t2;
  double t3;
  double t10;
  t1  = sqrt(0.154e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t7  = -t3;
  t10 = t7 * t7;
  t16 = t3 * t3;
  return (-t1 * t3 *
            (0.536481792000e12 * t5 + 0.321889075200e12 * t7 * t2 +
             0.13412044800e11 * t10) *
            x / 0.143061811200e12 +
          t1 * t16 *
            (0.2789705318400e13 * t2 * x + 0.697426329600e12 * t7 * x) /
            0.572247244800e12);
}



double
AssLegFunction::P_8_5(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t5;
  double t4;
  t1 = sqrt(0.2002e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (-t1 * t5 * t4 *
          (0.2789705318400e13 * t2 * x - 0.697426329600e12 * t3 * x) /
          0.14878428364800e14);
}



double
AssLegFunction::P_8_5_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t16;
  double t4;
  t1  = sqrt(0.2002e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t16 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.2789705318400e13 * t2 * x - 0.697426329600e12 * t3 * x) * x /
            0.2975685672960e13 -
          t1 * t4 * t16 * (0.10461394944000e14 * t2 - 0.697426329600e12) /
            0.14878428364800e14);
}



double
AssLegFunction::P_8_6(const double x) const
{
  double t4;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.429e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (t1 * t4 * t3 * (0.10461394944000e14 * t2 - 0.697426329600e12) /
          0.44635285094400e14);
}



double
AssLegFunction::P_8_6_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.429e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (-t1 * t4 * (0.10461394944000e14 * t2 - 0.697426329600e12) * x /
            0.7439214182400e13 +
          0.15e2 / 0.32e2 * t1 * t4 * t3 * x);
}



double
AssLegFunction::P_8_7(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  t1 = sqrt(0.1430e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = sqrt(t3);
  return (-0.3e1 / 0.64e2 * t1 * t6 * t4 * t3 * x);
}



double
AssLegFunction::P_8_7_Deriv(const double x) const
{
  double t3;
  double t4;
  double t5;
  double t1;
  double t2;
  t1 = sqrt(0.1430e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (0.21e2 / 0.64e2 * t1 * t5 * t4 * t2 -
          0.3e1 / 0.64e2 * t1 * t5 * t4 * t3);
}



double
AssLegFunction::P_8_8(const double x) const
{
  double t5;
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.1430e4);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (0.3e1 / 0.256e3 * t1 * t5);
}



double
AssLegFunction::P_8_8_Deriv(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  t1 = sqrt(0.1430e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (-0.3e1 / 0.32e2 * t1 * t4 * t3 * x);
}

double
AssLegFunction::P_9_0(const double x) const
{
  double t10;
  double t1;
  double t5;
  double t17;
  double t6;
  double t3;
  double t2;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t5  = t1 - 0.1e1;
  t6  = t1 * x;
  t10 = t5 * t5;
  t17 = t10 * t10;
  return (t3 * x + 0.18e2 * t5 * t2 * t6 + 0.189e3 / 0.4e1 * t10 * t2 * x +
          0.105e3 / 0.4e1 * t10 * t5 * t6 + 0.315e3 / 0.128e3 * t17 * x);
}

double
AssLegFunction::P_9_0_Deriv(const double x) const
{
  double t1;
  double t5;
  double t2;
  double t9;
  double t15;
  double t3;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t5  = t1 - 0.1e1;
  t9  = t5 * t5;
  t15 = t9 * t9;
  return (0.45e2 * t3 + 0.315e3 * t5 * t2 * t1 + 0.1575e4 / 0.4e1 * t9 * t2 +
          0.1575e4 / 0.16e2 * t9 * t5 * t1 + 0.315e3 / 0.128e3 * t15);
}



double
AssLegFunction::P_9_1(const double x) const
{
  double t19;
  double t1;
  double t3;
  double t2;
  double t13;
  double t4;
  double t6;
  double t7;
  double t9;
  t1  = sqrt(0.10e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t9  = -t3;
  t13 = t9 * t9;
  t19 = t13 * t13;
  return (-t1 * t4 *
          (0.8360755200e10 * t7 + 0.58525286400e11 * t9 * t6 * t2 +
           0.73156608000e11 * t13 * t6 + 0.18289152000e11 * t13 * t9 * t2 +
           0.457228800e9 * t19) /
          0.5573836800e10);
}



double
AssLegFunction::P_9_1_Deriv(const double x) const
{
  double t14;
  double t17;
  double t1;
  double t26;
  double t8;
  double t7;
  double t20;
  double t10;
  double t4;
  double t2;
  double t3;
  t1  = sqrt(0.10e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t20 = t14 * t14;
  t26 = t2 * x;
  return (t1 / t4 *
            (0.8360755200e10 * t8 + 0.58525286400e11 * t10 * t7 * t2 +
             0.73156608000e11 * t14 * t7 + 0.18289152000e11 * t17 * t2 +
             0.457228800e9 * t20) *
            x / 0.5573836800e10 -
          t1 * t4 *
            (0.183936614400e12 * t7 * t26 + 0.643778150400e12 * t10 * t7 * x +
             0.402361344000e12 * t14 * t26 + 0.40236134400e11 * t17 * x) /
            0.5573836800e10);
}



double
AssLegFunction::P_9_2(const double x) const
{
  double t2;
  double t9;
  double t3;
  double t13;
  double t5;
  double t1;
  double t6;
  t1  = sqrt(0.55e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  return (t1 * t3 *
          (0.183936614400e12 * t6 * t5 + 0.643778150400e12 * t9 * t6 * x +
           0.402361344000e12 * t13 * t5 + 0.40236134400e11 * t13 * t9 * x) /
          0.122624409600e12);
}



double
AssLegFunction::P_9_2_Deriv(const double x) const
{
  double t1;
  double t12;
  double t15;
  double t4;
  double t3;
  double t5;
  double t8;
  t1  = sqrt(0.55e2);
  t3  = x * x;
  t4  = t3 * x;
  t5  = t3 * t3;
  t8  = t3 - 0.1e1;
  t12 = t8 * t8;
  t15 = t12 * t8;
  return (-t1 * x *
            (0.183936614400e12 * t5 * t4 + 0.643778150400e12 * t8 * t5 * x +
             0.402361344000e12 * t12 * t4 + 0.40236134400e11 * t15 * x) /
            0.61312204800e11 -
          t1 * t8 *
            (0.2575112601600e13 * t5 * t3 + 0.4828336128000e13 * t8 * t5 +
             0.1448500838400e13 * t12 * t3 + 0.40236134400e11 * t15) /
            0.122624409600e12);
}



double
AssLegFunction::P_9_3(const double x) const
{
  double t10;
  double t2;
  double t3;
  double t13;
  double t4;
  double t7;
  double t1;
  t1  = sqrt(0.1155e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (-t1 * t4 * t3 *
          (0.2575112601600e13 * t7 * t2 + 0.4828336128000e13 * t10 * t7 +
           0.1448500838400e13 * t13 * t2 + 0.40236134400e11 * t13 * t10) /
          0.5150225203200e13);
}



double
AssLegFunction::P_9_3_Deriv(const double x) const
{
  double t9;
  double t6;
  double t2;
  double t3;
  double t4;
  double t1;
  double t12;
  t1  = sqrt(0.1155e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (t1 * t4 *
            (0.2575112601600e13 * t6 * t2 + 0.4828336128000e13 * t9 * t6 +
             0.1448500838400e13 * t12 * t2 + 0.40236134400e11 * t12 * t9) *
            x / 0.1716741734400e13 -
          t1 * t4 * t3 *
            (0.25107347865600e14 * t6 * x + 0.25107347865600e14 * t9 * t2 * x +
             0.3138418483200e13 * t12 * x) /
            0.5150225203200e13);
}



double
AssLegFunction::P_9_4(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t13;
  double t6;
  double t9;
  t1  = sqrt(0.10010e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  return (t1 * t4 *
          (0.25107347865600e14 * t6 * x + 0.25107347865600e14 * t9 * t2 * x +
           0.3138418483200e13 * t13 * x) /
          0.133905855283200e15);
}



double
AssLegFunction::P_9_4_Deriv(const double x) const
{
  double t3;
  double t12;
  double t5;
  double t19;
  double t1;
  double t8;
  double t2;
  t1  = sqrt(0.10010e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t8  = -t3;
  t12 = t8 * t8;
  t19 = t3 * t3;
  return (-t1 * t3 *
            (0.25107347865600e14 * t5 * x + 0.25107347865600e14 * t8 * t2 * x +
             0.3138418483200e13 * t12 * x) *
            x / 0.33476463820800e14 +
          t1 * t19 *
            (0.175751435059200e15 * t5 + 0.87875717529600e14 * t2 * t8 +
             0.3138418483200e13 * t12) /
            0.133905855283200e15);
}



double
AssLegFunction::P_9_5(const double x) const
{
  double t13;
  double t5;
  double t1;
  double t8;
  double t2;
  double t10;
  double t3;
  double t4;
  t1  = sqrt(0.143e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (-t1 * t5 * t4 *
          (0.175751435059200e15 * t8 + 0.87875717529600e14 * t10 * t2 +
           0.3138418483200e13 * t13) /
          0.133905855283200e15);
}



double
AssLegFunction::P_9_5_Deriv(const double x) const
{
  double t2;
  double t3;
  double t9;
  double t12;
  double t4;
  double t1;
  double t18;
  double t7;
  t1  = sqrt(0.143e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  t18 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.175751435059200e15 * t7 + 0.87875717529600e14 * t9 * t2 +
             0.3138418483200e13 * t12) *
            x / 0.26781171056640e14 -
          t1 * t4 * t18 *
            (0.878757175296000e15 * t2 * x + 0.188305108992000e15 * t9 * x) /
            0.133905855283200e15);
}



double
AssLegFunction::P_9_6(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.2145e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (t1 * t4 * t3 *
          (0.878757175296000e15 * t2 * x - 0.188305108992000e15 * t3 * x) /
          0.4017175658496000e16);
}



double
AssLegFunction::P_9_6_Deriv(const double x) const
{
  double t2;
  double t4;
  double t3;
  double t1;
  t1 = sqrt(0.2145e4);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  return (-t1 * t4 *
            (0.878757175296000e15 * t2 * x - 0.188305108992000e15 * t3 * x) *
            x / 0.669529276416000e15 +
          t1 * t4 * t3 * (0.3201186852864000e16 * t2 - 0.188305108992000e15) /
            0.4017175658496000e16);
}



double
AssLegFunction::P_9_7(const double x) const
{
  double t4;
  double t6;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.715e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = sqrt(t3);
  return (-t1 * t6 * t4 * t3 *
          (0.3201186852864000e16 * t2 - 0.188305108992000e15) /
          0.16068702633984000e17);
}



double
AssLegFunction::P_9_7_Deriv(const double x) const
{
  double t5;
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.715e3);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (t1 * t5 * t4 * (0.3201186852864000e16 * t2 - 0.188305108992000e15) *
            x / 0.2295528947712000e16 -
          0.51e2 / 0.128e3 * t1 * t5 * t4 * t3 * x);
}



double
AssLegFunction::P_9_8(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  t1 = sqrt(0.24310e5);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (0.3e1 / 0.256e3 * t1 * t5 * x);
}



double
AssLegFunction::P_9_8_Deriv(const double x) const
{
  double t2;
  double t3;
  double t9;
  double t4;
  double t1;
  t1 = sqrt(0.24310e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t9 = t4 * t4;
  return (-0.3e1 / 0.32e2 * t1 * t4 * t3 * t2 + 0.3e1 / 0.256e3 * t1 * t9);
}



double
AssLegFunction::P_9_9(const double x) const
{
  double t4;
  double t5;
  double t6;
  double t1;
  double t3;
  double t2;
  t1 = sqrt(0.12155e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (-t1 * t6 * t5 / 0.256e3);
}



double
AssLegFunction::P_9_9_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  t1 = sqrt(0.12155e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = sqrt(t3);
  return (0.9e1 / 0.256e3 * t1 * t6 * t4 * t3 * x);
}

double
AssLegFunction::P_10_0(const double x) const
{
  double t2;
  double t3;
  double t5;
  double t15;
  double t8;
  double t1;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t5  = t1 - 0.1e1;
  t8  = t5 * t5;
  t15 = t8 * t8;
  return (t3 * t1 + 0.45e2 / 0.2e1 * t3 * t5 + 0.315e3 / 0.4e1 * t8 * t2 * t1 +
          0.525e3 / 0.8e1 * t8 * t5 * t2 + 0.1575e4 / 0.128e3 * t15 * t1 +
          0.63e2 / 0.256e3 * t15 * t5);
}

double
AssLegFunction::P_10_0_Deriv(const double x) const
{
  double t18;
  double t2;
  double t1;
  double t3;
  double t7;
  double t11;
  double t6;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t6  = t1 - 0.1e1;
  t7  = t1 * x;
  t11 = t6 * t6;
  t18 = t11 * t11;
  return (0.55e2 * t3 * x + 0.495e3 * t6 * t2 * t7 +
          0.3465e4 / 0.4e1 * t11 * t2 * x + 0.5775e4 / 0.16e2 * t11 * t6 * t7 +
          0.3465e4 / 0.128e3 * t18 * x);
}



double
AssLegFunction::P_10_1(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t6;
  double t7;
  double t10;
  double t11;
  double t15;
  double t1;
  double t22;
  t1  = sqrt(0.110e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t22 = t15 * t15;
  return (-t1 * t4 *
          (0.204374016000e12 * t7 * x + 0.1839366144000e13 * t10 * t6 * t11 +
           0.3218890752000e13 * t15 * t6 * x +
           0.1341204480000e13 * t15 * t10 * t11 + 0.100590336000e12 * t22 * x) /
          0.408748032000e12);
}



double
AssLegFunction::P_10_1_Deriv(const double x) const
{
  double t16;
  double t7;
  double t20;
  double t8;
  double t1;
  double t23;
  double t2;
  double t11;
  double t3;
  double t4;
  double t12;
  t1  = sqrt(0.110e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t23 = t16 * t16;
  return (t1 / t4 *
            (0.204374016000e12 * t8 * x + 0.1839366144000e13 * t11 * t7 * t12 +
             0.3218890752000e13 * t16 * t7 * x +
             0.1341204480000e13 * t20 * t12 + 0.100590336000e12 * t23 * x) *
            x / 0.408748032000e12 -
          t1 * t4 *
            (0.5518098432000e13 * t8 + 0.25751126016000e14 * t11 * t7 * t2 +
             0.24141680640000e14 * t16 * t7 + 0.4828336128000e13 * t20 * t2 +
             0.100590336000e12 * t23) /
            0.408748032000e12);
}



double
AssLegFunction::P_10_2(const double x) const
{
  double t18;
  double t8;
  double t2;
  double t3;
  double t5;
  double t6;
  double t12;
  double t1;
  t1  = sqrt(0.330e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t8  = -t3;
  t12 = t8 * t8;
  t18 = t12 * t12;
  return (t1 * t3 *
          (0.5518098432000e13 * t6 + 0.25751126016000e14 * t8 * t5 * t2 +
           0.24141680640000e14 * t12 * t5 + 0.4828336128000e13 * t12 * t8 * t2 +
           0.100590336000e12 * t18) /
          0.7357464576000e13);
}



double
AssLegFunction::P_10_2_Deriv(const double x) const
{
  double t3;
  double t17;
  double t11;
  double t4;
  double t7;
  double t1;
  double t5;
  double t24;
  double t14;
  t1  = sqrt(0.330e3);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t3 - 0.1e1;
  t11 = t7 * t7;
  t14 = t11 * t7;
  t17 = t11 * t11;
  t24 = t3 * x;
  return (-t1 * x *
            (0.5518098432000e13 * t5 + 0.25751126016000e14 * t7 * t4 * t3 +
             0.24141680640000e14 * t11 * t4 + 0.4828336128000e13 * t14 * t3 +
             0.100590336000e12 * t17) /
            0.3678732288000e13 -
          t1 * t7 *
            (0.95647039488000e14 * t4 * t24 +
             0.251073478656000e15 * t7 * t4 * x +
             0.125536739328000e15 * t11 * t24 + 0.10461394944000e14 * t14 * x) /
            0.7357464576000e13);
}



double
AssLegFunction::P_10_3(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t3;
  double t15;
  double t7;
  double t8;
  double t11;
  t1  = sqrt(0.2145e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (-t1 * t4 * t3 *
          (0.95647039488000e14 * t8 * t7 + 0.251073478656000e15 * t11 * t8 * x +
           0.125536739328000e15 * t15 * t7 +
           0.10461394944000e14 * t15 * t11 * x) /
          0.191294078976000e15);
}



double
AssLegFunction::P_10_3_Deriv(const double x) const
{
  double t10;
  double t14;
  double t6;
  double t1;
  double t17;
  double t2;
  double t7;
  double t3;
  double t4;
  t1  = sqrt(0.2145e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  return (
    t1 * t4 *
      (0.95647039488000e14 * t7 * t6 + 0.251073478656000e15 * t10 * t7 * x +
       0.125536739328000e15 * t14 * t6 + 0.10461394944000e14 * t17 * x) *
      x / 0.63764692992000e14 -
    t1 * t4 * t3 *
      (0.1171676233728000e16 * t7 * t2 + 0.1757514350592000e16 * t10 * t7 +
       0.439378587648000e15 * t14 * t2 + 0.10461394944000e14 * t17) /
      0.191294078976000e15);
}



double
AssLegFunction::P_10_4(const double x) const
{
  double t2;
  double t3;
  double t12;
  double t4;
  double t9;
  double t6;
  double t1;
  t1  = sqrt(0.4290e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (t1 * t4 *
          (0.1171676233728000e16 * t6 * t2 + 0.1757514350592000e16 * t9 * t6 +
           0.439378587648000e15 * t12 * t2 + 0.10461394944000e14 * t12 * t9) /
          0.2678117105664000e16);
}



double
AssLegFunction::P_10_4_Deriv(const double x) const
{
  double t8;
  double t1;
  double t5;
  double t20;
  double t2;
  double t3;
  double t11;
  t1  = sqrt(0.4290e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  t20 = t3 * t3;
  return (-t1 * t3 *
            (0.1171676233728000e16 * t5 * t2 + 0.1757514350592000e16 * t5 * t8 +
             0.439378587648000e15 * t11 * t2 + 0.10461394944000e14 * t11 * t8) *
            x / 0.669529276416000e15 +
          t1 * t20 *
            (0.10545086103552000e17 * t5 * x +
             0.8787571752960000e16 * t8 * t2 * x +
             0.941525544960000e15 * t11 * x) /
            0.2678117105664000e16);
}



double
AssLegFunction::P_10_5(const double x) const
{
  double t11;
  double t5;
  double t4;
  double t2;
  double t3;
  double t8;
  double t15;
  double t1;
  t1  = sqrt(0.429e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (-t1 * t5 * t4 *
          (0.10545086103552000e17 * t8 * x +
           0.8787571752960000e16 * t11 * t2 * x +
           0.941525544960000e15 * t15 * x) /
          0.8034351316992000e16);
}



double
AssLegFunction::P_10_5_Deriv(const double x) const
{
  double t7;
  double t21;
  double t14;
  double t1;
  double t2;
  double t10;
  double t3;
  double t4;
  t1  = sqrt(0.429e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  t21 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.10545086103552000e17 * t7 * x +
             0.8787571752960000e16 * t10 * t2 * x +
             0.941525544960000e15 * t14 * x) *
            x / 0.1606870263398400e16 -
          t1 * t4 * t21 *
            (0.70300574023680000e17 * t7 + 0.30128817438720000e17 * t10 * t2 +
             0.941525544960000e15 * t14) /
            0.8034351316992000e16);
}



double
AssLegFunction::P_10_6(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t9;
  double t12;
  double t7;
  t1  = sqrt(0.2145e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (t1 * t4 * t3 *
          (0.70300574023680000e17 * t7 + 0.30128817438720000e17 * t9 * t2 +
           0.941525544960000e15 * t12) /
          0.160687026339840000e18);
}



double
AssLegFunction::P_10_6_Deriv(const double x) const
{
  double t6;
  double t2;
  double t8;
  double t3;
  double t11;
  double t4;
  double t1;
  t1  = sqrt(0.2145e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t8  = -t3;
  t11 = t8 * t8;
  return (
    -t1 * t4 *
      (0.70300574023680000e17 * t6 + 0.30128817438720000e17 * t2 * t8 +
       0.941525544960000e15 * t11) *
      x / 0.26781171056640000e17 +
    t1 * t4 * t3 *
      (0.341459930972160000e18 * t2 * x + 0.64023737057280000e17 * t8 * x) /
      0.160687026339840000e18);
}



double
AssLegFunction::P_10_7(const double x) const
{
  double t4;
  double t6;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.36465e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = sqrt(t3);
  return (-t1 * t6 * t4 * t3 *
          (0.341459930972160000e18 * t2 * x - 0.64023737057280000e17 * t3 * x) /
          0.5463358895554560000e19);
}



double
AssLegFunction::P_10_7_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.36465e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = sqrt(t3);
  return (
    t1 * t5 * t4 *
      (0.341459930972160000e18 * t2 * x - 0.64023737057280000e17 * t3 * x) * x /
      0.780479842222080000e18 -
    t1 * t5 * t4 * t3 *
      (-0.64023737057280000e17 + 0.1216451004088320000e19 * t2) /
      0.5463358895554560000e19);
}



double
AssLegFunction::P_10_8(const double x) const
{
  double t4;
  double t5;
  double t2;
  double t1;
  t1 = sqrt(0.24310e5);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (t1 * t5 * (-0.64023737057280000e17 + 0.1216451004088320000e19 * t2) /
          0.32780153373327360000e20);
}



double
AssLegFunction::P_10_8_Deriv(const double x) const
{
  double t1;
  double t2;
  double t12;
  double t4;
  double t3;
  t1  = sqrt(0.24310e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t12 = t4 * t4;
  return (-t1 * t4 * t3 *
            (-0.64023737057280000e17 + 0.1216451004088320000e19 * t2) * x /
            0.4097519171665920000e19 +
          0.19e2 / 0.256e3 * t1 * t12 * x);
}



double
AssLegFunction::P_10_9(const double x) const
{
  double t5;
  double t6;
  double t2;
  double t3;
  double t4;
  double t1;
  t1 = sqrt(0.230945e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (-t1 * t6 * t5 * x / 0.256e3);
}



double
AssLegFunction::P_10_9_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t6;
  double t1;
  double t11;
  t1  = sqrt(0.230945e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t11 = t4 * t4;
  return (0.9e1 / 0.256e3 * t1 * t6 * t4 * t3 * t2 - t1 * t6 * t11 / 0.256e3);
}



double
AssLegFunction::P_10_10(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t3;
  double t5;
  t1 = sqrt(0.46189e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (t1 * t5 * t3 / 0.512e3);
}



double
AssLegFunction::P_10_10_Deriv(const double x) const
{
  double t5;
  double t1;
  double t4;
  double t2;
  t1 = sqrt(0.46189e5);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (-0.5e1 / 0.256e3 * t1 * t5 * x);
}

double
AssLegFunction::P_11_0(const double x) const
{
  double t6;
  double t1;
  double t18;
  double t2;
  double t4;
  double t3;
  double t10;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t3;
  t6  = t1 - 0.1e1;
  t10 = t6 * t6;
  t18 = t10 * t10;
  return (t4 * t2 + 0.55e2 / 0.2e1 * t6 * t4 * x +
          0.495e3 / 0.4e1 * t10 * t3 * t2 +
          0.1155e4 / 0.8e1 * t10 * t6 * t3 * x + 0.5775e4 / 0.128e3 * t18 * t2 +
          0.693e3 / 0.256e3 * t18 * t6 * x);
}

double
AssLegFunction::P_11_0_Deriv(const double x) const
{
  double t16;
  double t1;
  double t6;
  double t9;
  double t2;
  double t3;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t6  = t1 - 0.1e1;
  t9  = t6 * t6;
  t16 = t9 * t9;
  return (0.66e2 * t3 * t1 + 0.1485e4 / 0.2e1 * t6 * t3 +
          0.3465e4 / 0.2e1 * t9 * t2 * t1 + 0.17325e5 / 0.16e2 * t9 * t6 * t2 +
          0.10395e5 / 0.64e2 * t16 * t1 + 0.693e3 / 0.256e3 * t16 * t6);
}



double
AssLegFunction::P_11_1(const double x) const
{
  double t6;
  double t20;
  double t1;
  double t10;
  double t7;
  double t2;
  double t3;
  double t4;
  double t13;
  t1  = sqrt(0.33e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t13 = t10 * t10;
  t20 = t13 * t13;
  return (-t1 * t4 *
          (0.5395474022400e13 * t7 * t2 + 0.60699082752000e14 * t10 * t7 +
           0.141631193088000e15 * t13 * t6 * t2 +
           0.88519495680000e14 * t13 * t10 * t6 +
           0.13277924352000e14 * t20 * t2 + 0.221298739200e12 * t20 * t10) /
          0.5395474022400e13);
}



double
AssLegFunction::P_11_1_Deriv(const double x) const
{
  double t11;
  double t14;
  double t7;
  double t8;
  double t32;
  double t1;
  double t18;
  double t21;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.33e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t18 = t14 * t11;
  t21 = t14 * t14;
  t32 = t2 * x;
  return (
    t1 / t4 *
      (0.5395474022400e13 * t8 * t2 + 0.60699082752000e14 * t11 * t8 +
       0.141631193088000e15 * t14 * t7 * t2 + 0.88519495680000e14 * t18 * t7 +
       0.13277924352000e14 * t21 * t2 + 0.221298739200e12 * t21 * t11) *
      x / 0.5395474022400e13 -
    t1 * t4 *
      (0.175352905728000e15 * t8 * x + 0.1052117434368000e16 * t11 * t7 * t32 +
       0.1380904132608000e16 * t14 * t7 * x + 0.460301377536000e15 * t18 * t32 +
       0.28768836096000e14 * t21 * x) /
      0.5395474022400e13);
}



double
AssLegFunction::P_11_2(const double x) const
{
  double t14;
  double t1;
  double t2;
  double t3;
  double t5;
  double t6;
  double t21;
  double t9;
  double t10;
  t1  = sqrt(0.4290e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t10 = t2 * x;
  t14 = t9 * t9;
  t21 = t14 * t14;
  return (
    t1 * t3 *
    (0.175352905728000e15 * t6 * x + 0.1052117434368000e16 * t9 * t5 * t10 +
     0.1380904132608000e16 * t14 * t5 * x +
     0.460301377536000e15 * t14 * t9 * t10 + 0.28768836096000e14 * t21 * x) /
    0.701411622912000e15);
}



double
AssLegFunction::P_11_2_Deriv(const double x) const
{
  double t13;
  double t17;
  double t1;
  double t20;
  double t9;
  double t8;
  double t3;
  double t4;
  double t5;
  t1  = sqrt(0.4290e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t3 - 0.1e1;
  t9  = t3 * x;
  t13 = t8 * t8;
  t17 = t13 * t8;
  t20 = t13 * t13;
  return (
    -t1 * x *
      (0.175352905728000e15 * t5 * x + 0.1052117434368000e16 * t8 * t4 * t9 +
       0.1380904132608000e16 * t13 * t4 * x + 0.460301377536000e15 * t17 * t9 +
       0.28768836096000e14 * t20 * x) /
      0.350705811456000e15 -
    t1 * t8 *
      (0.3682411020288000e16 * t5 + 0.12888438571008000e17 * t8 * t4 * t3 +
       0.9666328928256000e16 * t13 * t4 + 0.1611054821376000e16 * t17 * t3 +
       0.28768836096000e14 * t20) /
      0.701411622912000e15);
}



double
AssLegFunction::P_11_3(const double x) const
{
  double t3;
  double t4;
  double t7;
  double t2;
  double t8;
  double t14;
  double t20;
  double t10;
  double t1;
  t1  = sqrt(0.15015e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t20 = t14 * t14;
  return (-t1 * t4 * t3 *
          (0.3682411020288000e16 * t8 + 0.12888438571008000e17 * t10 * t7 * t2 +
           0.9666328928256000e16 * t14 * t7 +
           0.1611054821376000e16 * t14 * t10 * t2 + 0.28768836096000e14 * t20) /
          0.14729644081152000e17);
}



double
AssLegFunction::P_11_3_Deriv(const double x) const
{
  double t7;
  double t27;
  double t2;
  double t13;
  double t9;
  double t3;
  double t4;
  double t16;
  double t19;
  double t1;
  double t6;
  t1  = sqrt(0.15015e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t9  = -t3;
  t13 = t9 * t9;
  t16 = t13 * t9;
  t19 = t13 * t13;
  t27 = t2 * x;
  return (
    t1 * t4 *
      (0.3682411020288000e16 * t7 + 0.12888438571008000e17 * t9 * t6 * t2 +
       0.9666328928256000e16 * t13 * t6 + 0.1611054821376000e16 * t16 * t2 +
       0.28768836096000e14 * t19) *
      x / 0.4909881360384000e16 -
    t1 * t4 * t3 *
      (0.55236165304320000e17 * t6 * t27 +
       0.115995947139072000e18 * t9 * t6 * x +
       0.48331644641280000e17 * t13 * t27 + 0.3452260331520000e16 * t16 * x) /
      0.14729644081152000e17);
}



double
AssLegFunction::P_11_4(const double x) const
{
  double t1;
  double t10;
  double t14;
  double t2;
  double t3;
  double t4;
  double t6;
  double t7;
  t1  = sqrt(0.2002e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * x;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (t1 * t4 *
          (0.55236165304320000e17 * t7 * t6 +
           0.115995947139072000e18 * t10 * t7 * x +
           0.48331644641280000e17 * t14 * t6 +
           0.3452260331520000e16 * t14 * t10 * x) /
          0.58918576324608000e17);
}



double
AssLegFunction::P_11_4_Deriv(const double x) const
{
  double t1;
  double t5;
  double t2;
  double t3;
  double t13;
  double t16;
  double t9;
  double t6;
  double t23;
  t1  = sqrt(0.2002e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  t16 = t13 * t9;
  t23 = t3 * t3;
  return (
    -t1 * t3 *
      (0.55236165304320000e17 * t6 * t5 +
       0.115995947139072000e18 * t9 * t6 * x +
       0.48331644641280000e17 * t13 * t5 + 0.3452260331520000e16 * t16 * x) *
      x / 0.14729644081152000e17 +
    t1 * t23 *
      (0.618645051408384000e18 * t6 * t2 + 0.773306314260480000e18 * t6 * t9 +
       0.165708495912960000e18 * t13 * t2 + 0.3452260331520000e16 * t16) /
      0.58918576324608000e17);
}



double
AssLegFunction::P_11_5(const double x) const
{
  double t5;
  double t14;
  double t8;
  double t1;
  double t2;
  double t11;
  double t4;
  double t3;
  t1  = sqrt(0.286e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (
    -t1 * t5 * t4 *
    (0.618645051408384000e18 * t8 * t2 + 0.773306314260480000e18 * t11 * t8 +
     0.165708495912960000e18 * t14 * t2 + 0.3452260331520000e16 * t14 * t11) /
    0.235674305298432000e18);
}



double
AssLegFunction::P_11_5_Deriv(const double x) const
{
  double t3;
  double t4;
  double t10;
  double t22;
  double t13;
  double t7;
  double t1;
  double t2;
  t1  = sqrt(0.286e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  t22 = t3 * t3;
  return (
    t1 * t4 * t3 *
      (0.618645051408384000e18 * t7 * t2 + 0.773306314260480000e18 * t7 * t10 +
       0.165708495912960000e18 * t13 * t2 + 0.3452260331520000e16 * t10 * t13) *
      x / 0.47134861059686400e17 -
    t1 * t4 * t22 *
      (0.5258482936971264000e19 * t7 * x +
       0.3756059240693760000e19 * t10 * t2 * x +
       0.352130553815040000e18 * t13 * x) /
      0.235674305298432000e18);
}



double
AssLegFunction::P_11_6(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t7;
  double t10;
  double t14;
  double t1;
  t1  = sqrt(0.7293e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (t1 * t4 * t3 *
          (0.5258482936971264000e19 * t7 * x +
           0.3756059240693760000e19 * t10 * t2 * x +
           0.352130553815040000e18 * t14 * x) /
          0.12019389570220032000e20);
}



double
AssLegFunction::P_11_6_Deriv(const double x) const
{
  double t6;
  double t1;
  double t9;
  double t2;
  double t13;
  double t3;
  double t4;
  t1  = sqrt(0.7293e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t9  = -t3;
  t13 = t9 * t9;
  return (-t1 * t4 *
            (0.5258482936971264000e19 * t6 * x +
             0.3756059240693760000e19 * t9 * t2 * x +
             0.352130553815040000e18 * t13 * x) *
            x / 0.2003231595036672000e19 +
          t1 * t4 * t3 *
            (0.33804533166243840000e20 * t6 +
             0.12676699937341440000e20 * t9 * t2 +
             0.352130553815040000e18 * t13) /
            0.12019389570220032000e20);
}



double
AssLegFunction::P_11_7(const double x) const
{
  double t9;
  double t1;
  double t2;
  double t3;
  double t14;
  double t11;
  double t4;
  double t6;
  t1  = sqrt(0.72930e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (-t1 * t6 * t4 * t3 *
          (0.33804533166243840000e20 * t9 +
           0.12676699937341440000e20 * t11 * t2 +
           0.352130553815040000e18 * t14) /
          0.360581687106600960000e21);
}



double
AssLegFunction::P_11_7_Deriv(const double x) const
{
  double t2;
  double t4;
  double t3;
  double t8;
  double t10;
  double t1;
  double t13;
  double t5;
  t1  = sqrt(0.72930e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t5 * t4 *
            (0.33804533166243840000e20 * t8 +
             0.12676699937341440000e20 * t10 * t2 +
             0.352130553815040000e18 * t13) *
            x / 0.51511669586657280000e20 -
          t1 * t5 * t4 * t3 *
            (0.160571532539658240000e21 * t2 * x +
             0.26761922089943040000e20 * t10 * x) /
            0.360581687106600960000e21);
}



double
AssLegFunction::P_11_8(const double x) const
{
  double t5;
  double t1;
  double t3;
  double t2;
  double t4;
  t1 = sqrt(0.1385670e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (
    t1 * t5 *
    (0.160571532539658240000e21 * t2 * x - 0.26761922089943040000e20 * t3 * x) /
    0.13702104110050836480000e23);
}



double
AssLegFunction::P_11_8_Deriv(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t3;
  double t16;
  t1  = sqrt(0.1385670e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t16 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.160571532539658240000e21 * t2 * x -
             0.26761922089943040000e20 * t3 * x) *
            x / 0.1712763013756354560000e22 +
          t1 * t16 *
            (0.562000363888803840000e21 * t2 - 0.26761922089943040000e20) /
            0.13702104110050836480000e23);
}



double
AssLegFunction::P_11_9(const double x) const
{
  double t5;
  double t4;
  double t6;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.92378e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (-t1 * t6 * t5 *
          (0.562000363888803840000e21 * t2 - 0.26761922089943040000e20) /
          0.27404208220101672960000e23);
}



double
AssLegFunction::P_11_9_Deriv(const double x) const
{
  double t14;
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  t1  = sqrt(0.92378e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t14 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.562000363888803840000e21 * t2 - 0.26761922089943040000e20) * x /
            0.3044912024455741440000e22 -
          0.21e2 / 0.512e3 * t1 * t6 * t14 * x);
}



double
AssLegFunction::P_11_10(const double x) const
{
  double t3;
  double t4;
  double t1;
  double t5;
  double t2;
  t1 = sqrt(0.969969e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (t1 * t5 * t3 * x / 0.512e3);
}



double
AssLegFunction::P_11_10_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t5;
  double t4;
  t1 = sqrt(0.969969e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-0.5e1 / 0.256e3 * t1 * t5 * t2 + t1 * t5 * t3 / 0.512e3);
}



double
AssLegFunction::P_11_11(const double x) const
{
  double t7;
  double t1;
  double t3;
  double t2;
  double t4;
  double t5;
  t1 = sqrt(0.176358e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-t1 * t7 * t5 * t3 / 0.1024e4);
}



double
AssLegFunction::P_11_11_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  t1 = sqrt(0.176358e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (0.11e2 / 0.1024e4 * t1 * t6 * t5 * x);
}

double
AssLegFunction::P_12_0(const double x) const
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t5;
  double t9;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t5  = t1 - 0.1e1;
  t9  = t5 * t5;
  t16 = t9 * t9;
  return (t3 * t2 + 0.33e2 * t5 * t3 * t1 + 0.1485e4 / 0.8e1 * t9 * t3 +
          0.1155e4 / 0.4e1 * t9 * t5 * t2 * t1 +
          0.17325e5 / 0.128e3 * t16 * t2 + 0.2079e4 / 0.128e3 * t16 * t5 * t1 +
          0.231e3 / 0.1024e4 * t16 * t9);
}

double
AssLegFunction::P_12_0_Deriv(const double x) const
{
  double t11;
  double t19;
  double t3;
  double t4;
  double t1;
  double t2;
  double t7;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t3;
  t7  = t1 - 0.1e1;
  t11 = t7 * t7;
  t19 = t11 * t11;
  return (0.78e2 * t2 * t4 + 0.2145e4 / 0.2e1 * t7 * t4 * x +
          0.6435e4 / 0.2e1 * t11 * t3 * t2 +
          0.45045e5 / 0.16e2 * t11 * t7 * t3 * x +
          0.45045e5 / 0.64e2 * t19 * t2 + 0.9009e4 / 0.256e3 * t19 * t7 * x);
}



double
AssLegFunction::P_12_1(const double x) const
{
  double t2;
  double t3;
  double t23;
  double t4;
  double t1;
  double t8;
  double t7;
  double t6;
  double t15;
  double t11;
  t1  = sqrt(0.39e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t23 = t15 * t15;
  return (
    -t1 * t4 *
    (0.153035263180800e15 * t8 * t6 + 0.2104234868736000e16 * t11 * t8 * x +
     0.6312704606208000e16 * t15 * t7 * t6 +
     0.5523616530432000e16 * t15 * t11 * t7 * x +
     0.1380904132608000e16 * t23 * t6 + 0.69045206630400e14 * t23 * t11 * x) /
    0.153035263180800e15);
}



double
AssLegFunction::P_12_1_Deriv(const double x) const
{
  double t27;
  double t7;
  double t9;
  double t8;
  double t20;
  double t24;
  double t3;
  double t4;
  double t1;
  double t12;
  double t2;
  double t16;
  t1  = sqrt(0.39e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t24 = t16 * t16;
  t27 = t24 * t12;
  return (
    t1 / t4 *
      (0.153035263180800e15 * t7 * t9 + 0.2104234868736000e16 * t9 * t12 * x +
       0.6312704606208000e16 * t16 * t8 * t7 +
       0.5523616530432000e16 * t20 * t8 * x + 0.1380904132608000e16 * t24 * t7 +
       0.69045206630400e14 * t27 * x) *
      x / 0.153035263180800e15 -
    t1 * t4 *
      (0.5891857632460800e16 * t9 * t2 + 0.44188932243456000e17 * t12 * t9 +
       0.77330631426048000e17 * t16 * t8 * t2 +
       0.38665315713024000e17 * t20 * t8 + 0.4833164464128000e16 * t24 * t2 +
       0.69045206630400e14 * t27) /
      0.153035263180800e15);
}



double
AssLegFunction::P_12_2(const double x) const
{
  double t2;
  double t9;
  double t19;
  double t5;
  double t3;
  double t1;
  double t6;
  double t12;
  t1  = sqrt(0.6006e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t12 = t9 * t9;
  t19 = t12 * t12;
  return (t1 * t3 *
          (0.5891857632460800e16 * t6 * t2 + 0.44188932243456000e17 * t9 * t6 +
           0.77330631426048000e17 * t12 * t5 * t2 +
           0.38665315713024000e17 * t12 * t9 * t5 +
           0.4833164464128000e16 * t19 * t2 + 0.69045206630400e14 * t9 * t19) /
          0.23567430529843200e17);
}



double
AssLegFunction::P_12_2_Deriv(const double x) const
{
  double t11;
  double t15;
  double t30;
  double t1;
  double t8;
  double t4;
  double t5;
  double t18;
  double t3;
  t1  = sqrt(0.6006e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t3 - 0.1e1;
  t11 = t8 * t8;
  t15 = t11 * t8;
  t18 = t11 * t11;
  t30 = t3 * x;
  return (
    -t1 * x *
      (0.5891857632460800e16 * t5 * t3 + 0.44188932243456000e17 * t8 * t5 +
       0.77330631426048000e17 * t11 * t4 * t3 +
       0.38665315713024000e17 * t15 * t4 + 0.4833164464128000e16 * t18 * t3 +
       0.69045206630400e14 * t18 * t8) /
      0.11783715264921600e17 -
    t1 * t8 *
      (0.147296440811520000e18 * t5 * x +
       0.662833983651840000e18 * t8 * t4 * t30 +
       0.695975682834432000e18 * t11 * t4 * x +
       0.193326578565120000e18 * t15 * t30 + 0.10356780994560000e17 * t18 * x) /
      0.23567430529843200e17);
}



double
AssLegFunction::P_12_3(const double x) const
{
  double t23;
  double t7;
  double t8;
  double t11;
  double t12;
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.1001e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t23 = t16 * t16;
  return (-t1 * t4 * t3 *
          (0.147296440811520000e18 * t8 * x +
           0.662833983651840000e18 * t11 * t7 * t12 +
           0.695975682834432000e18 * t16 * t7 * x +
           0.193326578565120000e18 * t16 * t11 * t12 +
           0.10356780994560000e17 * t23 * x) /
          0.117837152649216000e18);
}



double
AssLegFunction::P_12_3_Deriv(const double x) const
{
  double t22;
  double t11;
  double t7;
  double t4;
  double t15;
  double t6;
  double t1;
  double t3;
  double t2;
  double t10;
  double t19;
  t1  = sqrt(0.1001e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t22 = t15 * t15;
  return (
    t1 * t4 *
      (0.147296440811520000e18 * t7 * x +
       0.662833983651840000e18 * t10 * t6 * t11 +
       0.695975682834432000e18 * t15 * t6 * x +
       0.193326578565120000e18 * t19 * t11 + 0.10356780994560000e17 * t22 * x) *
      x / 0.39279050883072000e17 -
    t1 * t4 * t3 *
      (0.2651335934607360000e19 * t7 +
       0.7423740616900608000e19 * t10 * t6 * t2 +
       0.4639837885562880000e19 * t15 * t6 +
       0.662833983651840000e18 * t19 * t2 + 0.10356780994560000e17 * t22) /
      0.117837152649216000e18);
}



double
AssLegFunction::P_12_4(const double x) const
{
  double t4;
  double t9;
  double t1;
  double t19;
  double t3;
  double t2;
  double t6;
  double t7;
  double t13;
  t1  = sqrt(0.1001e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t9  = -t3;
  t13 = t9 * t9;
  t19 = t13 * t13;
  return (
    t1 * t4 *
    (0.2651335934607360000e19 * t7 + 0.7423740616900608000e19 * t9 * t6 * t2 +
     0.4639837885562880000e19 * t13 * t6 +
     0.662833983651840000e18 * t13 * t9 * t2 + 0.10356780994560000e17 * t19) /
    0.1414045831790592000e19);
}



double
AssLegFunction::P_12_4_Deriv(const double x) const
{
  double t12;
  double t15;
  double t1;
  double t5;
  double t6;
  double t2;
  double t3;
  double t24;
  double t8;
  double t26;
  double t18;
  t1  = sqrt(0.1001e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t8  = -t3;
  t12 = t8 * t8;
  t15 = t12 * t8;
  t18 = t12 * t12;
  t24 = t3 * t3;
  t26 = t2 * x;
  return (
    -t1 * t3 *
      (0.2651335934607360000e19 * t6 + 0.7423740616900608000e19 * t8 * t5 * t2 +
       0.4639837885562880000e19 * t12 * t5 +
       0.662833983651840000e18 * t15 * t2 + 0.10356780994560000e17 * t18) *
      x / 0.353511457947648000e18 +
    t1 * t24 *
      (0.36058168710660096000e20 * t5 * t26 +
       0.63101795243655168000e20 * t8 * t5 * x +
       0.22536355444162560000e20 * t12 * t26 +
       0.1408522215260160000e19 * t15 * x) /
      0.1414045831790592000e19);
}



double
AssLegFunction::P_12_5(const double x) const
{
  double t8;
  double t5;
  double t9;
  double t1;
  double t16;
  double t12;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.34034e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * x;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (-t1 * t5 * t4 *
          (0.36058168710660096000e20 * t9 * t8 +
           0.63101795243655168000e20 * t12 * t9 * x +
           0.22536355444162560000e20 * t16 * t8 +
           0.1408522215260160000e19 * t16 * t12 * x) /
          0.96155116561760256000e20);
}



double
AssLegFunction::P_12_5_Deriv(const double x) const
{
  double t15;
  double t18;
  double t25;
  double t1;
  double t2;
  double t3;
  double t4;
  double t7;
  double t8;
  double t11;
  t1  = sqrt(0.34034e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t25 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.36058168710660096000e20 * t8 * t7 +
             0.63101795243655168000e20 * t11 * t8 * x +
             0.22536355444162560000e20 * t15 * t7 +
             0.1408522215260160000e19 * t18 * x) *
            x / 0.19231023312352051200e20 -
          t1 * t4 * t25 *
            (0.378610771461931008000e21 * t8 * t2 +
             0.405654397994926080000e21 * t11 * t8 +
             0.76060199624048640000e20 * t2 * t15 +
             0.1408522215260160000e19 * t18) /
            0.96155116561760256000e20);
}



double
AssLegFunction::P_12_6(const double x) const
{
  double t10;
  double t2;
  double t4;
  double t3;
  double t13;
  double t7;
  double t1;
  t1  = sqrt(0.2431e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t4 * t3 *
          (0.378610771461931008000e21 * t7 * t2 +
           0.405654397994926080000e21 * t10 * t7 +
           0.76060199624048640000e20 * t13 * t2 +
           0.1408522215260160000e19 * t13 * t10) /
          0.288465349685280768000e21);
}



double
AssLegFunction::P_12_6_Deriv(const double x) const
{
  double t1;
  double t9;
  double t3;
  double t2;
  double t12;
  double t4;
  double t6;
  t1  = sqrt(0.2431e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (-t1 * t4 *
            (0.378610771461931008000e21 * t2 * t6 +
             0.405654397994926080000e21 * t9 * t6 +
             0.76060199624048640000e20 * t12 * t2 +
             0.1408522215260160000e19 * t12 * t9) *
            x / 0.48077558280880128000e20 +
          t1 * t4 * t3 *
            (0.3082973424761438208000e22 * t6 * x +
             0.1926858390475898880000e22 * t9 * t2 * x +
             0.160571532539658240000e21 * t12 * x) /
            0.288465349685280768000e21);
}



double
AssLegFunction::P_12_7(const double x) const
{
  double t1;
  double t12;
  double t16;
  double t6;
  double t9;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.277134e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (-t1 * t6 * t4 * t3 *
          (0.3082973424761438208000e22 * t9 * x +
           0.1926858390475898880000e22 * t12 * t2 * x +
           0.160571532539658240000e21 * t16 * x) /
          0.32885049864122007552000e23);
}



double
AssLegFunction::P_12_7_Deriv(const double x) const
{
  double t11;
  double t4;
  double t15;
  double t8;
  double t5;
  double t2;
  double t3;
  double t1;
  t1  = sqrt(0.277134e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (t1 * t5 * t4 *
            (0.3082973424761438208000e22 * t8 * x +
             0.1926858390475898880000e22 * t11 * t2 * x +
             0.160571532539658240000e21 * t15 * x) *
            x / 0.4697864266303143936000e22 -
          t1 * t5 * t4 * t3 *
            (0.19268583904758988800000e23 * t8 +
             0.6422861301586329600000e22 * t11 * t2 +
             0.160571532539658240000e21 * t15) /
            0.32885049864122007552000e23);
}



double
AssLegFunction::P_12_8(const double x) const
{
  double t4;
  double t9;
  double t5;
  double t3;
  double t12;
  double t2;
  double t1;
  double t7;
  t1  = sqrt(0.277134e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (t1 * t5 *
          (0.19268583904758988800000e23 * t7 +
           0.6422861301586329600000e22 * t9 * t2 +
           0.160571532539658240000e21 * t12) /
          0.328850498641220075520000e24);
}



double
AssLegFunction::P_12_8_Deriv(const double x) const
{
  double t2;
  double t7;
  double t3;
  double t1;
  double t18;
  double t9;
  double t12;
  double t4;
  t1  = sqrt(0.277134e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  t18 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.19268583904758988800000e23 * t7 +
             0.6422861301586329600000e22 * t9 * t2 +
             0.160571532539658240000e21 * t12) *
            x / 0.41106312330152509440000e23 +
          t1 * t18 *
            (0.89920058222208614400000e23 * t2 * x +
             0.13488008733331292160000e23 * t9 * x) /
            0.328850498641220075520000e24);
}



double
AssLegFunction::P_12_9(const double x) const
{
  double t3;
  double t2;
  double t4;
  double t5;
  double t6;
  double t1;
  t1 = sqrt(0.646646e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (-t1 * t6 * t5 *
          (0.89920058222208614400000e23 * t2 * x -
           0.13488008733331292160000e23 * t3 * x) /
          0.4603906980977081057280000e25);
}



double
AssLegFunction::P_12_9_Deriv(const double x) const
{
  double t4;
  double t6;
  double t1;
  double t2;
  double t18;
  double t3;
  t1  = sqrt(0.646646e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t18 = t4 * t4;
  return (
    t1 * t6 * t4 * t3 *
      (0.89920058222208614400000e23 * t2 * x -
       0.13488008733331292160000e23 * t3 * x) *
      x / 0.511545220108564561920000e24 -
    t1 * t6 * t18 *
      (0.310224200866619719680000e24 * t2 - 0.13488008733331292160000e23) /
      0.4603906980977081057280000e25);
}



double
AssLegFunction::P_12_10(const double x) const
{
  double t5;
  double t1;
  double t2;
  double t4;
  double t3;
  t1 = sqrt(0.88179e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (t1 * t5 * t3 *
          (0.310224200866619719680000e24 * t2 - 0.13488008733331292160000e23) /
          0.13811720942931243171840000e26);
}



double
AssLegFunction::P_12_10_Deriv(const double x) const
{
  double t1;
  double t2;
  double t5;
  double t3;
  double t4;
  t1 = sqrt(0.88179e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (
    -t1 * t5 *
      (0.310224200866619719680000e24 * t2 - 0.13488008733331292160000e23) * x /
      0.1381172094293124317184000e25 +
    0.23e2 / 0.512e3 * t1 * t5 * t3 * x);
}



double
AssLegFunction::P_12_11(const double x) const
{
  double t5;
  double t4;
  double t7;
  double t2;
  double t1;
  double t3;
  t1 = sqrt(0.4056234e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-t1 * t7 * t5 * t3 * x / 0.1024e4);
}



double
AssLegFunction::P_12_11_Deriv(const double x) const
{
  double t6;
  double t1;
  double t2;
  double t4;
  double t3;
  double t5;
  t1 = sqrt(0.4056234e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (0.11e2 / 0.1024e4 * t1 * t6 * t5 * t2 - t1 * t6 * t5 * t3 / 0.1024e4);
}



double
AssLegFunction::P_12_12(const double x) const
{
  double t4;
  double t5;
  double t1;
  double t2;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (t1 * t5 * t4 / 0.2048e4);
}



double
AssLegFunction::P_12_12_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t1;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-0.3e1 / 0.512e3 * t1 * t5 * t3 * x);
}

double
AssLegFunction::P_13_0(const double x) const
{
  double t6;
  double t7;
  double t19;
  double t3;
  double t11;
  double t1;
  double t2;
  double t4;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * x;
  t4  = t2 * t2;
  t6  = t1 - 0.1e1;
  t7  = t1 * x;
  t11 = t6 * t6;
  t19 = t11 * t11;
  return (t4 * t3 + 0.39e2 * t6 * t4 * t7 + 0.2145e4 / 0.8e1 * t11 * t4 * x +
          0.2145e4 / 0.4e1 * t11 * t6 * t2 * t7 +
          0.45045e5 / 0.128e3 * t19 * t3 + 0.9009e4 / 0.128e3 * t19 * t6 * t7 +
          0.3003e4 / 0.1024e4 * t19 * t11 * x);
}

double
AssLegFunction::P_13_0_Deriv(const double x) const
{
  double t6;
  double t17;
  double t1;
  double t2;
  double t10;
  double t3;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t6  = t1 - 0.1e1;
  t10 = t6 * t6;
  t17 = t10 * t10;
  return (
    0.91e2 * t3 * t2 + 0.3003e4 / 0.2e1 * t6 * t3 * t1 +
    0.45045e5 / 0.8e1 * t10 * t3 + 0.105105e6 / 0.16e2 * t10 * t6 * t2 * t1 +
    0.315315e6 / 0.128e3 * t17 * t2 + 0.63063e5 / 0.256e3 * t1 * t6 * t17 +
    0.3003e4 / 0.1024e4 * t17 * t10);
}



double
AssLegFunction::P_13_1(const double x) const
{
  double t7;
  double t14;
  double t2;
  double t10;
  double t6;
  double t4;
  double t21;
  double t1;
  double t3;
  t1  = sqrt(0.182e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t21 = t14 * t14;
  return (-t1 * t4 *
          (0.4642069649817600e16 * t7 * t6 +
           0.76594149221990400e17 * t10 * t7 * t2 +
           0.287228059582464000e18 * t14 * t7 +
           0.335099402846208000e18 * t14 * t10 * t6 * t2 +
           0.125662276067328000e18 * t21 * t6 +
           0.12566227606732800e17 * t21 * t10 * t2 +
           0.149597947699200e15 * t21 * t14) /
          0.9284139299635200e16);
}



double
AssLegFunction::P_13_1_Deriv(const double x) const
{
  double t18;
  double t34;
  double t1;
  double t22;
  double t7;
  double t8;
  double t15;
  double t4;
  double t2;
  double t3;
  double t25;
  double t11;
  t1  = sqrt(0.182e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t22 = t15 * t15;
  t25 = t22 * t11;
  t34 = t2 * x;
  return (
    t1 / t4 *
      (0.4642069649817600e16 * t8 * t7 +
       0.76594149221990400e17 * t11 * t8 * t2 +
       0.287228059582464000e18 * t15 * t8 +
       0.335099402846208000e18 * t18 * t7 * t2 +
       0.125662276067328000e18 * t22 * t7 + 0.12566227606732800e17 * t25 * t2 +
       0.149597947699200e15 * t22 * t15) *
      x / 0.9284139299635200e16 -
    t1 * t4 *
      (0.208893134241792000e18 * t8 * t34 +
       0.1914853730549760000e19 * t11 * t8 * x +
       0.4308420893736960000e19 * t15 * t7 * t34 +
       0.3015894625615872000e19 * t18 * t7 * x +
       0.628311380336640000e18 * t22 * t34 + 0.26927630585856000e17 * t25 * x) /
      0.9284139299635200e16);
}



double
AssLegFunction::P_13_2(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t5;
  double t6;
  double t7;
  double t10;
  double t22;
  double t14;
  t1  = sqrt(0.910e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t22 = t14 * t14;
  return (t1 * t3 *
          (0.208893134241792000e18 * t7 * t5 +
           0.1914853730549760000e19 * t10 * t7 * x +
           0.4308420893736960000e19 * t14 * t6 * t5 +
           0.3015894625615872000e19 * t14 * t10 * t6 * x +
           0.628311380336640000e18 * t22 * t5 +
           0.26927630585856000e17 * t22 * t10 * x) /
          0.278524178989056000e18);
}



double
AssLegFunction::P_13_2_Deriv(const double x) const
{
  double t9;
  double t17;
  double t3;
  double t21;
  double t4;
  double t5;
  double t24;
  double t1;
  double t6;
  double t13;
  t1  = sqrt(0.910e3);
  t3  = x * x;
  t4  = t3 * x;
  t5  = t3 * t3;
  t6  = t5 * t5;
  t9  = t3 - 0.1e1;
  t13 = t9 * t9;
  t17 = t13 * t9;
  t21 = t13 * t13;
  t24 = t21 * t9;
  return (
    -t1 * x *
      (0.208893134241792000e18 * t6 * t4 +
       0.1914853730549760000e19 * t9 * t6 * x +
       0.4308420893736960000e19 * t13 * t5 * t4 +
       0.3015894625615872000e19 * t17 * t5 * x +
       0.628311380336640000e18 * t21 * t4 + 0.26927630585856000e17 * t24 * x) /
      0.139262089494528000e18 -
    t1 * t9 *
      (0.6127531937759232000e19 * t6 * t3 +
       0.34467367149895680000e20 * t9 * t6 +
       0.48254314009853952000e20 * t13 * t5 * t3 +
       0.20105964170772480000e20 * t17 * t5 +
       0.2154210446868480000e19 * t21 * t3 + 0.26927630585856000e17 * t24) /
      0.278524178989056000e18);
}



double
AssLegFunction::P_13_3(const double x) const
{
  double t21;
  double t14;
  double t8;
  double t7;
  double t1;
  double t2;
  double t3;
  double t11;
  double t4;
  t1  = sqrt(0.10010e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t21 = t14 * t14;
  return (
    -t1 * t4 * t3 *
    (0.6127531937759232000e19 * t8 * t2 + 0.34467367149895680000e20 * t11 * t8 +
     0.48254314009853952000e20 * t14 * t7 * t2 +
     0.20105964170772480000e20 * t14 * t11 * t7 +
     0.2154210446868480000e19 * t21 * t2 + 0.26927630585856000e17 * t21 * t11) /
    0.12255063875518464000e20);
}



double
AssLegFunction::P_13_3_Deriv(const double x) const
{
  double t17;
  double t4;
  double t1;
  double t20;
  double t13;
  double t6;
  double t7;
  double t3;
  double t2;
  double t10;
  double t33;
  t1  = sqrt(0.10010e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t13 = t10 * t10;
  t17 = t13 * t10;
  t20 = t13 * t13;
  t33 = t2 * x;
  return (t1 * t4 *
            (0.6127531937759232000e19 * t7 * t2 +
             0.34467367149895680000e20 * t10 * t7 +
             0.48254314009853952000e20 * t13 * t6 * t2 +
             0.20105964170772480000e20 * t17 * t6 +
             0.2154210446868480000e19 * t20 * t2 +
             0.26927630585856000e17 * t20 * t10) *
            x / 0.4085021291839488000e19 -
          t1 * t4 * t3 *
            (0.130210053677383680000e21 * t7 * x +
             0.468756193238581248000e21 * t10 * t6 * t33 +
             0.410161669083758592000e21 * t13 * t6 * x +
             0.97657540258037760000e20 * t17 * t33 +
             0.4577697199595520000e19 * t20 * x) /
            0.12255063875518464000e20);
}



double
AssLegFunction::P_13_4(const double x) const
{
  double t11;
  double t6;
  double t15;
  double t7;
  double t1;
  double t10;
  double t2;
  double t22;
  double t4;
  double t3;
  t1  = sqrt(0.17017e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t22 = t15 * t15;
  return (t1 * t4 *
          (0.130210053677383680000e21 * t7 * x +
           0.468756193238581248000e21 * t10 * t6 * t11 +
           0.410161669083758592000e21 * t15 * t6 * x +
           0.97657540258037760000e20 * t15 * t10 * t11 +
           0.4577697199595520000e19 * t22 * x) /
          0.208336085883813888000e21);
}



double
AssLegFunction::P_13_4_Deriv(const double x) const
{
  double t28;
  double t9;
  double t18;
  double t10;
  double t21;
  double t1;
  double t14;
  double t5;
  double t3;
  double t2;
  double t6;
  t1  = sqrt(0.17017e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t10 = t2 * x;
  t14 = t9 * t9;
  t18 = t14 * t9;
  t21 = t14 * t14;
  t28 = t3 * t3;
  return (-t1 * t3 *
            (0.130210053677383680000e21 * t6 * x +
             0.468756193238581248000e21 * t9 * t5 * t10 +
             0.410161669083758592000e21 * t14 * t5 * x +
             0.97657540258037760000e20 * t18 * t10 +
             0.4577697199595520000e19 * t21 * x) *
            x / 0.52084021470953472000e20 +
          t1 * t28 *
            (0.2109402869573615616000e22 * t6 +
             0.4921940029005103104000e22 * t9 * t5 * t2 +
             0.2636753586967019520000e22 * t14 * t5 +
             0.329594198370877440000e21 * t18 * t2 +
             0.4577697199595520000e19 * t21) /
            0.208336085883813888000e21);
}



double
AssLegFunction::P_13_5(const double x) const
{
  double t9;
  double t8;
  double t11;
  double t15;
  double t21;
  double t2;
  double t4;
  double t3;
  double t5;
  double t1;
  t1  = sqrt(0.34034e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t11 = -t3;
  t15 = t11 * t11;
  t21 = t15 * t15;
  return (-t1 * t5 * t4 *
          (0.2109402869573615616000e22 * t9 +
           0.4921940029005103104000e22 * t11 * t8 * t2 +
           0.2636753586967019520000e22 * t15 * t8 +
           0.329594198370877440000e21 * t15 * t11 * t2 +
           0.4577697199595520000e19 * t21) /
          0.3750049545908649984000e22);
}



double
AssLegFunction::P_13_5_Deriv(const double x) const
{
  double t4;
  double t14;
  double t1;
  double t7;
  double t8;
  double t26;
  double t10;
  double t3;
  double t2;
  double t20;
  double t17;
  double t29;
  t1  = sqrt(0.34034e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t20 = t14 * t14;
  t26 = t3 * t3;
  t29 = t2 * x;
  return (t1 * t4 * t3 *
            (0.2109402869573615616000e22 * t8 +
             0.4921940029005103104000e22 * t10 * t7 * t2 +
             0.2636753586967019520000e22 * t7 * t14 +
             0.329594198370877440000e21 * t17 * t2 +
             0.4577697199595520000e19 * t20) *
            x / 0.750009909181729996800e21 -
          t1 * t4 * t26 *
            (0.26719103014599131136000e23 * t7 * t29 +
             0.40078654521898696704000e23 * t10 * t7 * x +
             0.12524579538093342720000e23 * t14 * t29 +
             0.695809974338519040000e21 * t17 * x) /
            0.3750049545908649984000e22);
}



double
AssLegFunction::P_13_6(const double x) const
{
  double t11;
  double t7;
  double t15;
  double t1;
  double t2;
  double t8;
  double t3;
  double t4;
  t1  = sqrt(0.323323e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (t1 * t4 * t3 *
          (0.26719103014599131136000e23 * t8 * t7 +
           0.40078654521898696704000e23 * t11 * t8 * x +
           0.12524579538093342720000e23 * t7 * t15 +
           0.695809974338519040000e21 * t15 * t11 * x) /
          0.142501882744528699392000e24);
}



double
AssLegFunction::P_13_6_Deriv(const double x) const
{
  double t14;
  double t6;
  double t1;
  double t17;
  double t7;
  double t10;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.323323e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * x;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  return (-t1 * t4 *
            (0.26719103014599131136000e23 * t7 * t6 +
             0.40078654521898696704000e23 * t10 * t7 * x +
             0.12524579538093342720000e23 * t14 * t6 +
             0.695809974338519040000e21 * t17 * x) *
            x / 0.23750313790754783232000e23 +
          t1 * t4 * t3 *
            (0.267191030145991311360000e24 * t7 * t2 +
             0.250491590761866854400000e24 * t10 * t7 +
             0.41748598460311142400000e23 * t14 * t2 +
             0.695809974338519040000e21 * t17) /
            0.142501882744528699392000e24);
}



double
AssLegFunction::P_13_7(const double x) const
{
  double t9;
  double t12;
  double t15;
  double t2;
  double t3;
  double t4;
  double t6;
  double t1;
  t1  = sqrt(0.230945e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (-t1 * t6 * t4 * t3 *
          (0.267191030145991311360000e24 * t9 * t2 +
           0.250491590761866854400000e24 * t12 * t9 +
           0.41748598460311142400000e23 * t15 * t2 +
           0.695809974338519040000e21 * t12 * t15) /
          0.1425018827445286993920000e25);
}



double
AssLegFunction::P_13_7_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t8;
  double t11;
  double t14;
  t1  = sqrt(0.230945e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (t1 * t5 * t4 *
            (0.267191030145991311360000e24 * t8 * t2 +
             0.250491590761866854400000e24 * t11 * t8 +
             0.41748598460311142400000e23 * t14 * t2 +
             0.695809974338519040000e21 * t14 * t11) *
            x / 0.203574118206469570560000e24 -
          t1 * t5 * t4 * t3 *
            (0.2104129362399681576960000e25 * t8 * x +
             0.1168960756888711987200000e25 * t11 * t2 * x +
             0.87672056766653399040000e23 * t14 * x) /
            0.1425018827445286993920000e25);
}



double
AssLegFunction::P_13_8(const double x) const
{
  double t5;
  double t7;
  double t1;
  double t10;
  double t2;
  double t14;
  double t3;
  double t4;
  t1  = sqrt(0.3233230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (t1 * t5 *
          (0.2104129362399681576960000e25 * t7 * x +
           0.1168960756888711987200000e25 * t10 * t2 * x +
           0.87672056766653399040000e23 * t14 * x) /
          0.59850790752702053744640000e26);
}



double
AssLegFunction::P_13_8_Deriv(const double x) const
{
  double t7;
  double t14;
  double t2;
  double t10;
  double t3;
  double t4;
  double t21;
  double t1;
  t1  = sqrt(0.3233230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  t21 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.2104129362399681576960000e25 * t7 * x +
             0.1168960756888711987200000e25 * t10 * t2 * x +
             0.87672056766653399040000e23 * t14 * x) *
            x / 0.7481348844087756718080000e25 +
          t1 * t21 *
            (0.12858568325775831859200000e26 * t7 +
             0.3857570497732749557760000e25 * t10 * t2 +
             0.87672056766653399040000e23 * t14) /
            0.59850790752702053744640000e26);
}



double
AssLegFunction::P_13_9(const double x) const
{
  double t14;
  double t11;
  double t9;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  t1  = sqrt(0.29393e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (-t1 * t6 * t5 *
          (0.12858568325775831859200000e26 * t9 +
           0.3857570497732749557760000e25 * t11 * t2 +
           0.87672056766653399040000e23 * t14) /
          0.59850790752702053744640000e26);
}



double
AssLegFunction::P_13_9_Deriv(const double x) const
{
  double t20;
  double t6;
  double t1;
  double t9;
  double t2;
  double t3;
  double t14;
  double t4;
  double t11;
  t1  = sqrt(0.29393e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  t20 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.12858568325775831859200000e26 * t9 +
             0.3857570497732749557760000e25 * t11 * t2 +
             0.87672056766653399040000e23 * t14) *
            x / 0.6650087861411339304960000e25 -
          t1 * t6 * t20 *
            (0.59149414298568826552320000e26 * t2 * x +
             0.8065829222532112711680000e25 * t11 * x) /
            0.59850790752702053744640000e26);
}



double
AssLegFunction::P_13_10(const double x) const
{
  double t4;
  double t5;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (t1 * t5 * t3 *
          (0.59149414298568826552320000e26 * t2 * x -
           0.8065829222532112711680000e25 * t3 * x) /
          0.2753136374624294472253440000e28);
}



double
AssLegFunction::P_13_10_Deriv(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t3;
  double t5;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (
    -t1 * t5 *
      (0.59149414298568826552320000e26 * t2 * x -
       0.8065829222532112711680000e25 * t3 * x) *
      x / 0.275313637462429447225344000e27 +
    t1 * t5 * t3 *
      (0.201645730563302817792000000e27 * t2 - 0.8065829222532112711680000e25) /
      0.2753136374624294472253440000e28);
}



double
AssLegFunction::P_13_11(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t3;
  double t5;
  double t7;
  t1 = sqrt(0.1352078e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (
    -t1 * t7 * t5 * t3 *
    (0.201645730563302817792000000e27 * t2 - 0.8065829222532112711680000e25) /
    0.33037636495491533667041280000e29);
}



double
AssLegFunction::P_13_11_Deriv(const double x) const
{
  double t3;
  double t2;
  double t5;
  double t4;
  double t6;
  double t1;
  t1 = sqrt(0.1352078e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (
    t1 * t6 * t5 *
      (0.201645730563302817792000000e27 * t2 - 0.8065829222532112711680000e25) *
      x / 0.3003421499590139424276480000e28 -
    0.25e2 / 0.2048e4 * t1 * t6 * t5 * t3 * x);
}



double
AssLegFunction::P_13_12(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t5;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (0.5e1 / 0.2048e4 * t1 * t5 * t4 * x);
}



double
AssLegFunction::P_13_12_Deriv(const double x) const
{
  double t1;
  double t5;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.676039e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-0.15e2 / 0.512e3 * t1 * t5 * t3 * t2 +
          0.5e1 / 0.2048e4 * t1 * t5 * t4);
}



double
AssLegFunction::P_13_13(const double x) const
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t7;
  double t2;
  t1 = sqrt(0.104006e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-0.5e1 / 0.4096e4 * t1 * t7 * t5 * t4);
}



double
AssLegFunction::P_13_13_Deriv(const double x) const
{
  double t2;
  double t3;
  double t5;
  double t4;
  double t7;
  double t1;
  t1 = sqrt(0.104006e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (0.65e2 / 0.4096e4 * t1 * t7 * t5 * t3 * x);
}

double
AssLegFunction::P_14_0(const double x) const
{
  double t10;
  double t6;
  double t17;
  double t1;
  double t2;
  double t14;
  double t3;
  double t4;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t1;
  t4  = t2 * t2;
  t6  = t1 - 0.1e1;
  t10 = t6 * t6;
  t14 = t10 * t6;
  t17 = t10 * t10;
  return (
    t4 * t3 + 0.91e2 / 0.2e1 * t6 * t4 * t2 + 0.3003e4 / 0.8e1 * t10 * t4 * t1 +
    0.15015e5 / 0.16e2 * t14 * t4 + 0.105105e6 / 0.128e3 * t17 * t3 +
    0.63063e5 / 0.256e3 * t17 * t6 * t2 +
    0.21021e5 / 0.1024e4 * t17 * t10 * t1 + 0.429e3 / 0.2048e4 * t17 * t14);
}

double
AssLegFunction::P_14_0_Deriv(const double x) const
{
  double t20;
  double t4;
  double t1;
  double t12;
  double t7;
  double t2;
  double t8;
  double t3;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * x;
  t4  = t2 * t2;
  t7  = t1 - 0.1e1;
  t8  = t1 * x;
  t12 = t7 * t7;
  t20 = t12 * t12;
  return (0.105e3 * t4 * t3 + 0.4095e4 / 0.2e1 * t7 * t4 * t8 +
          0.75075e5 / 0.8e1 * t12 * t4 * x +
          0.225225e6 / 0.16e2 * t12 * t7 * t2 * t8 +
          0.945945e6 / 0.128e3 * t20 * t3 +
          0.315315e6 / 0.256e3 * t20 * t7 * t8 +
          0.45045e5 / 0.1024e4 * t20 * t12 * x);
}



double
AssLegFunction::P_14_1(const double x) const
{
  double t8;
  double t3;
  double t11;
  double t16;
  double t4;
  double t12;
  double t24;
  double t1;
  double t7;
  double t2;
  double t6;
  t1  = sqrt(0.210e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * x;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t24 = t16 * t16;
  return (-t1 * t4 *
          (0.149974557917184000e18 * t8 * t7 +
           0.2924503879385088000e19 * t11 * t8 * t12 +
           0.13403976113848320000e20 * t16 * t8 * x +
           0.20105964170772480000e20 * t16 * t11 * t6 * t12 +
           0.10555631189655552000e20 * t24 * t7 +
           0.1759271864942592000e19 * t24 * t11 * t12 +
           0.62831138033664000e17 * t24 * t16 * x) /
          0.299949115834368000e18);
}



double
AssLegFunction::P_14_1_Deriv(const double x) const
{
  double t31;
  double t7;
  double t3;
  double t8;
  double t17;
  double t4;
  double t28;
  double t9;
  double t12;
  double t21;
  double t1;
  double t25;
  double t13;
  double t2;
  t1  = sqrt(0.210e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * x;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t25 = t17 * t17;
  t28 = t25 * t12;
  t31 = t25 * t17;
  return (t1 / t4 *
            (0.149974557917184000e18 * t9 * t8 +
             0.2924503879385088000e19 * t12 * t9 * t13 +
             0.13403976113848320000e20 * t17 * t9 * x +
             0.20105964170772480000e20 * t21 * t7 * t13 +
             0.10555631189655552000e20 * t25 * t8 +
             0.1759271864942592000e19 * t28 * t13 +
             0.62831138033664000e17 * t31 * x) *
            x / 0.299949115834368000e18 -
          t1 * t4 *
            (0.7798677011693568000e19 * t9 * t7 +
             0.85785447128629248000e20 * t12 * t9 * t2 +
             0.241271570049269760000e21 * t17 * t9 +
             0.225186798712651776000e21 * t21 * t7 * t2 +
             0.70370874597703680000e20 * t25 * t7 +
             0.6031789251231744000e19 * t28 * t2 +
             0.62831138033664000e17 * t31) /
            0.299949115834368000e18);
}



double
AssLegFunction::P_14_2(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t20;
  double t13;
  double t5;
  double t6;
  double t9;
  t1  = sqrt(0.2730e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t13 = t9 * t9;
  t20 = t13 * t13;
  return (t1 * t3 *
          (0.7798677011693568000e19 * t6 * t5 +
           0.85785447128629248000e20 * t9 * t6 * t2 +
           0.241271570049269760000e21 * t13 * t6 +
           0.225186798712651776000e21 * t13 * t9 * t5 * t2 +
           0.70370874597703680000e20 * t20 * t5 +
           0.6031789251231744000e19 * t20 * t9 * t2 +
           0.62831138033664000e17 * t20 * t13) /
          0.15597354023387136000e20);
}



double
AssLegFunction::P_14_2_Deriv(const double x) const
{
  double t4;
  double t5;
  double t19;
  double t12;
  double t32;
  double t1;
  double t15;
  double t22;
  double t8;
  double t3;
  t1  = sqrt(0.2730e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t3 - 0.1e1;
  t12 = t8 * t8;
  t15 = t12 * t8;
  t19 = t12 * t12;
  t22 = t19 * t8;
  t32 = t3 * x;
  return (-t1 * x *
            (0.7798677011693568000e19 * t5 * t4 +
             0.85785447128629248000e20 * t8 * t5 * t3 +
             0.241271570049269760000e21 * t12 * t5 +
             0.225186798712651776000e21 * t15 * t4 * t3 +
             0.70370874597703680000e20 * t19 * t4 +
             0.6031789251231744000e19 * t22 * t3 +
             0.62831138033664000e17 * t19 * t12) /
            0.7798677011693568000e19 -
          t1 * t8 *
            (0.265155018397581312000e21 * t5 * t32 +
             0.1822940751483371520000e22 * t8 * t5 * x +
             0.3281293352670068736000e22 * t12 * t4 * t32 +
             0.1914087789057540096000e22 * t15 * t4 * x +
             0.341801390903132160000e21 * t19 * t32 +
             0.12817552158867456000e20 * t22 * x) /
            0.15597354023387136000e20);
}



double
AssLegFunction::P_14_3(const double x) const
{
  double t7;
  double t2;
  double t9;
  double t4;
  double t24;
  double t8;
  double t16;
  double t1;
  double t3;
  double t12;
  t1  = sqrt(0.15470e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t24 = t16 * t16;
  return (-t1 * t4 * t3 *
          (0.265155018397581312000e21 * t9 * t7 +
           0.1822940751483371520000e22 * t12 * t9 * x +
           0.3281293352670068736000e22 * t16 * t8 * t7 +
           0.1914087789057540096000e22 * t16 * t12 * t8 * x +
           0.341801390903132160000e21 * t24 * t7 +
           0.12817552158867456000e20 * t24 * t12 * x) /
          0.530310036795162624000e21);
}



double
AssLegFunction::P_14_3_Deriv(const double x) const
{
  double t23;
  double t11;
  double t1;
  double t26;
  double t2;
  double t15;
  double t3;
  double t4;
  double t19;
  double t7;
  double t6;
  double t8;
  t1  = sqrt(0.15470e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t23 = t15 * t15;
  t26 = t23 * t11;
  return (t1 * t4 *
            (0.265155018397581312000e21 * t6 * t8 +
             0.1822940751483371520000e22 * t11 * t8 * x +
             0.3281293352670068736000e22 * t15 * t7 * t6 +
             0.1914087789057540096000e22 * t19 * t7 * x +
             0.341801390903132160000e21 * t23 * t6 +
             0.12817552158867456000e20 * t26 * x) *
            x / 0.176770012265054208000e21 -
          t1 * t4 * t3 *
            (0.6562586705340137472000e22 * t8 * t2 +
             0.29531640174030618624000e23 * t11 * t8 +
             0.34453580203035721728000e23 * t15 * t7 * t2 +
             0.12304850072512757760000e23 * t19 * t7 +
             0.1153579694298071040000e22 * t2 * t23 +
             0.12817552158867456000e20 * t26) /
            0.530310036795162624000e21);
}



double
AssLegFunction::P_14_4(const double x) const
{
  double t13;
  double t20;
  double t3;
  double t2;
  double t4;
  double t10;
  double t6;
  double t7;
  double t1;
  t1  = sqrt(0.85085e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t13 = t10 * t10;
  t20 = t13 * t13;
  return (t1 * t4 *
          (0.6562586705340137472000e22 * t7 * t2 +
           0.29531640174030618624000e23 * t10 * t7 +
           0.34453580203035721728000e23 * t13 * t6 * t2 +
           0.12304850072512757760000e23 * t13 * t10 * t6 +
           0.1153579694298071040000e22 * t20 * t2 +
           0.12817552158867456000e20 * t20 * t10) /
          0.17500231214240366592000e23);
}



double
AssLegFunction::P_14_4_Deriv(const double x) const
{
  double t12;
  double t16;
  double t19;
  double t28;
  double t32;
  double t1;
  double t2;
  double t9;
  double t5;
  double t6;
  double t3;
  t1  = sqrt(0.85085e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t12 = t9 * t9;
  t16 = t12 * t9;
  t19 = t12 * t12;
  t28 = t3 * t3;
  t32 = t2 * x;
  return (-t1 * t3 *
            (0.6562586705340137472000e22 * t6 * t2 +
             0.29531640174030618624000e23 * t9 * t6 +
             0.34453580203035721728000e23 * t12 * t5 * t2 +
             0.12304850072512757760000e23 * t16 * t5 +
             0.1153579694298071040000e22 * t19 * t2 +
             0.12817552158867456000e20 * t19 * t9) *
            x / 0.4375057803560091648000e22 +
          t1 * t28 *
            (0.124689147401462611968000e24 * t6 * x +
             0.374067442204387835904000e24 * t9 * t5 * t32 +
             0.280550581653290876928000e24 * t12 * t5 * x +
             0.58448037844435599360000e23 * t16 * t32 +
             0.2435334910184816640000e22 * t19 * x) /
            0.17500231214240366592000e23);
}



double
AssLegFunction::P_14_5(const double x) const
{
  double t1;
  double t12;
  double t5;
  double t2;
  double t13;
  double t8;
  double t17;
  double t9;
  double t3;
  double t4;
  double t24;
  t1  = sqrt(0.646646e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t24 = t17 * t17;
  return (-t1 * t5 * t4 *
          (0.124689147401462611968000e24 * t9 * x +
           0.374067442204387835904000e24 * t12 * t8 * t13 +
           0.280550581653290876928000e24 * t17 * t8 * x +
           0.58448037844435599360000e23 * t17 * t12 * t13 +
           0.2435334910184816640000e22 * t24 * x) /
          0.665008786141133930496000e24);
}



double
AssLegFunction::P_14_5_Deriv(const double x) const
{
  double t20;
  double t7;
  double t8;
  double t1;
  double t23;
  double t11;
  double t2;
  double t3;
  double t12;
  double t4;
  double t16;
  double t30;
  t1  = sqrt(0.646646e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t23 = t16 * t16;
  t30 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.124689147401462611968000e24 * t8 * x +
             0.374067442204387835904000e24 * t11 * t7 * t12 +
             0.280550581653290876928000e24 * t16 * t7 * x +
             0.58448037844435599360000e23 * t20 * t12 +
             0.2435334910184816640000e22 * t23 * x) *
            x / 0.133001757228226786099200e24 -
          t1 * t4 * t30 *
            (0.1870337211021939179520000e25 * t8 +
             0.3740674422043878359040000e25 * t11 * t7 * t2 +
             0.1753441135333067980800000e25 * t16 * t7 +
             0.194826792814785331200000e24 * t20 * t2 +
             0.2435334910184816640000e22 * t23) /
            0.665008786141133930496000e24);
}



double
AssLegFunction::P_14_6(const double x) const
{
  double t1;
  double t20;
  double t3;
  double t2;
  double t7;
  double t8;
  double t4;
  double t10;
  double t14;
  t1  = sqrt(0.3233230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t20 = t14 * t14;
  return (t1 * t4 * t3 *
          (0.1870337211021939179520000e25 * t8 +
           0.3740674422043878359040000e25 * t10 * t7 * t2 +
           0.1753441135333067980800000e25 * t14 * t7 +
           0.194826792814785331200000e24 * t14 * t10 * t2 +
           0.2435334910184816640000e22 * t20) /
          0.19950263584234017914880000e26);
}



double
AssLegFunction::P_14_6_Deriv(const double x) const
{
  double t2;
  double t3;
  double t16;
  double t27;
  double t4;
  double t19;
  double t6;
  double t7;
  double t1;
  double t9;
  double t13;
  t1  = sqrt(0.3233230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t9  = -t3;
  t13 = t9 * t9;
  t16 = t9 * t13;
  t19 = t13 * t13;
  t27 = t2 * x;
  return (-t1 * t4 *
            (0.1870337211021939179520000e25 * t7 +
             0.3740674422043878359040000e25 * t9 * t6 * t2 +
             0.1753441135333067980800000e25 * t13 * t6 +
             0.194826792814785331200000e24 * t16 * t2 +
             0.2435334910184816640000e22 * t19) *
            x / 0.3325043930705669652480000e25 +
          t1 * t4 * t3 *
            (0.22444046532263270154240000e26 * t27 * t6 +
             0.29457811073595542077440000e26 * t9 * t6 * x +
             0.8182725298220983910400000e25 * t13 * t27 +
             0.409136264911049195520000e24 * t16 * x) /
            0.19950263584234017914880000e26);
}



double
AssLegFunction::P_14_7(const double x) const
{
  double t1;
  double t13;
  double t2;
  double t3;
  double t4;
  double t6;
  double t9;
  double t10;
  double t17;
  t1  = sqrt(0.692835e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (-t1 * t6 * t4 * t3 *
          (0.22444046532263270154240000e26 * t10 * t9 +
           0.29457811073595542077440000e26 * t13 * t10 * x +
           0.8182725298220983910400000e25 * t17 * t9 +
           0.409136264911049195520000e24 * t17 * t13 * x) /
          0.119701581505404107489280000e27);
}



double
AssLegFunction::P_14_7_Deriv(const double x) const
{
  double t4;
  double t16;
  double t9;
  double t12;
  double t5;
  double t19;
  double t1;
  double t2;
  double t8;
  double t3;
  t1  = sqrt(0.692835e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * x;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  t19 = t16 * t12;
  return (t1 * t5 * t4 *
            (0.22444046532263270154240000e26 * t9 * t8 +
             0.29457811073595542077440000e26 * t12 * t9 * x +
             0.8182725298220983910400000e25 * t16 * t8 +
             0.409136264911049195520000e24 * t19 * x) *
            x / 0.17100225929343443927040000e26 -
          t1 * t5 * t4 * t3 *
            (0.216023947873033975234560000e27 * t9 * t2 +
             0.180019956560861646028800000e27 * t12 * t9 +
             0.27002993484129246904320000e26 * t16 * t2 +
             0.409136264911049195520000e24 * t19) /
            0.119701581505404107489280000e27);
}



double
AssLegFunction::P_14_8(const double x) const
{
  double t13;
  double t1;
  double t7;
  double t2;
  double t4;
  double t3;
  double t5;
  double t10;
  t1  = sqrt(0.881790e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t5 *
          (0.216023947873033975234560000e27 * t7 * t2 +
           0.180019956560861646028800000e27 * t10 * t7 +
           0.27002993484129246904320000e26 * t13 * t2 +
           0.409136264911049195520000e24 * t13 * t10) /
          0.1675822141075657504849920000e28);
}



double
AssLegFunction::P_14_8_Deriv(const double x) const
{
  double t7;
  double t22;
  double t1;
  double t10;
  double t2;
  double t4;
  double t3;
  double t13;
  t1  = sqrt(0.881790e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  t22 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.216023947873033975234560000e27 * t7 * t2 +
             0.180019956560861646028800000e27 * t7 * t10 +
             0.27002993484129246904320000e26 * t13 * t2 +
             0.409136264911049195520000e24 * t13 * t10) *
            x / 0.209477767634457188106240000e27 +
          t1 * t22 *
            (0.1656183600359927143464960000e28 * t7 * x +
             0.828091800179963571732480000e27 * t10 * t2 * x +
             0.56460804557724788981760000e26 * t13 * x) /
            0.1675822141075657504849920000e28);
}



double
AssLegFunction::P_14_9(const double x) const
{
  double t1;
  double t12;
  double t9;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  t1  = sqrt(0.3380195e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (-t1 * t6 * t5 *
          (0.1656183600359927143464960000e28 * t9 * x +
           0.828091800179963571732480000e27 * t12 * t2 * x +
           0.56460804557724788981760000e26 * t16 * x) /
          0.38543909244740122611548160000e29);
}



double
AssLegFunction::P_14_9_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t12;
  double t16;
  double t6;
  double t1;
  double t9;
  double t23;
  t1  = sqrt(0.3380195e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  t23 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.1656183600359927143464960000e28 * t9 * x +
             0.828091800179963571732480000e27 * t12 * t2 * x +
             0.56460804557724788981760000e26 * t16 * x) *
            x / 0.4282656582748902512394240000e28 -
          t1 * t6 * t23 *
            (0.9937101602159562860789760000e28 * t9 +
             0.2710118618770789871124480000e28 * t12 * t2 +
             0.56460804557724788981760000e26 * t16) /
            0.38543909244740122611548160000e29);
}



double
AssLegFunction::P_14_10(const double x) const
{
  double t3;
  double t8;
  double t10;
  double t1;
  double t13;
  double t4;
  double t5;
  double t2;
  t1  = sqrt(0.4056234e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t5 * t3 *
          (0.9937101602159562860789760000e28 * t8 +
           0.2710118618770789871124480000e28 * t10 * t2 +
           0.56460804557724788981760000e26 * t13) /
          0.462526910936881471338577920000e30);
}



double
AssLegFunction::P_14_10_Deriv(const double x) const
{
  double t7;
  double t1;
  double t9;
  double t2;
  double t12;
  double t4;
  double t3;
  double t5;
  t1  = sqrt(0.4056234e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t9  = -t3;
  t12 = t9 * t9;
  return (-t1 * t5 *
            (0.9937101602159562860789760000e28 * t7 +
             0.2710118618770789871124480000e28 * t9 * t2 +
             0.56460804557724788981760000e26 * t12) *
            x / 0.46252691093688147133857792000e29 +
          t1 * t5 * t3 *
            (0.45168643646179831185408000000e29 * t2 * x +
             0.5646080455772478898176000000e28 * t9 * x) /
            0.462526910936881471338577920000e30);
}



double
AssLegFunction::P_14_11(const double x) const
{
  double t5;
  double t7;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.4056234e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-t1 * t7 * t5 * t3 *
          (0.45168643646179831185408000000e29 * t2 * x -
           0.5646080455772478898176000000e28 * t3 * x) /
          0.4625269109368814713385779200000e31);
}



double
AssLegFunction::P_14_11_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  t1 = sqrt(0.4056234e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = sqrt(t3);
  return (t1 * t6 * t5 *
            (0.45168643646179831185408000000e29 * t2 * x -
             0.5646080455772478898176000000e28 * t3 * x) *
            x / 0.420479009942619519398707200000e30 -
          t1 * t6 * t5 * t3 *
            (0.152444172305856930250752000000e30 * t2 -
             0.5646080455772478898176000000e28) /
            0.4625269109368814713385779200000e31);
}



double
AssLegFunction::P_14_12(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  t1 = sqrt(0.52003e5);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (t1 * t5 * t4 *
          (0.152444172305856930250752000000e30 * t2 -
           0.5646080455772478898176000000e28) /
          0.4625269109368814713385779200000e31);
}



double
AssLegFunction::P_14_12_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.52003e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-t1 * t5 * t3 *
            (0.152444172305856930250752000000e30 * t2 -
             0.5646080455772478898176000000e28) *
            x / 0.385439092447401226115481600000e30 +
          0.135e3 / 0.2048e4 * t1 * t5 * t4 * x);
}



double
AssLegFunction::P_14_13(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t1;
  t1 = sqrt(0.312018e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-0.15e2 / 0.4096e4 * t1 * t7 * t5 * t4 * x);
}



double
AssLegFunction::P_14_13_Deriv(const double x) const
{
  double t7;
  double t1;
  double t2;
  double t3;
  double t5;
  double t4;
  t1 = sqrt(0.312018e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (0.195e3 / 0.4096e4 * t1 * t7 * t5 * t3 * t2 -
          0.15e2 / 0.4096e4 * t1 * t7 * t5 * t4);
}



double
AssLegFunction::P_14_14(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  t1 = sqrt(0.44574e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  return (0.15e2 / 0.8192e4 * t1 * t6 * t4 * t3);
}



double
AssLegFunction::P_14_14_Deriv(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  t1 = sqrt(0.44574e5);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  return (-0.105e3 / 0.4096e4 * t1 * t5 * t4 * x);
}

double
AssLegFunction::P_15_0(const double x) const
{
  double t7;
  double t8;
  double t3;
  double t4;
  double t5;
  double t16;
  double t20;
  double t12;
  double t1;
  double t2;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t2;
  t5  = t3 * t3;
  t7  = t1 - 0.1e1;
  t8  = t3 * x;
  t12 = t7 * t7;
  t16 = t12 * t7;
  t20 = t12 * t12;
  return (t5 * t4 + 0.105e3 / 0.2e1 * t7 * t5 * t8 +
          0.4095e4 / 0.8e1 * t12 * t5 * t2 + 0.25025e5 / 0.16e2 * t16 * t5 * x +
          0.225225e6 / 0.128e3 * t4 * t20 +
          0.189189e6 / 0.256e3 * t20 * t7 * t8 +
          0.105105e6 / 0.1024e4 * t20 * t12 * t2 +
          0.6435e4 / 0.2048e4 * t20 * t16 * x);
}

double
AssLegFunction::P_15_0_Deriv(const double x) const
{
  double t3;
  double t11;
  double t4;
  double t15;
  double t7;
  double t1;
  double t18;
  double t2;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t1;
  t4  = t2 * t2;
  t7  = t1 - 0.1e1;
  t11 = t7 * t7;
  t15 = t11 * t7;
  t18 = t11 * t11;
  return (0.120e3 * t4 * t3 + 0.2730e4 * t7 * t4 * t2 +
          0.15015e5 * t11 * t4 * t1 + 0.225225e6 / 0.8e1 * t15 * t4 +
          0.315315e6 / 0.16e2 * t18 * t3 + 0.315315e6 / 0.64e2 * t18 * t7 * t2 +
          0.45045e5 / 0.128e3 * t18 * t11 * t1 +
          0.6435e4 / 0.2048e4 * t15 * t18);
}



double
AssLegFunction::P_15_1(const double x) const
{
  double t11;
  double t19;
  double t4;
  double t3;
  double t6;
  double t15;
  double t22;
  double t7;
  double t8;
  double t1;
  double t2;
  t1  = sqrt(0.15e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t2;
  t8  = t6 * t6;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t22 = t15 * t15;
  return (-t1 * t4 *
          (0.5141984842874880000e19 * t8 * t7 +
           0.116980155175403520000e21 * t11 * t8 * t6 +
           0.643390853464719360000e21 * t15 * t8 * t2 +
           0.1206357850246348800000e22 * t19 * t8 +
           0.844450495172444160000e21 * t22 * t7 +
           0.211112623793111040000e21 * t22 * t11 * t6 +
           0.15079473128079360000e20 * t22 * t15 * t2 +
           0.134638152929280000e18 * t22 * t19) /
          0.2570992421437440000e19);
}



double
AssLegFunction::P_15_1_Deriv(const double x) const
{
  double t7;
  double t8;
  double t9;
  double t3;
  double t4;
  double t38;
  double t12;
  double t23;
  double t41;
  double t1;
  double t16;
  double t26;
  double t20;
  double t29;
  double t2;
  t1  = sqrt(0.15e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t2;
  t9  = t7 * t7;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t23 = t16 * t16;
  t26 = t12 * t23;
  t29 = t23 * t16;
  t38 = t7 * x;
  t41 = t2 * x;
  return (t1 / t4 *
            (0.5141984842874880000e19 * t8 * t9 +
             0.116980155175403520000e21 * t12 * t9 * t7 +
             0.643390853464719360000e21 * t16 * t9 * t2 +
             0.1206357850246348800000e22 * t20 * t9 +
             0.844450495172444160000e21 * t23 * t8 +
             0.211112623793111040000e21 * t26 * t7 +
             0.15079473128079360000e20 * t29 * t2 +
             0.134638152929280000e18 * t23 * t20) *
            x / 0.2570992421437440000e19 -
          t1 * t4 *
            (0.305948098151055360000e21 * t9 * t38 +
             0.3977325275963719680000e22 * t12 * t9 * t41 +
             0.13672055636125286400000e23 * t16 * t9 * x +
             0.16406466763350343680000e23 * t20 * t7 * t41 +
             0.7177829208965775360000e22 * t23 * t38 +
             0.1025404172709396480000e22 * t26 * t41 +
             0.32043880397168640000e20 * t29 * x) /
            0.2570992421437440000e19);
}



double
AssLegFunction::P_15_2(const double x) const
{
  double t6;
  double t7;
  double t2;
  double t10;
  double t11;
  double t3;
  double t5;
  double t1;
  double t15;
  double t23;
  t1  = sqrt(0.3570e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * x;
  t7  = t5 * t5;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t23 = t15 * t15;
  return (t1 * t3 *
          (0.305948098151055360000e21 * t7 * t6 +
           0.3977325275963719680000e22 * t10 * t7 * t11 +
           0.13672055636125286400000e23 * t15 * t7 * x +
           0.16406466763350343680000e23 * t15 * t10 * t5 * t11 +
           0.7177829208965775360000e22 * t23 * t6 +
           0.1025404172709396480000e22 * t23 * t10 * t11 +
           0.32043880397168640000e20 * t23 * t15 * x) /
          0.611896196302110720000e21);
}



double
AssLegFunction::P_15_2_Deriv(const double x) const
{
  double t1;
  double t28;
  double t22;
  double t3;
  double t4;
  double t5;
  double t18;
  double t6;
  double t14;
  double t9;
  double t25;
  double t10;
  t1  = sqrt(0.3570e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * x;
  t6  = t4 * t4;
  t9  = t3 - 0.1e1;
  t10 = t3 * x;
  t14 = t9 * t9;
  t18 = t14 * t9;
  t22 = t14 * t14;
  t25 = t22 * t9;
  t28 = t22 * t14;
  return (-t1 * x *
            (0.305948098151055360000e21 * t6 * t5 +
             0.3977325275963719680000e22 * t9 * t6 * t10 +
             0.13672055636125286400000e23 * t14 * t6 * x +
             0.16406466763350343680000e23 * t18 * t4 * t10 +
             0.7177829208965775360000e22 * t22 * t5 +
             0.1025404172709396480000e22 * t25 * t10 +
             0.32043880397168640000e20 * t28 * x) /
            0.305948098151055360000e21 -
          t1 * t9 *
            (0.11931975827891159040000e23 * t6 * t4 +
             0.98438800580102062080000e23 * t9 * t6 * t3 +
             0.221487301305229639680000e24 * t14 * t6 +
             0.172267901015178608640000e24 * t18 * t4 * t3 +
             0.46143187771922841600000e23 * t22 * t4 +
             0.3460739082894213120000e22 * t25 * t3 +
             0.32043880397168640000e20 * t28) /
            0.611896196302110720000e21);
}



double
AssLegFunction::P_15_3(const double x) const
{
  double t11;
  double t4;
  double t2;
  double t15;
  double t22;
  double t8;
  double t1;
  double t3;
  double t7;
  t1  = sqrt(0.23205e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t22 = t15 * t15;
  return (-t1 * t4 * t3 *
          (0.11931975827891159040000e23 * t8 * t7 +
           0.98438800580102062080000e23 * t11 * t8 * t2 +
           0.221487301305229639680000e24 * t15 * t8 +
           0.172267901015178608640000e24 * t15 * t11 * t7 * t2 +
           0.46143187771922841600000e23 * t22 * t7 +
           0.3460739082894213120000e22 * t22 * t11 * t2 +
           0.32043880397168640000e20 * t22 * t15) /
          0.23863951655782318080000e23);
}



double
AssLegFunction::P_15_3_Deriv(const double x) const
{
  double t6;
  double t7;
  double t1;
  double t2;
  double t3;
  double t4;
  double t14;
  double t17;
  double t24;
  double t35;
  double t10;
  double t21;
  t1  = sqrt(0.23205e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t21 = t14 * t14;
  t24 = t10 * t21;
  t35 = t2 * x;
  return (t1 * t4 *
            (0.11931975827891159040000e23 * t7 * t6 +
             0.98438800580102062080000e23 * t10 * t7 * t2 +
             0.221487301305229639680000e24 * t14 * t7 +
             0.172267901015178608640000e24 * t17 * t6 * t2 +
             0.46143187771922841600000e23 * t21 * t6 +
             0.3460739082894213120000e22 * t24 * t2 +
             0.32043880397168640000e20 * t21 * t14) *
            x / 0.7954650551927439360000e22 -
          t1 * t4 * t3 *
            (0.340061311094898032640000e24 * t7 * t35 +
             0.1870337211021939179520000e25 * t10 * t7 * x +
             0.2805505816532908769280000e25 * t14 * t6 * t35 +
             0.1402752908266454384640000e25 * t17 * t6 * x +
             0.219180141916633497600000e24 * t21 * t35 +
             0.7306004730554449920000e22 * t24 * x) /
            0.23863951655782318080000e23);
}



double
AssLegFunction::P_15_4(const double x) const
{
  double t3;
  double t6;
  double t8;
  double t4;
  double t7;
  double t15;
  double t2;
  double t23;
  double t11;
  double t1;
  t1  = sqrt(0.146965e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t23 = t15 * t15;
  return (t1 * t4 *
          (0.340061311094898032640000e24 * t6 * t8 +
           0.1870337211021939179520000e25 * t11 * t8 * x +
           0.2805505816532908769280000e25 * t15 * t7 * t6 +
           0.1402752908266454384640000e25 * t15 * t11 * t7 * x +
           0.219180141916633497600000e24 * t23 * t6 +
           0.7306004730554449920000e22 * t23 * t11 * x) /
          0.906830162919728087040000e24);
}



double
AssLegFunction::P_15_4_Deriv(const double x) const
{
  double t25;
  double t18;
  double t1;
  double t32;
  double t2;
  double t3;
  double t14;
  double t5;
  double t10;
  double t6;
  double t7;
  double t22;
  t1  = sqrt(0.146965e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t18 = t14 * t10;
  t22 = t14 * t14;
  t25 = t22 * t10;
  t32 = t3 * t3;
  return (-t1 * t3 *
            (0.340061311094898032640000e24 * t7 * t5 +
             0.1870337211021939179520000e25 * t10 * t7 * x +
             0.2805505816532908769280000e25 * t14 * t6 * t5 +
             0.1402752908266454384640000e25 * t18 * t6 * x +
             0.219180141916633497600000e24 * t22 * t5 +
             0.7306004730554449920000e22 * t25 * x) *
            x / 0.226707540729932021760000e24 +
          t1 * t32 *
            (0.7481348844087756718080000e25 * t7 * t2 +
             0.28055058165329087692800000e26 * t10 * t7 +
             0.28055058165329087692800000e26 * t14 * t6 * t2 +
             0.8767205676665339904000000e25 * t18 * t6 +
             0.730600473055444992000000e24 * t22 * t2 +
             0.7306004730554449920000e22 * t25) /
            0.906830162919728087040000e24);
}



double
AssLegFunction::P_15_5(const double x) const
{
  double t15;
  double t12;
  double t8;
  double t9;
  double t5;
  double t1;
  double t2;
  double t4;
  double t22;
  double t3;
  t1  = sqrt(0.323323e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t15 = t12 * t12;
  t22 = t15 * t15;
  return (-t1 * t5 * t4 *
          (0.7481348844087756718080000e25 * t9 * t2 +
           0.28055058165329087692800000e26 * t12 * t9 +
           0.28055058165329087692800000e26 * t15 * t8 * t2 +
           0.8767205676665339904000000e25 * t15 * t12 * t8 +
           0.730600473055444992000000e24 * t22 * t2 +
           0.7306004730554449920000e22 * t22 * t12) /
          0.19950263584234017914880000e26);
}



double
AssLegFunction::P_15_5_Deriv(const double x) const
{
  double t7;
  double t8;
  double t11;
  double t30;
  double t14;
  double t18;
  double t21;
  double t1;
  double t35;
  double t3;
  double t2;
  double t4;
  t1  = sqrt(0.323323e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t18 = t14 * t11;
  t21 = t14 * t14;
  t30 = t3 * t3;
  t35 = t2 * x;
  return (t1 * t4 * t3 *
            (0.7481348844087756718080000e25 * t8 * t2 +
             0.28055058165329087692800000e26 * t8 * t11 +
             0.28055058165329087692800000e26 * t14 * t7 * t2 +
             0.8767205676665339904000000e25 * t18 * t7 +
             0.730600473055444992000000e24 * t21 * t2 +
             0.7306004730554449920000e22 * t21 * t11) *
            x / 0.3990052716846803582976000e25 -
          t1 * t4 * t30 *
            (0.130923604771535742566400000e27 * t8 * x +
             0.336660697983949052313600000e27 * t11 * t7 * t35 +
             0.220933583051966565580800000e27 * t14 * t7 * x +
             0.40913626491104919552000000e26 * t18 * t35 +
             0.1534260993416434483200000e25 * t21 * x) /
            0.19950263584234017914880000e26);
}



double
AssLegFunction::P_15_6(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t8;
  double t16;
  double t11;
  double t23;
  double t12;
  double t7;
  double t3;
  t1  = sqrt(0.1385670e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t23 = t16 * t16;
  return (t1 * t4 * t3 *
          (0.130923604771535742566400000e27 * t8 * x +
           0.336660697983949052313600000e27 * t11 * t7 * t12 +
           0.220933583051966565580800000e27 * t16 * t7 * x +
           0.40913626491104919552000000e26 * t16 * t11 * t12 +
           0.1534260993416434483200000e25 * t23 * x) /
          0.598507907527020537446400000e27);
}



double
AssLegFunction::P_15_6_Deriv(const double x) const
{
  double t6;
  double t7;
  double t1;
  double t15;
  double t4;
  double t2;
  double t3;
  double t19;
  double t22;
  double t10;
  double t11;
  t1  = sqrt(0.1385670e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t22 = t15 * t15;
  return (-t1 * t4 *
            (0.130923604771535742566400000e27 * t7 * x +
             0.336660697983949052313600000e27 * t10 * t6 * t11 +
             0.220933583051966565580800000e27 * t15 * t6 * x +
             0.40913626491104919552000000e26 * t19 * t11 +
             0.1534260993416434483200000e25 * t22 * x) *
            x / 0.99751317921170089574400000e26 +
          t1 * t4 * t3 *
            (0.1851633838911719787724800000e28 * t7 +
             0.3240359218095509628518400000e28 * t10 * t6 * t2 +
             0.1350149674206462345216000000e28 * t15 * t6 +
             0.135014967420646234521600000e27 * t19 * t2 +
             0.1534260993416434483200000e25 * t22) /
            0.598507907527020537446400000e27);
}



double
AssLegFunction::P_15_7(const double x) const
{
  double t22;
  double t3;
  double t4;
  double t16;
  double t12;
  double t6;
  double t9;
  double t2;
  double t1;
  double t10;
  t1  = sqrt(0.62985e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t12 = -t3;
  t16 = t12 * t12;
  t22 = t16 * t16;
  return (-t1 * t6 * t4 * t3 *
          (0.1851633838911719787724800000e28 * t10 +
           0.3240359218095509628518400000e28 * t12 * t9 * t2 +
           0.1350149674206462345216000000e28 * t16 * t9 +
           0.135014967420646234521600000e27 * t16 * t12 * t2 +
           0.1534260993416434483200000e25 * t22) /
          0.1795523722581061612339200000e28);
}



double
AssLegFunction::P_15_7_Deriv(const double x) const
{
  double t9;
  double t18;
  double t1;
  double t11;
  double t21;
  double t2;
  double t3;
  double t4;
  double t5;
  double t30;
  double t8;
  double t15;
  t1  = sqrt(0.62985e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t21 = t15 * t15;
  t30 = t2 * x;
  return (t1 * t5 * t4 *
            (0.1851633838911719787724800000e28 * t9 +
             0.3240359218095509628518400000e28 * t11 * t8 * t2 +
             0.1350149674206462345216000000e28 * t8 * t15 +
             0.135014967420646234521600000e27 * t18 * t2 +
             0.1534260993416434483200000e25 * t21) *
            x / 0.256503388940151658905600000e27 -
          t1 * t5 * t4 * t3 *
            (0.21293789147484777558835200000e29 * t30 * t8 +
             0.24842754005398907151974400000e29 * t11 * t8 * x +
             0.6210688501349726787993600000e28 * t15 * t30 +
             0.282304022788623944908800000e27 * t18 * x) /
            0.1795523722581061612339200000e28);
}



double
AssLegFunction::P_15_8(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t11;
  double t3;
  double t5;
  double t15;
  double t8;
  double t7;
  t1  = sqrt(0.2897310e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (t1 * t5 *
          (0.21293789147484777558835200000e29 * t8 * t7 +
           0.24842754005398907151974400000e29 * t11 * t8 * x +
           0.6210688501349726787993600000e28 * t15 * t7 +
           0.282304022788623944908800000e27 * t15 * t11 * x) /
          0.165188182477457668335206400000e30);
}



double
AssLegFunction::P_15_8_Deriv(const double x) const
{
  double t18;
  double t2;
  double t4;
  double t3;
  double t11;
  double t15;
  double t7;
  double t1;
  double t25;
  double t8;
  t1  = sqrt(0.2897310e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t25 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.21293789147484777558835200000e29 * t8 * t7 +
             0.24842754005398907151974400000e29 * t11 * t8 * x +
             0.6210688501349726787993600000e28 * t7 * t15 +
             0.282304022788623944908800000e27 * t18 * x) *
            x / 0.20648522809682208541900800000e29 +
          t1 * t25 *
            (0.198742032043191257215795200000e30 * t8 * t2 +
             0.149056524032393442911846400000e30 * t11 * t8 +
             0.20325889640780924033433600000e29 * t15 * t2 +
             0.282304022788623944908800000e27 * t18) /
            0.165188182477457668335206400000e30);
}



double
AssLegFunction::P_15_9(const double x) const
{
  double t12;
  double t6;
  double t9;
  double t1;
  double t15;
  double t2;
  double t3;
  double t5;
  double t4;
  t1  = sqrt(0.3380195e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (-t1 * t6 * t5 *
          (0.198742032043191257215795200000e30 * t9 * t2 +
           0.149056524032393442911846400000e30 * t12 * t9 +
           0.20325889640780924033433600000e29 * t15 * t2 +
           0.282304022788623944908800000e27 * t15 * t12) /
          0.2312634554684407356692889600000e31);
}



double
AssLegFunction::P_15_9_Deriv(const double x) const
{
  double t3;
  double t4;
  double t12;
  double t24;
  double t6;
  double t15;
  double t2;
  double t1;
  double t9;
  t1  = sqrt(0.3380195e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  t24 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.198742032043191257215795200000e30 * t9 * t2 +
             0.149056524032393442911846400000e30 * t12 * t9 +
             0.20325889640780924033433600000e29 * t2 * t15 +
             0.282304022788623944908800000e27 * t12 * t15) *
            x / 0.256959394964934150743654400000e30 -
          t1 * t6 * t24 *
            (0.1490565240323934429118464000000e31 * t9 * x +
             0.677529654692697467781120000000e30 * t12 * t2 * x +
             0.42345603418293591736320000000e29 * t15 * x) /
            0.2312634554684407356692889600000e31);
}



double
AssLegFunction::P_15_10(const double x) const
{
  double t11;
  double t8;
  double t15;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1  = sqrt(0.20281170e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (t1 * t5 * t3 *
          (0.1490565240323934429118464000000e31 * t8 * x +
           0.677529654692697467781120000000e30 * t11 * t2 * x +
           0.42345603418293591736320000000e29 * t15 * x) /
          0.69379036640532220700786688000000e32);
}



double
AssLegFunction::P_15_10_Deriv(const double x) const
{
  double t4;
  double t5;
  double t14;
  double t2;
  double t3;
  double t1;
  double t7;
  double t10;
  t1  = sqrt(0.20281170e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t10 = -t3;
  t14 = t10 * t10;
  return (-t1 * t5 *
            (0.1490565240323934429118464000000e31 * t7 * x +
             0.677529654692697467781120000000e30 * t10 * t2 * x +
             0.42345603418293591736320000000e29 * t14 * x) *
            x / 0.6937903664053222070078668800000e31 +
          t1 * t5 * t3 *
            (0.8807885511005067081154560000000e31 * t7 +
             0.2201971377751266770288640000000e31 * t10 * t2 +
             0.42345603418293591736320000000e29 * t14) /
            0.69379036640532220700786688000000e32);
}



double
AssLegFunction::P_15_11(const double x) const
{
  double t10;
  double t15;
  double t1;
  double t2;
  double t12;
  double t3;
  double t4;
  double t5;
  double t7;
  t1  = sqrt(0.156009e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (-t1 * t7 * t5 * t3 *
          (0.8807885511005067081154560000000e31 * t10 +
           0.2201971377751266770288640000000e31 * t12 * t2 +
           0.42345603418293591736320000000e29 * t15) /
          0.69379036640532220700786688000000e32);
}



double
AssLegFunction::P_15_11_Deriv(const double x) const
{
  double t1;
  double t2;
  double t9;
  double t3;
  double t6;
  double t4;
  double t5;
  double t11;
  double t14;
  t1  = sqrt(0.156009e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (t1 * t6 * t5 *
            (0.8807885511005067081154560000000e31 * t9 +
             0.2201971377751266770288640000000e31 * t11 * t2 +
             0.42345603418293591736320000000e29 * t14) *
            x / 0.6307185149139292790980608000000e31 -
          t1 * t6 * t5 * t3 *
            (0.39635484799522801865195520000000e32 * t2 * x +
             0.4573325169175707907522560000000e31 * t11 * x) /
            0.69379036640532220700786688000000e32);
}



double
AssLegFunction::P_15_12(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  double t5;
  t1 = sqrt(0.52003e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (t1 * t5 * t4 *
          (0.39635484799522801865195520000000e32 * t2 * x -
           0.4573325169175707907522560000000e31 * t3 * x) /
          0.416274219843193324204720128000000e33);
}



double
AssLegFunction::P_15_12_Deriv(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t5;
  double t4;
  t1 = sqrt(0.52003e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-t1 * t5 * t3 *
            (0.39635484799522801865195520000000e32 * t2 * x -
             0.4573325169175707907522560000000e31 * t3 * x) *
            x / 0.34689518320266110350393344000000e32 +
          t1 * t5 * t4 *
            (0.132626429906095529318154240000000e33 * t2 -
             0.4573325169175707907522560000000e31) /
            0.416274219843193324204720128000000e33);
}



double
AssLegFunction::P_15_13(const double x) const
{
  double t5;
  double t7;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.22287e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-t1 * t7 * t5 * t4 *
          (0.132626429906095529318154240000000e33 * t2 -
           0.4573325169175707907522560000000e31) /
          0.2497645319059159945228320768000000e34);
}



double
AssLegFunction::P_15_13_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t1;
  t1 = sqrt(0.22287e5);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (t1 * t7 * t5 * t3 *
            (0.132626429906095529318154240000000e33 * t2 -
             0.4573325169175707907522560000000e31) *
            x / 0.192126563004550765017563136000000e33 -
          0.435e3 / 0.4096e4 * t1 * t7 * t5 * t4 * x);
}



double
AssLegFunction::P_15_14(const double x) const
{
  double t4;
  double t6;
  double t1;
  double t2;
  double t3;
  t1 = sqrt(0.1292646e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  return (0.15e2 / 0.8192e4 * t1 * t6 * t4 * t3 * x);
}



double
AssLegFunction::P_15_14_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.1292646e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-0.105e3 / 0.4096e4 * t1 * t5 * t4 * t2 +
          0.15e2 / 0.8192e4 * t1 * t5 * t4 * t3);
}



double
AssLegFunction::P_15_15(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  t1 = sqrt(0.1077205e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  t8 = sqrt(t3);
  return (-0.3e1 / 0.8192e4 * t1 * t8 * t6 * t4 * t3);
}



double
AssLegFunction::P_15_15_Deriv(const double x) const
{
  double t7;
  double t2;
  double t3;
  double t4;
  double t5;
  double t1;
  t1 = sqrt(0.1077205e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (0.45e2 / 0.8192e4 * t1 * t7 * t5 * t4 * x);
}

double
AssLegFunction::P_16_0(const double x) const
{
  double t1;
  double t6;
  double t2;
  double t14;
  double t10;
  double t3;
  double t30;
  double t4;
  double t18;
  double t5;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t5  = t1 - 0.1e1;
  t6  = t2 * t1;
  t10 = t5 * t5;
  t14 = t10 * t5;
  t18 = t10 * t10;
  t30 = t18 * t18;
  return (t4 + 0.60e2 * t5 * t3 * t6 + 0.1365e4 / 0.2e1 * t10 * t3 * t2 +
          0.5005e4 / 0.2e1 * t14 * t3 * t1 + 0.225225e6 / 0.64e2 * t18 * t3 +
          0.63063e5 / 0.32e2 * t18 * t5 * t6 +
          0.105105e6 / 0.256e3 * t18 * t10 * t2 +
          0.6435e4 / 0.256e3 * t18 * t14 * t1 + 0.6435e4 / 0.32768e5 * t30);
}

double
AssLegFunction::P_16_0_Deriv(const double x) const
{
  double t17;
  double t2;
  double t8;
  double t9;
  double t3;
  double t4;
  double t5;
  double t21;
  double t13;
  double t1;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t2;
  t5  = t3 * t3;
  t8  = t1 - 0.1e1;
  t9  = t3 * x;
  t13 = t8 * t8;
  t17 = t13 * t8;
  t21 = t13 * t13;
  return (0.136e3 * t5 * t4 + 0.3570e4 * t8 * t5 * t9 +
          0.23205e5 * t13 * t5 * t2 + 0.425425e6 / 0.8e1 * t17 * t5 * x +
          0.765765e6 / 0.16e2 * t21 * t4 +
          0.1072071e7 / 0.64e2 * t21 * t8 * t9 +
          0.255255e6 / 0.128e3 * t21 * t13 * t2 +
          0.109395e6 / 0.2048e4 * t21 * t17 * x);
}



double
AssLegFunction::P_16_1(const double x) const
{
  double t17;
  double t6;
  double t21;
  double t13;
  double t7;
  double t12;
  double t8;
  double t1;
  double t9;
  double t2;
  double t3;
  double t4;
  double t25;
  t1  = sqrt(0.17e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t6;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t7 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t25 = t17 * t17;
  return (-t1 * t4 *
          (0.186482650301595648000e21 * t9 * t8 +
           0.4895169570416885760000e22 * t12 * t9 * t13 +
           0.31818602207709757440000e23 * t17 * t9 * t6 +
           0.72917630059334860800000e23 * t21 * t9 * x +
           0.65625867053401374720000e23 * t25 * t8 +
           0.22969053468690481152000e23 * t25 * t12 * t13 +
           0.2734411127225057280000e22 * t25 * t17 * t6 +
           0.73243155193528320000e20 * t25 * t21 * x) /
          0.93241325150797824000e20);
}



double
AssLegFunction::P_16_1_Deriv(const double x) const
{
  double t14;
  double t29;
  double t2;
  double t32;
  double t8;
  double t7;
  double t3;
  double t4;
  double t42;
  double t26;
  double t9;
  double t10;
  double t35;
  double t18;
  double t1;
  double t13;
  double t22;
  t1  = sqrt(0.17e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t7;
  t10 = t8 * t8;
  t13 = -t3;
  t14 = t8 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t26 = t18 * t18;
  t29 = t26 * t13;
  t32 = t26 * t18;
  t35 = t26 * t22;
  t42 = t8 * t2;
  return (t1 / t4 *
            (0.186482650301595648000e21 * t9 * t10 +
             0.4895169570416885760000e22 * t13 * t10 * t14 +
             0.31818602207709757440000e23 * t18 * t10 * t7 +
             0.72917630059334860800000e23 * t22 * t10 * x +
             0.65625867053401374720000e23 * t26 * t9 +
             0.22969053468690481152000e23 * t29 * t14 +
             0.2734411127225057280000e22 * t32 * t7 +
             0.73243155193528320000e20 * t35 * x) *
            x / 0.93241325150797824000e20 -
          t1 * t4 *
            (0.12587578895357706240000e23 * t10 * t42 +
             0.190911613246258544640000e24 * t13 * t10 * t8 +
             0.787510404640816496640000e24 * t18 * t10 * t2 +
             0.1181265606961224744960000e25 * t22 * t10 +
             0.689071604060714434560000e24 * t26 * t42 +
             0.147658200870153093120000e24 * t8 * t29 +
             0.9228637554384568320000e22 * t32 * t2 +
             0.73243155193528320000e20 * t35) /
            0.93241325150797824000e20);
}



double
AssLegFunction::P_16_2(const double x) const
{
  double t2;
  double t18;
  double t21;
  double t1;
  double t10;
  double t5;
  double t6;
  double t3;
  double t7;
  double t14;
  t1  = sqrt(0.510e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t2;
  t7  = t5 * t5;
  t10 = -t3;
  t14 = t10 * t10;
  t18 = t14 * t10;
  t21 = t14 * t14;
  return (t1 * t3 *
          (0.12587578895357706240000e23 * t7 * t6 +
           0.190911613246258544640000e24 * t10 * t7 * t5 +
           0.787510404640816496640000e24 * t14 * t7 * t2 +
           0.1181265606961224744960000e25 * t18 * t7 +
           0.689071604060714434560000e24 * t21 * t6 +
           0.147658200870153093120000e24 * t21 * t10 * t5 +
           0.9228637554384568320000e22 * t21 * t14 * t2 +
           0.73243155193528320000e20 * t18 * t21) /
          0.8391719263571804160000e22);
}



double
AssLegFunction::P_16_2_Deriv(const double x) const
{
  double t6;
  double t9;
  double t13;
  double t39;
  double t4;
  double t26;
  double t17;
  double t5;
  double t20;
  double t1;
  double t36;
  double t23;
  double t3;
  t1  = sqrt(0.510e3);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t3;
  t6  = t4 * t4;
  t9  = t3 - 0.1e1;
  t13 = t9 * t9;
  t17 = t9 * t13;
  t20 = t13 * t13;
  t23 = t20 * t9;
  t26 = t20 * t13;
  t36 = t4 * x;
  t39 = t3 * x;
  return (-t1 * x *
            (0.12587578895357706240000e23 * t6 * t5 +
             0.190911613246258544640000e24 * t9 * t6 * t4 +
             0.787510404640816496640000e24 * t13 * t6 * t3 +
             0.1181265606961224744960000e25 * t17 * t6 +
             0.689071604060714434560000e24 * t20 * t5 +
             0.147658200870153093120000e24 * t23 * t4 +
             0.9228637554384568320000e22 * t26 * t3 +
             0.73243155193528320000e20 * t20 * t17) /
            0.4195859631785902080000e22 -
          t1 * t9 *
            (0.558049331027524976640000e24 * t6 * t36 +
             0.5440980977518368522240000e25 * t9 * t6 * t39 +
             0.14962697688175513436160000e26 * t13 * t6 * x +
             0.14962697688175513436160000e26 * t17 * t4 * t39 +
             0.5611011633065817538560000e25 * t20 * t36 +
             0.701376454133227192320000e24 * t23 * t39 +
             0.19482679281478533120000e23 * t26 * x) /
            0.8391719263571804160000e22);
}



double
AssLegFunction::P_16_3(const double x) const
{
  double t3;
  double t2;
  double t4;
  double t1;
  double t7;
  double t17;
  double t8;
  double t12;
  double t25;
  double t13;
  double t9;
  t1  = sqrt(0.33915e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * x;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t25 = t17 * t17;
  return (-t1 * t4 * t3 *
          (0.558049331027524976640000e24 * t9 * t8 +
           0.5440980977518368522240000e25 * t12 * t9 * t13 +
           0.14962697688175513436160000e26 * t17 * t9 * x +
           0.14962697688175513436160000e26 * t17 * t12 * t7 * t13 +
           0.5611011633065817538560000e25 * t25 * t8 +
           0.701376454133227192320000e24 * t25 * t12 * t13 +
           0.19482679281478533120000e23 * t25 * t17 * x) /
          0.1116098662055049953280000e25);
}



double
AssLegFunction::P_16_3_Deriv(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t12;
  double t16;
  double t30;
  double t1;
  double t7;
  double t6;
  double t20;
  double t8;
  double t27;
  double t24;
  double t11;
  t1  = sqrt(0.33915e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * x;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t24 = t16 * t16;
  t27 = t24 * t11;
  t30 = t24 * t16;
  return (t1 * t4 *
            (0.558049331027524976640000e24 * t8 * t7 +
             0.5440980977518368522240000e25 * t11 * t8 * t12 +
             0.14962697688175513436160000e26 * t16 * t8 * x +
             0.14962697688175513436160000e26 * t20 * t6 * t12 +
             0.5611011633065817538560000e25 * t24 * t7 +
             0.701376454133227192320000e24 * t27 * t12 +
             0.19482679281478533120000e23 * t30 * x) *
            x / 0.372032887351683317760000e24 -
          t1 * t4 * t3 *
            (0.18136603258394561740800000e26 * t8 * t6 +
             0.119701581505404107489280000e27 * t11 * t8 * t2 +
             0.224440465322632701542400000e27 * t16 * t8 +
             0.149626976881755134361600000e27 * t20 * t6 * t2 +
             0.35068822706661359616000000e26 * t24 * t6 +
             0.2337921513777423974400000e25 * t27 * t2 +
             0.19482679281478533120000e23 * t30) /
            0.1116098662055049953280000e25);
}



double
AssLegFunction::P_16_4(const double x) const
{
  double t21;
  double t14;
  double t6;
  double t7;
  double t3;
  double t1;
  double t10;
  double t2;
  double t4;
  t1  = sqrt(0.88179e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t21 = t14 * t14;
  return (t1 * t4 *
          (0.18136603258394561740800000e26 * t7 * t6 +
           0.119701581505404107489280000e27 * t10 * t7 * t2 +
           0.224440465322632701542400000e27 * t14 * t7 +
           0.149626976881755134361600000e27 * t14 * t10 * t6 * t2 +
           0.35068822706661359616000000e26 * t21 * t6 +
           0.2337921513777423974400000e25 * t21 * t10 * t2 +
           0.19482679281478533120000e23 * t21 * t14) /
          0.29018565213431298785280000e26);
}



double
AssLegFunction::P_16_4_Deriv(const double x) const
{
  double t1;
  double t23;
  double t2;
  double t3;
  double t13;
  double t34;
  double t32;
  double t5;
  double t6;
  double t16;
  double t20;
  double t9;
  t1  = sqrt(0.88179e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t9  = -t3;
  t13 = t9 * t9;
  t16 = t13 * t9;
  t20 = t13 * t13;
  t23 = t20 * t9;
  t32 = t3 * t3;
  t34 = t2 * x;
  return (-t1 * t3 *
            (0.18136603258394561740800000e26 * t6 * t5 +
             0.119701581505404107489280000e27 * t9 * t6 * t2 +
             0.224440465322632701542400000e27 * t13 * t6 +
             0.149626976881755134361600000e27 * t16 * t5 * t2 +
             0.35068822706661359616000000e26 * t20 * t5 +
             0.2337921513777423974400000e25 * t23 * t2 +
             0.19482679281478533120000e23 * t20 * t13) *
            x / 0.7254641303357824696320000e25 +
          t1 * t32 *
            (0.457042402111542955868160000e27 * t6 * t34 +
             0.2094777676344571881062400000e28 * t9 * t6 * x +
             0.2693285583871592418508800000e28 * t13 * t5 * t34 +
             0.1178312442943821683097600000e28 * t16 * t5 * x +
             0.163654505964419678208000000e27 * t20 * t34 +
             0.4909635178932590346240000e25 * t23 * x) /
            0.29018565213431298785280000e26);
}



double
AssLegFunction::P_16_5(const double x) const
{
  double t9;
  double t8;
  double t25;
  double t17;
  double t1;
  double t10;
  double t13;
  double t2;
  double t3;
  double t4;
  double t5;
  t1  = sqrt(0.12597e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * x;
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t17 = t13 * t13;
  t25 = t17 * t17;
  return (-t1 * t5 * t4 *
          (0.457042402111542955868160000e27 * t10 * t8 +
           0.2094777676344571881062400000e28 * t13 * t10 * x +
           0.2693285583871592418508800000e28 * t17 * t9 * t8 +
           0.1178312442943821683097600000e28 * t17 * t13 * t9 * x +
           0.163654505964419678208000000e27 * t8 * t25 +
           0.4909635178932590346240000e25 * t25 * t13 * x) /
          0.174111391280587792711680000e27);
}



double
AssLegFunction::P_16_5_Deriv(const double x) const
{
  double t20;
  double t9;
  double t27;
  double t1;
  double t8;
  double t7;
  double t24;
  double t12;
  double t34;
  double t2;
  double t4;
  double t3;
  double t16;
  t1  = sqrt(0.12597e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t24 = t16 * t16;
  t27 = t24 * t12;
  t34 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.457042402111542955868160000e27 * t7 * t9 +
             0.2094777676344571881062400000e28 * t12 * t9 * x +
             0.2693285583871592418508800000e28 * t16 * t8 * t7 +
             0.1178312442943821683097600000e28 * t20 * t8 * x +
             0.163654505964419678208000000e27 * t24 * t7 +
             0.4909635178932590346240000e25 * t27 * x) *
            x / 0.34822278256117558542336000e26 -
          t1 * t4 * t34 *
            (0.9217021775916116276674560000e28 * t9 * t2 +
             0.29626141422587516603596800000e29 * t12 * t9 +
             0.25922873744764077028147200000e29 * t16 * t8 * t2 +
             0.7200798262434465841152000000e28 * t20 * t8 +
             0.540059869682584938086400000e27 * t24 * t2 +
             0.4909635178932590346240000e25 * t27) /
            0.174111391280587792711680000e27);
}



double
AssLegFunction::P_16_6(const double x) const
{
  double t8;
  double t21;
  double t2;
  double t11;
  double t4;
  double t7;
  double t1;
  double t14;
  double t3;
  t1  = sqrt(0.25194e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t21 = t14 * t14;
  return (t1 * t4 * t3 *
          (0.9217021775916116276674560000e28 * t8 * t2 +
           0.29626141422587516603596800000e29 * t11 * t8 +
           0.25922873744764077028147200000e29 * t14 * t7 * t2 +
           0.7200798262434465841152000000e28 * t14 * t11 * t7 +
           0.540059869682584938086400000e27 * t21 * t2 +
           0.4909635178932590346240000e25 * t21 * t11) /
          0.3830450608172931439656960000e28);
}



double
AssLegFunction::P_16_6_Deriv(const double x) const
{
  double t13;
  double t17;
  double t1;
  double t10;
  double t2;
  double t3;
  double t4;
  double t33;
  double t6;
  double t7;
  double t20;
  t1  = sqrt(0.25194e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t13 = t10 * t10;
  t17 = t10 * t13;
  t20 = t13 * t13;
  t33 = t2 * x;
  return (-t1 * t4 *
            (0.9217021775916116276674560000e28 * t7 * t2 +
             0.29626141422587516603596800000e29 * t10 * t7 +
             0.25922873744764077028147200000e29 * t13 * t6 * t2 +
             0.7200798262434465841152000000e28 * t17 * t6 +
             0.540059869682584938086400000e27 * t20 * t2 +
             0.4909635178932590346240000e25 * t20 * t10) *
            x / 0.638408434695488573276160000e27 +
          t1 * t4 * t3 *
            (0.151422500604336195973939200000e30 * t7 * x +
             0.340700626359756440941363200000e30 * t10 * t6 * t33 +
             0.198742032043191257215795200000e30 * t13 * t6 * x +
             0.33123672007198542869299200000e29 * t17 * t33 +
             0.1129216091154495779635200000e28 * t20 * x) /
            0.3830450608172931439656960000e28);
}



double
AssLegFunction::P_16_7(const double x) const
{
  double t6;
  double t18;
  double t1;
  double t25;
  double t2;
  double t13;
  double t3;
  double t9;
  double t10;
  double t4;
  double t14;
  t1  = sqrt(0.1448655e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t25 = t18 * t18;
  return (-t1 * t6 * t4 * t3 *
          (0.151422500604336195973939200000e30 * t10 * x +
           0.340700626359756440941363200000e30 * t13 * t9 * t14 +
           0.198742032043191257215795200000e30 * t18 * t9 * x +
           0.33123672007198542869299200000e29 * t18 * t13 * t14 +
           0.1129216091154495779635200000e28 * t25 * x) /
          0.440501819939887115560550400000e30);
}



double
AssLegFunction::P_16_7_Deriv(const double x) const
{
  double t21;
  double t8;
  double t9;
  double t24;
  double t12;
  double t2;
  double t13;
  double t4;
  double t3;
  double t5;
  double t17;
  double t1;
  t1  = sqrt(0.1448655e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t24 = t17 * t17;
  return (t1 * t5 * t4 *
            (0.151422500604336195973939200000e30 * t9 * x +
             0.340700626359756440941363200000e30 * t12 * t8 * t13 +
             0.198742032043191257215795200000e30 * t17 * t8 * x +
             0.33123672007198542869299200000e29 * t21 * t13 +
             0.1129216091154495779635200000e28 * t24 * x) *
            x / 0.62928831419983873651507200000e29 -
          t1 * t5 * t4 * t3 *
            (0.2044203758158538645648179200000e31 * t9 +
             0.3179872512691060115452723200000e31 * t12 * t8 * t2 +
             0.1192452192259147543294771200000e31 * t17 * t8 +
             0.108404744750831594844979200000e30 * t21 * t2 +
             0.1129216091154495779635200000e28 * t24) /
            0.440501819939887115560550400000e30);
}



double
AssLegFunction::P_16_8(const double x) const
{
  double t20;
  double t3;
  double t1;
  double t2;
  double t14;
  double t7;
  double t4;
  double t5;
  double t8;
  double t10;
  t1  = sqrt(0.965770e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t20 = t14 * t14;
  return (t1 * t5 *
          (0.2044203758158538645648179200000e31 * t8 +
           0.3179872512691060115452723200000e31 * t10 * t7 * t2 +
           0.1192452192259147543294771200000e31 * t14 * t7 +
           0.108404744750831594844979200000e30 * t14 * t10 * t2 +
           0.1129216091154495779635200000e28 * t20) /
          0.5286021839278645386726604800000e31);
}



double
AssLegFunction::P_16_8_Deriv(const double x) const
{
  double t26;
  double t7;
  double t17;
  double t28;
  double t8;
  double t20;
  double t1;
  double t2;
  double t3;
  double t10;
  double t4;
  double t14;
  t1  = sqrt(0.965770e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t20 = t14 * t14;
  t26 = t4 * t4;
  t28 = t2 * x;
  return (-t1 * t4 * t3 *
            (0.2044203758158538645648179200000e31 * t8 +
             0.3179872512691060115452723200000e31 * t10 * t7 * t2 +
             0.1192452192259147543294771200000e31 * t14 * t7 +
             0.108404744750831594844979200000e30 * t17 * t2 +
             0.1129216091154495779635200000e28 * t20) *
            x / 0.660752729909830673340825600000e30 +
          t1 * t26 *
            (0.22713375090650429396090880000000e32 * t7 * t28 +
             0.23849043845182950865895424000000e32 * t10 * t7 * x +
             0.5420237237541579742248960000000e31 * t14 * t28 +
             0.225843218230899155927040000000e30 * t17 * x) /
            0.5286021839278645386726604800000e31);
}



double
AssLegFunction::P_16_9(const double x) const
{
  double t17;
  double t9;
  double t10;
  double t2;
  double t1;
  double t13;
  double t3;
  double t4;
  double t5;
  double t6;
  t1  = sqrt(0.482885e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (-t1 * t6 * t5 *
          (0.22713375090650429396090880000000e32 * t9 * t10 +
           0.23849043845182950865895424000000e32 * t13 * t10 * x +
           0.5420237237541579742248960000000e31 * t17 * t9 +
           0.225843218230899155927040000000e30 * t17 * t13 * x) /
          0.52860218392786453867266048000000e32);
}



double
AssLegFunction::P_16_9_Deriv(const double x) const
{
  double t1;
  double t17;
  double t9;
  double t10;
  double t20;
  double t2;
  double t3;
  double t4;
  double t6;
  double t27;
  double t13;
  t1  = sqrt(0.482885e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  t20 = t17 * t13;
  t27 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.22713375090650429396090880000000e32 * t10 * t9 +
             0.23849043845182950865895424000000e32 * t13 * t10 * x +
             0.5420237237541579742248960000000e31 * t17 * t9 +
             0.225843218230899155927040000000e30 * t20 * x) *
            x / 0.5873357599198494874140672000000e31 -
          t1 * t6 * t27 *
            (0.206691713324918907504427008000000e33 * t10 * t2 +
             0.140926168176081073298472960000000e33 * t13 * t10 +
             0.17615771022010134162309120000000e32 * t17 * t2 +
             0.225843218230899155927040000000e30 * t20) /
            0.52860218392786453867266048000000e32);
}



double
AssLegFunction::P_16_10(const double x) const
{
  double t11;
  double t14;
  double t3;
  double t2;
  double t8;
  double t1;
  double t4;
  double t5;
  t1  = sqrt(0.520030e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (t1 * t5 * t3 *
          (0.206691713324918907504427008000000e33 * t8 * t2 +
           0.140926168176081073298472960000000e33 * t11 * t8 +
           0.17615771022010134162309120000000e32 * t14 * t2 +
           0.225843218230899155927040000000e30 * t14 * t11) /
          0.740043057499010354141724672000000e33);
}



double
AssLegFunction::P_16_10_Deriv(const double x) const
{
  double t2;
  double t7;
  double t3;
  double t13;
  double t4;
  double t5;
  double t1;
  double t10;
  t1  = sqrt(0.520030e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (-t1 * t5 *
            (0.206691713324918907504427008000000e33 * t7 * t2 +
             0.140926168176081073298472960000000e33 * t10 * t7 +
             0.17615771022010134162309120000000e32 * t13 * t2 +
             0.225843218230899155927040000000e30 * t13 * t10) *
            x / 0.74004305749901035414172467200000e32 +
          t1 * t5 * t3 *
            (0.1522002616301675591623507968000000e34 * t7 * x +
             0.634167756792364829843128320000000e33 * t10 * t2 * x +
             0.36586601353405663260180480000000e32 * t13 * x) /
            0.740043057499010354141724672000000e33);
}



double
AssLegFunction::P_16_11(const double x) const
{
  double t3;
  double t2;
  double t7;
  double t4;
  double t10;
  double t13;
  double t1;
  double t5;
  double t17;
  t1  = sqrt(0.260015e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (-t1 * t7 * t5 * t3 *
          (0.1522002616301675591623507968000000e34 * t10 * x +
           0.634167756792364829843128320000000e33 * t13 * t2 * x +
           0.36586601353405663260180480000000e32 * t17 * x) /
          0.6660387517491093187275522048000000e34);
}



double
AssLegFunction::P_16_11_Deriv(const double x) const
{
  double t9;
  double t16;
  double t12;
  double t5;
  double t6;
  double t1;
  double t4;
  double t2;
  double t3;
  t1  = sqrt(0.260015e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (t1 * t6 * t5 *
            (0.1522002616301675591623507968000000e34 * t9 * x +
             0.634167756792364829843128320000000e33 * t12 * t2 * x +
             0.36586601353405663260180480000000e32 * t16 * x) *
            x / 0.605489774317372107934138368000000e33 -
          t1 * t6 * t5 * t3 *
            (0.8878348595093107617803796480000000e34 * t9 +
             0.2048849675790717142570106880000000e34 * t12 * t2 +
             0.36586601353405663260180480000000e32 * t16) /
            0.6660387517491093187275522048000000e34);
}



double
AssLegFunction::P_16_12(const double x) const
{
  double t8;
  double t10;
  double t3;
  double t2;
  double t13;
  double t4;
  double t5;
  double t1;
  t1  = sqrt(0.7429e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t5 * t4 *
          (0.8878348595093107617803796480000000e34 * t8 +
           0.2048849675790717142570106880000000e34 * t10 * t2 +
           0.36586601353405663260180480000000e32 * t13) /
          0.13320775034982186374551044096000000e35);
}



double
AssLegFunction::P_16_12_Deriv(const double x) const
{
  double t1;
  double t5;
  double t2;
  double t3;
  double t4;
  double t13;
  double t8;
  double t10;
  t1  = sqrt(0.7429e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (-t1 * t5 * t3 *
            (0.8878348595093107617803796480000000e34 * t8 +
             0.2048849675790717142570106880000000e34 * t10 * t2 +
             0.36586601353405663260180480000000e32 * t13) *
            x / 0.1110064586248515531212587008000000e34 +
          t1 * t5 * t4 *
            (0.39611093731953864756355399680000000e35 * t2 * x +
             0.4244045756995056938180935680000000e34 * t10 * x) /
            0.13320775034982186374551044096000000e35);
}



double
AssLegFunction::P_16_13(const double x) const
{
  double t5;
  double t7;
  double t2;
  double t1;
  double t3;
  double t4;
  t1 = sqrt(0.215441e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (-t1 * t7 * t5 * t4 *
          (0.39611093731953864756355399680000000e35 * t2 * x -
           0.4244045756995056938180935680000000e34 * t3 * x) /
          0.772604952028966809723960557568000000e36);
}



double
AssLegFunction::P_16_13_Deriv(const double x) const
{
  double t7;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.215441e6);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (t1 * t7 * t5 * t3 *
            (0.39611093731953864756355399680000000e35 * t2 * x -
             0.4244045756995056938180935680000000e34 * t3 * x) *
            x / 0.59431150156074369978766196736000000e35 -
          t1 * t7 * t5 * t4 *
            (0.131565418466846765083609006080000000e36 * t2 -
             0.4244045756995056938180935680000000e34) /
            0.772604952028966809723960557568000000e36);
}



double
AssLegFunction::P_16_14(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  double t6;
  t1 = sqrt(0.2154410e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
          (0.131565418466846765083609006080000000e36 * t2 -
           0.4244045756995056938180935680000000e34) /
          0.23178148560869004291718816727040000000e38);
}



double
AssLegFunction::P_16_14_Deriv(const double x) const
{
  double t1;
  double t4;
  double t5;
  double t2;
  double t3;
  t1 = sqrt(0.2154410e7);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-t1 * t5 * t4 *
            (0.131565418466846765083609006080000000e36 * t2 -
             0.4244045756995056938180935680000000e34) *
            x / 0.1655582040062071735122772623360000000e37 +
          0.93e2 / 0.8192e4 * t1 * t5 * t4 * t3 * x);
}



double
AssLegFunction::P_16_15(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  double t1;
  t1 = sqrt(0.33393355e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  t8 = sqrt(t3);
  return (-0.3e1 / 0.8192e4 * t1 * t8 * t6 * t4 * t3 * x);
}



double
AssLegFunction::P_16_15_Deriv(const double x) const
{
  double t7;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.33393355e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (0.45e2 / 0.8192e4 * t1 * t7 * t5 * t4 * t2 -
          0.3e1 / 0.8192e4 * t1 * t7 * t5 * t4 * t3);
}



double
AssLegFunction::P_16_16(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t5;
  double t6;
  t1 = sqrt(0.66786710e8);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (0.3e1 / 0.65536e5 * t1 * t6);
}



double
AssLegFunction::P_16_16_Deriv(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  double t6;
  t1 = sqrt(0.66786710e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  return (-0.3e1 / 0.4096e4 * t1 * t6 * t4 * t3 * x);
}

double
AssLegFunction::P_17_0(const double x) const
{
  double t6;
  double t17;
  double t12;
  double t1;
  double t8;
  double t7;
  double t21;
  double t13;
  double t2;
  double t3;
  double t4;
  double t34;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t6  = t1 - 0.1e1;
  t7  = t1 * x;
  t8  = t2 * t7;
  t12 = t6 * t6;
  t13 = t2 * x;
  t17 = t12 * t6;
  t21 = t12 * t12;
  t34 = t21 * t21;
  return (
    t4 * x + 0.68e2 * t6 * t3 * t8 + 0.1785e4 / 0.2e1 * t12 * t3 * t13 +
    0.7735e4 / 0.2e1 * t17 * t3 * t7 + 0.425425e6 / 0.64e2 * t21 * t3 * x +
    0.153153e6 / 0.32e2 * t21 * t6 * t8 +
    0.357357e6 / 0.256e3 * t21 * t12 * t13 +
    0.36465e5 / 0.256e3 * t21 * t17 * t7 + 0.109395e6 / 0.32768e5 * t34 * x);
}

double
AssLegFunction::P_17_0_Deriv(const double x) const
{
  double t11;
  double t31;
  double t7;
  double t6;
  double t15;
  double t19;
  double t4;
  double t1;
  double t2;
  double t3;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t6  = t1 - 0.1e1;
  t7  = t2 * t1;
  t11 = t6 * t6;
  t15 = t11 * t6;
  t19 = t11 * t11;
  t31 = t19 * t19;
  return (
    0.153e3 * t4 + 0.4590e4 * t6 * t3 * t7 + 0.69615e5 / 0.2e1 * t11 * t3 * t2 +
    0.765765e6 / 0.8e1 * t15 * t3 * t1 + 0.6891885e7 / 0.64e2 * t19 * t3 +
    0.3216213e7 / 0.64e2 * t19 * t6 * t7 +
    0.2297295e7 / 0.256e3 * t19 * t11 * t2 +
    0.984555e6 / 0.2048e4 * t19 * t15 * t1 + 0.109395e6 / 0.32768e5 * t31);
}



double
AssLegFunction::P_17_1(const double x) const
{
  double t23;
  double t19;
  double t4;
  double t6;
  double t1;
  double t7;
  double t15;
  double t10;
  double t2;
  double t3;
  double t35;
  double t11;
  double t8;
  t1  = sqrt(0.34e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t10 = -t3;
  t11 = t6 * t2;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t23 = t15 * t15;
  t35 = t23 * t23;
  return (-t1 * t4 *
          (0.7132961374036033536000e22 * t8 +
           0.213988841221081006080000e24 * t10 * t7 * t11 +
           0.1622748712593197629440000e25 * t15 * t7 * t6 +
           0.4462558959631293480960000e25 * t19 * t7 * t2 +
           0.5020378829585205166080000e25 * t23 * t7 +
           0.2342843453806429077504000e25 * t23 * t10 * t11 +
           0.418364902465433763840000e24 * t23 * t15 * t6 +
           0.22412405489219665920000e23 * t23 * t19 * t2 +
           0.155641704786247680000e21 * t35) /
          0.4755307582690689024000e22);
}



double
AssLegFunction::P_17_1_Deriv(const double x) const
{
  double t11;
  double t42;
  double t43;
  double t12;
  double t24;
  double t33;
  double t46;
  double t27;
  double t16;
  double t1;
  double t2;
  double t4;
  double t3;
  double t30;
  double t20;
  double t8;
  double t7;
  double t36;
  double t9;
  t1  = sqrt(0.34e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t11 = -t3;
  t12 = t7 * t2;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t24 = t16 * t16;
  t27 = t24 * t11;
  t30 = t24 * t16;
  t33 = t24 * t20;
  t36 = t24 * t24;
  t42 = t2 * x;
  t43 = t7 * t42;
  t46 = t7 * x;
  return (t1 / t4 *
            (0.7132961374036033536000e22 * t9 +
             0.213988841221081006080000e24 * t11 * t8 * t12 +
             0.1622748712593197629440000e25 * t16 * t8 * t7 +
             0.4462558959631293480960000e25 * t20 * t8 * t2 +
             0.5020378829585205166080000e25 * t24 * t8 +
             0.2342843453806429077504000e25 * t27 * t12 +
             0.418364902465433763840000e24 * t30 * t7 +
             0.22412405489219665920000e23 * t33 * t2 +
             0.155641704786247680000e21 * t36) *
            x / 0.4755307582690689024000e22 -
          t1 * t4 *
            (0.542105064426738548736000e24 * t8 * t43 +
             0.9486838627467924602880000e25 * t11 * t8 * t46 +
             0.46248338308906132439040000e26 * t16 * t8 * t42 +
             0.84788620232994576138240000e26 * t20 * t8 * x +
             0.63591465174745932103680000e26 * t24 * t43 +
             0.19077439552423779631104000e26 * t46 * t27 +
             0.1987233286710810378240000e25 * t30 * t42 +
             0.47315078255019294720000e23 * t33 * x) /
            0.4755307582690689024000e22);
}



double
AssLegFunction::P_17_2(const double x) const
{
  double t7;
  double t8;
  double t12;
  double t1;
  double t3;
  double t6;
  double t20;
  double t2;
  double t5;
  double t11;
  double t16;
  double t24;
  t1  = sqrt(0.646e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t7  = t6 * t5;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t6 * x;
  t16 = t11 * t11;
  t20 = t11 * t16;
  t24 = t16 * t16;
  return (t1 * t3 *
          (0.542105064426738548736000e24 * t8 * t7 +
           0.9486838627467924602880000e25 * t11 * t8 * t12 +
           0.46248338308906132439040000e26 * t16 * t8 * t5 +
           0.84788620232994576138240000e26 * t20 * t8 * x +
           0.63591465174745932103680000e26 * t24 * t7 +
           0.19077439552423779631104000e26 * t24 * t11 * t12 +
           0.1987233286710810378240000e25 * t24 * t16 * t5 +
           0.47315078255019294720000e23 * t24 * t20 * x) /
          0.361403376284492365824000e24);
}



double
AssLegFunction::P_17_2_Deriv(const double x) const
{
  double t29;
  double t11;
  double t15;
  double t7;
  double t10;
  double t3;
  double t6;
  double t4;
  double t5;
  double t40;
  double t32;
  double t23;
  double t26;
  double t19;
  double t1;
  t1  = sqrt(0.646e3);
  t3  = x * x;
  t4  = t3 * x;
  t5  = t3 * t3;
  t6  = t5 * t4;
  t7  = t5 * t5;
  t10 = t3 - 0.1e1;
  t11 = t5 * x;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t23 = t15 * t15;
  t26 = t23 * t10;
  t29 = t23 * t15;
  t32 = t23 * t19;
  t40 = t3 * t5;
  return (-t1 * x *
            (0.542105064426738548736000e24 * t7 * t6 +
             0.9486838627467924602880000e25 * t10 * t7 * t11 +
             0.46248338308906132439040000e26 * t15 * t7 * t4 +
             0.84788620232994576138240000e26 * t19 * t7 * x +
             0.63591465174745932103680000e26 * t23 * t6 +
             0.19077439552423779631104000e26 * t26 * t11 +
             0.1987233286710810378240000e25 * t29 * t4 +
             0.47315078255019294720000e23 * t32 * x) /
            0.180701688142246182912000e24 -
          t1 * t10 *
            (0.27105253221336927436800000e26 * t7 * t40 +
             0.308322255392707549593600000e27 * t10 * t7 * t5 +
             0.1017463442795934913658880000e28 * t15 * t7 * t3 +
             0.1271829303494918642073600000e28 * t19 * t7 +
             0.635914651747459321036800000e27 * t23 * t40 +
             0.119233997202648622694400000e27 * t26 * t5 +
             0.6624110955702701260800000e25 * t29 * t3 +
             0.47315078255019294720000e23 * t32) /
            0.361403376284492365824000e24);
}



double
AssLegFunction::P_17_3(const double x) const
{
  double t2;
  double t3;
  double t16;
  double t20;
  double t7;
  double t12;
  double t4;
  double t8;
  double t1;
  double t23;
  double t9;
  t1  = sqrt(0.1938e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t2;
  t9  = t7 * t7;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t23 = t16 * t16;
  return (-t1 * t4 * t3 *
          (0.27105253221336927436800000e26 * t9 * t8 +
           0.308322255392707549593600000e27 * t12 * t9 * t7 +
           0.1017463442795934913658880000e28 * t16 * t9 * t2 +
           0.1271829303494918642073600000e28 * t20 * t9 +
           0.635914651747459321036800000e27 * t23 * t8 +
           0.119233997202648622694400000e27 * t23 * t12 * t7 +
           0.6624110955702701260800000e25 * t23 * t16 * t2 +
           0.47315078255019294720000e23 * t23 * t20) /
          0.10842101288534770974720000e26);
}



double
AssLegFunction::P_17_3_Deriv(const double x) const
{
  double t6;
  double t1;
  double t19;
  double t7;
  double t25;
  double t8;
  double t28;
  double t2;
  double t15;
  double t3;
  double t4;
  double t42;
  double t22;
  double t11;
  double t39;
  t1  = sqrt(0.1938e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t2;
  t8  = t6 * t6;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t22 = t15 * t15;
  t25 = t22 * t11;
  t28 = t22 * t15;
  t39 = t6 * x;
  t42 = t2 * x;
  return (t1 * t4 *
            (0.27105253221336927436800000e26 * t8 * t7 +
             0.308322255392707549593600000e27 * t11 * t8 * t6 +
             0.1017463442795934913658880000e28 * t15 * t8 * t2 +
             0.1271829303494918642073600000e28 * t19 * t8 +
             0.635914651747459321036800000e27 * t22 * t7 +
             0.119233997202648622694400000e27 * t25 * t6 +
             0.6624110955702701260800000e25 * t28 * t2 +
             0.47315078255019294720000e23 * t22 * t19) *
            x / 0.3614033762844923658240000e25 -
          t1 * t4 * t3 *
            (0.996118055884132083302400000e27 * t8 * t39 +
             0.7769720835896230249758720000e28 * t11 * t8 * t42 +
             0.17805610248928860989030400000e29 * t15 * t8 * x +
             0.15261951641939023704883200000e29 * t19 * t6 * t42 +
             0.5007827882511242153164800000e28 * t22 * t39 +
             0.556425320279026905907200000e27 * t25 * t42 +
             0.13910633006975672647680000e26 * t28 * x) /
            0.10842101288534770974720000e26);
}



double
AssLegFunction::P_17_4(const double x) const
{
  double t24;
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t16;
  double t12;
  double t6;
  double t7;
  double t8;
  t1  = sqrt(0.323e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * x;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t24 = t16 * t16;
  return (t1 * t4 *
          (0.996118055884132083302400000e27 * t7 * t8 +
           0.7769720835896230249758720000e28 * t11 * t8 * t12 +
           0.17805610248928860989030400000e29 * t16 * t8 * x +
           0.15261951641939023704883200000e29 * t16 * t11 * t6 * t12 +
           0.5007827882511242153164800000e28 * t24 * t7 +
           0.556425320279026905907200000e27 * t24 * t11 * t12 +
           0.13910633006975672647680000e26 * t24 * t16 * x) /
          0.75894709019743396823040000e26);
}



double
AssLegFunction::P_17_4_Deriv(const double x) const
{
  double t5;
  double t23;
  double t15;
  double t6;
  double t7;
  double t1;
  double t26;
  double t36;
  double t10;
  double t11;
  double t3;
  double t2;
  double t29;
  double t19;
  t1  = sqrt(0.323e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * x;
  t7  = t5 * t5;
  t10 = -t3;
  t11 = t2 * x;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t23 = t15 * t15;
  t26 = t23 * t10;
  t29 = t23 * t15;
  t36 = t3 * t3;
  return (-t1 * t3 *
            (0.996118055884132083302400000e27 * t7 * t6 +
             0.7769720835896230249758720000e28 * t10 * t7 * t11 +
             0.17805610248928860989030400000e29 * t15 * t7 * x +
             0.15261951641939023704883200000e29 * t19 * t5 * t11 +
             0.5007827882511242153164800000e28 * t23 * t6 +
             0.556425320279026905907200000e27 * t26 * t11 +
             0.13910633006975672647680000e26 * t29 * x) *
            x / 0.18973677254935849205760000e26 +
          t1 * t36 *
            (0.28488976398286177582448640000e29 * t7 * t5 +
             0.156689370190573976703467520000e30 * t10 * t7 * t2 +
             0.251822202091993891130572800000e30 * t15 * t7 +
             0.146896284553663103159500800000e30 * t19 * t5 * t2 +
             0.30603392615346479824896000000e29 * t23 * t5 +
             0.1836203556920788789493760000e28 * t26 * t2 +
             0.13910633006975672647680000e26 * t29) /
            0.75894709019743396823040000e26);
}



double
AssLegFunction::P_17_5(const double x) const
{
  double t9;
  double t5;
  double t3;
  double t23;
  double t12;
  double t4;
  double t2;
  double t16;
  double t8;
  double t1;
  t1  = sqrt(0.92378e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t23 = t16 * t16;
  return (-t1 * t5 * t4 *
          (0.28488976398286177582448640000e29 * t9 * t8 +
           0.156689370190573976703467520000e30 * t12 * t9 * t2 +
           0.251822202091993891130572800000e30 * t16 * t9 +
           0.146896284553663103159500800000e30 * t16 * t12 * t8 * t2 +
           0.30603392615346479824896000000e29 * t23 * t8 +
           0.1836203556920788789493760000e28 * t23 * t12 * t2 +
           0.13910633006975672647680000e26 * t23 * t16) /
          0.21705886779646611491389440000e29);
}



double
AssLegFunction::P_17_5_Deriv(const double x) const
{
  double t2;
  double t4;
  double t3;
  double t11;
  double t22;
  double t37;
  double t34;
  double t25;
  double t15;
  double t1;
  double t7;
  double t8;
  double t18;
  t1  = sqrt(0.92378e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t11 * t15;
  t22 = t15 * t15;
  t25 = t22 * t11;
  t34 = t3 * t3;
  t37 = t2 * x;
  return (t1 * t4 * t3 *
            (0.28488976398286177582448640000e29 * t7 * t8 +
             0.156689370190573976703467520000e30 * t11 * t8 * t2 +
             0.251822202091993891130572800000e30 * t15 * t8 +
             0.146896284553663103159500800000e30 * t18 * t7 * t2 +
             0.30603392615346479824896000000e29 * t22 * t7 +
             0.1836203556920788789493760000e28 * t25 * t2 +
             0.13910633006975672647680000e26 * t22 * t15) *
            x / 0.4341177355929322298277888000e28 -
          t1 * t4 * t34 *
            (0.655246457160582084396318720000e30 * t8 * t37 +
             0.2574182510273715331556966400000e31 * t11 * t8 * x +
             0.2895955324057929748001587200000e31 * t15 * t7 * t37 +
             0.1126204848244750457556172800000e31 * t18 * t7 * x +
             0.140775606030593807194521600000e30 * t22 * t37 +
             0.3839334709925285650759680000e28 * t25 * x) /
            0.21705886779646611491389440000e29);
}



double
AssLegFunction::P_17_6(const double x) const
{
  double t4;
  double t16;
  double t24;
  double t1;
  double t8;
  double t9;
  double t7;
  double t3;
  double t12;
  double t2;
  t1  = sqrt(0.6374082e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t24 = t16 * t16;
  return (t1 * t4 * t3 *
          (0.655246457160582084396318720000e30 * t7 * t9 +
           0.2574182510273715331556966400000e31 * t12 * t9 * x +
           0.2895955324057929748001587200000e31 * t16 * t8 * t7 +
           0.1126204848244750457556172800000e31 * t16 * t12 * t8 * x +
           0.140775606030593807194521600000e30 * t24 * t7 +
           0.3839334709925285650759680000e28 * t24 * t12 * x) /
          0.2995412375591232385811742720000e31);
}



double
AssLegFunction::P_17_6_Deriv(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t4;
  double t6;
  double t7;
  double t8;
  double t11;
  double t15;
  double t19;
  double t23;
  double t26;
  t1  = sqrt(0.6374082e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t23 = t15 * t15;
  t26 = t23 * t11;
  return (-t1 * t4 *
            (0.655246457160582084396318720000e30 * t8 * t6 +
             0.2574182510273715331556966400000e31 * t11 * t8 * x +
             0.2895955324057929748001587200000e31 * t15 * t7 * t6 +
             0.1126204848244750457556172800000e31 * t19 * t7 * x +
             0.140775606030593807194521600000e30 * t23 * t6 +
             0.3839334709925285650759680000e28 * t26 * x) *
            x / 0.499235395931872064301957120000e30 +
          t1 * t4 * t3 *
            (0.12356076049313833591473438720000e32 * t2 * t8 +
             0.34751463888695156976019046400000e32 * t11 * t8 +
             0.27028916357874010981348147200000e32 * t15 * t7 * t2 +
             0.6757229089468502745337036800000e31 * t19 * t7 +
             0.460720165191034278091161600000e30 * t23 * t2 +
             0.3839334709925285650759680000e28 * t26) /
            0.2995412375591232385811742720000e31);
}



double
AssLegFunction::P_17_7(const double x) const
{
  double t6;
  double t13;
  double t16;
  double t10;
  double t9;
  double t23;
  double t1;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.96577e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t16 = t13 * t13;
  t23 = t16 * t16;
  return (-t1 * t6 * t4 * t3 *
          (0.12356076049313833591473438720000e32 * t10 * t2 +
           0.34751463888695156976019046400000e32 * t13 * t10 +
           0.27028916357874010981348147200000e32 * t16 * t9 * t2 +
           0.6757229089468502745337036800000e31 * t16 * t13 * t9 +
           0.460720165191034278091161600000e30 * t23 * t2 +
           0.3839334709925285650759680000e28 * t23 * t13) /
          0.5990824751182464771623485440000e31);
}



double
AssLegFunction::P_17_7_Deriv(const double x) const
{
  double t36;
  double t1;
  double t2;
  double t5;
  double t3;
  double t8;
  double t4;
  double t15;
  double t12;
  double t9;
  double t19;
  double t22;
  t1  = sqrt(0.96577e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t15 = t12 * t12;
  t19 = t15 * t12;
  t22 = t15 * t15;
  t36 = t2 * x;
  return (t1 * t5 * t4 *
            (0.12356076049313833591473438720000e32 * t9 * t2 +
             0.34751463888695156976019046400000e32 * t12 * t9 +
             0.27028916357874010981348147200000e32 * t15 * t8 * t2 +
             0.6757229089468502745337036800000e31 * t19 * t8 +
             0.460720165191034278091161600000e30 * t22 * t2 +
             0.3839334709925285650759680000e28 * t22 * t12) *
            x / 0.855832107311780681660497920000e30 -
          t1 * t5 * t4 * t3 *
            (0.193063688270528649866772480000000e33 * t9 * x +
             0.386127376541057299733544960000000e33 * t12 * t8 * t36 +
             0.202716872684055082360111104000000e33 * t15 * t8 * x +
             0.30714677679402285206077440000000e32 * t19 * t36 +
             0.959833677481321412689920000000e30 * t22 * x) /
            0.5990824751182464771623485440000e31);
}



double
AssLegFunction::P_17_8(const double x) const
{
  double t16;
  double t1;
  double t12;
  double t2;
  double t23;
  double t3;
  double t4;
  double t5;
  double t7;
  double t8;
  double t11;
  t1  = sqrt(0.965770e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t23 = t16 * t16;
  return (t1 * t5 *
          (0.193063688270528649866772480000000e33 * t8 * x +
           0.386127376541057299733544960000000e33 * t11 * t7 * t12 +
           0.202716872684055082360111104000000e33 * t16 * t7 * x +
           0.30714677679402285206077440000000e32 * t16 * t11 * t12 +
           0.959833677481321412689920000000e30 * t23 * x) /
          0.299541237559123238581174272000000e33);
}



double
AssLegFunction::P_17_8_Deriv(const double x) const
{
  double t20;
  double t11;
  double t12;
  double t1;
  double t30;
  double t7;
  double t3;
  double t2;
  double t23;
  double t4;
  double t16;
  double t8;
  t1  = sqrt(0.965770e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t23 = t16 * t16;
  t30 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.193063688270528649866772480000000e33 * t8 * x +
             0.386127376541057299733544960000000e33 * t11 * t7 * t12 +
             0.202716872684055082360111104000000e33 * t16 * t7 * x +
             0.30714677679402285206077440000000e32 * t20 * t12 +
             0.959833677481321412689920000000e30 * t23 * x) *
            x / 0.37442654694890404822646784000000e32 +
          t1 * t30 *
            (0.2509827947516872448268042240000000e34 * t8 +
             0.3513759126523621427575259136000000e34 * t11 * t7 * t2 +
             0.1197872429496689123037020160000000e34 * t16 * t7 +
             0.99822702458057426919751680000000e32 * t20 * t2 +
             0.959833677481321412689920000000e30 * t23) /
            0.299541237559123238581174272000000e33);
}



double
AssLegFunction::P_17_9(const double x) const
{
  double t12;
  double t4;
  double t10;
  double t9;
  double t6;
  double t3;
  double t5;
  double t22;
  double t1;
  double t2;
  double t16;
  t1  = sqrt(0.37145e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t12 = -t3;
  t16 = t12 * t12;
  t22 = t16 * t16;
  return (-t1 * t6 * t5 *
          (0.2509827947516872448268042240000000e34 * t10 +
           0.3513759126523621427575259136000000e34 * t12 * t9 * t2 +
           0.1197872429496689123037020160000000e34 * t16 * t9 +
           0.99822702458057426919751680000000e32 * t16 * t12 * t2 +
           0.959833677481321412689920000000e30 * t22) /
          0.898623712677369715743522816000000e33);
}



double
AssLegFunction::P_17_9_Deriv(const double x) const
{
  double t22;
  double t1;
  double t9;
  double t10;
  double t28;
  double t12;
  double t3;
  double t2;
  double t4;
  double t19;
  double t16;
  double t6;
  double t31;
  t1  = sqrt(0.37145e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t12 = -t3;
  t16 = t12 * t12;
  t19 = t16 * t12;
  t22 = t16 * t16;
  t28 = t4 * t4;
  t31 = t2 * x;
  return (t1 * t6 * t4 * t3 *
            (0.2509827947516872448268042240000000e34 * t10 +
             0.3513759126523621427575259136000000e34 * t12 * t9 * t2 +
             0.1197872429496689123037020160000000e34 * t16 * t9 +
             0.99822702458057426919751680000000e32 * t19 * t2 +
             0.959833677481321412689920000000e30 * t22) *
            x / 0.99847079186374412860391424000000e32 -
          t1 * t6 * t28 *
            (0.27106141833182222441294856192000000e35 * t9 * t31 +
             0.25874044477128485057599635456000000e35 * t12 * t9 * x +
             0.5390425932735101053666590720000000e34 * t16 * t31 +
             0.207324074335965425141022720000000e33 * t19 * x) /
            0.898623712677369715743522816000000e33);
}



double
AssLegFunction::P_17_10(const double x) const
{
  double t3;
  double t9;
  double t5;
  double t4;
  double t12;
  double t1;
  double t16;
  double t8;
  double t2;
  t1  = sqrt(0.222870e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * x;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (t1 * t5 * t3 *
          (0.27106141833182222441294856192000000e35 * t9 * t8 +
           0.25874044477128485057599635456000000e35 * t12 * t9 * x +
           0.5390425932735101053666590720000000e34 * t16 * t8 +
           0.207324074335965425141022720000000e33 * t16 * t12 * x) /
          0.32350453656385309766766821376000000e35);
}



double
AssLegFunction::P_17_10_Deriv(const double x) const
{
  double t3;
  double t2;
  double t15;
  double t11;
  double t4;
  double t5;
  double t1;
  double t8;
  double t7;
  double t18;
  t1  = sqrt(0.222870e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * x;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  return (-t1 * t5 *
            (0.27106141833182222441294856192000000e35 * t8 * t7 +
             0.25874044477128485057599635456000000e35 * t11 * t8 * x +
             0.5390425932735101053666590720000000e34 * t15 * t7 +
             0.207324074335965425141022720000000e33 * t18 * x) *
            x / 0.3235045365638530976676682137600000e34 +
          t1 * t5 * t3 *
            (0.241491081786532527204263264256000000e36 * t2 * t8 +
             0.150931926116582829502664540160000000e36 * t11 * t8 +
             0.17415222244221095711845908480000000e35 * t2 * t15 +
             0.207324074335965425141022720000000e33 * t18) /
            0.32350453656385309766766821376000000e35);
}



double
AssLegFunction::P_17_11(const double x) const
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t10;
  double t5;
  double t13;
  double t7;
  t1  = sqrt(0.222870e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t16 = t13 * t13;
  return (-t1 * t7 * t5 * t3 *
          (0.241491081786532527204263264256000000e36 * t10 * t2 +
           0.150931926116582829502664540160000000e36 * t13 * t10 +
           0.17415222244221095711845908480000000e35 * t16 * t2 +
           0.207324074335965425141022720000000e33 * t16 * t13) /
          0.452906351189394336734735499264000000e36);
}



double
AssLegFunction::P_17_11_Deriv(const double x) const
{
  double t1;
  double t9;
  double t2;
  double t12;
  double t4;
  double t3;
  double t5;
  double t6;
  double t15;
  t1  = sqrt(0.222870e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (t1 * t6 * t5 *
            (0.241491081786532527204263264256000000e36 * t9 * t2 +
             0.150931926116582829502664540160000000e36 * t12 * t9 +
             0.17415222244221095711845908480000000e35 * t15 * t2 +
             0.207324074335965425141022720000000e33 * t15 * t12) *
            x / 0.41173304653581303339521409024000000e35 -
          t1 * t6 * t5 * t3 *
            (0.1750810342952360822230908665856000000e37 * t9 * x +
             0.673388593443215700858041794560000000e36 * t12 * t2 * x +
             0.36074388934457983974537953280000000e35 * t15 * x) /
            0.452906351189394336734735499264000000e36);
}



double
AssLegFunction::P_17_12(const double x) const
{
  double t15;
  double t3;
  double t4;
  double t11;
  double t5;
  double t2;
  double t1;
  double t8;
  t1  = sqrt(0.1077205e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (t1 * t5 * t4 *
          (0.1750810342952360822230908665856000000e37 * t8 * x +
           0.673388593443215700858041794560000000e36 * t11 * t2 * x +
           0.36074388934457983974537953280000000e35 * t15 * x) /
          0.13134284184492435765307329478656000000e38);
}



double
AssLegFunction::P_17_12_Deriv(const double x) const
{
  double t11;
  double t1;
  double t15;
  double t2;
  double t3;
  double t4;
  double t5;
  double t8;
  t1  = sqrt(0.1077205e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (-t1 * t5 * t3 *
            (0.1750810342952360822230908665856000000e37 * t8 * x +
             0.673388593443215700858041794560000000e36 * t11 * t2 * x +
             0.36074388934457983974537953280000000e35 * t15 * x) *
            x / 0.1094523682041036313775610789888000000e37 +
          t1 * t5 * t4 *
            (0.10100828901648235512870626918400000000e38 * t8 +
             0.2164463336067479038472277196800000000e37 * t11 * t2 +
             0.36074388934457983974537953280000000e35 * t15) /
            0.13134284184492435765307329478656000000e38);
}



double
AssLegFunction::P_17_13(const double x) const
{
  double t12;
  double t7;
  double t10;
  double t15;
  double t2;
  double t3;
  double t4;
  double t5;
  double t1;
  t1  = sqrt(0.6463230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (-t1 * t7 * t5 * t4 *
          (0.10100828901648235512870626918400000000e38 * t10 +
           0.2164463336067479038472277196800000000e37 * t12 * t2 +
           0.36074388934457983974537953280000000e35 * t15) /
          0.394028525534773072959219884359680000000e39);
}



double
AssLegFunction::P_17_13_Deriv(const double x) const
{
  double t10;
  double t12;
  double t15;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  t1  = sqrt(0.6463230e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (t1 * t7 * t5 * t3 *
            (0.10100828901648235512870626918400000000e38 * t10 +
             0.2164463336067479038472277196800000000e37 * t12 * t2 +
             0.36074388934457983974537953280000000e35 * t15) *
            x / 0.30309886579597928689170760335360000000e38 -
          t1 * t7 * t5 * t4 *
            (0.44732242278727900128427062067200000000e38 * t2 * x +
             0.4473224227872790012842706206720000000e37 * t12 * x) /
            0.394028525534773072959219884359680000000e39);
}



double
AssLegFunction::P_17_14(const double x) const
{
  double t6;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.200360130e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
          (0.44732242278727900128427062067200000000e38 * t2 * x -
           0.4473224227872790012842706206720000000e37 * t3 * x) /
          0.24429768583155930523471632830300160000000e41);
}



double
AssLegFunction::P_17_14_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.200360130e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  return (-t1 * t5 * t4 *
            (0.44732242278727900128427062067200000000e38 * t2 * x -
             0.4473224227872790012842706206720000000e37 * t3 * x) *
            x / 0.1744983470225423608819402345021440000000e40 +
          t1 * t5 * t4 * t3 *
            (0.147616399519802070423809304821760000000e39 * t2 -
             0.4473224227872790012842706206720000000e37) /
            0.24429768583155930523471632830300160000000e41);
}



double
AssLegFunction::P_17_15(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t6;
  double t8;
  double t4;
  t1 = sqrt(0.33393355e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  t8 = sqrt(t3);
  return (-t1 * t8 * t6 * t4 * t3 *
          (0.147616399519802070423809304821760000000e39 * t2 -
           0.4473224227872790012842706206720000000e37) /
          0.97719074332623722093886531321200640000000e41);
}



double
AssLegFunction::P_17_15_Deriv(const double x) const
{
  double t2;
  double t3;
  double t1;
  double t4;
  double t5;
  double t7;
  t1 = sqrt(0.33393355e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (t1 * t7 * t5 * t4 *
            (0.147616399519802070423809304821760000000e39 * t2 -
             0.4473224227872790012842706206720000000e37) *
            x / 0.6514604955508248139592435421413376000000e40 -
          0.99e2 / 0.32768e5 * t1 * t7 * t5 * t4 * t3 * x);
}



double
AssLegFunction::P_17_16(const double x) const
{
  double t6;
  double t5;
  double t1;
  double t2;
  double t4;
  t1 = sqrt(0.2203961430e10);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (0.3e1 / 0.65536e5 * t1 * t6 * x);
}



double
AssLegFunction::P_17_16_Deriv(const double x) const
{
  double t3;
  double t4;
  double t1;
  double t6;
  double t11;
  double t2;
  t1  = sqrt(0.2203961430e10);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t11 = t6 * t6;
  return (-0.3e1 / 0.4096e4 * t1 * t6 * t4 * t3 * t2 +
          0.3e1 / 0.65536e5 * t1 * t11);
}



double
AssLegFunction::P_17_17(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  t1 = sqrt(0.64822395e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (-0.3e1 / 0.65536e5 * t1 * t7 * t6);
}



double
AssLegFunction::P_17_17_Deriv(const double x) const
{
  double t4;
  double t6;
  double t1;
  double t2;
  double t3;
  double t8;
  t1 = sqrt(0.64822395e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  t8 = sqrt(t3);
  return (0.51e2 / 0.65536e5 * t1 * t8 * t6 * t4 * t3 * x);
}

double
AssLegFunction::P_18_0(const double x) const
{
  double t4;
  double t6;
  double t9;
  double t31;
  double t2;
  double t10;
  double t14;
  double t3;
  double t1;
  double t18;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t6  = t1 - 0.1e1;
  t9  = t6 * t6;
  t10 = t2 * t1;
  t14 = t6 * t9;
  t18 = t9 * t9;
  t31 = t18 * t18;
  return (t4 * t1 + 0.153e3 / 0.2e1 * t6 * t4 +
          0.2295e4 / 0.2e1 * t9 * t3 * t10 + 0.23205e5 / 0.4e1 * t14 * t3 * t2 +
          0.765765e6 / 0.64e2 * t18 * t3 * t1 +
          0.1378377e7 / 0.128e3 * t18 * t6 * t3 +
          0.1072071e7 / 0.256e3 * t18 * t9 * t10 +
          0.328185e6 / 0.512e3 * t18 * t14 * t2 +
          0.984555e6 / 0.32768e5 * t31 * t1 + 0.12155e5 / 0.65536e5 * t31 * t6);
}

double
AssLegFunction::P_18_0_Deriv(const double x) const
{
  double t22;
  double t1;
  double t2;
  double t7;
  double t18;
  double t3;
  double t4;
  double t8;
  double t35;
  double t13;
  double t14;
  double t9;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t7  = t1 - 0.1e1;
  t8  = t1 * x;
  t9  = t2 * t8;
  t13 = t7 * t7;
  t14 = t2 * x;
  t18 = t13 * t7;
  t22 = t13 * t13;
  t35 = t22 * t22;
  return (0.171e3 * t4 * x + 0.5814e4 * t7 * t3 * t9 +
          0.101745e6 / 0.2e1 * t13 * t3 * t14 +
          0.1322685e7 / 0.8e1 * t18 * t3 * t8 +
          0.14549535e8 / 0.64e2 * t22 * t3 * x +
          0.8729721e7 / 0.64e2 * t22 * t7 * t9 +
          0.8729721e7 / 0.256e3 * t22 * t13 * t14 +
          0.6235515e7 / 0.2048e4 * t22 * t18 * t8 +
          0.2078505e7 / 0.32768e5 * t35 * x);
}



double
AssLegFunction::P_18_1(const double x) const
{
  double t7;
  double t1;
  double t8;
  double t22;
  double t2;
  double t39;
  double t6;
  double t3;
  double t4;
  double t26;
  double t13;
  double t17;
  double t11;
  double t18;
  double t12;
  t1  = sqrt(0.38e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t13 = t6 * t12;
  t17 = t11 * t11;
  t18 = t6 * x;
  t22 = t17 * t11;
  t26 = t17 * t17;
  t39 = t26 * t26;
  return (-t1 * t4 *
          (0.286996798814155702272000e24 * t8 * x +
           0.9757891159681293877248000e25 * t11 * t7 * t13 +
           0.85381547647211321425920000e26 * t17 * t7 * t18 +
           0.277490029853436794634240000e27 * t22 * t7 * t12 +
           0.381548791048475592622080000e27 * t26 * t7 * x +
           0.228929274629085355573248000e27 * t26 * t11 * t13 +
           0.57232318657271338893312000e26 * t26 * t17 * t18 +
           0.5110028451542083829760000e25 * t26 * t22 * t12 +
           0.106458926073793413120000e24 * t39 * x) /
          0.191331199209437134848000e24);
}



double
AssLegFunction::P_18_1_Deriv(const double x) const
{
  double t34;
  double t13;
  double t37;
  double t23;
  double t14;
  double t7;
  double t9;
  double t8;
  double t40;
  double t27;
  double t48;
  double t1;
  double t18;
  double t3;
  double t2;
  double t19;
  double t31;
  double t4;
  double t12;
  t1  = sqrt(0.38e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t12 = -t3;
  t13 = t2 * x;
  t14 = t7 * t13;
  t18 = t12 * t12;
  t19 = t7 * x;
  t23 = t18 * t12;
  t27 = t18 * t18;
  t31 = t27 * t12;
  t34 = t27 * t18;
  t37 = t27 * t23;
  t40 = t27 * t27;
  t48 = t7 * t2;
  return (t1 / t4 *
            (0.286996798814155702272000e24 * t9 * x +
             0.9757891159681293877248000e25 * t12 * t8 * t14 +
             0.85381547647211321425920000e26 * t18 * t8 * t19 +
             0.277490029853436794634240000e27 * t23 * t8 * t13 +
             0.381548791048475592622080000e27 * t27 * t8 * x +
             0.228929274629085355573248000e27 * t31 * t14 +
             0.57232318657271338893312000e26 * t34 * t19 +
             0.5110028451542083829760000e25 * t13 * t37 +
             0.106458926073793413120000e24 * t40 * x) *
            x / 0.191331199209437134848000e24 -
          t1 * t4 *
            (0.24394727899203234693120000e26 * t9 +
             0.487894557984064693862400000e27 * t12 * t8 * t48 +
             0.2774900298534367946342400000e28 * t18 * t8 * t7 +
             0.6104780656775609481953280000e28 * t23 * t8 * t2 +
             0.5723231865727133889331200000e28 * t27 * t8 +
             0.2289292746290853555732480000e28 * t31 * t48 +
             0.357701991607945868083200000e27 * t34 * t7 +
             0.17033428171806946099200000e26 * t37 * t2 +
             0.106458926073793413120000e24 * t40) /
            0.191331199209437134848000e24);
}



double
AssLegFunction::P_18_2(const double x) const
{
  double t7;
  double t14;
  double t9;
  double t22;
  double t5;
  double t1;
  double t34;
  double t2;
  double t3;
  double t6;
  double t18;
  double t10;
  t1  = sqrt(0.3230e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t7  = t6 * t6;
  t9  = -t3;
  t10 = t5 * t2;
  t14 = t9 * t9;
  t18 = t14 * t9;
  t22 = t14 * t14;
  t34 = t22 * t22;
  return (t1 * t3 *
          (0.24394727899203234693120000e26 * t7 +
           0.487894557984064693862400000e27 * t9 * t6 * t10 +
           0.2774900298534367946342400000e28 * t14 * t6 * t5 +
           0.6104780656775609481953280000e28 * t18 * t6 * t2 +
           0.5723231865727133889331200000e28 * t22 * t6 +
           0.2289292746290853555732480000e28 * t22 * t9 * t10 +
           0.357701991607945868083200000e27 * t22 * t14 * t5 +
           0.17033428171806946099200000e26 * t22 * t18 * t2 +
           0.106458926073793413120000e24 * t34) /
          0.32526303865604312924160000e26);
}



double
AssLegFunction::P_18_2_Deriv(const double x) const
{
  double t24;
  double t27;
  double t13;
  double t41;
  double t17;
  double t21;
  double t30;
  double t3;
  double t8;
  double t9;
  double t44;
  double t4;
  double t6;
  double t5;
  double t33;
  double t40;
  double t1;
  t1  = sqrt(0.3230e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = t5 * t5;
  t8  = t3 - 0.1e1;
  t9  = t4 * t3;
  t13 = t8 * t8;
  t17 = t13 * t8;
  t21 = t13 * t13;
  t24 = t21 * t8;
  t27 = t21 * t13;
  t30 = t21 * t17;
  t33 = t21 * t21;
  t40 = t3 * x;
  t41 = t4 * t40;
  t44 = t4 * x;
  return (-t1 * x *
            (0.24394727899203234693120000e26 * t6 +
             0.487894557984064693862400000e27 * t8 * t5 * t9 +
             0.2774900298534367946342400000e28 * t13 * t5 * t4 +
             0.6104780656775609481953280000e28 * t17 * t5 * t3 +
             0.5723231865727133889331200000e28 * t21 * t5 +
             0.2289292746290853555732480000e28 * t24 * t9 +
             0.357701991607945868083200000e27 * t27 * t4 +
             0.17033428171806946099200000e26 * t30 * t3 +
             0.106458926073793413120000e24 * t33) /
            0.16263151932802156462080000e26 -
          t1 * t8 *
            (0.1366104762355381142814720000e28 * t5 * t41 +
             0.17930125005914377499443200000e29 * t8 * t5 * t44 +
             0.69927487523066072247828480000e29 * t13 * t5 * t40 +
             0.106833661493573165934182400000e30 * t17 * t5 * x +
             0.68678782388725606671974400000e29 * t21 * t41 +
             0.18028180377040471751393280000e29 * t24 * t44 +
             0.1669275960837080717721600000e28 * t27 * t40 +
             0.35770199160794586808320000e26 * t30 * x) /
            0.32526303865604312924160000e26);
}



double
AssLegFunction::P_18_3(const double x) const
{
  double t13;
  double t8;
  double t2;
  double t14;
  double t4;
  double t9;
  double t7;
  double t3;
  double t26;
  double t1;
  double t18;
  double t22;
  double t10;
  t1  = sqrt(0.67830e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t7;
  t10 = t8 * t8;
  t13 = -t3;
  t14 = t8 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t26 = t18 * t18;
  return (-t1 * t4 * t3 *
          (0.1366104762355381142814720000e28 * t10 * t9 +
           0.17930125005914377499443200000e29 * t13 * t10 * t14 +
           0.69927487523066072247828480000e29 * t18 * t10 * t7 +
           0.106833661493573165934182400000e30 * t22 * t10 * x +
           0.68678782388725606671974400000e29 * t26 * t9 +
           0.18028180377040471751393280000e29 * t26 * t13 * t14 +
           0.1669275960837080717721600000e28 * t26 * t18 * t7 +
           0.35770199160794586808320000e26 * t26 * t22 * x) /
          0.2732209524710762285629440000e28);
}



double
AssLegFunction::P_18_3_Deriv(const double x) const
{
  double t17;
  double t21;
  double t2;
  double t3;
  double t4;
  double t9;
  double t25;
  double t12;
  double t13;
  double t43;
  double t28;
  double t7;
  double t6;
  double t8;
  double t1;
  double t34;
  double t31;
  t1  = sqrt(0.67830e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t6;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t7 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t25 = t17 * t17;
  t28 = t25 * t12;
  t31 = t25 * t17;
  t34 = t21 * t25;
  t43 = t7 * t2;
  return (t1 * t4 *
            (0.1366104762355381142814720000e28 * t9 * t8 +
             0.17930125005914377499443200000e29 * t12 * t9 * t13 +
             0.69927487523066072247828480000e29 * t17 * t9 * t6 +
             0.106833661493573165934182400000e30 * t21 * t9 * x +
             0.68678782388725606671974400000e29 * t25 * t8 +
             0.18028180377040471751393280000e29 * t28 * t13 +
             0.1669275960837080717721600000e28 * t31 * t6 +
             0.35770199160794586808320000e26 * t34 * x) *
            x / 0.910736508236920761876480000e27 -
          t1 * t4 * t3 *
            (0.56351821447159472141107200000e29 * t9 * t43 +
             0.512801575169151196484075520000e30 * t12 * t9 * t7 +
             0.1410204331715165790331207680000e31 * t17 * t9 * t2 +
             0.1510933212551963346783436800000e31 * t21 * t9 +
             0.661033280491483964217753600000e30 * t25 * t43 +
             0.110172213415247327369625600000e30 * t28 * t7 +
             0.5508610670762366368481280000e28 * t31 * t2 +
             0.35770199160794586808320000e26 * t34) /
            0.2732209524710762285629440000e28);
}



double
AssLegFunction::P_18_4(const double x) const
{
  double t6;
  double t11;
  double t22;
  double t19;
  double t8;
  double t7;
  double t3;
  double t1;
  double t2;
  double t15;
  double t4;
  t1  = sqrt(0.24871e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t2;
  t8  = t6 * t6;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t22 = t15 * t15;
  return (t1 * t4 *
          (0.56351821447159472141107200000e29 * t8 * t7 +
           0.512801575169151196484075520000e30 * t11 * t8 * t6 +
           0.1410204331715165790331207680000e31 * t15 * t8 * t2 +
           0.1510933212551963346783436800000e31 * t19 * t8 +
           0.661033280491483964217753600000e30 * t22 * t7 +
           0.110172213415247327369625600000e30 * t22 * t11 * t6 +
           0.5508610670762366368481280000e28 * t22 * t15 * t2 +
           0.35770199160794586808320000e26 * t22 * t19) /
          0.30054304771818385141923840000e29);
}



double
AssLegFunction::P_18_4_Deriv(const double x) const
{
  double t6;
  double t10;
  double t27;
  double t36;
  double t24;
  double t38;
  double t7;
  double t1;
  double t14;
  double t2;
  double t3;
  double t21;
  double t5;
  double t41;
  double t18;
  t1  = sqrt(0.24871e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t2;
  t7  = t5 * t5;
  t10 = -t3;
  t14 = t10 * t10;
  t18 = t14 * t10;
  t21 = t14 * t14;
  t24 = t21 * t10;
  t27 = t21 * t14;
  t36 = t3 * t3;
  t38 = t5 * x;
  t41 = t2 * x;
  return (-t1 * t3 *
            (0.56351821447159472141107200000e29 * t7 * t6 +
             0.512801575169151196484075520000e30 * t10 * t7 * t5 +
             0.1410204331715165790331207680000e31 * t14 * t7 * t2 +
             0.1510933212551963346783436800000e31 * t18 * t7 +
             0.661033280491483964217753600000e30 * t21 * t6 +
             0.110172213415247327369625600000e30 * t24 * t5 +
             0.5508610670762366368481280000e28 * t27 * t2 +
             0.35770199160794586808320000e26 * t21 * t18) *
            x / 0.7513576192954596285480960000e28 +
          t1 * t36 *
            (0.1814528650598535002943651840000e31 * t7 * t38 +
             0.11794436228890477519133736960000e32 * t10 * t7 * t41 +
             0.23167642592463437984012697600000e32 * t14 * t7 * x +
             0.17375731944347578488009523200000e32 * t18 * t5 * t41 +
             0.5067921817101377059002777600000e31 * t21 * t38 +
             0.506792181710137705900277760000e30 * t24 * t41 +
             0.11518004129775856952279040000e29 * t27 * x) /
            0.30054304771818385141923840000e29);
}



double
AssLegFunction::P_18_5(const double x) const
{
  double t18;
  double t1;
  double t13;
  double t14;
  double t2;
  double t3;
  double t26;
  double t9;
  double t8;
  double t4;
  double t5;
  double t10;
  t1  = sqrt(0.163438e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * x;
  t10 = t8 * t8;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t26 = t18 * t18;
  return (-t1 * t5 * t4 *
          (0.1814528650598535002943651840000e31 * t10 * t9 +
           0.11794436228890477519133736960000e32 * t13 * t10 * t14 +
           0.23167642592463437984012697600000e32 * t18 * t10 * x +
           0.17375731944347578488009523200000e32 * t18 * t13 * t8 * t14 +
           0.5067921817101377059002777600000e31 * t26 * t9 +
           0.506792181710137705900277760000e30 * t26 * t13 * t14 +
           0.11518004129775856952279040000e29 * t26 * t18 * x) /
          0.1382498019503645716528496640000e31);
}



double
AssLegFunction::P_18_5_Deriv(const double x) const
{
  double t28;
  double t21;
  double t31;
  double t7;
  double t1;
  double t12;
  double t8;
  double t2;
  double t17;
  double t38;
  double t3;
  double t9;
  double t4;
  double t13;
  double t25;
  t1  = sqrt(0.163438e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * x;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t25 = t17 * t17;
  t28 = t25 * t12;
  t31 = t25 * t17;
  t38 = t3 * t3;
  return (t1 * t4 * t3 *
            (0.1814528650598535002943651840000e31 * t9 * t8 +
             0.11794436228890477519133736960000e32 * t12 * t9 * t13 +
             0.23167642592463437984012697600000e32 * t17 * t9 * x +
             0.17375731944347578488009523200000e32 * t21 * t7 * t13 +
             0.5067921817101377059002777600000e31 * t25 * t8 +
             0.506792181710137705900277760000e30 * t28 * t13 +
             0.11518004129775856952279040000e29 * t31 * x) *
            x / 0.276499603900729143305699328000e30 -
          t1 * t4 * t38 *
            (0.47177744915561910076534947840000e32 * t9 * t7 +
             0.222409368887649004646521896960000e33 * t12 * t9 * t2 +
             0.312763174998256412784171417600000e33 * t17 * t9 +
             0.162173498147244065888088883200000e33 * t21 * t7 * t2 +
             0.30407530902608262354016665600000e32 * t25 * t7 +
             0.1658592594687723401128181760000e31 * t28 * t2 +
             0.11518004129775856952279040000e29 * t31) /
            0.1382498019503645716528496640000e31);
}



double
AssLegFunction::P_18_6(const double x) const
{
  double t22;
  double t1;
  double t11;
  double t4;
  double t15;
  double t3;
  double t2;
  double t8;
  double t7;
  t1  = sqrt(0.3187041e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t22 = t15 * t15;
  return (t1 * t4 * t3 *
          (0.47177744915561910076534947840000e32 * t8 * t7 +
           0.222409368887649004646521896960000e33 * t11 * t8 * t2 +
           0.312763174998256412784171417600000e33 * t15 * t8 +
           0.162173498147244065888088883200000e33 * t15 * t11 * t7 * t2 +
           0.30407530902608262354016665600000e32 * t22 * t7 +
           0.1658592594687723401128181760000e31 * t22 * t11 * t2 +
           0.11518004129775856952279040000e29 * t22 * t15) /
          0.107834845521284365889222737920000e33);
}



double
AssLegFunction::P_18_6_Deriv(const double x) const
{
  double t21;
  double t2;
  double t4;
  double t3;
  double t10;
  double t24;
  double t35;
  double t1;
  double t14;
  double t6;
  double t17;
  double t7;
  t1  = sqrt(0.3187041e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t21 = t14 * t14;
  t24 = t21 * t10;
  t35 = t2 * x;
  return (-t1 * t4 *
            (0.47177744915561910076534947840000e32 * t7 * t6 +
             0.222409368887649004646521896960000e33 * t10 * t7 * t2 +
             0.312763174998256412784171417600000e33 * t14 * t7 +
             0.162173498147244065888088883200000e33 * t17 * t6 * t2 +
             0.30407530902608262354016665600000e32 * t21 * t6 +
             0.1658592594687723401128181760000e31 * t24 * t2 +
             0.11518004129775856952279040000e29 * t21 * t14) *
            x / 0.17972474253547394314870456320000e32 +
          t1 * t4 * t3 *
            (0.1010951676762040930211463168000000e34 * t7 * t35 +
             0.3475146388869515697601904640000000e34 * t10 * t7 * x +
             0.3475146388869515697601904640000000e34 * t14 * t6 * t35 +
             0.1216301236104330494160666624000000e34 * t17 * t6 * x +
             0.138216049557310283427348480000000e33 * t21 * t35 +
             0.3455401238932757085683712000000e31 * t24 * x) /
            0.107834845521284365889222737920000e33);
}



double
AssLegFunction::P_18_7(const double x) const
{
  double t2;
  double t3;
  double t26;
  double t14;
  double t4;
  double t10;
  double t11;
  double t1;
  double t9;
  double t6;
  double t18;
  t1  = sqrt(0.1062347e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t11 = t10 * t10;
  t14 = -t3;
  t18 = t14 * t14;
  t26 = t18 * t18;
  return (-t1 * t6 * t4 * t3 *
          (0.1010951676762040930211463168000000e34 * t11 * t9 +
           0.3475146388869515697601904640000000e34 * t14 * t11 * x +
           0.3475146388869515697601904640000000e34 * t18 * t10 * t9 +
           0.1216301236104330494160666624000000e34 * t18 * t14 * t10 * x +
           0.138216049557310283427348480000000e33 * t26 * t9 +
           0.3455401238932757085683712000000e31 * t26 * t14 * x) /
          0.1078348455212843658892227379200000e34);
}



double
AssLegFunction::P_18_7_Deriv(const double x) const
{
  double t1;
  double t10;
  double t3;
  double t2;
  double t25;
  double t13;
  double t4;
  double t5;
  double t17;
  double t28;
  double t21;
  double t8;
  double t9;
  t1  = sqrt(0.1062347e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * x;
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t17 = t13 * t13;
  t21 = t17 * t13;
  t25 = t17 * t17;
  t28 = t25 * t13;
  return (t1 * t5 * t4 *
            (0.1010951676762040930211463168000000e34 * t10 * t8 +
             0.3475146388869515697601904640000000e34 * t13 * t10 * x +
             0.3475146388869515697601904640000000e34 * t17 * t9 * t8 +
             0.1216301236104330494160666624000000e34 * t21 * t9 * x +
             0.138216049557310283427348480000000e33 * t25 * t8 +
             0.3455401238932757085683712000000e31 * t28 * x) *
            x / 0.154049779316120522698889625600000e33 -
          t1 * t5 * t4 * t3 *
            (0.18070761222121481627529904128000000e35 * t10 * t2 +
             0.45176903055303704068824760320000000e35 * t13 * t10 +
             0.31623832138712592848177332224000000e35 * t17 * t9 * t2 +
             0.7187234576980134738222120960000000e34 * t21 * t9 +
             0.449202161061258421138882560000000e33 * t25 * t2 +
             0.3455401238932757085683712000000e31 * t28) /
            0.1078348455212843658892227379200000e34);
}



double
AssLegFunction::P_18_8(const double x) const
{
  double t3;
  double t11;
  double t4;
  double t1;
  double t7;
  double t21;
  double t14;
  double t5;
  double t8;
  double t2;
  t1  = sqrt(0.14858e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t21 = t14 * t14;
  return (t1 * t5 *
          (0.18070761222121481627529904128000000e35 * t8 * t2 +
           0.45176903055303704068824760320000000e35 * t11 * t8 +
           0.31623832138712592848177332224000000e35 * t14 * t7 * t2 +
           0.7187234576980134738222120960000000e34 * t14 * t11 * t7 +
           0.449202161061258421138882560000000e33 * t21 * t2 +
           0.3455401238932757085683712000000e31 * t21 * t11) /
          0.2156696910425687317784454758400000e34);
}



double
AssLegFunction::P_18_8_Deriv(const double x) const
{
  double t1;
  double t34;
  double t3;
  double t2;
  double t4;
  double t18;
  double t7;
  double t8;
  double t11;
  double t21;
  double t14;
  double t30;
  t1  = sqrt(0.14858e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t18 = t14 * t11;
  t21 = t14 * t14;
  t30 = t4 * t4;
  t34 = t2 * x;
  return (-t1 * t4 * t3 *
            (0.18070761222121481627529904128000000e35 * t8 * t2 +
             0.45176903055303704068824760320000000e35 * t11 * t8 +
             0.31623832138712592848177332224000000e35 * t14 * t7 * t2 +
             0.7187234576980134738222120960000000e34 * t7 * t18 +
             0.449202161061258421138882560000000e33 * t21 * t2 +
             0.3455401238932757085683712000000e31 * t11 * t21) *
            x / 0.269587113803210914723056844800000e33 +
          t1 * t30 *
            (0.271061418331822224412948561920000000e36 * t8 * x +
             0.487910552997280003943307411456000000e36 * t11 * t7 * t34 +
             0.232866400294156365518396719104000000e36 * t14 * t7 * x +
             0.32342555596410606321999544320000000e35 * t18 * t34 +
             0.932958334511844413134602240000000e33 * t21 * x) /
            0.2156696910425687317784454758400000e34);
}



double
AssLegFunction::P_18_9(const double x) const
{
  double t2;
  double t1;
  double t3;
  double t4;
  double t5;
  double t6;
  double t9;
  double t10;
  double t13;
  double t14;
  double t25;
  double t18;
  t1  = sqrt(0.111435e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t25 = t18 * t18;
  return (-t1 * t6 * t5 *
          (0.271061418331822224412948561920000000e36 * t10 * x +
           0.487910552997280003943307411456000000e36 * t13 * t9 * t14 +
           0.232866400294156365518396719104000000e36 * t18 * t9 * x +
           0.32342555596410606321999544320000000e35 * t18 * t13 * t14 +
           0.932958334511844413134602240000000e33 * t25 * x) /
          0.97051360969155929300300464128000000e35);
}



double
AssLegFunction::P_18_9_Deriv(const double x) const
{
  double t18;
  double t13;
  double t4;
  double t14;
  double t22;
  double t9;
  double t1;
  double t25;
  double t10;
  double t32;
  double t6;
  double t2;
  double t3;
  t1  = sqrt(0.111435e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t25 = t18 * t18;
  t32 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.271061418331822224412948561920000000e36 * t10 * x +
             0.487910552997280003943307411456000000e36 * t13 * t9 * t14 +
             0.232866400294156365518396719104000000e36 * t18 * t9 * x +
             0.32342555596410606321999544320000000e35 * t22 * t14 +
             0.932958334511844413134602240000000e33 * t25 * x) *
            x / 0.10783484552128436588922273792000000e35 -
          t1 * t6 * t32 *
            (0.3415373870980960027603151880192000000e37 * t10 +
             0.4346839472157585489676738756608000000e37 * t13 * t9 * t2 +
             0.1358387335049245465523980861440000000e37 * t18 * t9 +
             0.104491333465326574271075450880000000e36 * t22 * t2 +
             0.932958334511844413134602240000000e33 * t25) /
            0.97051360969155929300300464128000000e35);
}



double
AssLegFunction::P_18_10(const double x) const
{
  double t2;
  double t3;
  double t11;
  double t4;
  double t5;
  double t9;
  double t21;
  double t15;
  double t1;
  double t8;
  t1  = sqrt(0.780045e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t11 = -t3;
  t15 = t11 * t11;
  t21 = t15 * t15;
  return (t1 * t5 * t3 *
          (0.3415373870980960027603151880192000000e37 * t9 +
           0.4346839472157585489676738756608000000e37 * t11 * t8 * t2 +
           0.1358387335049245465523980861440000000e37 * t15 * t8 +
           0.104491333465326574271075450880000000e36 * t15 * t11 * t2 +
           0.932958334511844413134602240000000e33 * t21) /
          0.4076157160704549030612619493376000000e37);
}



double
AssLegFunction::P_18_10_Deriv(const double x) const
{
  double t3;
  double t4;
  double t5;
  double t7;
  double t8;
  double t10;
  double t28;
  double t14;
  double t17;
  double t1;
  double t20;
  double t2;
  t1  = sqrt(0.780045e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t10 = -t3;
  t14 = t10 * t10;
  t17 = t14 * t10;
  t20 = t14 * t14;
  t28 = t2 * x;
  return (-t1 * t5 *
            (0.3415373870980960027603151880192000000e37 * t8 +
             0.4346839472157585489676738756608000000e37 * t10 * t7 * t2 +
             0.1358387335049245465523980861440000000e37 * t14 * t7 +
             0.104491333465326574271075450880000000e36 * t17 * t2 +
             0.932958334511844413134602240000000e33 * t20) *
            x / 0.407615716070454903061261949337600000e36 +
          t1 * t5 * t3 *
            (0.36016669912162851200178692554752000000e38 * t7 * t28 +
             0.31514586173142494800156355985408000000e38 * t10 * t7 * x +
             0.6060497340988941307722376151040000000e37 * t14 * t28 +
             0.216446333606747903847227719680000000e36 * t17 * x) /
            0.4076157160704549030612619493376000000e37);
}



double
AssLegFunction::P_18_11(const double x) const
{
  double t11;
  double t14;
  double t1;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t10;
  t1  = sqrt(0.45242610e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * x;
  t11 = t2 * t2;
  t14 = -t3;
  t18 = t14 * t14;
  return (-t1 * t7 * t5 * t3 *
          (0.36016669912162851200178692554752000000e38 * t11 * t10 +
           0.31514586173142494800156355985408000000e38 * t14 * t11 * x +
           0.6060497340988941307722376151040000000e37 * t18 * t10 +
           0.216446333606747903847227719680000000e36 * t18 * t14 * x) /
          0.472834230641727687551063861231616000000e39);
}



double
AssLegFunction::P_18_11_Deriv(const double x) const
{
  double t6;
  double t13;
  double t1;
  double t17;
  double t9;
  double t20;
  double t2;
  double t10;
  double t5;
  double t3;
  double t4;
  t1  = sqrt(0.45242610e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  t20 = t17 * t13;
  return (t1 * t6 * t5 *
            (0.36016669912162851200178692554752000000e38 * t10 * t9 +
             0.31514586173142494800156355985408000000e38 * t13 * t10 * x +
             0.6060497340988941307722376151040000000e37 * t17 * t9 +
             0.216446333606747903847227719680000000e36 * t20 * x) *
            x / 0.42984930058338880686460351021056000000e38 -
          t1 * t6 * t5 * t3 *
            (0.315145861731424948001563559854080000000e39 * t10 * t2 +
             0.181814920229668239231671284531200000000e39 * t13 * t10 +
             0.19480170024607311346250494771200000000e38 * t17 * t2 +
             0.216446333606747903847227719680000000e36 * t20) /
            0.472834230641727687551063861231616000000e39);
}



double
AssLegFunction::P_18_12(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t8;
  double t11;
  double t14;
  double t1;
  t1  = sqrt(0.215441e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (t1 * t5 * t4 *
          (0.315145861731424948001563559854080000000e39 * t8 * t2 +
           0.181814920229668239231671284531200000000e39 * t11 * t8 +
           0.19480170024607311346250494771200000000e38 * t14 * t2 +
           0.216446333606747903847227719680000000e36 * t14 * t11) /
          0.472834230641727687551063861231616000000e39);
}



double
AssLegFunction::P_18_12_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t14;
  double t5;
  double t8;
  double t11;
  t1  = sqrt(0.215441e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (-t1 * t5 * t3 *
            (0.315145861731424948001563559854080000000e39 * t8 * t2 +
             0.181814920229668239231671284531200000000e39 * t11 * t8 +
             0.19480170024607311346250494771200000000e38 * t14 * t2 +
             0.216446333606747903847227719680000000e36 * t14 * t11) *
            x / 0.39402852553477307295921988435968000000e38 +
          t1 * t5 * t4 *
            (0.2254505010847886166472723928186880000000e40 * t8 * x +
             0.805180361017102202311687117209600000000e39 * t11 * t2 * x +
             0.40259018050855110115584355860480000000e38 * t14 * x) /
            0.472834230641727687551063861231616000000e39);
}



double
AssLegFunction::P_18_13(const double x) const
{
  double t10;
  double t4;
  double t2;
  double t3;
  double t5;
  double t17;
  double t13;
  double t1;
  double t7;
  t1  = sqrt(0.40072026e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (-t1 * t7 * t5 * t4 *
          (0.2254505010847886166472723928186880000000e40 * t10 * x +
           0.805180361017102202311687117209600000000e39 * t13 * t2 * x +
           0.40259018050855110115584355860480000000e38 * t17 * x) /
          0.87947166899361349884497878189080576000000e41);
}



double
AssLegFunction::P_18_13_Deriv(const double x) const
{
  double t17;
  double t2;
  double t4;
  double t3;
  double t10;
  double t5;
  double t7;
  double t1;
  double t13;
  t1  = sqrt(0.40072026e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (t1 * t7 * t5 * t3 *
            (0.2254505010847886166472723928186880000000e40 * t10 * x +
             0.805180361017102202311687117209600000000e39 * t13 * t2 * x +
             0.40259018050855110115584355860480000000e38 * t17 * x) *
            x / 0.6765166684566257683422913706852352000000e40 -
          t1 * t7 * t5 * t4 *
            (0.12882885776273635236986993875353600000000e41 * t10 +
             0.2576577155254727047397398775070720000000e40 * t13 * t2 +
             0.40259018050855110115584355860480000000e38 * t17) /
            0.87947166899361349884497878189080576000000e41);
}



double
AssLegFunction::P_18_14(const double x) const
{
  double t9;
  double t11;
  double t14;
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  t1  = sqrt(0.100180065e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (t1 * t6 * t4 * t3 *
          (0.12882885776273635236986993875353600000000e41 * t9 +
           0.2576577155254727047397398775070720000000e40 * t11 * t2 +
           0.40259018050855110115584355860480000000e38 * t14) /
          0.1758943337987226997689957563781611520000000e43);
}



double
AssLegFunction::P_18_14_Deriv(const double x) const
{
  double t10;
  double t13;
  double t1;
  double t2;
  double t4;
  double t3;
  double t5;
  double t8;
  t1  = sqrt(0.100180065e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (-t1 * t5 * t4 *
            (0.12882885776273635236986993875353600000000e41 * t8 +
             0.2576577155254727047397398775070720000000e40 * t10 * t2 +
             0.40259018050855110115584355860480000000e38 * t13) *
            x / 0.125638809856230499834996968841543680000000e42 +
          t1 * t5 * t4 * t3 *
            (0.56684697415603995042742773051555840000000e41 * t2 * x +
             0.5314190382712874535257134973583360000000e40 * t10 * x) /
            0.1758943337987226997689957563781611520000000e43);
}



double
AssLegFunction::P_18_15(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  t1 = sqrt(0.367326905e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t6 = t4 * t4;
  t8 = sqrt(t3);
  return (-t1 * t8 * t6 * t4 * t3 *
          (0.56684697415603995042742773051555840000000e41 * t2 * x -
           0.5314190382712874535257134973583360000000e40 * t3 * x) /
          0.38696753435718993949179066403195453440000000e44);
}



double
AssLegFunction::P_18_15_Deriv(const double x) const
{
  double t7;
  double t2;
  double t3;
  double t5;
  double t4;
  double t1;
  t1 = sqrt(0.367326905e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t7 = sqrt(t3);
  return (t1 * t7 * t5 * t4 *
            (0.56684697415603995042742773051555840000000e41 * t2 * x -
             0.5314190382712874535257134973583360000000e40 * t3 * x) *
            x / 0.2579783562381266263278604426879696896000000e43 -
          t1 * t7 * t5 * t4 * t3 *
            (0.185996663394950608733999724075417600000000e42 * t2 -
             0.5314190382712874535257134973583360000000e40) /
            0.38696753435718993949179066403195453440000000e44);
}



double
AssLegFunction::P_18_16(const double x) const
{
  double t1;
  double t2;
  double t5;
  double t4;
  double t6;
  t1 = sqrt(0.129644790e9);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (t1 * t6 *
          (0.185996663394950608733999724075417600000000e42 * t2 -
           0.5314190382712874535257134973583360000000e40) /
          0.232180520614313963695074398419172720640000000e45);
}



double
AssLegFunction::P_18_16_Deriv(const double x) const
{
  double t6;
  double t14;
  double t1;
  double t3;
  double t2;
  double t4;
  t1  = sqrt(0.129644790e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t14 = t6 * t6;
  return (-t1 * t6 * t4 * t3 *
            (0.185996663394950608733999724075417600000000e42 * t2 -
             0.5314190382712874535257134973583360000000e40) *
            x / 0.14511282538394622730942149901198295040000000e44 +
          0.105e3 / 0.65536e5 * t1 * t14 * x);
}



double
AssLegFunction::P_18_17(const double x) const
{
  double t4;
  double t5;
  double t6;
  double t1;
  double t2;
  double t7;
  double t3;
  t1 = sqrt(0.90751353e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (-0.15e2 / 0.65536e5 * t1 * t7 * t6 * x);
}



double
AssLegFunction::P_18_17_Deriv(const double x) const
{
  double t3;
  double t4;
  double t1;
  double t6;
  double t13;
  double t8;
  double t2;
  t1  = sqrt(0.90751353e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t8  = sqrt(t3);
  t13 = t6 * t6;
  return (0.255e3 / 0.65536e5 * t1 * t8 * t6 * t4 * t3 * t2 -
          0.15e2 / 0.65536e5 * t1 * t8 * t13);
}



double
AssLegFunction::P_18_18(const double x) const
{
  double t5;
  double t6;
  double t2;
  double t1;
  double t4;
  double t3;
  t1 = sqrt(0.90751353e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (0.5e1 / 0.131072e6 * t1 * t6 * t3);
}



double
AssLegFunction::P_18_18_Deriv(const double x) const
{
  double t4;
  double t6;
  double t5;
  double t2;
  double t1;
  t1 = sqrt(0.90751353e8);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (-0.45e2 / 0.65536e5 * t1 * t6 * x);
}

double
AssLegFunction::P_19_0(const double x) const
{
  double t11;
  double t16;
  double t1;
  double t17;
  double t2;
  double t12;
  double t3;
  double t4;
  double t5;
  double t21;
  double t7;
  double t35;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t1 - 0.1e1;
  t11 = t7 * t7;
  t12 = t3 * t2;
  t16 = t11 * t7;
  t17 = t3 * x;
  t21 = t11 * t11;
  t35 = t21 * t21;
  return (
    t5 * t2 + 0.171e3 / 0.2e1 * t7 * t5 * x +
    0.2907e4 / 0.2e1 * t11 * t4 * t12 + 0.33915e5 / 0.4e1 * t16 * t4 * t17 +
    0.1322685e7 / 0.64e2 * t21 * t4 * t2 +
    0.2909907e7 / 0.128e3 * t21 * t7 * t4 * x +
    0.2909907e7 / 0.256e3 * t21 * t11 * t12 +
    0.1247103e7 / 0.512e3 * t21 * t16 * t17 +
    0.6235515e7 / 0.32768e5 * t35 * t2 + 0.230945e6 / 0.65536e5 * t35 * t7 * x);
}

double
AssLegFunction::P_19_0_Deriv(const double x) const
{
  double t15;
  double t19;
  double t10;
  double t1;
  double t32;
  double t11;
  double t2;
  double t3;
  double t4;
  double t7;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t7  = t1 - 0.1e1;
  t10 = t7 * t7;
  t11 = t1 * t2;
  t15 = t10 * t7;
  t19 = t10 * t10;
  t32 = t19 * t19;
  return (0.190e3 * t4 * t1 + 0.14535e5 / 0.2e1 * t7 * t4 +
          0.72675e5 * t10 * t3 * t11 + 0.2204475e7 / 0.8e1 * t15 * t3 * t2 +
          0.14549535e8 / 0.32e2 * t19 * t3 * t1 +
          0.43648605e8 / 0.128e3 * t19 * t7 * t3 +
          0.14549535e8 / 0.128e3 * t19 * t10 * t11 +
          0.31177575e8 / 0.2048e4 * t19 * t15 * t2 +
          0.10392525e8 / 0.16384e5 * t32 * t1 +
          0.230945e6 / 0.65536e5 * t32 * t7);
}



double
AssLegFunction::P_19_1(const double x) const
{
  double t11;
  double t14;
  double t15;
  double t2;
  double t6;
  double t3;
  double t19;
  double t7;
  double t4;
  double t8;
  double t23;
  double t36;
  double t1;
  t1  = sqrt(0.95e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t15 = t6 * t2;
  t19 = t14 * t11;
  t23 = t14 * t14;
  t36 = t23 * t23;
  return (-t1 * t4 *
          (0.12117642616597685207040000e26 * t8 * t2 +
           0.463499830084861459169280000e27 * t11 * t8 +
           0.4634998300848614591692800000e28 * t14 * t7 * t15 +
           0.17574368557384330326835200000e29 * t19 * t7 * t6 +
           0.28997708119684145039278080000e29 * t23 * t7 * t2 +
           0.21748281089763108779458560000e29 * t23 * t11 * t7 +
           0.7249427029921036259819520000e28 * t23 * t14 * t15 +
           0.970905405792995927654400000e27 * t23 * t19 * t6 +
           0.40454391908041496985600000e26 * t36 * t2 +
           0.224746621711341649920000e24 * t36 * t11) /
          0.12117642616597685207040000e26);
}



double
AssLegFunction::P_19_1_Deriv(const double x) const
{
  double t1;
  double t20;
  double t28;
  double t8;
  double t7;
  double t9;
  double t2;
  double t3;
  double t31;
  double t49;
  double t48;
  double t34;
  double t4;
  double t12;
  double t37;
  double t24;
  double t15;
  double t16;
  double t53;
  t1  = sqrt(0.95e2);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t12 = -t3;
  t15 = t12 * t12;
  t16 = t7 * t2;
  t20 = t15 * t12;
  t24 = t15 * t15;
  t28 = t24 * t12;
  t31 = t24 * t15;
  t34 = t24 * t20;
  t37 = t24 * t24;
  t48 = t2 * x;
  t49 = t7 * t48;
  t53 = t7 * x;
  return (t1 / t4 *
            (0.12117642616597685207040000e26 * t2 * t9 +
             0.463499830084861459169280000e27 * t12 * t9 +
             0.4634998300848614591692800000e28 * t15 * t8 * t16 +
             0.17574368557384330326835200000e29 * t20 * t8 * t7 +
             0.28997708119684145039278080000e29 * t24 * t8 * t2 +
             0.21748281089763108779458560000e29 * t28 * t8 +
             0.7249427029921036259819520000e28 * t31 * t16 +
             0.970905405792995927654400000e27 * t34 * t7 +
             0.40454391908041496985600000e26 * t37 * t2 +
             0.224746621711341649920000e24 * t37 * t12) *
            x / 0.12117642616597685207040000e26 -
          t1 * t4 *
            (0.1145117227268481252065280000e28 * t9 * x +
             0.25955990484752241713479680000e29 * t12 * t8 * t49 +
             0.170336187556186586244710400000e30 * t15 * t8 * t53 +
             0.442874087646085124236247040000e30 * t20 * t8 * t48 +
             0.507459892094472538187366400000e30 * t24 * t8 * x +
             0.260979373077157305353502720000e30 * t28 * t49 +
             0.57089237860628160546078720000e29 * t31 * t53 +
             0.4530891893700647662387200000e28 * t48 * t34 +
             0.84954223006887143669760000e26 * t37 * x) /
            0.12117642616597685207040000e26);
}



double
AssLegFunction::P_19_2(const double x) const
{
  double t16;
  double t10;
  double t1;
  double t11;
  double t6;
  double t2;
  double t38;
  double t25;
  double t5;
  double t12;
  double t3;
  double t21;
  double t17;
  double t7;
  t1  = sqrt(0.3990e4);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t7  = t6 * t6;
  t10 = -t3;
  t11 = t2 * x;
  t12 = t5 * t11;
  t16 = t10 * t10;
  t17 = t5 * x;
  t21 = t16 * t10;
  t25 = t16 * t16;
  t38 = t25 * t25;
  return (t1 * t3 *
          (0.1145117227268481252065280000e28 * t7 * x +
           0.25955990484752241713479680000e29 * t10 * t6 * t12 +
           0.170336187556186586244710400000e30 * t16 * t6 * t17 +
           0.442874087646085124236247040000e30 * t21 * t6 * t11 +
           0.507459892094472538187366400000e30 * t25 * t6 * x +
           0.260979373077157305353502720000e30 * t25 * t10 * t12 +
           0.57089237860628160546078720000e29 * t25 * t16 * t17 +
           0.4530891893700647662387200000e28 * t25 * t21 * t11 +
           0.84954223006887143669760000e26 * t38 * x) /
          0.1526822969691308336087040000e28);
}



double
AssLegFunction::P_19_2_Deriv(const double x) const
{
  double t37;
  double t1;
  double t34;
  double t28;
  double t15;
  double t3;
  double t16;
  double t9;
  double t4;
  double t10;
  double t5;
  double t6;
  double t24;
  double t46;
  double t11;
  double t31;
  double t20;
  t1  = sqrt(0.3990e4);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = t5 * t5;
  t9  = t3 - 0.1e1;
  t10 = t3 * x;
  t11 = t4 * t10;
  t15 = t9 * t9;
  t16 = t4 * x;
  t20 = t15 * t9;
  t24 = t15 * t15;
  t28 = t24 * t9;
  t31 = t24 * t15;
  t34 = t24 * t20;
  t37 = t24 * t24;
  t46 = t4 * t3;
  return (-t1 * x *
            (0.1145117227268481252065280000e28 * t6 * x +
             0.25955990484752241713479680000e29 * t9 * t5 * t11 +
             0.170336187556186586244710400000e30 * t15 * t5 * t16 +
             0.442874087646085124236247040000e30 * t20 * t5 * t10 +
             0.507459892094472538187366400000e30 * t24 * t5 * x +
             0.260979373077157305353502720000e30 * t28 * t11 +
             0.57089237860628160546078720000e29 * t31 * t16 +
             0.4530891893700647662387200000e28 * t34 * t10 +
             0.84954223006887143669760000e26 * t37 * x) /
            0.763411484845654168043520000e27 -
          t1 * t9 *
            (0.71378973833068664712069120000e29 * t6 +
             0.1070684607496029970681036800000e31 * t9 * t5 * t46 +
             0.4871614964106936366598717440000e31 * t15 * t5 * t4 +
             0.8931294100862716672097648640000e31 * t20 * t5 * t3 +
             0.7176932759621825897221324800000e31 * t24 * t5 +
             0.2511926465867639064027463680000e31 * t28 * t46 +
             0.348878675814949870003814400000e30 * t31 * t4 +
             0.14951943249212137285877760000e29 * t34 * t3 +
             0.84954223006887143669760000e26 * t37) /
            0.1526822969691308336087040000e28);
}



double
AssLegFunction::P_19_3(const double x) const
{
  double t24;
  double t8;
  double t20;
  double t11;
  double t36;
  double t7;
  double t9;
  double t2;
  double t1;
  double t12;
  double t3;
  double t4;
  double t16;
  t1  = sqrt(0.373065e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t11 = -t3;
  t12 = t7 * t2;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t24 = t16 * t16;
  t36 = t24 * t24;
  return (-t1 * t4 * t3 *
          (0.71378973833068664712069120000e29 * t9 +
           0.1070684607496029970681036800000e31 * t11 * t8 * t12 +
           0.4871614964106936366598717440000e31 * t16 * t8 * t7 +
           0.8931294100862716672097648640000e31 * t20 * t8 * t2 +
           0.7176932759621825897221324800000e31 * t24 * t8 +
           0.2511926465867639064027463680000e31 * t24 * t11 * t12 +
           0.348878675814949870003814400000e30 * t24 * t16 * t7 +
           0.14951943249212137285877760000e29 * t24 * t20 * t2 +
           0.84954223006887143669760000e26 * t36) /
          0.285515895332274658848276480000e30);
}



double
AssLegFunction::P_19_3_Deriv(const double x) const
{
  double t4;
  double t19;
  double t11;
  double t35;
  double t32;
  double t10;
  double t6;
  double t7;
  double t8;
  double t29;
  double t47;
  double t1;
  double t3;
  double t2;
  double t43;
  double t44;
  double t15;
  double t26;
  double t23;
  t1  = sqrt(0.373065e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t10 = -t3;
  t11 = t6 * t2;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t23 = t15 * t15;
  t26 = t23 * t10;
  t29 = t23 * t15;
  t32 = t23 * t19;
  t35 = t23 * t23;
  t43 = t2 * x;
  t44 = t6 * t43;
  t47 = t6 * x;
  return (t1 * t4 *
            (0.71378973833068664712069120000e29 * t8 +
             0.1070684607496029970681036800000e31 * t10 * t7 * t11 +
             0.4871614964106936366598717440000e31 * t15 * t7 * t6 +
             0.8931294100862716672097648640000e31 * t19 * t7 * t2 +
             0.7176932759621825897221324800000e31 * t23 * t7 +
             0.2511926465867639064027463680000e31 * t26 * t11 +
             0.348878675814949870003814400000e30 * t29 * t6 +
             0.14951943249212137285877760000e29 * t32 * t2 +
             0.84954223006887143669760000e26 * t35) *
            x / 0.95171965110758219616092160000e29 -
          t1 * t4 * t3 *
            (0.3283432796321158576755179520000e31 * t7 * t44 +
             0.34476044361372165055929384960000e32 * t10 * t7 * t47 +
             0.112047144174459536431770501120000e33 * t15 * t7 * t43 +
             0.146728403085601773898747084800000e33 * t19 * t7 * x +
             0.82534726735650997818045235200000e32 * t23 * t44 +
             0.19258102904985232824210554880000e32 * t26 * t47 +
             0.1604841908748769402017546240000e31 * t29 * t43 +
             0.31263154066534468870471680000e29 * t32 * x) /
            0.285515895332274658848276480000e30);
}



double
AssLegFunction::P_19_4(const double x) const
{
  double t12;
  double t1;
  double t2;
  double t9;
  double t3;
  double t21;
  double t4;
  double t13;
  double t17;
  double t6;
  double t25;
  double t7;
  double t8;
  t1  = sqrt(0.8580495e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t6;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t7 * x;
  t17 = t12 * t12;
  t21 = t17 * t12;
  t25 = t17 * t17;
  return (t1 * t4 *
          (0.3283432796321158576755179520000e31 * t9 * t8 +
           0.34476044361372165055929384960000e32 * t12 * t9 * t13 +
           0.112047144174459536431770501120000e33 * t17 * t9 * t6 +
           0.146728403085601773898747084800000e33 * t21 * t9 * x +
           0.82534726735650997818045235200000e32 * t25 * t8 +
           0.19258102904985232824210554880000e32 * t25 * t12 * t13 +
           0.1604841908748769402017546240000e31 * t25 * t17 * t6 +
           0.31263154066534468870471680000e29 * t25 * t21 * x) /
          0.26267462370569268614041436160000e32);
}



double
AssLegFunction::P_19_4_Deriv(const double x) const
{
  double t20;
  double t40;
  double t42;
  double t24;
  double t1;
  double t27;
  double t3;
  double t2;
  double t6;
  double t5;
  double t7;
  double t30;
  double t8;
  double t33;
  double t11;
  double t12;
  double t16;
  t1  = sqrt(0.8580495e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * x;
  t6  = t2 * t2;
  t7  = t6 * t5;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t6 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t24 = t16 * t16;
  t27 = t24 * t11;
  t30 = t24 * t16;
  t33 = t24 * t20;
  t40 = t3 * t3;
  t42 = t6 * t2;
  return (-t1 * t3 *
            (0.3283432796321158576755179520000e31 * t7 * t8 +
             0.34476044361372165055929384960000e32 * t11 * t8 * t12 +
             0.112047144174459536431770501120000e33 * t16 * t8 * t5 +
             0.146728403085601773898747084800000e33 * t20 * t8 * x +
             0.82534726735650997818045235200000e32 * t24 * t7 +
             0.19258102904985232824210554880000e32 * t27 * t12 +
             0.1604841908748769402017546240000e31 * t30 * t5 +
             0.31263154066534468870471680000e29 * t33 * x) *
            x / 0.6566865592642317153510359040000e31 +
          t1 * t40 *
            (0.118203580667561708763186462720000e33 * t8 * t42 +
             0.896377153395676291454164008960000e33 * t11 * t8 * t6 +
             0.2112889004432665544141958021120000e34 * t16 * t8 * t2 +
             0.1980833441655623947633085644800000e34 * t20 * t8 +
             0.770324116199409312968422195200000e33 * t24 * t42 +
             0.115548617429911396945263329280000e33 * t27 * t6 +
             0.5252209883177790770239242240000e31 * t30 * t2 +
             0.31263154066534468870471680000e29 * t33) /
            0.26267462370569268614041436160000e32);
}



double
AssLegFunction::P_19_5(const double x) const
{
  double t2;
  double t3;
  double t1;
  double t4;
  double t5;
  double t17;
  double t9;
  double t13;
  double t24;
  double t21;
  double t10;
  double t8;
  t1  = sqrt(0.3432198e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t2;
  t10 = t8 * t8;
  t13 = -t3;
  t17 = t13 * t13;
  t21 = t17 * t13;
  t24 = t17 * t17;
  return (-t1 * t5 * t4 *
          (0.118203580667561708763186462720000e33 * t10 * t9 +
           0.896377153395676291454164008960000e33 * t13 * t10 * t8 +
           0.2112889004432665544141958021120000e34 * t17 * t10 * t2 +
           0.1980833441655623947633085644800000e34 * t21 * t10 +
           0.770324116199409312968422195200000e33 * t24 * t9 +
           0.115548617429911396945263329280000e33 * t24 * t13 * t8 +
           0.5252209883177790770239242240000e31 * t24 * t17 * t2 +
           0.31263154066534468870471680000e29 * t24 * t21) /
          0.315209548446831223368497233920000e33);
}



double
AssLegFunction::P_19_5_Deriv(const double x) const
{
  double t7;
  double t8;
  double t9;
  double t12;
  double t38;
  double t41;
  double t44;
  double t1;
  double t2;
  double t3;
  double t4;
  double t20;
  double t23;
  double t26;
  double t16;
  double t29;
  t1  = sqrt(0.3432198e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t2;
  t9  = t7 * t7;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t23 = t16 * t16;
  t26 = t23 * t12;
  t29 = t23 * t16;
  t38 = t3 * t3;
  t41 = t7 * x;
  t44 = t2 * x;
  return (t1 * t4 * t3 *
            (0.118203580667561708763186462720000e33 * t8 * t9 +
             0.896377153395676291454164008960000e33 * t12 * t9 * t7 +
             0.2112889004432665544141958021120000e34 * t16 * t9 * t2 +
             0.1980833441655623947633085644800000e34 * t20 * t9 +
             0.770324116199409312968422195200000e33 * t23 * t8 +
             0.115548617429911396945263329280000e33 * t26 * t7 +
             0.5252209883177790770239242240000e31 * t29 * t2 +
             0.31263154066534468870471680000e29 * t23 * t20) *
            x / 0.63041909689366244673699446784000e32 -
          t1 * t4 * t38 *
            (0.3447604436137216505592938496000000e34 * t9 * t41 +
             0.19208081858478777674017800192000000e35 * t12 * t9 * t44 +
             0.33013890694260399127218094080000000e35 * t16 * t9 * x +
             0.22009260462840266084812062720000000e35 * t20 * t7 * t44 +
             0.5777430871495569847263166464000000e34 * t23 * t41 +
             0.525220988317779077023924224000000e33 * t26 * t44 +
             0.10942103923287064104665088000000e32 * t29 * x) /
            0.315209548446831223368497233920000e33);
}



double
AssLegFunction::P_19_6(const double x) const
{
  double t7;
  double t8;
  double t9;
  double t2;
  double t12;
  double t13;
  double t3;
  double t1;
  double t17;
  double t4;
  double t25;
  t1  = sqrt(0.245157e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * x;
  t9  = t7 * t7;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t25 = t17 * t17;
  return (t1 * t4 * t3 *
          (0.3447604436137216505592938496000000e34 * t8 * t9 +
           0.19208081858478777674017800192000000e35 * t12 * t9 * t13 +
           0.33013890694260399127218094080000000e35 * t17 * t9 * x +
           0.22009260462840266084812062720000000e35 * t17 * t12 * t7 * t13 +
           0.5777430871495569847263166464000000e34 * t25 * t8 +
           0.525220988317779077023924224000000e33 * t25 * t12 * t13 +
           0.10942103923287064104665088000000e32 * t25 * t17 * x) /
          0.1576047742234156116842486169600000e34);
}



double
AssLegFunction::P_19_6_Deriv(const double x) const
{
  double t11;
  double t7;
  double t12;
  double t2;
  double t30;
  double t20;
  double t3;
  double t4;
  double t8;
  double t24;
  double t16;
  double t27;
  double t1;
  double t6;
  t1  = sqrt(0.245157e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * x;
  t8  = t6 * t6;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t24 = t16 * t16;
  t27 = t24 * t11;
  t30 = t24 * t16;
  return (-t1 * t4 *
            (0.3447604436137216505592938496000000e34 * t8 * t7 +
             0.19208081858478777674017800192000000e35 * t11 * t8 * t12 +
             0.33013890694260399127218094080000000e35 * t16 * t8 * x +
             0.22009260462840266084812062720000000e35 * t20 * t6 * t12 +
             0.5777430871495569847263166464000000e34 * t24 * t7 +
             0.525220988317779077023924224000000e33 * t27 * t12 +
             0.10942103923287064104665088000000e32 * t30 * x) *
            x / 0.262674623705692686140414361600000e33 +
          t1 * t4 * t3 *
            (0.83235021386741369920743800832000000e35 * t8 * t6 +
             0.343344463220308150923068178432000000e36 * t11 * t8 * t2 +
             0.429180579025385188653835223040000000e36 * t16 * t8 +
             0.200284270211846421371789770752000000e36 * t20 * t6 * t2 +
             0.34139364240655640006555074560000000e35 * t24 * t6 +
             0.1706968212032782000327753728000000e34 * t27 * t2 +
             0.10942103923287064104665088000000e32 * t30) /
            0.1576047742234156116842486169600000e34);
}



double
AssLegFunction::P_19_7(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t9;
  double t17;
  double t10;
  double t13;
  double t24;
  t1  = sqrt(0.490314e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t17 = t13 * t13;
  t24 = t17 * t17;
  return (-t1 * t6 * t4 * t3 *
          (0.83235021386741369920743800832000000e35 * t10 * t9 +
           0.343344463220308150923068178432000000e36 * t13 * t10 * t2 +
           0.429180579025385188653835223040000000e36 * t17 * t10 +
           0.200284270211846421371789770752000000e36 * t17 * t13 * t9 * t2 +
           0.34139364240655640006555074560000000e35 * t24 * t9 +
           0.1706968212032782000327753728000000e34 * t24 * t13 * t2 +
           0.10942103923287064104665088000000e32 * t24 * t17) /
          0.40977241298088059037904640409600000e35);
}



double
AssLegFunction::P_19_7_Deriv(const double x) const
{
  double t9;
  double t2;
  double t3;
  double t12;
  double t23;
  double t4;
  double t5;
  double t16;
  double t26;
  double t1;
  double t19;
  double t38;
  double t8;
  t1  = sqrt(0.490314e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t19 = t16 * t12;
  t23 = t16 * t16;
  t26 = t23 * t12;
  t38 = t2 * x;
  return (t1 * t5 * t4 *
            (0.83235021386741369920743800832000000e35 * t9 * t8 +
             0.343344463220308150923068178432000000e36 * t12 * t9 * t2 +
             0.429180579025385188653835223040000000e36 * t16 * t9 +
             0.200284270211846421371789770752000000e36 * t19 * t8 * t2 +
             0.34139364240655640006555074560000000e35 * t23 * t8 +
             0.1706968212032782000327753728000000e34 * t26 * t2 +
             0.10942103923287064104665088000000e32 * t23 * t16) *
            x / 0.5853891614012579862557805772800000e34 -
          t1 * t5 * t4 * t3 *
            (0.1685509183081512740895061966848000000e37 * t9 * t38 +
             0.5150166948304622263846022676480000000e37 * t12 * t9 * x +
             0.4635150253474160037461420408832000000e37 * t16 * t8 * t38 +
             0.1474820535196323648283179220992000000e37 * t19 * t8 * x +
             0.153627139082950380029497835520000000e36 * t23 * t38 +
             0.3545241671145008769911488512000000e34 * t26 * x) /
            0.40977241298088059037904640409600000e35);
}



double
AssLegFunction::P_19_8(const double x) const
{
  double t12;
  double t4;
  double t16;
  double t24;
  double t1;
  double t5;
  double t8;
  double t9;
  double t2;
  double t3;
  double t7;
  t1  = sqrt(0.490314e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t24 = t16 * t16;
  return (t1 * t5 *
          (0.1685509183081512740895061966848000000e37 * t9 * t7 +
           0.5150166948304622263846022676480000000e37 * t12 * t9 * x +
           0.4635150253474160037461420408832000000e37 * t16 * t8 * t7 +
           0.1474820535196323648283179220992000000e37 * t16 * t12 * t8 * x +
           0.153627139082950380029497835520000000e36 * t24 * t7 +
           0.3545241671145008769911488512000000e34 * t24 * t12 * x) /
          0.737590343365585062682283527372800000e36);
}



double
AssLegFunction::P_19_8_Deriv(const double x) const
{
  double t1;
  double t8;
  double t16;
  double t9;
  double t27;
  double t2;
  double t20;
  double t3;
  double t7;
  double t12;
  double t34;
  double t4;
  double t24;
  t1  = sqrt(0.490314e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t24 = t16 * t16;
  t27 = t24 * t12;
  t34 = t4 * t4;
  return (-t1 * t4 * t3 *
            (0.1685509183081512740895061966848000000e37 * t9 * t7 +
             0.5150166948304622263846022676480000000e37 * t12 * t9 * x +
             0.4635150253474160037461420408832000000e37 * t16 * t8 * t7 +
             0.1474820535196323648283179220992000000e37 * t20 * t8 * x +
             0.153627139082950380029497835520000000e36 * t24 * t7 +
             0.3545241671145008769911488512000000e34 * t27 * x) *
            x / 0.92198792920698132835285440921600000e35 +
          t1 * t34 *
            (0.28840934910505884677537726988288000000e38 * t9 * t2 +
             0.64892103548638240524459885723648000000e38 * t12 * t9 +
             0.41294974985497062151929018187776000000e38 * t16 * t8 * t2 +
             0.8603119788645221281651878789120000000e37 * t20 * t8 +
             0.496333833960301227787608391680000000e36 * t24 * t2 +
             0.3545241671145008769911488512000000e34 * t27) /
            0.737590343365585062682283527372800000e36);
}



double
AssLegFunction::P_19_9(const double x) const
{
  double t9;
  double t10;
  double t4;
  double t5;
  double t13;
  double t3;
  double t16;
  double t6;
  double t2;
  double t1;
  double t23;
  t1  = sqrt(0.312018e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t16 = t13 * t13;
  t23 = t16 * t16;
  return (-t1 * t6 * t5 *
          (0.28840934910505884677537726988288000000e38 * t10 * t2 +
           0.64892103548638240524459885723648000000e38 * t13 * t10 +
           0.41294974985497062151929018187776000000e38 * t16 * t9 * t2 +
           0.8603119788645221281651878789120000000e37 * t16 * t13 * t9 +
           0.496333833960301227787608391680000000e36 * t23 * t2 +
           0.3545241671145008769911488512000000e34 * t23 * t13) /
          0.10326264807118190877551969383219200000e38);
}



double
AssLegFunction::P_19_9_Deriv(const double x) const
{
  double t4;
  double t13;
  double t37;
  double t23;
  double t32;
  double t6;
  double t16;
  double t9;
  double t10;
  double t1;
  double t20;
  double t3;
  double t2;
  t1  = sqrt(0.312018e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t16 = t13 * t13;
  t20 = t16 * t13;
  t23 = t16 * t16;
  t32 = t4 * t4;
  t37 = t2 * x;
  return (t1 * t6 * t4 * t3 *
            (0.28840934910505884677537726988288000000e38 * t10 * t2 +
             0.64892103548638240524459885723648000000e38 * t13 * t10 +
             0.41294974985497062151929018187776000000e38 * t16 * t9 * t2 +
             0.8603119788645221281651878789120000000e37 * t20 * t9 +
             0.496333833960301227787608391680000000e36 * t23 * t2 +
             0.3545241671145008769911488512000000e34 * t23 * t13) *
            x / 0.1147362756346465653061329931468800000e37 -
          t1 * t6 * t32 *
            (0.418193556202335327824297041330176000000e39 * t10 * x +
             0.684316728331094172803395158540288000000e39 * t13 * t9 * t37 +
             0.299388568644853700601485381861376000000e39 * t16 * t9 * x +
             0.38383149826263294948908382289920000000e38 * t20 * t37 +
             0.1028120084632052543274331668480000000e37 * t23 * x) /
            0.10326264807118190877551969383219200000e38);
}



double
AssLegFunction::P_19_10(const double x) const
{
  double t17;
  double t4;
  double t5;
  double t13;
  double t24;
  double t1;
  double t8;
  double t9;
  double t2;
  double t3;
  double t12;
  t1  = sqrt(0.22621305e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t13 = t2 * x;
  t17 = t12 * t12;
  t24 = t17 * t17;
  return (t1 * t5 * t3 *
          (0.418193556202335327824297041330176000000e39 * t9 * x +
           0.684316728331094172803395158540288000000e39 * t12 * t8 * t13 +
           0.299388568644853700601485381861376000000e39 * t17 * t8 * x +
           0.38383149826263294948908382289920000000e38 * t17 * t12 * t13 +
           0.1028120084632052543274331668480000000e37 * t24 * x) /
          0.1497308397032137677245035560566784000000e40);
}



double
AssLegFunction::P_19_10_Deriv(const double x) const
{
  double t4;
  double t5;
  double t16;
  double t1;
  double t11;
  double t8;
  double t12;
  double t2;
  double t3;
  double t20;
  double t7;
  double t23;
  t1  = sqrt(0.22621305e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t16 = t11 * t11;
  t20 = t16 * t11;
  t23 = t16 * t16;
  return (-t1 * t5 *
            (0.418193556202335327824297041330176000000e39 * t8 * x +
             0.684316728331094172803395158540288000000e39 * t11 * t7 * t12 +
             0.299388568644853700601485381861376000000e39 * t16 * t7 * x +
             0.38383149826263294948908382289920000000e38 * t20 * t12 +
             0.1028120084632052543274331668480000000e37 * t23 * x) *
            x / 0.149730839703213767724503556056678400000e39 +
          t1 * t5 * t3 *
            (0.5132375462483206296025463689052160000000e40 * t8 +
             0.5987771372897074012029707637227520000000e40 * t11 * t7 * t2 +
             0.1727241742181848272700877203046400000000e40 * t16 * t7 +
             0.123374410155846305192919800217600000000e39 * t20 * t2 +
             0.1028120084632052543274331668480000000e37 * t23) /
            0.1497308397032137677245035560566784000000e40);
}



double
AssLegFunction::P_19_11(const double x) const
{
  double t13;
  double t1;
  double t2;
  double t3;
  double t23;
  double t5;
  double t4;
  double t17;
  double t7;
  double t11;
  double t10;
  t1  = sqrt(0.3016174e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t11 = t10 * t10;
  t13 = -t3;
  t17 = t13 * t13;
  t23 = t17 * t17;
  return (-t1 * t7 * t5 * t3 *
          (0.5132375462483206296025463689052160000000e40 * t11 +
           0.5987771372897074012029707637227520000000e40 * t13 * t10 * t2 +
           0.1727241742181848272700877203046400000000e40 * t17 * t10 +
           0.123374410155846305192919800217600000000e39 * t17 * t13 * t2 +
           0.1028120084632052543274331668480000000e37 * t23) /
          0.8983850382192826063470213363400704000000e40);
}



double
AssLegFunction::P_19_11_Deriv(const double x) const
{
  double t16;
  double t4;
  double t2;
  double t3;
  double t5;
  double t6;
  double t22;
  double t9;
  double t10;
  double t19;
  double t31;
  double t1;
  double t12;
  t1  = sqrt(0.3016174e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t12 = -t3;
  t16 = t12 * t12;
  t19 = t16 * t12;
  t22 = t16 * t16;
  t31 = t2 * x;
  return (t1 * t6 * t5 *
            (0.5132375462483206296025463689052160000000e40 * t10 +
             0.5987771372897074012029707637227520000000e40 * t12 * t9 * t2 +
             0.1727241742181848272700877203046400000000e40 * t16 * t9 +
             0.123374410155846305192919800217600000000e39 * t19 * t2 +
             0.1028120084632052543274331668480000000e37 * t22) *
            x / 0.816713671108438733042746669400064000000e39 -
          t1 * t6 * t5 * t3 *
            (0.53034546445659798392263124786872320000000e41 * t9 * t31 +
             0.42835595206109837162981754635550720000000e41 * t12 * t9 * x +
             0.7649213429662470921961027613491200000000e40 * t16 * t31 +
             0.254973780988749030732034253783040000000e39 * t19 * x) /
            0.8983850382192826063470213363400704000000e40);
}



double
AssLegFunction::P_19_12(const double x) const
{
  double t8;
  double t9;
  double t16;
  double t1;
  double t2;
  double t3;
  double t12;
  double t4;
  double t5;
  t1  = sqrt(0.46750697e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * x;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (t1 * t5 * t4 *
          (0.53034546445659798392263124786872320000000e41 * t9 * t8 +
           0.42835595206109837162981754635550720000000e41 * t12 * t9 * x +
           0.7649213429662470921961027613491200000000e40 * t16 * t8 +
           0.254973780988749030732034253783040000000e39 * t16 * t12 * x) /
          0.556998723695955215935153228530843648000000e42);
}



double
AssLegFunction::P_19_12_Deriv(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t5;
  double t19;
  double t16;
  double t8;
  double t9;
  double t12;
  double t1;
  t1  = sqrt(0.46750697e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * x;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  t19 = t16 * t12;
  return (-t1 * t5 * t3 *
            (0.53034546445659798392263124786872320000000e41 * t9 * t8 +
             0.42835595206109837162981754635550720000000e41 * t12 * t9 * x +
             0.7649213429662470921961027613491200000000e40 * t8 * t16 +
             0.254973780988749030732034253783040000000e39 * t19 * x) *
            x / 0.46416560307996267994596102377570304000000e41 +
          t1 * t5 * t4 *
            (0.456913015531838263071805382779207680000000e42 * t9 * t2 +
             0.244774829749199069502752883631718400000000e42 * t12 * t9 +
             0.24477482974919906950275288363171840000000e41 * t16 * t2 +
             0.254973780988749030732034253783040000000e39 * t19) /
            0.556998723695955215935153228530843648000000e42);
}



double
AssLegFunction::P_19_13(const double x) const
{
  double t13;
  double t2;
  double t16;
  double t3;
  double t1;
  double t4;
  double t10;
  double t7;
  double t5;
  t1  = sqrt(0.13357342e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t16 = t13 * t13;
  return (-t1 * t7 * t5 * t4 *
          (0.456913015531838263071805382779207680000000e42 * t10 * t2 +
           0.244774829749199069502752883631718400000000e42 * t13 * t10 +
           0.24477482974919906950275288363171840000000e41 * t2 * t16 +
           0.254973780988749030732034253783040000000e39 * t16 * t13) /
          0.4455989789567641727481225828246749184000000e43);
}



double
AssLegFunction::P_19_13_Deriv(const double x) const
{
  double t2;
  double t13;
  double t3;
  double t4;
  double t5;
  double t16;
  double t7;
  double t10;
  double t1;
  t1  = sqrt(0.13357342e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t16 = t13 * t13;
  return (t1 * t7 * t5 * t3 *
            (0.456913015531838263071805382779207680000000e42 * t10 * t2 +
             0.244774829749199069502752883631718400000000e42 * t13 * t10 +
             0.24477482974919906950275288363171840000000e41 * t16 * t2 +
             0.254973780988749030732034253783040000000e39 * t16 * t13) *
            x / 0.342768445351357055960094294480519168000000e42 -
          t1 * t7 * t5 * t4 *
            (0.3231027752689427717436338063938682880000000e43 * t10 * x +
             0.1077009250896475905812112687979560960000000e43 * t13 * t2 * x +
             0.50484808635772308084942782249041920000000e41 * t16 * x) /
            0.4455989789567641727481225828246749184000000e43);
}



double
AssLegFunction::P_19_14(const double x) const
{
  double t12;
  double t1;
  double t16;
  double t9;
  double t6;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.73465381e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t9  = t2 * t2;
  t12 = -t3;
  t16 = t12 * t12;
  return (t1 * t6 * t4 * t3 *
          (0.3231027752689427717436338063938682880000000e43 * t9 * x +
           0.1077009250896475905812112687979560960000000e43 * t12 * t2 * x +
           0.50484808635772308084942782249041920000000e41 * t16 * x) /
          0.147047663055732177006880452332142723072000000e45);
}



double
AssLegFunction::P_19_14_Deriv(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t8;
  double t11;
  double t15;
  double t1;
  t1  = sqrt(0.73465381e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t15 = t11 * t11;
  return (-t1 * t5 * t4 *
            (0.3231027752689427717436338063938682880000000e43 * t8 * x +
             0.1077009250896475905812112687979560960000000e43 * t2 * t11 * x +
             0.50484808635772308084942782249041920000000e41 * t15 * x) *
            x / 0.10503404503980869786205746595153051648000000e44 +
          t1 * t5 * t4 * t3 *
            (0.18309157265240090398805915695652536320000000e44 * t8 +
             0.3432966987232516949776109192934850560000000e43 * t11 * t2 +
             0.50484808635772308084942782249041920000000e41 * t15) /
            0.147047663055732177006880452332142723072000000e45);
}



double
AssLegFunction::P_19_15(const double x) const
{
  double t16;
  double t4;
  double t3;
  double t6;
  double t8;
  double t1;
  double t11;
  double t2;
  double t13;
  t1  = sqrt(0.43214930e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t8  = sqrt(t3);
  t11 = t2 * t2;
  t13 = -t3;
  t16 = t13 * t13;
  return (-t1 * t8 * t6 * t4 * t3 *
          (0.18309157265240090398805915695652536320000000e44 * t11 +
           0.3432966987232516949776109192934850560000000e43 * t13 * t2 +
           0.50484808635772308084942782249041920000000e41 * t16) /
          0.1470476630557321770068804523321427230720000000e46);
}



double
AssLegFunction::P_19_15_Deriv(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t12;
  double t15;
  double t5;
  double t1;
  double t7;
  double t10;
  t1  = sqrt(0.43214930e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (t1 * t7 * t5 * t4 *
            (0.18309157265240090398805915695652536320000000e44 * t10 +
             0.3432966987232516949776109192934850560000000e43 * t12 * t2 +
             0.50484808635772308084942782249041920000000e41 * t15) *
            x / 0.98031775370488118004586968221428482048000000e44 -
          t1 * t7 * t5 * t4 * t3 *
            (0.80102563035425395494775881168479846400000000e44 * t2 * x +
             0.7067873209008123131891989514865868800000000e43 * t12 * x) /
            0.1470476630557321770068804523321427230720000000e46);
}



double
AssLegFunction::P_19_16(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  t1 = sqrt(0.60500902e8);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (t1 * t6 *
          (0.80102563035425395494775881168479846400000000e44 * t2 * x -
           0.7067873209008123131891989514865868800000000e43 * t3 * x) /
          0.20586672827802504780963263326499981230080000000e47);
}



double
AssLegFunction::P_19_16_Deriv(const double x) const
{
  double t18;
  double t1;
  double t6;
  double t2;
  double t3;
  double t4;
  t1  = sqrt(0.60500902e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t18 = t6 * t6;
  return (-t1 * t6 * t4 * t3 *
            (0.80102563035425395494775881168479846400000000e44 * t2 * x -
             0.7067873209008123131891989514865868800000000e43 * t3 * x) *
            x / 0.1286667051737656548810203957906248826880000000e46 +
          t1 * t18 *
            (0.261511308733300555880003612050037145600000000e45 * t2 -
             0.7067873209008123131891989514865868800000000e43) /
            0.20586672827802504780963263326499981230080000000e47);
}



double
AssLegFunction::P_19_17(const double x) const
{
  double t1;
  double t3;
  double t2;
  double t5;
  double t4;
  double t6;
  double t7;
  t1 = sqrt(0.181502706e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (-t1 * t7 * t6 *
          (0.261511308733300555880003612050037145600000000e45 * t2 -
           0.7067873209008123131891989514865868800000000e43) /
          0.370560110900445086057338739876999662141440000000e48);
}



double
AssLegFunction::P_19_17_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t8;
  double t16;
  t1  = sqrt(0.181502706e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t8  = sqrt(t3);
  t16 = t6 * t6;
  return (t1 * t8 * t6 * t4 * t3 *
            (0.261511308733300555880003612050037145600000000e45 * t2 -
             0.7067873209008123131891989514865868800000000e43) *
            x / 0.21797653582379122709255219992764686008320000000e47 -
          0.185e3 / 0.131072e6 * t1 * t8 * t16 * x);
}



double
AssLegFunction::P_19_18(const double x) const
{
  double t1;
  double t6;
  double t2;
  double t3;
  double t4;
  double t5;
  t1 = sqrt(0.3357800061e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (0.5e1 / 0.131072e6 * t1 * t6 * t3 * x);
}



double
AssLegFunction::P_19_18_Deriv(const double x) const
{
  double t1;
  double t4;
  double t2;
  double t3;
  double t5;
  double t6;
  t1 = sqrt(0.3357800061e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (-0.45e2 / 0.65536e5 * t1 * t6 * t2 +
          0.5e1 / 0.131072e6 * t1 * t6 * t3);
}



double
AssLegFunction::P_19_19(const double x) const
{
  double t3;
  double t4;
  double t5;
  double t6;
  double t8;
  double t1;
  double t2;
  t1 = sqrt(0.353452638e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t8 = sqrt(t3);
  return (-0.5e1 / 0.262144e6 * t1 * t8 * t6 * t3);
}



double
AssLegFunction::P_19_19_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t6;
  double t5;
  double t7;
  t1 = sqrt(0.353452638e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (0.95e2 / 0.262144e6 * t1 * t7 * t6 * x);
}

double
AssLegFunction::P_20_0(const double x) const
{
  double t13;
  double t18;
  double t4;
  double t14;
  double t1;
  double t2;
  double t3;
  double t32;
  double t10;
  double t6;
  t1  = x * x;
  t2  = t1 * t1;
  t3  = t2 * t2;
  t4  = t3 * t3;
  t6  = t1 - 0.1e1;
  t10 = t6 * t6;
  t13 = t10 * t6;
  t14 = t2 * t1;
  t18 = t10 * t10;
  t32 = t18 * t18;
  return (t4 * t2 + 0.95e2 * t6 * t4 * t1 + 0.14535e5 / 0.8e1 * t10 * t4 +
          0.24225e5 / 0.2e1 * t13 * t3 * t14 +
          0.2204475e7 / 0.64e2 * t18 * t3 * t2 +
          0.2909907e7 / 0.64e2 * t18 * t6 * t3 * t1 +
          0.14549535e8 / 0.512e3 * t18 * t10 * t3 +
          0.2078505e7 / 0.256e3 * t18 * t13 * t14 +
          0.31177575e8 / 0.32768e5 * t32 * t2 +
          0.1154725e7 / 0.32768e5 * t32 * t6 * t1 +
          0.46189e5 / 0.262144e6 * t32 * t10);
}

double
AssLegFunction::P_20_0_Deriv(const double x) const
{
  double t13;
  double t17;
  double t36;
  double t1;
  double t2;
  double t18;
  double t8;
  double t3;
  double t4;
  double t5;
  double t12;
  double t22;
  t1  = x * x;
  t2  = t1 * x;
  t3  = t1 * t1;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t1 - 0.1e1;
  t12 = t8 * t8;
  t13 = t3 * t2;
  t17 = t12 * t8;
  t18 = t3 * x;
  t22 = t12 * t12;
  t36 = t22 * t22;
  return (0.210e3 * t5 * t2 + 0.17955e5 / 0.2e1 * t8 * t5 * x +
          0.101745e6 * t12 * t4 * t13 + 0.3561075e7 / 0.8e1 * t17 * t4 * t18 +
          0.27776385e8 / 0.32e2 * t22 * t4 * t2 +
          0.101846745e9 / 0.128e3 * t22 * t8 * t4 * x +
          0.43648605e8 / 0.128e3 * t22 * t12 * t13 +
          0.130945815e9 / 0.2048e4 * t22 * t17 * t18 +
          0.72747675e8 / 0.16384e5 * t36 * t2 +
          0.4849845e7 / 0.65536e5 * t36 * t8 * x);
}



double
AssLegFunction::P_20_1(const double x) const
{
  double t12;
  double t40;
  double t16;
  double t1;
  double t9;
  double t22;
  double t26;
  double t7;
  double t21;
  double t17;
  double t2;
  double t3;
  double t6;
  double t8;
  double t4;
  t1  = sqrt(0.105e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * x;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t12 = -t3;
  t16 = t12 * t12;
  t17 = t7 * t6;
  t21 = t16 * t12;
  t22 = t7 * x;
  t26 = t16 * t16;
  t40 = t26 * t26;
  return (-t1 * t4 *
          (0.535727357786423977574400000e27 * t9 * t6 +
           0.22902344545369625041305600000e29 * t12 * t9 * x +
           0.259559904847522417134796800000e30 * t16 * t8 * t17 +
           0.1135574583707910574964736000000e31 * t21 * t8 * t22 +
           0.2214370438230425621181235200000e31 * t26 * t8 * t6 +
           0.2029839568377890152749465600000e31 * t26 * t12 * t8 * x +
           0.869931243590524351178342400000e30 * t26 * t16 * t17 +
           0.163112108173223315845939200000e30 * t26 * t21 * t22 +
           0.11327229734251619155968000000e29 * t40 * t6 +
           0.188787162237526985932800000e27 * t40 * t12 * x) /
          0.535727357786423977574400000e27);
}



double
AssLegFunction::P_20_1_Deriv(const double x) const
{
  double t10;
  double t1;
  double t2;
  double t3;
  double t17;
  double t18;
  double t23;
  double t22;
  double t4;
  double t13;
  double t55;
  double t31;
  double t9;
  double t7;
  double t8;
  double t27;
  double t44;
  double t41;
  double t38;
  double t35;
  t1  = sqrt(0.105e3);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t10 = t9 * t9;
  t13 = -t3;
  t17 = t13 * t13;
  t18 = t8 * t7;
  t22 = t17 * t13;
  t23 = t8 * x;
  t27 = t17 * t17;
  t31 = t27 * t13;
  t35 = t27 * t17;
  t38 = t27 * t22;
  t41 = t27 * t27;
  t44 = t41 * t13;
  t55 = t8 * t2;
  return (t1 / t4 *
            (0.535727357786423977574400000e27 * t10 * t7 +
             0.22902344545369625041305600000e29 * t13 * t10 * x +
             0.259559904847522417134796800000e30 * t17 * t9 * t18 +
             0.1135574583707910574964736000000e31 * t22 * t9 * t23 +
             0.2214370438230425621181235200000e31 * t27 * t9 * t7 +
             0.2029839568377890152749465600000e31 * t31 * t9 * x +
             0.869931243590524351178342400000e30 * t35 * t18 +
             0.163112108173223315845939200000e30 * t38 * t23 +
             0.11327229734251619155968000000e29 * t41 * t7 +
             0.188787162237526985932800000e27 * t44 * x) *
            x / 0.535727357786423977574400000e27 -
          t1 * t4 *
            (0.55983508888681305656524800000e29 * t2 * t10 +
             0.1427579476661373294241382400000e31 * t13 * t10 +
             0.10706846074960299706810368000000e32 * t17 * t9 * t55 +
             0.32477433094046242443991449600000e32 * t22 * t9 * t8 +
             0.44656470504313583360488243200000e32 * t27 * t9 * t2 +
             0.28707731038487303588885299200000e32 * t31 * t9 +
             0.8373088219558796880091545600000e31 * t35 * t55 +
             0.996796216614142485725184000000e30 * t38 * t8 +
             0.37379858123030343214694400000e29 * t41 * t2 +
             0.188787162237526985932800000e27 * t44) /
            0.535727357786423977574400000e27);
}



double
AssLegFunction::P_20_2(const double x) const
{
  double t6;
  double t13;
  double t14;
  double t5;
  double t7;
  double t18;
  double t35;
  double t22;
  double t1;
  double t2;
  double t10;
  double t3;
  t1  = sqrt(0.43890e5);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t7  = t6 * t6;
  t10 = -t3;
  t13 = t10 * t10;
  t14 = t5 * t2;
  t18 = t13 * t10;
  t22 = t13 * t13;
  t35 = t22 * t22;
  return (t1 * t3 *
          (0.55983508888681305656524800000e29 * t7 * t2 +
           0.1427579476661373294241382400000e31 * t10 * t7 +
           0.10706846074960299706810368000000e32 * t13 * t6 * t14 +
           0.32477433094046242443991449600000e32 * t18 * t6 * t5 +
           0.44656470504313583360488243200000e32 * t22 * t6 * t2 +
           0.28707731038487303588885299200000e32 * t22 * t10 * t6 +
           0.8373088219558796880091545600000e31 * t22 * t13 * t14 +
           0.996796216614142485725184000000e30 * t22 * t18 * t5 +
           0.37379858123030343214694400000e29 * t35 * t2 +
           0.188787162237526985932800000e27 * t35 * t10) /
          0.223934035554725222626099200000e30);
}



double
AssLegFunction::P_20_2_Deriv(const double x) const
{
  double t1;
  double t21;
  double t31;
  double t9;
  double t12;
  double t25;
  double t13;
  double t28;
  double t34;
  double t3;
  double t46;
  double t47;
  double t17;
  double t4;
  double t5;
  double t51;
  double t6;
  t1  = sqrt(0.43890e5);
  t3  = x * x;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = t5 * t5;
  t9  = t3 - 0.1e1;
  t12 = t9 * t9;
  t13 = t4 * t3;
  t17 = t12 * t9;
  t21 = t12 * t12;
  t25 = t21 * t9;
  t28 = t21 * t12;
  t31 = t21 * t17;
  t34 = t21 * t21;
  t46 = t3 * x;
  t47 = t4 * t46;
  t51 = t4 * x;
  return (-t1 * x *
            (0.55983508888681305656524800000e29 * t6 * t3 +
             0.1427579476661373294241382400000e31 * t9 * t6 +
             0.10706846074960299706810368000000e32 * t12 * t5 * t13 +
             0.32477433094046242443991449600000e32 * t17 * t5 * t4 +
             0.44656470504313583360488243200000e32 * t21 * t5 * t3 +
             0.28707731038487303588885299200000e32 * t25 * t5 +
             0.8373088219558796880091545600000e31 * t28 * t13 +
             0.996796216614142485725184000000e30 * t31 * t4 +
             0.37379858123030343214694400000e29 * t34 * t3 +
             0.188787162237526985932800000e27 * t34 * t9) /
            0.111967017777362611313049600000e30 -
          t1 * t9 *
            (0.3862862113319010090300211200000e31 * t6 * x +
             0.65668655926423171535103590400000e32 * t9 * t5 * t47 +
             0.344760443613721650559293849600000e33 * t12 * t5 * t51 +
             0.746980961163063576211803340800000e33 * t17 * t5 * t46 +
             0.733642015428008869493735424000000e33 * t21 * t5 * x +
             0.330138906942603991272180940800000e33 * t25 * t47 +
             0.64193676349950776080701849600000e32 * t28 * t51 +
             0.4585262596425055434335846400000e31 * t31 * t46 +
             0.78157885166336172176179200000e29 * t34 * x) /
            0.223934035554725222626099200000e30);
}



double
AssLegFunction::P_20_3(const double x) const
{
  double t27;
  double t8;
  double t18;
  double t13;
  double t3;
  double t7;
  double t1;
  double t19;
  double t23;
  double t40;
  double t14;
  double t4;
  double t9;
  double t2;
  double t12;
  t1  = sqrt(0.504735e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * t2;
  t8  = t7 * t7;
  t9  = t8 * t8;
  t12 = -t3;
  t13 = t2 * x;
  t14 = t7 * t13;
  t18 = t12 * t12;
  t19 = t7 * x;
  t23 = t18 * t12;
  t27 = t18 * t18;
  t40 = t27 * t27;
  return (-t1 * t4 * t3 *
          (0.3862862113319010090300211200000e31 * t9 * x +
           0.65668655926423171535103590400000e32 * t12 * t8 * t14 +
           0.344760443613721650559293849600000e33 * t18 * t8 * t19 +
           0.746980961163063576211803340800000e33 * t23 * t8 * t13 +
           0.733642015428008869493735424000000e33 * t27 * t8 * x +
           0.330138906942603991272180940800000e33 * t27 * t12 * t14 +
           0.64193676349950776080701849600000e32 * t27 * t18 * t19 +
           0.4585262596425055434335846400000e31 * t27 * t23 * t13 +
           0.78157885166336172176179200000e29 * t40 * x) /
          0.15451448453276040361200844800000e32);
}



double
AssLegFunction::P_20_3_Deriv(const double x) const
{
  double t22;
  double t49;
  double t12;
  double t39;
  double t13;
  double t6;
  double t26;
  double t7;
  double t33;
  double t1;
  double t17;
  double t18;
  double t8;
  double t2;
  double t3;
  double t11;
  double t4;
  double t36;
  double t30;
  t1  = sqrt(0.504735e6);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t11 = -t3;
  t12 = t2 * x;
  t13 = t6 * t12;
  t17 = t11 * t11;
  t18 = t6 * x;
  t22 = t17 * t11;
  t26 = t17 * t17;
  t30 = t26 * t11;
  t33 = t26 * t17;
  t36 = t26 * t22;
  t39 = t26 * t26;
  t49 = t6 * t2;
  return (t1 * t4 *
            (0.3862862113319010090300211200000e31 * t8 * x +
             0.65668655926423171535103590400000e32 * t11 * t7 * t13 +
             0.344760443613721650559293849600000e33 * t17 * t7 * t18 +
             0.746980961163063576211803340800000e33 * t22 * t12 * t7 +
             0.733642015428008869493735424000000e33 * t26 * t7 * x +
             0.330138906942603991272180940800000e33 * t30 * t13 +
             0.64193676349950776080701849600000e32 * t33 * t18 +
             0.4585262596425055434335846400000e31 * t36 * t12 +
             0.78157885166336172176179200000e29 * t39 * x) *
            x / 0.5150482817758680120400281600000e31 -
          t1 * t4 * t3 *
            (0.197005967779269514605310771200000e33 * t8 +
             0.2364071613351234175263729254400000e34 * t11 * t7 * t49 +
             0.8963771533956762914541640089600000e34 * t17 * t7 * t6 +
             0.14085926696217770294279720140800000e35 * t22 * t7 * t2 +
             0.9904167208278119738165428224000000e34 * t26 * t7 +
             0.3081296464797637251873688780800000e34 * t30 * t49 +
             0.385162058099704656484211097600000e33 * t33 * t6 +
             0.15006313951936545057826406400000e32 * t36 * t2 +
             0.78157885166336172176179200000e29 * t39) /
            0.15451448453276040361200844800000e32);
}



double
AssLegFunction::P_20_4(const double x) const
{
  double t4;
  double t8;
  double t23;
  double t10;
  double t1;
  double t35;
  double t11;
  double t19;
  double t6;
  double t15;
  double t7;
  double t2;
  double t3;
  t1  = sqrt(0.5720330e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t6;
  t8  = t7 * t7;
  t10 = -t3;
  t11 = t6 * t2;
  t15 = t10 * t10;
  t19 = t15 * t10;
  t23 = t15 * t15;
  t35 = t23 * t23;
  return (t1 * t4 *
          (0.197005967779269514605310771200000e33 * t8 +
           0.2364071613351234175263729254400000e34 * t10 * t7 * t11 +
           0.8963771533956762914541640089600000e34 * t15 * t7 * t6 +
           0.14085926696217770294279720140800000e35 * t19 * t7 * t2 +
           0.9904167208278119738165428224000000e34 * t23 * t7 +
           0.3081296464797637251873688780800000e34 * t23 * t10 * t11 +
           0.385162058099704656484211097600000e33 * t23 * t15 * t6 +
           0.15006313951936545057826406400000e32 * t23 * t19 * t2 +
           0.78157885166336172176179200000e29 * t35) /
          0.1050698494822770744561657446400000e34);
}



double
AssLegFunction::P_20_4_Deriv(const double x) const
{
  double t18;
  double t14;
  double t9;
  double t2;
  double t31;
  double t22;
  double t10;
  double t40;
  double t3;
  double t43;
  double t42;
  double t6;
  double t5;
  double t25;
  double t7;
  double t28;
  double t34;
  double t46;
  double t1;
  t1  = sqrt(0.5720330e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t5  = t2 * t2;
  t6  = t5 * t5;
  t7  = t6 * t6;
  t9  = -t3;
  t10 = t5 * t2;
  t14 = t9 * t9;
  t18 = t14 * t9;
  t22 = t14 * t14;
  t25 = t22 * t9;
  t28 = t22 * t14;
  t31 = t22 * t18;
  t34 = t22 * t22;
  t40 = t3 * t3;
  t42 = t2 * x;
  t43 = t5 * t42;
  t46 = t5 * x;
  return (-t1 * t3 *
            (0.197005967779269514605310771200000e33 * t7 +
             0.2364071613351234175263729254400000e34 * t9 * t6 * t10 +
             0.8963771533956762914541640089600000e34 * t14 * t6 * t5 +
             0.14085926696217770294279720140800000e35 * t18 * t6 * t2 +
             0.9904167208278119738165428224000000e34 * t22 * t6 +
             0.3081296464797637251873688780800000e34 * t25 * t10 +
             0.385162058099704656484211097600000e33 * t28 * t5 +
             0.15006313951936545057826406400000e32 * t31 * t2 +
             0.78157885166336172176179200000e29 * t34) *
            x / 0.262674623705692686140414361600000e33 +
          t1 * t40 *
            (0.7880238711170780584212430848000000e34 * t6 * t43 +
             0.68952088722744330111858769920000000e35 * t9 * t6 * t46 +
             0.192080818584787776740178001920000000e36 * t14 * t6 * t42 +
             0.220092604628402660848120627200000000e36 * t18 * t6 * x +
             0.110046302314201330424060313600000000e36 * t22 * t43 +
             0.23109723485982279389052665856000000e35 * t25 * t46 +
             0.1750736627725930256746414080000000e34 * t28 * t42 +
             0.31263154066534468870471680000000e32 * t31 * x) /
            0.1050698494822770744561657446400000e34);
}



double
AssLegFunction::P_20_5(const double x) const
{
  double t9;
  double t1;
  double t5;
  double t15;
  double t14;
  double t11;
  double t27;
  double t4;
  double t23;
  double t19;
  double t3;
  double t2;
  double t8;
  double t10;
  t1  = sqrt(0.5720330e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * x;
  t9  = t2 * t2;
  t10 = t9 * t8;
  t11 = t9 * t9;
  t14 = -t3;
  t15 = t9 * x;
  t19 = t14 * t14;
  t23 = t19 * t14;
  t27 = t19 * t19;
  return (-t1 * t5 * t4 *
          (0.7880238711170780584212430848000000e34 * t11 * t10 +
           0.68952088722744330111858769920000000e35 * t14 * t11 * t15 +
           0.192080818584787776740178001920000000e36 * t19 * t11 * t8 +
           0.220092604628402660848120627200000000e36 * t23 * t11 * x +
           0.110046302314201330424060313600000000e36 * t27 * t10 +
           0.23109723485982279389052665856000000e35 * t27 * t14 * t15 +
           0.1750736627725930256746414080000000e34 * t27 * t19 * t8 +
           0.31263154066534468870471680000000e32 * t27 * t23 * x) /
          0.21013969896455414891233148928000000e35);
}



double
AssLegFunction::P_20_5_Deriv(const double x) const
{
  double t22;
  double t32;
  double t35;
  double t7;
  double t8;
  double t10;
  double t9;
  double t13;
  double t14;
  double t4;
  double t42;
  double t45;
  double t1;
  double t3;
  double t2;
  double t26;
  double t29;
  double t18;
  t1  = sqrt(0.5720330e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = sqrt(t3);
  t7  = t2 * x;
  t8  = t2 * t2;
  t9  = t8 * t7;
  t10 = t8 * t8;
  t13 = -t3;
  t14 = t8 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t26 = t18 * t18;
  t29 = t26 * t13;
  t32 = t26 * t18;
  t35 = t26 * t22;
  t42 = t3 * t3;
  t45 = t8 * t2;
  return (t1 * t4 * t3 *
            (0.7880238711170780584212430848000000e34 * t9 * t10 +
             0.68952088722744330111858769920000000e35 * t13 * t10 * t14 +
             0.192080818584787776740178001920000000e36 * t18 * t10 * t7 +
             0.220092604628402660848120627200000000e36 * t22 * t10 * x +
             0.110046302314201330424060313600000000e36 * t26 * t9 +
             0.23109723485982279389052665856000000e35 * t29 * t14 +
             0.1750736627725930256746414080000000e34 * t32 * t7 +
             0.31263154066534468870471680000000e32 * t35 * x) *
            x / 0.4202793979291082978246629785600000e34 -
          t1 * t4 * t42 *
            (0.256107758113050368986904002560000000e36 * t10 * t45 +
             0.1664700427734827398414876016640000000e37 * t13 * t10 * t8 +
             0.3433444632203081509230681784320000000e37 * t18 * t10 * t2 +
             0.2861203860169234591025568153600000000e37 * t22 * t10 +
             0.1001421351059232106858948853760000000e37 * t26 * t45 +
             0.136557456962622560026220298240000000e36 * t29 * t8 +
             0.5689894040109273334425845760000000e34 * t32 * t2 +
             0.31263154066534468870471680000000e32 * t35) /
            0.21013969896455414891233148928000000e35);
}



double
AssLegFunction::P_20_6(const double x) const
{
  double t4;
  double t16;
  double t1;
  double t8;
  double t23;
  double t9;
  double t2;
  double t12;
  double t7;
  double t20;
  double t3;
  t1  = sqrt(0.22309287e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t2;
  t9  = t7 * t7;
  t12 = -t3;
  t16 = t12 * t12;
  t20 = t16 * t12;
  t23 = t16 * t16;
  return (t1 * t4 * t3 *
          (0.256107758113050368986904002560000000e36 * t9 * t8 +
           0.1664700427734827398414876016640000000e37 * t12 * t9 * t7 +
           0.3433444632203081509230681784320000000e37 * t16 * t9 * t2 +
           0.2861203860169234591025568153600000000e37 * t20 * t9 +
           0.1001421351059232106858948853760000000e37 * t23 * t8 +
           0.136557456962622560026220298240000000e36 * t23 * t12 * t7 +
           0.5689894040109273334425845760000000e34 * t23 * t16 * t2 +
           0.31263154066534468870471680000000e32 * t23 * t20) /
          0.819544825961761180758092808192000000e36);
}



double
AssLegFunction::P_20_6_Deriv(const double x) const
{
  double t6;
  double t7;
  double t8;
  double t28;
  double t42;
  double t1;
  double t11;
  double t22;
  double t39;
  double t2;
  double t19;
  double t25;
  double t3;
  double t4;
  double t15;
  t1  = sqrt(0.22309287e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t2 * t2;
  t7  = t6 * t2;
  t8  = t6 * t6;
  t11 = -t3;
  t15 = t11 * t11;
  t19 = t15 * t11;
  t22 = t15 * t15;
  t25 = t22 * t11;
  t28 = t22 * t15;
  t39 = t6 * x;
  t42 = t2 * x;
  return (-t1 * t4 *
            (0.256107758113050368986904002560000000e36 * t7 * t8 +
             0.1664700427734827398414876016640000000e37 * t11 * t8 * t6 +
             0.3433444632203081509230681784320000000e37 * t15 * t8 * t2 +
             0.2861203860169234591025568153600000000e37 * t19 * t8 +
             0.1001421351059232106858948853760000000e37 * t7 * t22 +
             0.136557456962622560026220298240000000e36 * t25 * t6 +
             0.5689894040109273334425845760000000e34 * t28 * t2 +
             0.31263154066534468870471680000000e32 * t22 * t19) *
            x / 0.136590804326960196793015468032000000e36 +
          t1 * t4 * t3 *
            (0.6914909469052359962646408069120000000e37 * t8 * t39 +
             0.33710183661630254817901239336960000000e38 * t11 * t8 * t42 +
             0.51501669483046222638460226764800000000e38 * t15 * t8 * x +
             0.30901001689827733583076136058880000000e38 * t19 * t6 * t42 +
             0.7374102675981618241415896104960000000e37 * t22 * t39 +
             0.614508556331801520117991342080000000e36 * t25 * t42 +
             0.11817472237150029233038295040000000e35 * t28 * x) /
            0.819544825961761180758092808192000000e36);
}



double
AssLegFunction::P_20_7(const double x) const
{
  double t6;
  double t2;
  double t10;
  double t1;
  double t14;
  double t11;
  double t19;
  double t27;
  double t3;
  double t9;
  double t15;
  double t4;
  t1  = sqrt(0.2124694e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * x;
  t11 = t9 * t9;
  t14 = -t3;
  t15 = t2 * x;
  t19 = t14 * t14;
  t27 = t19 * t19;
  return (-t1 * t6 * t4 * t3 *
          (0.6914909469052359962646408069120000000e37 * t11 * t10 +
           0.33710183661630254817901239336960000000e38 * t14 * t11 * t15 +
           0.51501669483046222638460226764800000000e38 * t19 * t11 * x +
           0.30901001689827733583076136058880000000e38 * t19 * t14 * t9 * t15 +
           0.7374102675981618241415896104960000000e37 * t27 * t10 +
           0.614508556331801520117991342080000000e36 * t27 * t14 * t15 +
           0.11817472237150029233038295040000000e35 * t27 * t19 * x) /
          0.4917268955770567084548556849152000000e37);
}



double
AssLegFunction::P_20_7_Deriv(const double x) const
{
  double t3;
  double t4;
  double t5;
  double t1;
  double t2;
  double t14;
  double t13;
  double t29;
  double t8;
  double t9;
  double t32;
  double t26;
  double t10;
  double t18;
  double t22;
  t1  = sqrt(0.2124694e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = sqrt(t3);
  t8  = t2 * t2;
  t9  = t8 * x;
  t10 = t8 * t8;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t26 = t18 * t18;
  t29 = t26 * t13;
  t32 = t26 * t18;
  return (t1 * t5 * t4 *
            (0.6914909469052359962646408069120000000e37 * t10 * t9 +
             0.33710183661630254817901239336960000000e38 * t13 * t10 * t14 +
             0.51501669483046222638460226764800000000e38 * t18 * t10 * x +
             0.30901001689827733583076136058880000000e38 * t22 * t8 * t14 +
             0.7374102675981618241415896104960000000e37 * t26 * t9 +
             0.614508556331801520117991342080000000e36 * t29 * t14 +
             0.11817472237150029233038295040000000e35 * t32 * x) *
            x / 0.702466993681509583506936692736000000e36 -
          t1 * t5 * t4 * t3 *
            (0.157314190420941189150205783572480000000e39 * t10 * t8 +
             0.576818698210117693550754539765760000000e39 * t13 * t10 * t2 +
             0.648921035486382405244598857236480000000e39 * t18 * t10 +
             0.275299833236647081012860121251840000000e39 * t22 * t8 * t2 +
             0.43015598943226106408259393945600000000e38 * t26 * t8 +
             0.1985335335841204911150433566720000000e37 * t29 * t2 +
             0.11817472237150029233038295040000000e35 * t32) /
            0.4917268955770567084548556849152000000e37);
}



double
AssLegFunction::P_20_8(const double x) const
{
  double t5;
  double t15;
  double t4;
  double t7;
  double t11;
  double t1;
  double t2;
  double t8;
  double t22;
  double t3;
  t1  = sqrt(0.1144066e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t22 = t15 * t15;
  return (t1 * t5 *
          (0.157314190420941189150205783572480000000e39 * t8 * t7 +
           0.576818698210117693550754539765760000000e39 * t11 * t8 * t2 +
           0.648921035486382405244598857236480000000e39 * t15 * t8 +
           0.275299833236647081012860121251840000000e39 * t15 * t11 * t7 * t2 +
           0.43015598943226106408259393945600000000e38 * t22 * t7 +
           0.1985335335841204911150433566720000000e37 * t22 * t11 * t2 +
           0.11817472237150029233038295040000000e35 * t22 * t15) /
          0.68841765380787939183679795888128000000e38);
}



double
AssLegFunction::P_20_8_Deriv(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t15;
  double t22;
  double t7;
  double t8;
  double t36;
  double t34;
  double t25;
  double t11;
  double t1;
  double t18;
  t1  = sqrt(0.1144066e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t22 = t15 * t15;
  t25 = t22 * t11;
  t34 = t4 * t4;
  t36 = t2 * x;
  return (-t1 * t4 * t3 *
            (0.157314190420941189150205783572480000000e39 * t8 * t7 +
             0.576818698210117693550754539765760000000e39 * t11 * t8 * t2 +
             0.648921035486382405244598857236480000000e39 * t15 * t8 +
             0.275299833236647081012860121251840000000e39 * t18 * t7 * t2 +
             0.43015598943226106408259393945600000000e38 * t22 * t7 +
             0.1985335335841204911150433566720000000e37 * t25 * t2 +
             0.11817472237150029233038295040000000e35 * t22 * t15) *
            x / 0.8605220672598492397959974486016000000e37 +
          t1 * t34 *
            (0.3041407681471529656903978482401280000000e40 * t8 * t36 +
             0.8363871124046706556485940826603520000000e40 * t11 * t8 * x +
             0.6843167283310941728033951585402880000000e40 * t15 * t7 * t36 +
             0.1995923790965691337343235879075840000000e40 * t18 * t7 * x +
             0.191915749131316474744541911449600000000e39 * t22 * t36 +
             0.4112480338528210173097326673920000000e37 * t25 * x) /
            0.68841765380787939183679795888128000000e38);
}



double
AssLegFunction::P_20_9(const double x) const
{
  double t6;
  double t5;
  double t2;
  double t3;
  double t4;
  double t9;
  double t10;
  double t11;
  double t1;
  double t26;
  double t18;
  double t14;
  t1  = sqrt(0.99533742e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t11 = t10 * t10;
  t14 = -t3;
  t18 = t14 * t14;
  t26 = t18 * t18;
  return (-t1 * t6 * t5 *
          (0.3041407681471529656903978482401280000000e40 * t11 * t9 +
           0.8363871124046706556485940826603520000000e40 * t14 * t11 * x +
           0.6843167283310941728033951585402880000000e40 * t18 * t10 * t9 +
           0.1995923790965691337343235879075840000000e40 * t18 * t14 * t10 * x +
           0.191915749131316474744541911449600000000e39 * t26 * t9 +
           0.4112480338528210173097326673920000000e37 * t26 * t14 * x) /
          0.11978467176257101417960284484534272000000e41);
}



double
AssLegFunction::P_20_9_Deriv(const double x) const
{
  double t14;
  double t1;
  double t18;
  double t2;
  double t3;
  double t4;
  double t22;
  double t6;
  double t36;
  double t9;
  double t26;
  double t10;
  double t29;
  double t11;
  t1  = sqrt(0.99533742e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = sqrt(t3);
  t9  = t2 * x;
  t10 = t2 * t2;
  t11 = t10 * t10;
  t14 = -t3;
  t18 = t14 * t14;
  t22 = t18 * t14;
  t26 = t18 * t18;
  t29 = t26 * t14;
  t36 = t4 * t4;
  return (t1 * t6 * t4 * t3 *
            (0.3041407681471529656903978482401280000000e40 * t11 * t9 +
             0.8363871124046706556485940826603520000000e40 * t14 * t11 * x +
             0.6843167283310941728033951585402880000000e40 * t18 * t10 * t9 +
             0.1995923790965691337343235879075840000000e40 * t22 * t10 * x +
             0.191915749131316474744541911449600000000e39 * t9 * t26 +
             0.4112480338528210173097326673920000000e37 * t29 * x) *
            x / 0.1330940797361900157551142720503808000000e40 -
          t1 * t6 * t36 *
            (0.50183226744280239338915644959621120000000e41 * t11 * t2 +
             0.102647509249664125920509273781043200000000e42 * t14 * t11 +
             0.59877713728970740120297076372275200000000e41 * t18 * t10 * t2 +
             0.11514944947878988484672514686976000000000e41 * t22 * t10 +
             0.616872050779231525964599001088000000000e39 * t26 * t2 +
             0.4112480338528210173097326673920000000e37 * t29) /
            0.11978467176257101417960284484534272000000e41);
}



double
AssLegFunction::P_20_10(const double x) const
{
  double t1;
  double t15;
  double t12;
  double t22;
  double t5;
  double t8;
  double t9;
  double t2;
  double t4;
  double t3;
  t1  = sqrt(0.7540435e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t12 = -t3;
  t15 = t12 * t12;
  t22 = t15 * t15;
  return (t1 * t5 * t3 *
          (0.50183226744280239338915644959621120000000e41 * t9 * t2 +
           0.102647509249664125920509273781043200000000e42 * t12 * t9 +
           0.59877713728970740120297076372275200000000e41 * t15 * t8 * t2 +
           0.11514944947878988484672514686976000000000e41 * t15 * t12 * t8 +
           0.616872050779231525964599001088000000000e39 * t22 * t2 +
           0.4112480338528210173097326673920000000e37 * t22 * t12) /
          0.59892335881285507089801422422671360000000e41);
}



double
AssLegFunction::P_20_10_Deriv(const double x) const
{
  double t1;
  double t14;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t21;
  double t8;
  double t34;
  double t11;
  t1  = sqrt(0.7540435e7);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = t2 * t2;
  t8  = t7 * t7;
  t11 = -t3;
  t14 = t11 * t11;
  t18 = t14 * t11;
  t21 = t14 * t14;
  t34 = t2 * x;
  return (-t1 * t5 *
            (0.50183226744280239338915644959621120000000e41 * t8 * t2 +
             0.102647509249664125920509273781043200000000e42 * t11 * t8 +
             0.59877713728970740120297076372275200000000e41 * t14 * t7 * t2 +
             0.11514944947878988484672514686976000000000e41 * t18 * t7 +
             0.616872050779231525964599001088000000000e39 * t21 * t2 +
             0.4112480338528210173097326673920000000e37 * t21 * t11) *
            x / 0.5989233588128550708980142242267136000000e40 +
          t1 * t5 * t3 *
            (0.707127285942130645230174997158297600000000e42 * t8 * x +
             0.1060690928913195967845262495737446400000000e43 * t11 * t7 * t34 +
             0.428355952061098371629817546355507200000000e42 * t14 * t7 * x +
             0.50994756197749806146406850756608000000000e41 * t18 * t34 +
             0.1274868904943745153660171268915200000000e40 * t21 * x) /
            0.59892335881285507089801422422671360000000e41);
}



double
AssLegFunction::P_20_11(const double x) const
{
  double t14;
  double t15;
  double t10;
  double t11;
  double t1;
  double t19;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t26;
  t1  = sqrt(0.93501394e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t11 = t10 * t10;
  t14 = -t3;
  t15 = t2 * x;
  t19 = t14 * t14;
  t26 = t19 * t19;
  return (-t1 * t7 * t5 * t3 *
          (0.707127285942130645230174997158297600000000e42 * t11 * x +
           0.1060690928913195967845262495737446400000000e43 * t14 * t10 * t15 +
           0.428355952061098371629817546355507200000000e42 * t19 * t10 * x +
           0.50994756197749806146406850756608000000000e41 * t19 * t14 * t15 +
           0.1274868904943745153660171268915200000000e40 * t26 * x) /
          0.3713324824639701439567688190205624320000000e43);
}



double
AssLegFunction::P_20_11_Deriv(const double x) const
{
  double t2;
  double t13;
  double t14;
  double t3;
  double t25;
  double t1;
  double t4;
  double t9;
  double t5;
  double t6;
  double t22;
  double t10;
  double t18;
  t1  = sqrt(0.93501394e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = sqrt(t3);
  t9  = t2 * t2;
  t10 = t9 * t9;
  t13 = -t3;
  t14 = t2 * x;
  t18 = t13 * t13;
  t22 = t18 * t13;
  t25 = t18 * t18;
  return (t1 * t6 * t5 *
            (0.707127285942130645230174997158297600000000e42 * t10 * x +
             0.1060690928913195967845262495737446400000000e43 * t13 * t9 * t14 +
             0.428355952061098371629817546355507200000000e42 * t18 * t9 * x +
             0.50994756197749806146406850756608000000000e41 * t22 * t14 +
             0.1274868904943745153660171268915200000000e40 * t25 * x) *
            x / 0.337574984058154676324335290018693120000000e42 -
          t1 * t6 * t5 * t3 *
            (0.8485527431305567742762099965899571200000000e43 * t10 +
             0.9138260310636765261436107655584153600000000e43 * t13 * t9 * t2 +
             0.2447748297491990695027528836317184000000000e43 * t18 * t9 +
             0.163183219832799379668501922421145600000000e42 * t22 * t2 +
             0.1274868904943745153660171268915200000000e40 * t25) /
            0.3713324824639701439567688190205624320000000e43);
}



double
AssLegFunction::P_20_12(const double x) const
{
  double t21;
  double t8;
  double t9;
  double t11;
  double t1;
  double t15;
  double t3;
  double t2;
  double t4;
  double t5;
  t1  = sqrt(0.46750697e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t11 = -t3;
  t15 = t11 * t11;
  t21 = t15 * t15;
  return (t1 * t5 * t4 *
          (0.8485527431305567742762099965899571200000000e43 * t9 +
           0.9138260310636765261436107655584153600000000e43 * t11 * t8 * t2 +
           0.2447748297491990695027528836317184000000000e43 * t15 * t8 +
           0.163183219832799379668501922421145600000000e42 * t15 * t11 * t2 +
           0.1274868904943745153660171268915200000000e40 * t21) /
          0.44559897895676417274812258282467491840000000e44);
}



double
AssLegFunction::P_20_12_Deriv(const double x) const
{
  double t11;
  double t15;
  double t1;
  double t18;
  double t2;
  double t29;
  double t21;
  double t8;
  double t9;
  double t3;
  double t5;
  double t4;
  t1  = sqrt(0.46750697e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t9  = t8 * t8;
  t11 = -t3;
  t15 = t11 * t11;
  t18 = t15 * t11;
  t21 = t15 * t15;
  t29 = t2 * x;
  return (-t1 * t5 * t3 *
            (0.8485527431305567742762099965899571200000000e43 * t9 +
             0.9138260310636765261436107655584153600000000e43 * t11 * t8 * t2 +
             0.2447748297491990695027528836317184000000000e43 * t15 * t8 +
             0.163183219832799379668501922421145600000000e42 * t18 * t2 +
             0.1274868904943745153660171268915200000000e40 * t21) *
            x / 0.3713324824639701439567688190205624320000000e43 +
          t1 * t5 * t4 *
            (0.86160740071718072464969015038364876800000000e44 * t8 * t29 +
             0.64620555053788554348726761278773657600000000e44 * t11 * t8 * x +
             0.10770092508964759058121126879795609600000000e44 * t29 * t15 +
             0.336565390905148720566285214993612800000000e42 * t18 * x) /
            0.44559897895676417274812258282467491840000000e44);
}



double
AssLegFunction::P_20_13(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t10;
  double t11;
  double t14;
  double t18;
  t1  = sqrt(0.3085546002e10);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * x;
  t11 = t2 * t2;
  t14 = -t3;
  t18 = t14 * t14;
  return (-t1 * t7 * t5 * t4 *
          (0.86160740071718072464969015038364876800000000e44 * t11 * t10 +
           0.64620555053788554348726761278773657600000000e44 * t14 * t11 * x +
           0.10770092508964759058121126879795609600000000e44 * t18 * t10 +
           0.336565390905148720566285214993612800000000e42 * t18 * t14 * x) /
          0.5881906522229287080275218093285708922880000000e46);
}



double
AssLegFunction::P_20_13_Deriv(const double x) const
{
  double t4;
  double t5;
  double t1;
  double t14;
  double t18;
  double t7;
  double t21;
  double t2;
  double t3;
  double t11;
  double t10;
  t1  = sqrt(0.3085546002e10);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * x;
  t11 = t2 * t2;
  t14 = -t3;
  t18 = t14 * t14;
  t21 = t18 * t14;
  return (t1 * t7 * t5 * t3 *
            (0.86160740071718072464969015038364876800000000e44 * t11 * t10 +
             0.64620555053788554348726761278773657600000000e44 * t14 * t11 * x +
             0.10770092508964759058121126879795609600000000e44 * t18 * t10 +
             0.336565390905148720566285214993612800000000e42 * t21 * x) *
            x / 0.452454347863791313867324468714285301760000000e45 -
          t1 * t7 * t5 * t4 *
            (0.732366290609603615952236627826101452800000000e45 * t11 * t2 +
             0.366183145304801807976118313913050726400000000e45 * t14 * t11 +
             0.34329669872325169497761091929348505600000000e44 * t18 * t2 +
             0.336565390905148720566285214993612800000000e42 * t21) /
            0.5881906522229287080275218093285708922880000000e46);
}



double
AssLegFunction::P_20_14(const double x) const
{
  double t2;
  double t3;
  double t4;
  double t9;
  double t6;
  double t15;
  double t12;
  double t1;
  t1  = sqrt(0.12964479e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t9  = t2 * t2;
  t12 = -t3;
  t15 = t12 * t12;
  return (t1 * t6 * t4 * t3 *
          (0.732366290609603615952236627826101452800000000e45 * t9 * t2 +
           0.366183145304801807976118313913050726400000000e45 * t9 * t12 +
           0.34329669872325169497761091929348505600000000e44 * t15 * t2 +
           0.336565390905148720566285214993612800000000e42 * t12 * t15) /
          0.5881906522229287080275218093285708922880000000e46);
}



double
AssLegFunction::P_20_14_Deriv(const double x) const
{
  double t1;
  double t14;
  double t8;
  double t11;
  double t2;
  double t3;
  double t5;
  double t4;
  t1  = sqrt(0.12964479e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t8  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  return (
    -t1 * t5 * t4 *
      (0.732366290609603615952236627826101452800000000e45 * t8 * t2 +
       0.366183145304801807976118313913050726400000000e45 * t11 * t8 +
       0.34329669872325169497761091929348505600000000e44 * t14 * t2 +
       0.336565390905148720566285214993612800000000e42 * t14 * t11) *
      x / 0.420136180159234791448229863806122065920000000e45 +
    t1 * t5 * t4 * t3 *
      (0.5126564034267225311665656394782710169600000000e46 * t8 * x +
       0.1602051260708507909895517623369596928000000000e46 * t11 * t2 * x +
       0.70678732090081231318919895148658688000000000e44 * t14 * x) /
      0.5881906522229287080275218093285708922880000000e46);
}



double
AssLegFunction::P_20_15(const double x) const
{
  double t8;
  double t11;
  double t18;
  double t1;
  double t2;
  double t14;
  double t3;
  double t4;
  double t6;
  t1  = sqrt(0.302504510e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t8  = sqrt(t3);
  t11 = t2 * t2;
  t14 = -t3;
  t18 = t14 * t14;
  return (-t1 * t8 * t6 * t4 * t3 *
          (0.5126564034267225311665656394782710169600000000e46 * t11 * x +
           0.1602051260708507909895517623369596928000000000e46 * t14 * t2 * x +
           0.70678732090081231318919895148658688000000000e44 * t18 * x) /
          0.411733456556050095619265266529999624601600000000e48);
}



double
AssLegFunction::P_20_15_Deriv(const double x) const
{
  double t4;
  double t13;
  double t5;
  double t7;
  double t17;
  double t1;
  double t2;
  double t3;
  double t10;
  t1  = sqrt(0.302504510e9);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t7  = sqrt(t3);
  t10 = t2 * t2;
  t13 = -t3;
  t17 = t13 * t13;
  return (
    t1 * t7 * t5 * t4 *
      (0.5126564034267225311665656394782710169600000000e46 * t10 * x +
       0.1602051260708507909895517623369596928000000000e46 * t13 * t2 * x +
       0.70678732090081231318919895148658688000000000e44 * t17 * x) *
      x / 0.27448897103736673041284351101999974973440000000e47 -
    t1 * t7 * t5 * t4 * t3 *
      (0.28836922692753142378119317220652744704000000000e47 * t10 +
       0.5088868710485848654962232450703425536000000000e46 * t13 * t2 +
       0.70678732090081231318919895148658688000000000e44 * t17) /
      0.411733456556050095619265266529999624601600000000e48);
}



double
AssLegFunction::P_20_16(const double x) const
{
  double t8;
  double t10;
  double t5;
  double t2;
  double t3;
  double t4;
  double t1;
  double t6;
  double t13;
  t1  = sqrt(0.60500902e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t5  = t4 * t4;
  t6  = t5 * t5;
  t8  = t2 * t2;
  t10 = -t3;
  t13 = t10 * t10;
  return (t1 * t6 *
          (0.28836922692753142378119317220652744704000000000e47 * t8 +
           0.5088868710485848654962232450703425536000000000e46 * t10 * t2 +
           0.70678732090081231318919895148658688000000000e44 * t13) /
          0.2470400739336300573715591599179997747609600000000e49);
}



double
AssLegFunction::P_20_16_Deriv(const double x) const
{
  double t1;
  double t14;
  double t20;
  double t2;
  double t3;
  double t4;
  double t6;
  double t9;
  double t11;
  t1  = sqrt(0.60500902e8);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t9  = t2 * t2;
  t11 = -t3;
  t14 = t11 * t11;
  t20 = t6 * t6;
  return (-t1 * t6 * t4 * t3 *
            (0.28836922692753142378119317220652744704000000000e47 * t9 +
             0.5088868710485848654962232450703425536000000000e46 * t11 * t2 +
             0.70678732090081231318919895148658688000000000e44 * t14) *
            x / 0.154400046208518785857224474948749859225600000000e48 +
          t1 * t20 *
            (0.125525428191984266822401733784017829888000000000e48 * t2 * x +
             0.10460452349332022235200144482001485824000000000e47 * t11 * x) /
            0.2470400739336300573715591599179997747609600000000e49);
}



double
AssLegFunction::P_20_17(const double x) const
{
  double t7;
  double t1;
  double t6;
  double t2;
  double t5;
  double t3;
  double t4;
  t1 = sqrt(0.2238533374e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (-t1 * t7 * t6 *
          (0.125525428191984266822401733784017829888000000000e48 * t2 * x -
           0.10460452349332022235200144482001485824000000000e47 * t3 * x) /
          0.182809654710886242454953778339319833323110400000000e51);
}



double
AssLegFunction::P_20_17_Deriv(const double x) const
{
  double t3;
  double t4;
  double t6;
  double t1;
  double t20;
  double t8;
  double t2;
  t1  = sqrt(0.2238533374e10);
  t2  = x * x;
  t3  = 0.1e1 - t2;
  t4  = t3 * t3;
  t6  = t4 * t4;
  t8  = sqrt(t3);
  t20 = t6 * t6;
  return (t1 * t8 * t6 * t4 * t3 *
            (0.125525428191984266822401733784017829888000000000e48 * t2 * x -
             0.10460452349332022235200144482001485824000000000e47 * t3 * x) *
            x / 0.10753509100640367203232575196430578430771200000000e50 -
          t1 * t8 * t20 *
            (0.407957641623948867172805634798057947136000000000e48 * t2 -
             0.10460452349332022235200144482001485824000000000e47) /
            0.182809654710886242454953778339319833323110400000000e51);
}



double
AssLegFunction::P_20_18(const double x) const
{
  double t2;
  double t4;
  double t3;
  double t1;
  double t5;
  double t6;
  t1 = sqrt(0.176726319e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (t1 * t6 * t3 *
          (0.407957641623948867172805634798057947136000000000e48 * t2 -
           0.10460452349332022235200144482001485824000000000e47) /
          0.548428964132658727364861335017959499969331200000000e51);
}



double
AssLegFunction::P_20_18_Deriv(const double x) const
{
  double t4;
  double t2;
  double t3;
  double t5;
  double t6;
  double t1;
  t1 = sqrt(0.176726319e9);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (-t1 * t6 *
            (0.407957641623948867172805634798057947136000000000e48 * t2 -
             0.10460452349332022235200144482001485824000000000e47) *
            x / 0.30468275785147707075825629723219972220518400000000e50 +
          0.195e3 / 0.131072e6 * t1 * t6 * t3 * x);
}



double
AssLegFunction::P_20_19(const double x) const
{
  double t5;
  double t6;
  double t8;
  double t1;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.1531628098e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t8 = sqrt(t3);
  return (-0.15e2 / 0.262144e6 * t1 * t8 * t6 * t3 * x);
}



double
AssLegFunction::P_20_19_Deriv(const double x) const
{
  double t1;
  double t5;
  double t7;
  double t6;
  double t2;
  double t3;
  double t4;
  t1 = sqrt(0.1531628098e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  t7 = sqrt(t3);
  return (0.285e3 / 0.262144e6 * t1 * t7 * t6 * t2 -
          0.15e2 / 0.262144e6 * t1 * t7 * t6 * t3);
}



double
AssLegFunction::P_20_20(const double x) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t6;
  t1 = sqrt(0.3829070245e10);
  t2 = x * x;
  t4 = pow(0.1e1 - t2, 0.2e1);
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (0.3e1 / 0.524288e6 * t1 * t6 * t4);
}



double
AssLegFunction::P_20_20_Deriv(const double x) const
{
  double t1;
  double t2;
  double t3;
  double t6;
  double t4;
  double t5;
  t1 = sqrt(0.3829070245e10);
  t2 = x * x;
  t3 = 0.1e1 - t2;
  t4 = t3 * t3;
  t5 = t4 * t4;
  t6 = t5 * t5;
  return (-0.15e2 / 0.131072e6 * t1 * t6 * t3 * x);
}
