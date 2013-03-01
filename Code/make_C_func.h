void make_C_block1(unsigned int r,
		      const DiffusionEncodingDirection &n,
		      vnl_matrix<double> &C)
{
    C[r][0] = -1*n.gx*n.gx;
    C[r][1] = -2*n.gx*n.gy;
    C[r][2] = -2*n.gx*n.gz;
    C[r][3] = -1*n.gy*n.gy;
    C[r][4] = -2*n.gy*n.gz;
    C[r][5] = -1*n.gz*n.gz;
    C[r][6] = 0;
    C[r][7] = 0;
    C[r][8] = 0;
    C[r][9] = 0;
    C[r][10] = 0;
    C[r][11] = 0;
    C[r][12] = 0;
    C[r][13] = 0;
    C[r][14] = 0;
    C[r][15] = 0;
    C[r][16] = 0;
    C[r][17] = 0;
    C[r][18] = 0;
    C[r][19] = 0;
    C[r][20] = 0;
}
void make_C_block2(unsigned int r,
		      const DiffusionEncodingDirection &n,
		      vnl_matrix<double> &C)
{
    C[r][0] = 0;
    C[r][1] = 0;
    C[r][2] = 0;
    C[r][3] = 0;
    C[r][4] = 0;
    C[r][5] = 0;
    C[r][6] = -1*n.gx*n.gx*n.gx*n.gx;
    C[r][7] = -4*n.gx*n.gx*n.gx*n.gy;
    C[r][8] = -4*n.gx*n.gx*n.gx*n.gz;
    C[r][9] = -6*n.gx*n.gx*n.gy*n.gy;
    C[r][10] = -12*n.gx*n.gx*n.gy*n.gz;
    C[r][11] = -6*n.gx*n.gx*n.gz*n.gz;
    C[r][12] = -4*n.gx*n.gy*n.gy*n.gy;
    C[r][13] = -12*n.gx*n.gy*n.gy*n.gz;
    C[r][14] = -12*n.gx*n.gy*n.gz*n.gz;
    C[r][15] = -4*n.gx*n.gz*n.gz*n.gz;
    C[r][16] = -1*n.gy*n.gy*n.gy*n.gy;
    C[r][17] = -4*n.gy*n.gy*n.gy*n.gz;
    C[r][18] = -6*n.gy*n.gy*n.gz*n.gz;
    C[r][19] = -4*n.gy*n.gz*n.gz*n.gz;
    C[r][20] = -1*n.gz*n.gz*n.gz*n.gz;
}
void make_C_block3(unsigned int r,
				const DiffusionEncodingDirection &n,
				const double bmax,
				vnl_matrix<double> &C)
{
    C[r][0] = (-3/bmax)*1*n.gx*n.gx;
    C[r][1] = (-3/bmax)*2*n.gx*n.gy;
    C[r][2] = (-3/bmax)*2*n.gx*n.gz;
    C[r][3] = (-3/bmax)*1*n.gy*n.gy;
    C[r][4] = (-3/bmax)*2*n.gy*n.gz;
    C[r][5] = (-3/bmax)*1*n.gz*n.gz;
    C[r][6] = 1*n.gx*n.gx*n.gx*n.gx;
    C[r][7] = 4*n.gx*n.gx*n.gx*n.gy;
    C[r][8] = 4*n.gx*n.gx*n.gx*n.gz;
    C[r][9] = 6*n.gx*n.gx*n.gy*n.gy;
    C[r][10] = 12*n.gx*n.gx*n.gy*n.gz;
    C[r][11] = 6*n.gx*n.gx*n.gz*n.gz;
    C[r][12] = 4*n.gx*n.gy*n.gy*n.gy;
    C[r][13] = 12*n.gx*n.gy*n.gy*n.gz;
    C[r][14] = 12*n.gx*n.gy*n.gz*n.gz;
    C[r][15] = 4*n.gx*n.gz*n.gz*n.gz;
    C[r][16] = 1*n.gy*n.gy*n.gy*n.gy;
    C[r][17] = 4*n.gy*n.gy*n.gy*n.gz;
    C[r][18] = 6*n.gy*n.gy*n.gz*n.gz;
    C[r][19] = 4*n.gy*n.gz*n.gz*n.gz;
    C[r][20] = 1*n.gz*n.gz*n.gz*n.gz;
}
