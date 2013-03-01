void make_A(unsigned int r, const DiffusionEncodingDirection &n, vnl_matrix<double> &A)
{
    A[r][0] = -n.bvalue*1*n.gx*n.gx;
    A[r][1] = -n.bvalue*2*n.gx*n.gy;
    A[r][2] = -n.bvalue*2*n.gx*n.gz;
    A[r][3] = -n.bvalue*1*n.gy*n.gy;
    A[r][4] = -n.bvalue*2*n.gy*n.gz;
    A[r][5] = -n.bvalue*1*n.gz*n.gz;
    A[r][6] = (n.bvalue*n.bvalue/6)*1*n.gx*n.gx*n.gx*n.gx;
    A[r][7] = (n.bvalue*n.bvalue/6)*4*n.gx*n.gx*n.gx*n.gy;
    A[r][8] = (n.bvalue*n.bvalue/6)*4*n.gx*n.gx*n.gx*n.gz;
    A[r][9] = (n.bvalue*n.bvalue/6)*6*n.gx*n.gx*n.gy*n.gy;
    A[r][10] = (n.bvalue*n.bvalue/6)*12*n.gx*n.gx*n.gy*n.gz;
    A[r][11] = (n.bvalue*n.bvalue/6)*6*n.gx*n.gx*n.gz*n.gz;
    A[r][12] = (n.bvalue*n.bvalue/6)*4*n.gx*n.gy*n.gy*n.gy;
    A[r][13] = (n.bvalue*n.bvalue/6)*12*n.gx*n.gy*n.gy*n.gz;
    A[r][14] = (n.bvalue*n.bvalue/6)*12*n.gx*n.gy*n.gz*n.gz;
    A[r][15] = (n.bvalue*n.bvalue/6)*4*n.gx*n.gz*n.gz*n.gz;
    A[r][16] = (n.bvalue*n.bvalue/6)*1*n.gy*n.gy*n.gy*n.gy;
    A[r][17] = (n.bvalue*n.bvalue/6)*4*n.gy*n.gy*n.gy*n.gz;
    A[r][18] = (n.bvalue*n.bvalue/6)*6*n.gy*n.gy*n.gz*n.gz;
    A[r][19] = (n.bvalue*n.bvalue/6)*4*n.gy*n.gz*n.gz*n.gz;
    A[r][20] = (n.bvalue*n.bvalue/6)*1*n.gz*n.gz*n.gz*n.gz;
}
