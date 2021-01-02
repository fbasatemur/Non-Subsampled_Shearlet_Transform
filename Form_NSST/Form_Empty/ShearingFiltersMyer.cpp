#include "ShearingFiltersMyer.h"
#include "MatlabFuncs.h"
#include <math.h>

Matrix* ShearingFiltersMyer(int n, int level)
{
    //[x11, y11, x12, y12, F1] = gen_x_y_cordinates(n1);

    Matrix* wf;// = Windowing(ones(2 * n, 1), pow(2, level));

    Matrix* w_s = new Matrix;
    w_s->CreateMatrix(n, n, pow(2, level));
    w_s->mat = zeros(n, n, pow(2, level));

    double *one = ones(n, 1);
    Matrix* temp;
    for (int i = 0; i < pow(2, level); i++)
    {
        temp = MatrixCut(wf->mat, wf->height, wf->width, 0, wf->height, 0, i);
        for (int i = 0; i < n; i++)
            temp->mat[i] *= one[i];

        //w_s(:, : , k) = rec_from_pol(temp, n1, x11, y11, x12, y12, F1);% convert window array into Cartesian coord.
        //w_s(:, : , k) = real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);
    }

    return w_s;
}
