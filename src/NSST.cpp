#include <iostream>
#include <memory>
#include "NSST.h"
#include "fftwf_wrapper.h"
#include "Circ_Shift.h"

/// <summary>
///		ATROUSFILTERS	Generate pyramid 2D filters
///		'maxflat':		Filters derived from 1-D using maximally flat mapping function with 4 vanishing moments
///		calculate h0, h1, g0, g1:	pyramid filters for 2-D nonsubsampled filter bank (lowpass and highpass)
/// </summary>
void NSST::AtrousFilters() {

    constexpr int h0Width = 7, h0Height = 7;
    constexpr int h1Width = 10, h1Height = 10;
    constexpr int g0Width = 10, g0Height = 10;
    constexpr int g1Width = 7, g1Height = 7;

    if (!strcmp(Conf::lpfilt, "maxflat")) {

        const Eigen::Matrix<float, h0Height, h0Width> h0{
                {-7.900496718847182e-07,   0., 0.000014220894093924927, 0.000025281589500310983, -0.000049773129328737247, -0.00022753430550279883, -0.00033182086219158167},
                {0,                        0,  0,                       0,                       0,                        0,                       0},
                {0.000014220894093924927,  0., -0.0002559760936906487,  -0.00045506861100559767, 0.0008959163279172705,    0.004095617499050379,    0.00597277551944847},
                {0.000025281589500310983,  0., -0.00045506861100559767, 0.0009765625,            0.0015927401385195919,    -0.0087890625,           -0.01795090623402861},
                {-0.000049773129328737247, 0., 0.0008959163279172705,   0.0015927401385195919,   -0.0031357071477104465,   -0.014334661246676327,   -0.020904714318069645},
                {-0.00022753430550279883,  0., 0.004095617499050379,    -0.0087890625,           -0.014334661246676327,    0.0791015625,            0.16155815610625748},
                {-0.00033182086219158167,  0., 0.00597277551944847,     -0.01795090623402861,    -0.020904714318069645,    0.16155815610625748,     0.3177420190660832}
        };

        const Eigen::Matrix<float, g0Height, g0Width> g0{
                {-6.391587676622346e-010, 0., 1.7257286726880333e-08,   3.067962084778726e-08,   -1.3805829381504267e-07,  -5.522331752601707e-07,
                                                                                                                               -3.3747582932565985e-07, 1.9328161134105974e-06,   5.6949046198705095e-06,  7.649452131381623e-06},
                {0.,                      0., 0.,                       0.,                      0.,                       0., 0.,                      0.,                       0.,                      0.},
                {1.7257286726880333e-08,  0., -4.65946741625769e-07,    -8.283497628902559e-07,  3.727573933006152e-06,    0.000014910295732024608,
                                                                                                                               9.111847391792816e-06,   -0.000052186035062086126, -0.00015376242473650378, -0.00020653520754730382},
                {3.067962084778726e-08,   0., -8.283497628902559e-07,   -1.2809236054493144e-06, 6.6267981031220475e-06,   0.00002305662489808766,
                                                                                                                               0.000010064497559808503, -0.0000806981871433068,   -0.00021814634152337594, -0.00028666046030363884},
                {-1.3805829381504267e-07, 0., 3.727573933006152e-06,    6.6267981031220475e-06,  -0.000029820591464049215, -0.00011928236585619686,
                                                                                                                               -0.00007289477913434253, 0.000417488280496689,     0.0012300993978920302,   0.0016522816603784306},
                {-5.522331752601707e-07,  0., 0.000014910295732024608,  0.00002305662489808766,  -0.00011928236585619686,  -0.00041501924816557786,
                                                                                                                               -0.00018116095607655303, 0.0014525673685795225,    0.0039266341474207675,   0.005159888285465499},
                {-3.3747582932565985e-07, 0., 9.111847391792816e-06,    0.000010064497559808503, -0.00007289477913434253,  -0.00018116095607655303,
                                                                                                                               0.001468581806076247,    0.0006340633462679356,    -0.01181401175635013,    -0.021745034491193898},
                {1.9328161134105974e-06,  0., -0.000052186035062086126, -0.0000806981871433068,  0.000417488280496689,     0.0014525673685795225,
                                                                                                                               0.0006340633462679356,   -0.005083985790028328,    -0.013743219515972684,   -0.018059608999129246},
                {5.6949046198705095e-06,  0., -0.00015376242473650378,  -0.00021814634152337594, 0.0012300993978920302,    0.0039266341474207675,
                                                                                                                               -0.01181401175635013,    -0.013743219515972684,    0.0826466923977296,      0.1638988884584603},
                {7.649452131381623e-06,   0., -0.00020653520754730382,  -0.00028666046030363884, 0.0016522816603784306,    0.005159888285465499,
                                                                                                                               -0.021745034491193898,   -0.018059608999129246,    0.1638988884584603,      0.31358726209239235}
        };

        const Eigen::Matrix<float, g1Height, g1Width> g1{
                {-7.900496718847182e-07,   0., 0.000014220894093924927, 0.000025281589500310983, -0.000049773129328737247, -0.00022753430550279883, -0.00033182086219158167},
                {0,                        0,  0,                       0,                       0,                        0,                       0},
                {0.000014220894093924927,  0., -0.0002559760936906487,  -0.00045506861100559767, 0.0008959163279172705,    0.004095617499050379,    0.00597277551944847},
                {0.000025281589500310983,  0., -0.00045506861100559767, -0.0009765625,           0.0015927401385195919,    0.0087890625,            0.01329909376597139},
                {-0.000049773129328737247, 0., 0.0008959163279172705,   0.0015927401385195919,   -0.0031357071477104465,   -0.014334661246676327,   -0.020904714318069645},
                {-0.00022753430550279883,  0., 0.004095617499050379,    0.0087890625,            -0.014334661246676327,    -0.0791015625,           -0.1196918438937425},
                {-0.00033182086219158167,  0., 0.00597277551944847,     0.01329909376597139,     -0.020904714318069645,    -0.1196918438937425,     0.8177420190660831}
        };

        const Eigen::Matrix<float, h1Height, h1Width> h1{
                {6.391587676622346e-010,  0.0, -1.7257286726880333e-08,  -3.067962084778726e-08,  1.3805829381504267e-07,  5.522331752601707e-07,
                                                                                                                                3.3747582932565985e-07, -1.9328161134105974e-06, -5.6949046198705095e-06, -7.649452131381623e-06},
                {0.0,                     0.0, 0.0,                      0.0,                     0.0,                     0.0, 0.0,                    0.0,                     0.0,                     0.0},
                {-1.7257286726880333e-08, 0.0, 4.65946741625769e-07,     8.283497628902559e-07,   -3.727573933006152e-06,  -0.000014910295732024608,
                                                                                                                                -9.111847391792816e-06, 0.000052186035062086126, 0.00015376242473650378,  0.00020653520754730382},
                {-3.067962084778726e-08,  0.0, 8.283497628902559e-07,    -2.9917573832012203e-07, -6.6267981031220475e-06, 5.3851632897621965e-06,
                                                                                                                                0.00004049868144081346, -0.00001884807151416769, -0.00023692226948222173, -0.0003769812640795245},
                {1.3805829381504267e-07,  0.0, -3.727573933006152e-06,   -6.6267981031220475e-06, 0.000029820591464049215, 0.00011928236585619686,
                                                                                                                                0.00007289477913434253, -0.000417488280496689,   -0.0012300993978920302,  -0.0016522816603784306},
                {5.522331752601707e-07,   0.0, -0.000014910295732024608, 5.3851632897621965e-06,  0.00011928236585619686,  -0.00009693293921571956,
                                                                                                                                -0.0007289762659346422, 0.00033926528725501844,  0.004264600850679991,    0.006785662753431441},
                {3.3747582932565985e-07,  0.0, -9.111847391792816e-06,   0.00004049868144081346,  0.00007289477913434253,  -0.0007289762659346422,
                                                                                                                                -0.001468581806076247,  0.002551416930771248,    0.01181401175635013,     0.017093222023136675},
                {-1.9328161134105974e-06, 0.0, 0.000052186035062086126,  -0.00001884807151416769, -0.000417488280496689,   0.00033926528725501844,
                                                                                                                                0.002551416930771248,   -0.0011874285053925643,  -0.01492610297737997,    -0.023749819637010044},
                {-5.6949046198705095e-06, 0.0, 0.00015376242473650378,   -0.00023692226948222173, -0.0012300993978920302,  0.004264600850679991,
                                                                                                                                0.01181401175635013,    -0.01492610297737997,    -0.0826466923977296,     -0.12203257624594532},
                {-7.649452131381623e-06,  0.0, 0.00020653520754730382,   -0.0003769812640795245,  -0.0016522816603784306,  0.006785662753431441,
                                                                                                                                0.017093222023136675,   -0.023749819637010044,   -0.12203257624594532,    0.821896776039774}
        };


        /*
            g0 = [g0   fliplr(g0(:, 1 : end - 1))];
            g0 = [g0; flipud(g0(1:end - 1, : ))];
            h0 = [h0   fliplr(h0(:, 1 : end - 1))];
            h0 = [h0; flipud(h0(1:end - 1, : ))];

            g1 = [g1   fliplr(g1(:, 1 : end - 1))];
            g1 = [g1; flipud(g1(1:end - 1, : ))];
            h1 = [h1   fliplr(h1(:, 1 : end - 1))];
            h1 = [h1; flipud(h1(1:end - 1, : ))];
        */

        atrous_filters = new Cont(4);
        atrous_filters->Create_Cell(0, 1)->Set(2 * g0Height - 1, 2 * g0Width - 1).Create::Default();
        atrous_filters->Create_Cell(1, 1)->Set(2 * h0Height - 1, 2 * h0Width - 1).Create::Default();
        atrous_filters->Create_Cell(2, 1)->Set(2 * g1Height - 1, 2 * g1Width - 1).Create::Default();
        atrous_filters->Create_Cell(3, 1)->Set(2 * h1Height - 1, 2 * h1Width - 1).Create::Default();

        // lr -> left to right flip, ud -> up to down flip, h_ -> horizontal, v_ -> vertical

        // g0
        Eigen::Matrix<float, g0Height, 2 * g0Width - 1> h_extended_g0;
        atrous_filters->_conts[0]->Set(2 * g0Height - 1, 2 * g0Width - 1);
        h_extended_g0 << g0, g0.block<g0Height, g0Width - 1>(0, 0).rowwise().reverse();
        atrous_filters->_conts[0]->_mat << h_extended_g0,
                h_extended_g0.block<g0Height - 1, h_extended_g0.cols()>(0, 0).colwise().reverse();

        // h0
        Eigen::Matrix<float, h0Height, 2 * h0Width - 1> h_extended_h0;
        atrous_filters->_conts[1]->Set(2 * h0Height - 1, 2 * h0Width - 1);
        h_extended_h0 << h0, h0.block<h0Height, h0Width - 1>(0, 0).rowwise().reverse();
        atrous_filters->_conts[1]->_mat << h_extended_h0,
                h_extended_h0.block<h0Height - 1, h_extended_h0.cols()>(0, 0).colwise().reverse();

        // g1
        Eigen::Matrix<float, g1Height, 2 * g1Width - 1> h_extended_g1;
        atrous_filters->_conts[2]->Set(2 * g1Height - 1, 2 * g1Width - 1);
        h_extended_g1 << g1, g1.block<g1Height, g1Width - 1>(0, 0).rowwise().reverse();
        atrous_filters->_conts[2]->_mat << h_extended_g1,
                h_extended_g1.block<g1Height - 1, h_extended_g1.cols()>(0, 0).colwise().reverse();

        // h1
        Eigen::Matrix<float, h1Height, 2 * h1Width - 1> h_extended_h1;
        atrous_filters->_conts[3]->Set(2 * h1Height - 1, 2 * h1Width - 1);
        h_extended_h1 << h1, h1.block<h1Height, h1Width - 1>(0, 0).rowwise().reverse();
        atrous_filters->_conts[3]->_mat << h_extended_h1,
                h_extended_h1.block<h1Height - 1, h_extended_h1.cols()>(0, 0).colwise().reverse();

    } else exit(1);
}


/// <summary>
///		This function generates the matrix that contains the number of times the polar grid points go through the rectangular grid
/// </summary>
/// <param name="L : ">
///		L is the order of the block matrix
/// </param>
/// <returns>
///		D is the common grid point values
/// </returns>
Tensor *NSST::AvgPol(int L, const ArrayXXf &x1, const ArrayXXf &y1, const ArrayXXf &x2, const ArrayXXf &y2) {
    Tensor *d = new Tensor;
    d->Set(L, L).Create::Zero();

    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j) {
            d->_mat(static_cast<int>((y1(i, j) - 1)), static_cast<int>((x1(i, j) - 1)))++;
            d->_mat(static_cast<int>((y2(i, j) - 1)), static_cast<int>((x2(i, j) - 1)))++;
        }

    return d;
}

/// <summary>
///		This function generates the x and y vectors that contain the i, j coordinates to extract radial slices
/// </summary>
/// <param name="n : ">
///		n is the order of the block to be used
/// </param>
/// <returns>
///		** return [x1n, y1n, x2n, y2n, D] **
///		x1, y1 are the i, j values that correspond to the radial slices from the endpoints 1, 1 to n, 1 through the origin <para></para>
///		x2, y2 are the i, j values that correspond to the radial slices from the endpoints 1, 1 to 1, n through the origin <para></para>
///		D is the matrix that contains the number of times the polar grid points go through the rectangular grid
/// </returns>
Tensor *NSST::GenXYCoordinates(int n) {
    ++n;
    ArrayXXf x1 = ArrayXXf::Zero(n, n);
    ArrayXXf y1 = ArrayXXf::Zero(n, n);
    ArrayXXf x2 = ArrayXXf::Zero(n, n);
    ArrayXXf y2 = ArrayXXf::Zero(n, n);

    float *xt = nullptr;
    std::unique_ptr<float> m(ZERO(n));

    constexpr int y0 = 1;
    int flag = 0, x0 = 0, xN = 0, yN = n;

    for (int i = 0; i < n; ++i) {
        x0 = i + 1;
        xN = n - x0 + 1;

        if (xN == x0) flag = 1;
        else {
            m.get()[i] = static_cast<float>(yN - y0) / static_cast<float>(xN - x0);
            flag = 0;
        }

        xt = Linspace(x0, xN, n);

        for (int j = 0; j < n; ++j) {
            if (flag == 0) {
                x1(i, j) = y2(i, j) = roundf(xt[j]);
                x2(i, j) = y1(i, j) = roundf(m.get()[i] * static_cast<float>(xt[j] - x0) + y0);
            } else {
                x1(i, j) = y2(i, j) = static_cast<int>((n - 1) / 2) + 1;
                x2(i, j) = y1(i, j) = j + 1;
            }
        }

        delete xt;
    }

    --n;
    // correct for portion outside boundry
    x1 = x1.block(0, 0, n, n).colwise().reverse().eval();
    y1 = y1.block(0, 0, n, n).eval();
    x2 = x2.block(1, 0, n, n).eval();
    y2 = y2.block(1, 0, n, n).eval();
    y2(n - 1, 0) = n;

    // return [x1n,y1n,x2n,y2n,D]
    Tensor *ret = new Tensor[5];
    ret[0].Set(n, n);
    ret[0] = x1;
    ret[1].Set(n, n);
    ret[1] = y1;
    ret[2].Set(n, n);
    ret[2] = x2;
    ret[3].Set(n, n);
    ret[3] = y2;

    ret[4].Set(n, n);
    ret[4] = *AvgPol(n, x1, y1, x2, y2);

    return ret;
}


/// <summary>
///		This function computes the Meyer window function of signal x
/// </summary>
/// <param name="x : ">
///		x - signal to be windowed
/// </param>
/// <returns>
///		y - windowed signal
/// </returns>
inline float MeyerWind(const float &x) {

    float y = 0.0F;
    constexpr float x_less = -1.0F / 3.0F + 1.0F / 2.0F, x_greater = 1.0F / 3.0F + 1.0F / 2.0F;
    constexpr float x_less_equal1 = 1.0F / 3.0F + 1.0F / 2.0F, x_greater_equal1 = 2.0F / 3.0F + 1.0F / 2.0F;
    constexpr float x_less_equal2 = -2.0F / 3.0F + 1.0F / 2.0F, x_greater_equal2 = 1.0F / 3.0F + 1.0F / 2.0F;

    if (x_less < x && x < x_greater)
        y = 1.0F;

    else if ((x_less_equal1 <= x && x <= x_greater_equal1) || (x_less_equal2 <= x && x <= x_greater_equal2)) {

        float w = 3.0F * abs(x - 1.0F / 2.0F) - 1.0F;
        float z = powf(w, 4.0F) * (35.0F - 84.0F * w + 70.0F * powf(w, 2.0F) - 20.0F * powf(w, 3.0F));
        y = powf(cosf(PI / 2.0F * (z)), 2.0F);
    }

    return y;
}

/// <summary>
///		This function computes L band sequence of a Meyer windowed signal.
/// </summary>
/// <param name="x : ">
///		x - signal to be decomposed
/// </param>
/// <param name="N : ">
///     N - len(x)
/// </param>
/// <param name="L : ">
///		L - number of bandpass filters
/// </param>
/// <returns>
///		y - windowed signal of x where second index indicates which bandpass index
/// </returns>
Tensor *NSST::Windowing(const float *x, int N, int L) {

    Tensor *y = new Tensor;
    y->Set(N, L).Create::Zero();

    float T = N / L;
    std::unique_ptr<float> g(ZERO(2 * T));

    float n = 0;
    for (int j = 0; j < 2 * T; ++j) {
        n = -1.0F * T / 2.0F + j;
        g.get()[j] = MeyerWind(n / T);
    }

    int index = 0, in_sig = 0;
    for (int j = 0; j < L; ++j) {
        index = 0;
        for (int k = -1 * T / 2; k <= 1.5 * T - 1; ++k) {
            in_sig = floor(Realmod(static_cast<int>(k + j * T), N));
            y->_mat(in_sig, j) = g.get()[index] * x[in_sig];
            index++;
        }
    }

    return y;
}

/// <summary>
///		This funcion re - assembles the radial slice into block.
/// </summary>
/// <param name="l : ">
///		l is the radial slice matrix
/// </param>
/// <param name="n : ">
///		n is the order of the matrix that is to be re - assembled
///		x1, y1, x2, y2 are the polar coordinates generated from function
/// </param>
/// <param name="gen : ">
///		D is the matrix containing common polar grid points
/// </param>
/// <returns>
///		C is the re - assembled block matrix
/// </returns>
Tensor &NSST::RecFromPol(const ArrayXXf &l, int n, const Tensor *gen) {
    Tensor matC;
    matC.Set(n, n).Create::Zero();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matC._mat(static_cast<int>(gen[1]._mat(i, j) - 1), static_cast<int>(gen[0]._mat(i, j) - 1)) += l(i, j);
            matC._mat(static_cast<int>(gen[3]._mat(i, j) - 1), static_cast<int>(gen[2]._mat(i, j) - 1)) += l((i + n),j);
        }
    }

    return matC / gen[4];
}

/// <summary>
///		This function computes the directional / shearing filters using the Meyer window function.
/// </summary>
/// <param name="n : ">
///		n indicates the supports size of the directional filter is n1xn1 level indicates that the number of directions is 2 ^ level
/// </param>
/// <param name="level : ">
///		level indicates
/// </param>
/// <returns>
///		a sequence of 2D directional/shearing filters w_s where the third index determines which directional filter is used
/// </returns>
void NSST::ShearingFiltersMyer() {
    // Pseudo-polar koordinat sistemi icin indexler uretilir
    // gen : [x11, y11, x12, y12, F1]

    int size = 0;
    ArrayXXf temp;
    shearing_filters = new Tensor *[Conf::sp.dcompSize];
    FFTWF<ForwEigen, BackEigen> *ifft2d = nullptr;

    for (int t = 0; t < Conf::sp.dcompSize; ++t) {

        Tensor *gen = GenXYCoordinates(Conf::sp.dsize[t]);
        std::unique_ptr<float> ones_(ONE(2 * Conf::sp.dsize[t]));
        std::unique_ptr<Tensor> wf(Windowing(ones_.get(), 2 * Conf::sp.dsize[t], pow(2, Conf::sp.dcomp[t])));

        size = static_cast<int>(pow(2, Conf::sp.dcomp[t]));
        shearing_filters[t] = new Tensor[size];
        Eigen::MatrixXf one = Eigen::MatrixXf::Constant(1, Conf::sp.dsize[t], 1);
        ifft2d = FFTWF_Wrapper::Eigen::_2D::Get(Conf::sp.dsize[t], Conf::sp.dsize[t]);

        for (int i = 0; i < size; ++i) {
            // w_s(:, : , k) = real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);
            temp = wf.get()->_mat.block(0, i, wf.get()->_mat.rows(), 1).matrix() * one;

            shearing_filters[t][i].Set(Conf::sp.dsize[t], Conf::sp.dsize[t]);
            shearing_filters[t][i] = RecFromPol(temp, Conf::sp.dsize[t], gen);
            shearing_filters[t][i]._mat = fftshift(shearing_filters[t][i]._mat).eval();  // fftshift(w_s(:, : , k))

            ifft2d->Backward()->SetIn(shearing_filters[t][i]._mat);
            ifft2d->Backward()->Execute();                                  // ifft2(fftshift(w_s(:, : , k)))
            ifft2d->Backward()->GetOut(shearing_filters[t][i]._mat);

            // window dizisini kartezyen koordinat sistemine donusturur.
            shearing_filters[t][i]._mat = fftshift(shearing_filters[t][i]._mat);        // real(fftshift(ifft2(fftshift(w_s(:, : , k)))))

            shearing_filters[t][i] /= sqrtf(Conf::sp.dsize[t]);             // real(fftshift(ifft2(fftshift(w_s(:, : , k))))). / sqrt(n1);
        }

        delete[] gen;
    }

    delete ifft2d;
}


/// <summary>
///		This function computes the(local) nonsubsampled shearlet transform as given
///
///     shear_parameters.dcomp - a vector such that.dcomp(i) indicates that the
/// 	ith decomposition level has 2 ^ decomp(i)
/// 	directions.The length of the vector plus 1 is
/// 	total the number of decompostions.
/// </summary>
/// <param name="image">
///		input image
/// <returns></returns>
Cont *NSST::Dec(Tensor *image) {
    // Laplacian Pyramid decomposition	//NSLP //cok olceklilik
    // Katsayilara gore alt-goruntuler elde edilir.

    int level = Conf::sp.dcompSize;
    Cont *y = AtrousDec(image, level);

    Cont *dst = new Cont(level + 1);
    dst->_conts[0] = new Tensor(y->_conts[0]);
    ArrayXXf shear_f;
    Tensor* temp = nullptr;

    int size = 0;
    for (int i = 0; i < level; ++i) {
        size = (int) pow(2, Conf::sp.dcomp[i]);    // goruntuye uygulanacak yon sayisi belirlenir.
        dst->Create_Cell(i + 1, size);

        for (int k = 0; k < size; ++k) {
            shear_f = shearing_filters[i][k]._mat * sqrt(Conf::sp.dsize[i]);
            temp = Conv2D(y->_conts[i + 1][0]._mat, shear_f, conv_type::same);     // Cok yonluluk
            temp->CopyTo(&(dst->_conts[i + 1][k]));
            delete temp;
        }
        dst->_conts[i + 1]->_d = size;
    }

    delete y;
    return dst;
}


/// <summary>
///		ATROUSDEC - computes the 2 - D atrous decomposition using symmetric extension.
/// </summary>
/// <param name="image : ">
///		input image
/// </param>
/// <param name="level : ">
///		N levels - number of decomposition levels
/// </param>
/// <returns>
///		y->[0] => LFC (Low frequency coefficients)     <=> (AFK)
///		<para></para>
///		y->[1..4] => HFC (High frequency coefficients) <=> (YFK)
/// </returns>
Cont *NSST::AtrousDec(Tensor *image, int level) {
    Tensor *y0, *y1, *sym, *upSample;
    float shift[2] = {1.0F, 1.0F};

    // NSLP - Non Subsampled Laplacian Pyramid
    // Orijinal goruntuye alcak ve yuksek katsayilarin uygulanmasi
    // Symext kullanilarak katsayilara gore genisletilir, ardindan conv2 ile kucultulurek ayni boyuta getirilir.
    sym = Symext(image, atrous_filters->_conts[1], shift);
    y0 = Conv2D(sym->_mat, atrous_filters->_conts[1]->_mat, conv_type::valid);    // Alcak
    delete sym;

    sym = Symext(image, atrous_filters->_conts[3], shift);
    y1 = Conv2D(sym->_mat, atrous_filters->_conts[3]->_mat, conv_type::valid);    // Yuksek
    delete sym;

    Cont *y = new Cont(level + 1);
    y->_conts[level] = y1;
    image = y0;

    Eigen::Array<float, 2, 2> I2, I2L;
    I2 << 1, 1, 1, 1;
    int L;

    // Seviye 1'e inene kadar her alcak frekans katsayilarinin alcak ve yuksek katsayisi hesaplanir.
    // Yuksek frekansa katsayilari y[4..1] olarak tutulur.
    for (int i = 1; i <= level - 1; ++i) {
        shift[0] = shift[1] = 2.0F - powf(2.0F, float(i - 1));
        L = static_cast<int>(pow(2, i));
        I2L = I2 * L;

        upSample = Upsample2df(atrous_filters->_conts[1], i);
        sym = Symext(image, upSample, shift);
        y0 = Atrousc(sym, atrous_filters->_conts[1], I2L);
        delete upSample;
        delete sym;

        upSample = Upsample2df(atrous_filters->_conts[3], i);
        sym = Symext(image, upSample, shift);
        y1 = Atrousc(sym, atrous_filters->_conts[3], I2L);
        delete upSample;
        delete sym;
        delete image;

        y->_conts[level - i] = y1;
        image = y0;
    }
    y->_conts[0] = image;        // Alcak frekans katsayisi atanir.

    return y;
}

/// <summary>
///		Performs symmetric extension for image x, filter h. <para></para>
///		The filter h is assumed have odd dimensions. <para></para>
///		If the filter has horizontaland vertical symmetry, then the nonsymmetric part of conv2(h, x) has the same size of x.
/// </summary>
/// <param name="x : ">
///		mxn image
/// </param>
/// <param name="h : ">
///		2 - D filter coefficients
/// </param>
/// <param name="shift : ">
///		optional shift
/// </param>
/// <returns>
///		yT image symetrically extended(H / V symmetry)
/// </returns>
Tensor *NSST::Symext(const Tensor *x, const Tensor *h, const float *shift) {
    int m = x->_h;
    int n = x->_w;
    int p = h->_h;
    int q = h->_w;

    int s1 = static_cast<int>(shift[0]);
    int s2 = static_cast<int>(shift[1]);

    int ss = int(floorf(p / 2.0F)) - s1 + 1;
    int rr = int(floorf(q / 2.0F)) - s2 + 1;

    ArrayXXf extended_fliplr(x->_h, x->_w + ss);
    ArrayXXf yT(x->_h, extended_fliplr.cols() + p + s1);

    // [fliplr(x(:,1:ss)) x  x(:,n  :-1: n-p-s1+1)]
    extended_fliplr << x->_mat.block(0, 0, x->_h, ss).rowwise().reverse(), x->_mat;         // fliplr(x(:,1:ss))
    yT << extended_fliplr, x->_mat.block(0, n - (p + s1), x->_h,
                                         p + s1).rowwise().reverse();  // x(:, n : -1 : n - p - s1 + 1)

    ArrayXXf extended_flipud(yT.rows() + rr, yT.cols());
    ArrayXXf yT2(extended_flipud.rows() + q + s2, extended_flipud.cols());

    // [flipud(yT(1:rr, : )); yT;  yT(m  :-1 : m - q - s2 + 1, : )]
    extended_flipud << yT.block(0, 0, rr, yT.cols()).colwise().reverse(), yT;       // flipud(yT(1:rr, : ))
    yT2 << extended_flipud, yT.block(m - (q + s2), 0, q + s2, yT.cols());           // yT(m  :-1 : m - q - s2 + 1, : )

    // yT(1:m+p-1 ,1:n+q-1)
    Tensor *ret_yT = new Tensor();
    ret_yT->Set(m + p - 1, n + q - 1);
    ret_yT->_mat = yT2.block(0, 0, m + p - 1, n + q - 1);

    return ret_yT;
}


/// <summary>
///		This function does not actually upsample the filter,
///		it computes the convolution as if the filter had been upsampled.
///		This is the ultimate optimized version. Further optimized for separable(diagonal) upsampling matrices.
/// </summary>
/// <param name="signal : ">
///		A 2D Signal
/// </param>
/// <param name="filter : ">
///		2D filter
/// </param>
/// <param name="upMatrix : ">
///		separable upsampling matrix
/// </param>
/// <returns>
///		2D result of convolution with filter upsampled by a m, only the 'valid' part is returned. <para></para>
///		Similar to conv2(x, h, 'valid'), where h is the upsampled filter.
/// </returns>
Tensor *NSST::Atrousc(const Tensor *signal, const Tensor *filter, const Eigen::Array22f &upMatrix) {
    /*
        FArray   - Filter coefficients
        SArray   - Signal coefficients
        outArray - Output coefficients
        M        - upsampling matrix
    */
    int SColLength, SRowLength, FColLength, FRowLength, O_SColLength, O_SRowLength;
    int SFColLength, SFRowLength;
    int n1, n2, k1, k2, f1, f2, kk2, kk1;
    float sum, *M;
    int M0, M3, sM0, sM3;

    SColLength = signal->_w;
    SRowLength = signal->_h;
    FColLength = filter->_w;
    FRowLength = filter->_h;

    SFColLength = FColLength - 1;
    SFRowLength = FRowLength - 1;

    M0 = static_cast<int>(upMatrix(0));
    M3 = static_cast<int>(upMatrix(3));
    sM0 = M0 - 1;
    sM3 = M3 - 1;

    O_SColLength = SColLength - M0 * FColLength + 1;
    O_SRowLength = SRowLength - M3 * FRowLength + 1;

    Tensor *outArray = new Tensor();
    outArray->Set(O_SRowLength, O_SColLength).Create::Default();
    float *out_p = outArray->_mat.data();
    const float *signal_p = signal->_mat.data(), *filter_p = filter->_mat.data();

    /* Convolution loop */
    for (n1 = 0; n1 < O_SRowLength; n1++) {
        for (n2 = 0; n2 < O_SColLength; n2++) {
            sum = 0;
            kk1 = n1 + sM0;
            for (k1 = 0; k1 < FRowLength; k1++) {
                kk2 = n2 + sM3;
                for (k2 = 0; k2 < FColLength; k2++) {
                    f1 = SFRowLength - k1;    /* flipped index */
                    f2 = SFColLength - k2;
                    sum += filter_p[f1 * FColLength + f2] * signal_p[kk1 * SColLength + kk2];
                    kk2 += M3;
                }
                kk1 += M0;
            }
            out_p[n1 * O_SColLength + n2] = sum;
        }
    }

    return outArray;
}


/// <summary>
///		This function performs the inverse(local) nonsubsampled shearlet transform as given
/// </summary>
/// <param name="dst">
///		dst - the nonsubsampled shearlet coefficients
/// </param>
/// <param name="lpfilt">
///		lpfilt - the filter to be used for the Laplacian Pyramid / ATrous decomposition using the codes
/// </param>
/// <returns>
///		the reconstructed image
/// </returns>
Tensor *NSST::Rec(const Cont *__restrict dst) {

    // , const Cont* filters -> atrousfilter
    int level = dst->_cont_num - 1;
    Cont *y = new Cont(level + 1);
    y->_conts[0] = new Tensor(dst->_conts[0]);

    for (int i = 1; i <= level; ++i)
        y->_conts[i] = Sum_Axis3(dst->_conts[i]);

    Tensor *ret = AtrousRec(y);

    delete y;
    return ret;
}


/// <summary>
///		SATROUSREC - computes the inverse of 2 - D atrous decomposition computed with ATROUSDEC
///		y = satrousdrc(x, fname)
///     lpfilt: can be any filter available in the function atrousfilters
/// </summary>
/// <param name="y : ">
///		image
/// </param>
/// <param name="lpfilt : ">
///
///	</param>
/// <returns>
///		x : reconstructed image
/// </returns>
Tensor *NSST::AtrousRec(const Cont *y) {

    int NLevels = (y->_cont_num - 1), L;
    Tensor *x = new Tensor, *temp_p, temp;
    y->_conts[0]->CopyTo(x);

    Eigen::Array<float, 2, 2, Eigen::RowMajor> I2 {{1.0F, 1.0F},{1.0F, 1.0F}}, I2L;

    float shift[2] = {1.0F, 1.0F};
    Tensor *up, *sym, *atrousc1, *atrousc2;
    temp_p = x;

    for (int i = NLevels - 1; i >= 1; i--) {

        shift[0] = shift[1] = 2.0F - powf(2.0F, static_cast<float>(i - 1));

        L = static_cast<int>(pow(2, i));
        I2L = I2 * L;

        // x = atrousc(Symext(x, upsample2df(g0, i), shift), g0, L * I2) + atrousc(Symext(y1, upsample2df(g1, i), shift), g1, L * I2);

        up = Upsample2df(atrous_filters->_conts[0], i);
        sym = Symext(temp_p, up, shift);
        atrousc1 = Atrousc(sym, atrous_filters->_conts[0], I2L);
        delete up;
        delete sym;

        up = Upsample2df(atrous_filters->_conts[2], i);
        sym = Symext(y->_conts[NLevels - i], up, shift);
        atrousc2 = Atrousc(sym, atrous_filters->_conts[2], I2L);
        delete up;
        delete sym;

        temp = atrousc1->_mat + atrousc2->_mat;
        temp_p = &temp;

        delete atrousc1;
        delete atrousc2;
    }

    shift[0] = shift[1] = 1.0F;

    std::unique_ptr<Tensor> symetxX(Symext(temp_p, atrous_filters->_conts[0], shift));
    std::unique_ptr<Tensor> symetxY(Symext(y->_conts[NLevels], atrous_filters->_conts[2], shift));

    std::unique_ptr<Tensor> conv2X(Conv2D(symetxX->_mat, atrous_filters->_conts[0]->_mat, conv_type::valid));
    std::unique_ptr<Tensor> conv2Y(Conv2D(symetxY->_mat, atrous_filters->_conts[2]->_mat, conv_type::valid));

    x->Set(conv2X->_h, conv2X->_w);
    x->_mat = conv2X->_mat + conv2Y->_mat;

    return x;
}