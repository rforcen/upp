#pragma once

#include "Thread.h"
#include <complex>

using std::complex;

typedef uint32_t u32;

template <typename T> class Mandelbrot
{
  private:
    vector<u32> firePallete = {
        0,        0,        4,        12,       16,       24,       32,       36,       44,       48,       56,
        64,       68,       76,       80,       88,       96,       100,      108,      116,      120,      128,
        132,      140,      148,      152,      160,      164,      172,      180,      184,      192,      200,
        1224,     3272,     4300,     6348,     7376,     9424,     10448,    12500,    14548,    15576,    17624,
        18648,    20700,    21724,    23776,    25824,    26848,    28900,    29924,    31976,    33000,    35048,
        36076,    38124,    40176,    41200,    43248,    44276,    46324,    47352,    49400,    51452,    313596,
        837884,   1363196,  1887484,  2412796,  2937084,  3461372,  3986684,  4510972,  5036284,  5560572,  6084860,
        6610172,  7134460,  7659772,  8184060,  8708348,  9233660,  9757948,  10283260, 10807548, 11331836, 11857148,
        12381436, 12906748, 13431036, 13955324, 14480636, 15004924, 15530236, 16054524, 16579836, 16317692, 16055548,
        15793404, 15269116, 15006972, 14744828, 14220540, 13958396, 13696252, 13171964, 12909820, 12647676, 12123388,
        11861244, 11599100, 11074812, 10812668, 10550524, 10288380, 9764092,  9501948,  9239804,  8715516,  8453372,
        8191228,  7666940,  7404796,  7142652,  6618364,  6356220,  6094076,  5569788,  5307644,  5045500,  4783356,
        4259068,  3996924,  3734780,  3210492,  2948348,  2686204,  2161916,  1899772,  1637628,  1113340,  851196,
        589052,   64764,    63740,    62716,    61692,    59644,    58620,    57596,    55548,    54524,    53500,
        51452,    50428,    49404,    47356,    46332,    45308,    43260,    42236,    41212,    40188,    38140,
        37116,    36092,    34044,    33020,    31996,    29948,    28924,    27900,    25852,    24828,    23804,
        21756,    20732,    19708,    18684,    16636,    15612,    14588,    12540,    11516,    10492,    8444,
        7420,     6396,     4348,     3324,     2300,     252,      248,      244,      240,      236,      232,
        228,      224,      220,      216,      212,      208,      204,      200,      196,      192,      188,
        184,      180,      176,      172,      168,      164,      160,      156,      152,      148,      144,
        140,      136,      132,      128,      124,      120,      116,      112,      108,      104,      100,
        96,       92,       88,       84,       80,       76,       72,       68,       64,       60,       56,
        52,       48,       44,       40,       36,       32,       28,       24,       20,       16,       12,
        8,        0,        0};

  public:
    typedef complex<T> Cmplx;

    int w = 0, h = 0, iters = 200;
    Cmplx center = Cmplx(0.5, 0.0), range = Cmplx(-2.0, 2.0), cr;
    T rir, scale;
    u32 *image = nullptr;

  private:
    inline Cmplx doScale(T iw, T jh)
    {
        return cr + rir * Cmplx(iw, jh);
    }

    void swapBR() // firePallete blue <-> red
    {
        auto _swapBR = [=](u32 &value) {
            uint32_t byte0 = (value & 0x000000FF);
            uint32_t byte2 = (value & 0x00FF0000) >> 16;
            value &= 0xFF00FF00; // Clear the 0th and 2nd bytes
            value |= (byte0 << 16) | byte2;
        };
        for (auto &f : firePallete)
            _swapBR(f);
    }

  public:
    // works from external image array
    Mandelbrot(u32 *image, u32 w, u32 h, u32 iters, Cmplx center, Cmplx range)
        : w(w), h(h), iters(iters), center(center), range(range), image(image), cr(Cmplx(range.real(), range.real())),
          rir((range.imag() - range.real())), scale(0.8 * T(w) / h)
    {
        swapBR();
    }

    int size_bytes()
    {
        return w * h * sizeof(*image);
    }
    int image_size()
    {
        return w * h;
    }

    void genPixel(int index)
    {
        const Cmplx c0 = scale * doScale(T(index % w) / w, T(index / w) / h) - center;
        Cmplx z = c0;

        int ix = iters;
        for (int it = 0; it < iters; it++)
        {
            z = z * z + c0;
            if (norm(z) > 4.0)
            {
                ix = it;
                break;
            }
        }
        image[index] = 0xff000000 | ((ix == iters) ? 0
                                                   //: firePallete[(256 * ix) / 50]);
                                                   : firePallete[(ix << 2) % firePallete.size()]);
    }

    // single thread
    void maneldebrotST()
    {
        for (auto index = 0; index < image_size(); index++)
            genPixel(index);
    }

    // multithread version
    void maneldebrotMT()
    {
        Thread(image_size()).run([this](int index) { genPixel(index); });
    }

	// helper func
    static vector<u32> genMandelbrot(u32 w, u32 h, u32 iters, Cmplx center, Cmplx range)
    {
        vector<u32>vi(w*h);
        Mandelbrot<T>(vi.data(), w, h, iters, center, range).maneldebrotMT();
        return vi;
    }
};
