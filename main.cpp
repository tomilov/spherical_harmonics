#include "Math/Vector.h"

//#include <boost/math/special_functions.hpp>

#include <random>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <limits>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ostream>

#include <cmath>
#include <cassert>

struct spherical_harmonics
{

    using size_type = std::size_t;

    float const eps = std::numeric_limits< float >::epsilon();
    float const one = 1.0f;
    float const zero = 0.0f;

    float SH_Cartesian(int l, int m, float3 n)
    {
        //n /= length(n);
        float value = 1000.0;
        switch (l) {
        case 0 : {
            value = 1.0;
            break;
        }
        case 1 : {
            switch (m) {
            case -1 : {
                value = -n.y;
                break;
            }
            case 0 : {
                value = n.z;
                break;
            }
            case 1 : {
                value = -n.x;
                break;
            }
            }
            break;
        }
        case 2 : {
            switch (m) {
            case -2 : {
                value = sqrt(3.0) * n.x * n.y;
                break;
            }
            case -1 : {
                value = -sqrt(3.0) * n.y * n.z;
                break;
            }
            case 0 : {
                value = n.z * n.z - 0.5 * (n.x * n.x + n.y * n.y);
                break;
            }
            case 1 : {
                value = -sqrt(3.0) * n.z * n.x;
                break;
            }
            case 2 : {
                value = sqrt(3.0 / 4.0) * (n.x * n.x - n.y * n.y);
                break;
            }
            }
            break;
        }
        case 3 : {
            switch (m) {
            case -3 : {
                value = -sqrt(5.0 / 8.0) * (3.0 * n.x * n.x - n.y * n.y) * n.y;
                break;
            }
            case -2 : {
                value = sqrt(15.0) * n.x * n.y * n.z;
                break;
            }
            case -1 : {
                value = -sqrt(3.0 / 8.0) * n.y * (4.0 * n.z * n.z - (n.x * n.x + n.y * n.y));
                break;
            }
            case 0 : {
                value = n.z * (n.z * n.z - 1.5 * (n.x * n.x + n.y * n.y));
                break;
            }
            case 1 : {
                value = -sqrt(3.0 / 8.0) * n.x * (4.0 * n.z * n.z - (n.x * n.x + n.y * n.y));
                break;
            }
            case 2 : {
                value = sqrt(15.0 / 4.0) * (n.x * n.x - n.y * n.y) * n.z;
                break;
            }
            case 3 : {
                value = -sqrt(5.0 / 8.0) * (n.x * n.x - 3.0 * n.y * n.y) * n.x;
                break;
            }
            }
            break;
        }
        case 4 : {
            switch (m) {
            case -4 : {
                value = sqrt(35.0 / 4.0) * n.x * n.y * (n.x * n.x - n.y * n.y);
                break;
            }
            case -3 : {
                value = -sqrt(35.0 / 8.0) * (3.0 * n.x * n.x - n.y * n.y) * n.y * n.z;
                break;
            }
            case -2 : {
                value = sqrt(5.0 / 4.0) * n.x * n.y * (7.0 * n.z * n.z - 1.0);
                break;
            }
            case -1 : {
                value = -sqrt(5.0 / 8.0) * n.y * n.z * (7.0 * n.z * n.z - 3.0);
                break;
            }
            case 0 : {
                n.z *= n.z;
                value = 0.125 * ((35.0 * n.z - 30.0) * n.z + 3.0);
                break;
            }
            case 1 : {
                value = -sqrt(5.0 / 8.0) * n.x * n.z * (7.0 * n.z * n.z - 3.0);
                break;
            }
            case 2 : {
                value = sqrt(5.0 / 16.0) * (n.x * n.x - n.y * n.y) * (7.0 * n.z * n.z - 1.0);
                break;
            }
            case 3 : {
                value = -sqrt(35.0 / 8.0) * (n.x * n.x - 3.0 * n.y * n.y) * n.x * n.z;
                break;
            }
            case 4 : {
                n.x *= n.x;
                n.y *= n.y;
                value = sqrt(35.0 / 64.0) * (n.x * n.x + n.y * n.y - 6.0 * n.x * n.y);
                break;
            }
            }
            break;
        }
        }
        return value;
    }

    float factorial(int x)
    {
        float f = 1.0;
        for (int i = 2; i <= x; i++)
        {
            f *= i;
        }

        return f;
    }

    // Evaluate an Associated Legendre Polynomial P(l, m, x) at x
    float P(int l, int m, float x)
    {
        float pmm = 1.0;
        if (m > 0)
        {
            float somx2 = sqrt((1.0 - x) * (1.0 + x));

            float fact = 1.0;
            for (int i = 1; i <= m; i++)
            {
                pmm *= (-fact) * somx2;
                fact += 2.0;
            }
        }
        if (l == m) return pmm;

        float pmmp1 = x * (2.0 * m + 1.0) * pmm;
        if (l == m + 1) return pmmp1;

        float pll = 0.0;
        for (int ll = m + 2; ll <= l; ++ll)
        {
            pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pll;
    }

    // Normalization constant
    float K(int l, int m)
    {
        if (m == 0) {
            return sqrt((2.0 * l + 1.0) / (4.0 * PI));
        } else {
            return sqrt(((2.0 * l + 1.0) * factorial(l - m)) / (4.0 * PI * factorial(l + m)));
        }
    }

    float SH(const int l, const int m, float3 pos)
    {
        float x;
        float len = length(pos);

        float phi = atan2(pos.y, pos.x);
        float t = pos.z / len;

        float sqrt2 = 1.4142135623731;

        if (m == 0)
            x = K(l, 0) * P(l, m, t);
        else if (m > 0)
            x = sqrt2 * K(l, m) * cos(m * phi) * P(l, m, t);
        else
            x = sqrt2 * K(l, -m) * sin(-m * phi) * P(l, -m, t);

        float y = sqrt((2.0 * l + 1.0) / (4.0 * PI)) * SH_Cartesian(l, m, pos);
        if (std::abs(x - y) > 1E-5) {
            std::cerr << x << ' ' << y << std::endl;
        }
        return y;
    }

    using seed_type = typename std::mt19937_64::result_type;
    seed_type seed_;
    std::mt19937_64 random_;
    std::normal_distribution< float > normal_distribution_; // standard normal distribution

    static constexpr size_type BANDS = 5;
    static constexpr size_type NVERTICES = 20;
    static constexpr size_type NSAMPLES = 10000000;
    static_assert((NSAMPLES % NVERTICES) == 0, "!");

    std::vector< float3 > uniform_sphere[NVERTICES]; // uniformely distributed samples near the dodecahedron vertices on unit sphere

    static constexpr size_type max_size = NSAMPLES / NVERTICES;

    size_type min_size() const
    {
        size_type size = max_size;
        for (auto const & pyramid : uniform_sphere) {
            if (pyramid.size() < size) {
                size = pyramid.size();
            }
        }
        return size;
    }

    std::vector< float > sh;

    float const phi = (one + std::sqrt(5.0f)) / 2.0f; // golden ratio
    float const rphi = one / phi;
    float3 dodecahedron[NVERTICES] = {{one, one, one}, {one, -one, one}, {-one, one, one}, {-one, -one, one},
                                      {zero, rphi, phi}, {zero, -rphi, phi},
                                      {phi, zero, rphi}, {-phi, zero, rphi},
                                      {rphi, phi, zero}, {rphi, -phi, zero}, {-rphi, phi, zero}, {-rphi, -phi, zero},
                                      {phi, zero, -rphi}, {-phi, zero, -rphi},
                                      {zero, rphi, -phi}, {zero, -rphi, -phi},
                                      {one, one, -one}, {one, -one, -one}, {-one, one, -one}, {-one, -one, -one}}; // circumsphere radius = sqrt(3)
    //float const cone_cos = 0.9f;

    float dot_products[NVERTICES];

    size_type cone_index(float3 direction)
    {
        assert(!(length(direction) < -eps) && !(one + eps < length(direction)));
        for (size_type i = 0; i < NVERTICES; ++i) {
            dot_products[i] = dot(dodecahedron[i], direction);
        }
        auto const beg = std::cbegin(dot_products);
        auto const m = std::max_element(beg, std::cend(dot_products));
        return static_cast< size_type >(std::distance(beg, m));
    }

    double cosine[BANDS];
    double mean[NVERTICES][BANDS * BANDS];

    void operator () (std::ostream & gnuplot)
    {
        float const sqrt3 = std::sqrt(3.0f);
        for (float3 & vertex : dodecahedron) {
            vertex /= sqrt3;
        }
        for (auto & pyramid : uniform_sphere) {
            pyramid.reserve(max_size);
        }
        while (min_size() < max_size) {
            float3 point{normal_distribution_(random_), normal_distribution_(random_), normal_distribution_(random_)};
            float len = length(point);
            if (eps < len) {
                point /= len;
                auto & pyramid = uniform_sphere[cone_index(point)];
                if (pyramid.size() < max_size) {
                    pyramid.push_back(std::move(point));
                }
            }
        }
        for (auto & pyramid : uniform_sphere) { // truncate
            assert(!(pyramid.size() < max_size));
            pyramid.resize(max_size);
        }
        sh.reserve(BANDS * BANDS * NSAMPLES);
        std::cout << std::setprecision(10);
        for (size_type v = 0; v < NVERTICES; ++v) {
            //std::cout << '{';
            auto const & pyramid = uniform_sphere[v];
            auto & vmean = mean[v];
            size_type j = 0;
            for (int l = 0; l < BANDS; ++l) {
                for (int m = -l; m <= l; ++m) {
                    auto & cm = vmean[j];
                    cm = zero;
                    for (size_type s = 0; s < max_size; ++s) {
                        sh.push_back(SH(l, m, pyramid[s]));
                        cm += sh.back();
                    }
                    cm *= 1.0 / (5.0 * max_size);
                    //std::cout << cm << ", ";
                    ++j;
                }
            }
            //std::cout << "},\n";
        }
        assert(sh.size() == (BANDS * BANDS * NSAMPLES));
#if 1
        {
            for (auto & c : cosine) {
                c = zero;
            }
            size_type i = 0;
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                for (int l = 0; l < BANDS; ++l) {
                    auto & c = cosine[l];
                    i += l * max_size;
                    for (size_type s = 0; s < max_size; ++s) {
                        float dot_product = pyramid[s].z;
                        if (zero < dot_product) {
                            c += dot_product * sh[i];
                        }
                        ++i;
                    }
                    i += l * max_size;
                }
            }
            assert(i == sh.size());
            size_type k = 0;
            for (auto & c : cosine) {
                c *= ((4 * PI) / NSAMPLES);
                std::cout << c << ' ' << k++ << std::endl;
            }
        }
        for (auto const & vmean : mean) {
            std::cout << '{';
            size_type j = 0;
            for (int l = 0; l < BANDS; ++l) {
                auto & c = cosine[l];
                for (int m = -l; m <= l; ++m) {
                    std::cout << vmean[j++] * c << ", ";
                }
            }
            std::cout << "},\n";
        }
        std::cout << std::flush;
        // rotation:
        float3 direction{1.0f, 2.0f, 3.0f};
        direction /= length(direction);
        float rcosine[BANDS * BANDS];
        {
            size_type j = 0;
            for (int l = 0; l < BANDS; ++l) {
                auto const & c = cosine[l];
                for (int m = -l; m <= l; ++m) {
                    rcosine[j] = c * std::sqrt(4 * PI / (2 * l + 1)) * SH(l, m, direction);
                    //rcosine[j] = c * SH_(l, m, direction);
                    //std::cerr << std::abs(c * std::sqrt(4 * PI / (2 * l + 1)) * SH(l, m, direction) - rcosine[j]) << std::endl;
                    ++j;
                }
            }
        }
#endif
#if 1
        auto const print = [&] (float3 const & p) { gnuplot << p.x << ' ' << p.y << ' ' << p.z << '\n'; };

        auto const rot = rotateZ3(std::atan2(-direction.x, direction.y)) * rotateX3(std::atan2(-std::hypot(direction.x, direction.y), direction.z));
        gnuplot << "$cosine <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float3 const & point = pyramid[s];
                float c = zero;
                for (int l = 0; l < BANDS; ++l) {
                    c += cosine[l] * SH(l, 0, point);
                }
                print(rot * point * c);
            }
            gnuplot << '\n';
        }
        gnuplot << "EOD\n";
        gnuplot << "$rcosine <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float3 const & point = pyramid[s];
                float c = zero;
                size_type j = 0;
                for (int l = 0; l < BANDS; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        c += rcosine[j] * SH(l, m, point);
                        ++j;
                    }
                }
                print(point * c);
            }
            gnuplot << '\n';
        }
        gnuplot << "EOD\n";
#if 0
        gnuplot << "$mean <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            float const m = mean[v][2];
            for (size_type s = 0; s < max_size; ++s) {
                float3 point = pyramid[s] * std::abs(m);
                // l * (l + 1) + m
                print(point);
            }
            gnuplot << '\n';
        }
        gnuplot << "EOD\n";
        gnuplot << "$sh <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float const scale = SH(1, 0, pyramid[s]);
                float3 const & point = pyramid[s] * std::abs(scale);
                print(point);
            }
            gnuplot << '\n';
        }
        gnuplot << "EOD\n";
#endif
        gnuplot << "set term wxt\n"
                   "set view equal xyz\n"
                   "set autoscale\n"
                   "set key left\n"
                   "set xrange [-1:1]\n"
                   "set yrange [-1:1]\n"
                   "set zrange [-1:1]\n"
                   "set xyplane at -1\n";
        gnuplot << "set arrow 1 from 0,0,0 to 1,0,0 linecolor rgbcolor 'red'\n";
        gnuplot << "set arrow 2 from 0,0,0 to 0,1,0 linecolor rgbcolor 'green'\n";
        gnuplot << "set arrow 3 from 0,0,0 to 0,0,1 linecolor rgbcolor 'blue'\n";
        gnuplot << "set arrow 4 from 0,0,0 to "
                << direction.x << ',' << direction.y << ',' << direction.z
                << " linecolor rgbcolor 'black'\n";
        gnuplot << "splot '$cosine' with points pointtype 1"
                   ", '$rcosine' with points pointtype 1\n";
#endif
        gnuplot << std::flush;
    }

};

#include <fstream>

#include <cstdlib>

int main(int argc, char * argv [])
{
    (void(argc), void(argv));
    spherical_harmonics spherical_harmonics_;
    std::ofstream of("sh.plt");
    spherical_harmonics_(of);
    std::system("gnuplot -p sh.plt");
    return EXIT_SUCCESS;
}
