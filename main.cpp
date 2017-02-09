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
#include <fstream>

#include <cmath>
#include <cassert>

struct spherical_harmonics
{

    using size_type = std::size_t;

    float const eps = std::numeric_limits< float >::epsilon();
    float const one = 1.0f;
    float const zero = 0.0f;

    float factorial(const int x)
    {
        float f = 1.0f;
        for (int i = 2; i <= x; i++)
        {
            f *= i;
        }

        return f;
    }

    // Evaluate an Associated Legendre Polynomial P(l, m, x) at x
    float P(const int l, const int m, const float x)
    {
        float pmm = 1.0f;
        if (m > 0)
        {
            float somx2 = sqrtf((1.0f - x) * (1.0f + x));

            float fact = 1.0;
            for (int i = 1; i <= m; i++)
            {
                pmm *= (-fact) * somx2;
                fact += 2.0;
            }
        }
        if (l == m) return pmm;

        float pmmp1 = x * (2.0f * m + 1.0f) * pmm;
        if (l == m + 1) return pmmp1;

        float pll = 0.0;
        for (int ll = m + 2; ll <= l; ++ll)
        {
            pll = ((2.0f * ll - 1.0f) * x * pmmp1 - (ll + m - 1.0f) * pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pll;
    }

    // Normalization constant
    float K(const int l, const int m)
    {
        return sqrtf(((2.0f * l + 1.0f) * factorial(l - m)) / (4.0f * PI * factorial(l + m)));
    }

    // SH coefficient computation
    float SH(const int l, const int m, const float theta, const float phi)
    {
        const float sqrt2 = 1.4142135623731f;

        if (m == 0)
            return K(l, 0) * P(l, m, cosf(theta));
        else if (m > 0)
            return sqrt2 * K(l, m) * cosf(m * phi) * P(l, m, cosf(theta));
        else
            return sqrt2 * K(l, -m) * sinf(-m * phi) * P(l, -m, cosf(theta));
    }

    float SH(const int l, const int m, const float3 &pos)
    {
        float len = length(pos);

        float phi = atan2f(pos.y, pos.x);
        float t = pos.z / len;

        const float sqrt2 = 1.4142135623731f;

        if (m == 0)
            return K(l, 0) * P(l, m, t);
        else if (m > 0)
            return sqrt2 * K(l, m) * cosf(m * phi) * P(l, m, t);
        else
            return sqrt2 * K(l, -m) * sinf(-m * phi) * P(l, -m, t);
    }

    float SH_A(const int l, const int m, const float3 &pos)
    {
        float d = dot(pos, pos);
        float len = sqrtf(d);

        float p = atan2f(pos.y, pos.x);
        float t = acosf(pos.z / len);

        return SH(l, m, t, p) * powf(d, -1.5f);
    }

    using seed_type = typename std::mt19937_64::result_type;
    seed_type seed_;
    std::mt19937_64 random_;
    std::normal_distribution< float > normal_distribution_; // standard normal distribution

    static constexpr size_type BANDS = 6;
    static constexpr size_type NVERTICES = 20;
    static constexpr size_type NSAMPLES = 1000;
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
    float3 dodecahedron_[NVERTICES] = {{one, one, one}, {one, -one, one}, {-one, one, one}, {-one, -one, one},
                                       {zero, rphi, phi}, {zero, -rphi, phi},
                                       {phi, zero, rphi}, {-phi, zero, rphi},
                                       {rphi, phi, zero}, {rphi, -phi, zero}, {-rphi, phi, zero}, {-rphi, -phi, zero},
                                       {phi, zero, -rphi}, {-phi, zero, -rphi},
                                       {zero, rphi, -phi}, {zero, -rphi, -phi},
                                       {one, one, -one}, {one, -one, -one}, {-one, one, -one}, {-one, -one, -one}}; // circumsphere radius = sqrt(3)
    //float const cone_cos = 0.9f;

    float dot_products[NVERTICES];

    size_type d_index(float3 direction)
    {
        assert(!(length(direction) < -eps) && !(one + eps < length(direction)));
        for (size_type i = 0; i < NVERTICES; ++i) {
            dot_products[i] = dot(dodecahedron_[i], direction);
        }
        auto const beg = std::cbegin(dot_products);
        auto const m = std::max_element(beg, std::cend(dot_products));
        return static_cast< size_type >(std::distance(beg, m));
    }

    float cosine[BANDS];
    float mean[NVERTICES][BANDS * BANDS];

    void operator () ()
    {
        float const sqrt3 = std::sqrt(3.0f);
        for (float3 & vertex : dodecahedron_) {
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
                uniform_sphere[d_index(point)].push_back(std::move(point));
            }
        }
        for (auto & pyramid : uniform_sphere) { // truncate
            assert(!(pyramid.size() < max_size));
            pyramid.resize(max_size);
        }
        sh.reserve(BANDS * BANDS * NSAMPLES);
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            auto & vmean = mean[v];
            size_type j = 0;
            for (int l = 0; l < BANDS; ++l) {
                for (int m = -l; m <= l; ++m) {
                    float & cm = (vmean[j] = zero);
                    for (size_type s = 0; s < max_size; ++s) {
                        sh.push_back(SH(l, m, pyramid[s]));
                        cm += sh.back();
                    }
                    cm /= float(max_size);
                    //std::cout << std::setw(14) << cm << " - " << v << ' ' << l << ' ' << m << std::endl;
                    ++j;
                }
            }
        }
        assert(sh.size() == (BANDS * BANDS * NSAMPLES));
        {
            for (float & c : cosine) {
                c = zero;
            }
            size_type i = 0;
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                for (int l = 0; l < BANDS; ++l) {
                    float & c = cosine[l];
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
            size_type k = 0;
            for (float & c : cosine) {
                c *= ((4 * PI) / NSAMPLES);
                std::cout << c << ' ' << k++ << std::endl;
            }
        }
        // rotation:
        float3 direction{1.0f, 2.0f, 5.0f};
        direction /= length(direction);
        float rcosine[BANDS * BANDS];
        {
            size_type j = 0;
            for (int l = 0; l < BANDS; ++l) {
                float const & c = cosine[l];
                for (int m = -l; m <= l; ++m) {
                    rcosine[j] = c * std::sqrt(4 * PI / (2 * l + 1)) * SH(l, m, direction);
                    ++j;
                }
            }
        }
#if 1
        std::ofstream of("sh.plt");
        std::ostream & gnuplot = of;
        auto const print = [&] (float3 const & p) { gnuplot << p.x << ' ' << p.y << ' ' << p.z << '\n'; };

        auto const rot = rotateZ3(std::atan2(-direction.x, direction.y)) * rotateX3(std::atan2(-std::hypot(direction.x, direction.y), direction.z));
        gnuplot << "$cosine <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float3 point = pyramid[s];
                float c = zero;
                for (int l = 0; l < BANDS; ++l) {
                    c += cosine[l] * SH(l, 0, point);
                }
                point *= c;
                point = rot * point;
                print(point);
            }
            gnuplot << '\n';
        }
        gnuplot << "EOD\n";
        gnuplot << "$rcosine <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float3 point = pyramid[s];
                float c = zero;
                size_type j = 0;
                for (int l = 0; l < BANDS; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        c += rcosine[j] * SH(l, m, point);
                        ++j;
                    }
                }
                point *= c;
                print(point);
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
        gnuplot << "set view equal xyz\n"
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
    }

};

#include <cstdlib>

int main(int argc, char * argv [])
{
    (void(argc), void(argv));
    spherical_harmonics spherical_harmonics_;
    spherical_harmonics_();
    return EXIT_SUCCESS;
}
