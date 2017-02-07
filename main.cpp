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

        float p = atan2f(pos.z, pos.x);
        float t = acosf(pos.y / len);

        return SH(l, m, t, p);
    }

    float SH_A(const int l, const int m, const float3 &pos)
    {
        float d = dot(pos, pos);
        float len = sqrtf(d);

        float p = atan2f(pos.z, pos.x);
        float t = acosf(pos.y / len);

        return SH(l, m, t, p) * powf(d, -1.5f);
    }

    using seed_type = typename std::mt19937_64::result_type;
    seed_type seed_;
    std::mt19937_64 random_;
    std::normal_distribution< float > normal_distribution_; // standard normal distribution

    static constexpr size_type BANDS = 6;
    static constexpr size_type NVERTICES = 20;
    static constexpr size_type NSAMPLES = 1000000;

    std::vector< float3 > uniform_sphere[NVERTICES]; // uniformely distributed samples near the dodecahedron vertices on unit sphere

    static constexpr size_type max_size = NSAMPLES / NVERTICES;

    size_type min_size() const
    {
        size_type min_size_ = max_size;
        for (auto const & pyramid_ : uniform_sphere) {
            if (pyramid_.size() < min_size_) {
                min_size_ = pyramid_.size();
            }
        }
        return min_size_;
    }

    std::vector< float > sh_;

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

    size_type d_index(float3 direction) const
    {
        assert(!(length(direction) < -eps) && !(one + eps < length(direction)));
        float cosines[NVERTICES];
        for (size_type i = 0; i < NVERTICES; ++i) {
            cosines[i] = dot(dodecahedron_[i], direction);
        }
        auto const beg = std::cbegin(cosines);
        auto const mm = std::minmax_element(beg, std::cend(cosines));
        return static_cast< size_type >(std::distance(beg, mm.second));
    }

    float const solid_angle = (4 * std::acos(-one)) / NVERTICES; // pi / 5

    float cosine[BANDS * BANDS];
    float mean[NVERTICES][BANDS * BANDS];

    void init()
    {
        {
            float const sqrt3 = std::sqrt(3.0f);
            for (float3 & vertex : dodecahedron_) {
                vertex /= sqrt3;
            }
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
        sh_.reserve(BANDS * BANDS * NSAMPLES);
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (int l = 0; l < BANDS; ++l) {
                for (int m = -l; m <= l; ++m) {
                    for (size_type s = 0; s < max_size; ++s) {
                        sh_.push_back(SH(l, m, pyramid[s]));
                    }
                }
            }
        }
        assert(sh_.size() == (BANDS * BANDS * NSAMPLES));
        {
            for (float & c : cosine) {
                c = zero;
            }
            size_type i = 0;
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                auto & vmean = mean[v];
                size_type j = 0;
                for (int l = 0; l < BANDS; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        float & c = cosine[j];
                        float & cm = (vmean[j] = zero);
                        for (size_type s = 0; s < max_size; ++s) {
                            if (zero < pyramid[s].z) {
                                c += pyramid[s].z * sh_[i];
                            }
                            cm += sh_[i];
                            ++i;
                        }
                        c *= (solid_angle / float(max_size));
                        std::cout << std::setw(14) << cm << " - " << v << ' ' << l << ' ' << m << std::endl;
                        ++j;
                    }
                }
            }
            float const _4pi = 4.0f * std::acos(-one);
            size_type k = 0;
            for (float & c : cosine) {
                c *= (_4pi / float(NSAMPLES));
                std::cout << c << ' ' << k++ << std::endl;
            }
        }
    }

};

#include <cstdlib>

int main(int argc, char * argv [])
{
    (void(argc), void(argv));
    spherical_harmonics spherical_harmonics_;
    spherical_harmonics_.init();
    return EXIT_SUCCESS;
}
