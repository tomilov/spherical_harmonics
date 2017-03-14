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

        if (4 < l) {
            return x;
        } else {
            float y = sqrt((2.0 * l + 1.0) / (4.0 * PI)) * SH_Cartesian(l, m, pos);
            if (std::abs(x - y) > 1E-5) {
                std::cerr << x << ' ' << y << std::endl;
            }
            return y;
        }
    }

    using seed_type = typename std::mt19937_64::result_type;
    seed_type seed_;
    std::mt19937_64 random_;
    std::normal_distribution< float > normal_distribution_; // standard normal distribution

    static constexpr size_type BANDS = 5;
    static constexpr size_type NVERTICES = 20;
    static constexpr size_type NFACETS = 12;
    static constexpr size_type NSAMPLES = NVERTICES * 1000;
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

    std::vector< float > sh_;

    float const phi = (one + std::sqrt(5.0f)) / 2.0f; // golden ratio
    float const rphi = one / phi;
    float3 dodecahedron[NVERTICES] = {
        {one, one, one}, {one, -one, one}, {-one, one, one}, {-one, -one, one},
        {zero, rphi, phi}, {zero, -rphi, phi},
        {phi, zero, rphi}, {-phi, zero, rphi},
        {rphi, phi, zero}, {rphi, -phi, zero}, {-rphi, phi, zero}, {-rphi, -phi, zero},
        {phi, zero, -rphi}, {-phi, zero, -rphi},
        {zero, rphi, -phi}, {zero, -rphi, -phi},
        {one, one, -one}, {one, -one, -one}, {-one, one, -one}, {-one, -one, -one}
    }; // circumsphere radius = sqrt(3)
    //float const cone_cos = 0.9f;

    float3 icosahedron[NFACETS] = {
        {-1.0, 0.0, phi}, { 1.0, 0.0,  phi}, {-1.0,  0.0, -phi}, { 1.0,  0.0, -phi},
        { 0.0, phi, 1.0}, { 0.0, phi, -1.0}, { 0.0, -phi,  1.0}, { 0.0, -phi, -1.0},
        { phi, 1.0, 0.0}, {-phi, 1.0,  0.0}, { phi, -1.0,  0.0}, {-phi, -1.0,  0.0}
    };

    size_type facets[NVERTICES][3] = {
    #if 0
        {1,  4, 0}, { 4, 9, 0}, {4,  5, 9}, {8, 5,  4}, { 1, 8, 4},
        {1, 10, 8}, {10, 3, 8}, {8,  3, 5}, {3, 2,  5}, { 3, 7, 2},
        {3, 10, 7}, {10, 6, 7}, {6, 11, 7}, {6, 0, 11}, { 6, 1, 0},
        {10, 1, 6}, {11, 0, 9}, {2, 11, 9}, {5, 2,  9}, {11, 2, 7},
    #else
        {1, 8, 4}, {10,  1, 6}, { 4,  9, 0}, {6,  0, 11}, { 1, 4, 0},
        {6, 1, 0}, { 1, 10, 8}, {11,  0, 9}, {8,  5,  4}, {10, 6, 7},
        {4, 5, 9}, { 6, 11, 7}, {10,  3, 8}, {2, 11,  9}, { 3, 2, 5},
        {3, 7, 2}, { 8,  3, 5}, { 3, 10, 7}, {5,  2,  9}, {11, 2, 7}
    #endif
    };

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

    double mean[NVERTICES][BANDS * BANDS];

    void operator () (std::ostream & sh, std::ostream & ps)
    {
        {
            float const sqrt3 = length(dodecahedron[0]);
            for (float3 & vertex : dodecahedron) {
                vertex /= sqrt3;
            }
        }
        {
            float const l = length(icosahedron[0]);
            for (float3 & vertex : icosahedron) {
                vertex /= l;
            }
        }
        {
            auto const print = [&] (float3 const & p, size_type const i = 0)
            {
                ps << p.x << ' ' << p.y << ' ' << p.z;
                if (0 < i) {
                    ps << ' ' << (i - 1);
                }
                ps << '\n';
            };
            ps << "set term wxt\n"
                  "set view equal xyz\n"
                  "set autoscale\n"
                  "set key left\n"
                  "set xrange [*:*]\n"
                  "set yrange [*:*]\n"
                  "set zrange [*:*]\n"
                  "set xyplane at -1\n";
            ps << "$icosahedron << EOD\n";
            for (size_type i = 0; i < NFACETS; ++i) {
                print(icosahedron[i], i + 1);
            }
            ps << "EOD\n";
            ps << "$dodecahedron << EOD\n";
            for (size_type i = 0; i < NVERTICES; ++i) {
                print(dodecahedron[i], i + 1);
            }
            ps << "EOD\n";
            ps << "$ifaces << EOD\n";
            for (size_type i = 0; i < NVERTICES; ++i) {
                auto const & facet = facets[i];
                size_type const u = facet[0];
                print(icosahedron[u], u + 1);
                size_type const v = facet[1];
                print(icosahedron[v], v + 1);
                size_type const w = facet[2];
                print(icosahedron[w], w + 1);
                print(icosahedron[u], 0);
                ps << "\n\n";
            }
            ps << "EOD\n";
            float3 cfaces[NVERTICES];
            ps << "$cfaces << EOD\n";
            for (size_type i = 0; i < NVERTICES; ++i) {
                auto const & facet = facets[i];
                float3 const center = (icosahedron[facet[0]] + icosahedron[facet[1]] + icosahedron[facet[2]]) / 3.0;
                print(center, i + 1);
                cfaces[i] = center;
            }
            ps << "EOD\n";
            constexpr size_type const strip[NVERTICES] = {0,  4,  5,  3, 11, 19, 13,  7,  2, 10, 8, 16, 12,  6,  1, 9, 17, 15, 14, 18};
            ps << "$strip << EOD\n";
            for (auto const i : strip) {
                print(cfaces[i], 0);
            }
            ps << "EOD\n";
            ps << "splot '$icosahedron' with points pointtype 1"
                  //<< ", '' with labels offset character 0, character 1 notitle"
               << ", '$ifaces' with lines"
                  //<< ", '' with labels offset character 0, character -1 notitle"
               << ", '$cfaces' with points pointtype 1"
               //<< ", '' with labels offset character 0, character 1 notitle"
               << ", '$dodecahedron' with points pointtype 1"
               << ", '' with labels offset character 0, character 1 notitle"
               << ", '$strip' with lines"
               << std::endl;
            {
                std::set< size_type > prev, next;
                {
                    auto const & facet = facets[strip[0]];
                    prev.insert({facet[0], facet[1], facet[2]});
                    std::cout << "first triple: "
                              << facet[0] << ' '
                              << facet[1] << ' '
                              << facet[2] << std::endl;
                }
                for (size_type i = 1; i < NVERTICES; ++i) {
                    auto const & facet = facets[strip[i]];
                    next.insert({facet[0], facet[1], facet[2]});
                    size_type next_vertex = NVERTICES;
                    std::set_difference(std::cbegin(prev), std::cend(prev),
                                        std::cbegin(next), std::cend(next),
                                        &next_vertex);
                    std::cout << next_vertex;
                    std::set_difference(std::cbegin(next), std::cend(next),
                                        std::cbegin(prev), std::cend(prev),
                                        &next_vertex);
                    std::cout << " replaced with " << next_vertex << std::endl;
                    prev = std::move(next);
                }
                // 1, 4, 0, 6, 11, 7, 2, 9, 0, 4, 5, 8, 3, 10, 1, 6, 7, 3, 2, 5, 9
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
        sh_.reserve(BANDS * BANDS * NSAMPLES);
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
                        sh_.push_back(SH(l, m, pyramid[s]));
                        cm += sh_.back();
                    }
                    cm *= 1.0 / (5.0 * max_size);
                    //std::cout << cm << ", ";
                    ++j;
                }
            }
            //std::cout << "},\n";
        }
        assert(sh_.size() == (BANDS * BANDS * NSAMPLES));

        auto const print = [&] (float3 const & p) { sh << p.x << ' ' << p.y << ' ' << p.z << '\n'; };
        float3 direction{1.0f, 2.0f, 3.0f};
        direction /= length(direction);
#if 1
        sh << "set term wxt\n"
              "set view equal xyz\n"
              "set autoscale\n"
              "set key left\n"
              "set xrange [-1:1]\n"
              "set yrange [-1:1]\n"
              "set zrange [-1:1]\n"
              "set xyplane at -1\n";
        sh << "set arrow 1 from 0,0,0 to 1,0,0 linecolor rgbcolor 'red'\n";
        sh << "set arrow 2 from 0,0,0 to 0,1,0 linecolor rgbcolor 'green'\n";
        sh << "set arrow 3 from 0,0,0 to 0,0,1 linecolor rgbcolor 'blue'\n";
        sh << "set arrow 4 from 0,0,0 to "
           << direction.x << ',' << direction.y << ',' << direction.z
           << " linecolor rgbcolor 'black'\n";

        {
            constexpr size_type i = 2;
            sh << "$mean <<EOD\n";
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                float const m = mean[v][i];
                for (size_type s = 0; s < max_size; ++s) {
                    float3 point = pyramid[s] * std::abs(m);
                    // l * l + m
                    print(point);
                }
                sh << '\n';
            }
            sh << "EOD\n";
        }
        {
            sh << "$sh <<EOD\n";
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                for (size_type s = 0; s < max_size; ++s) {
                    float const scale = SH(1, 0, pyramid[s]);
                    float3 const & point = pyramid[s] * std::abs(scale);
                    print(point);
                }
                sh << '\n';
            }
            sh << "EOD\n";
        }

#if 0
        float cone[BANDS * BANDS];
        std::deque< float3 > true_cone;
#if 0
        {
            for (auto & c : cone) {
                c = zero;
            }
            size_type k = 0;
            for (size_type v = 0; v < NVERTICES; ++v) {
                assert(k == BANDS * BANDS * max_size * v);
                auto const & pyramid = uniform_sphere[v];
                size_type i = 0;
                for (int l = 0; l < BANDS; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        auto & c = cone[i];
                        for (size_type s = 0; s < max_size; ++s) {
                            float dot_product = dot(pyramid[s], direction);
                            if (0.9 <= dot_product) {
                                c += sh[k];
                                true_cone.push_back(pyramid[s]);
                            }
                            ++k;
                        }
                        ++i;
                    }
                }
                assert(i == BANDS * BANDS);
            }
            assert(k == sh.size());
            for (auto & c : cone) {
                c *= ((4 * PI) / NSAMPLES);
            }
        }
#else
        sh << "unset arrow 4\n";
        {
            for (auto & c : cone) {
                c = zero;
            }
            constexpr size_type v = 0;
            assert(v < NVERTICES);
            size_type k = BANDS * BANDS * max_size * v;
            auto const & pyramid = uniform_sphere[v];
            size_type i = 0;
            for (int l = 0; l < BANDS; ++l) {
                for (int m = -l; m <= l; ++m) {
                    auto & c = cone[i];
                    for (size_type s = 0; s < max_size; ++s) {
                        c += sh[k];
                        true_cone.push_back(pyramid[s]);
                        ++k;
                    }
                    ++i;
                }
            }
            assert(i == BANDS * BANDS);
            assert(k == BANDS * BANDS * max_size * (v + 1));
            for (auto & c : cone) {
                c *= ((4 * PI) / NSAMPLES);
            }
        }
#endif
        {
            std::deque< float3 > cone_neg;
            sh << "$cone <<EOD\n";
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                for (size_type s = 0; s < max_size; ++s) {
                    float3 point = pyramid[s];
                    float c = zero;
                    size_type i = 0;
                    for (int l = 0; l < BANDS; ++l) {
                        for (int m = -l; m <= l; ++m) {
                            c += cone[i] * SH(l, m, point);
                            ++i;
                        }
                    }
                    point *= c;
                    if (zero < c) {
                        print(point);
                    } else {
                        cone_neg.push_back(point);
                    }
                }
                sh << '\n';
            }
            sh << "EOD\n";

            sh << "$cone_neg <<EOD\n";
            for (float3 const & point : cone_neg) {
                print(point);
                sh << '\n';
            }
            sh << "EOD\n";
        }
        {
            sh << "$true_cone <<EOD\n";
            for (float3 const & point : true_cone) {
                print(point);
                sh << '\n';
            }
            sh << "EOD\n";
        }
        sh << "splot '$cone' with points pointtype 1"
                   ", '$cone_neg' with points pointtype 1"
                   ", '$true_cone' with points pointtype 1\n";
#else
        double cosine[BANDS];
        {
            for (auto & c : cosine) {
                c = zero;
            }
            size_type k = 0;
            for (size_type v = 0; v < NVERTICES; ++v) {
                auto const & pyramid = uniform_sphere[v];
                for (int l = 0; l < BANDS; ++l) {
                    auto & c = cosine[l];
                    k += l * max_size;
                    for (size_type s = 0; s < max_size; ++s) {
                        float dot_product = pyramid[s].z;
                        if (zero < dot_product) {
                            c += dot_product * sh_[k];
                        }
                        ++k;
                    }
                    k += l * max_size;
                }
            }
            assert(k == sh_.size());
            size_type l = 0;
            for (auto & c : cosine) {
                c *= ((4 * PI) / NSAMPLES);
                std::cout << c << ' ' << l++ << std::endl;
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

        auto const rot = rotateZ3(std::atan2(-direction.x, direction.y)) * rotateX3(std::atan2(-std::hypot(direction.x, direction.y), direction.z));
        sh << "$cosine <<EOD\n";
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
            sh << '\n';
        }
        sh << "EOD\n";
        sh << "$rcosine <<EOD\n";
        for (size_type v = 0; v < NVERTICES; ++v) {
            auto const & pyramid = uniform_sphere[v];
            for (size_type s = 0; s < max_size; ++s) {
                float3 const & point = pyramid[s];
                float c = zero;
                size_type i = 0;
                for (int l = 0; l < BANDS; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        c += rcosine[i] * SH(l, m, point);
                        ++i;
                    }
                }
                print(point * c);
            }
            sh << '\n';
        }
        sh << "EOD\n";
        sh << "splot '$cosine' with points pointtype 1"
              ", '$rcosine' with points pointtype 1\n";
#endif
#endif
        sh << std::flush;
    }

};

#include <fstream>

#include <cstdlib>

int main(int argc, char * argv [])
{
    (void(argc), void(argv));
    {
        spherical_harmonics spherical_harmonics_;
        std::ofstream sh("sh.plt");
        std::ofstream ps("ps.plt");
        spherical_harmonics_(sh, ps);
    } // flush files before using them
    if (std::system("gnuplot -p ps.plt") != 0) {
        return EXIT_FAILURE;
    }
    if (std::system("gnuplot -p sh.plt") != 0) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
