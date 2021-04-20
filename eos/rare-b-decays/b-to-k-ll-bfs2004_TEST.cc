/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
 * Copyright (c) 2021 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/observable.hh>
#include <eos/rare-b-decays/b-to-k-ll.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BToKDileptonBFS2004BobethCompatibilityTest :
    public TestCase
{
    public:
        BToKDileptonBFS2004BobethCompatibilityTest() :
            TestCase("b_to_k_dilepton_BFS2004_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            // Christoph uses \Delta C instead of C for C9, C10
            // important to agree on alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["b->s::c1"] = -0.3231323312;
            p["b->s::c2"] = 1.009301831;
            p["b->s::c3"] = -0.005233499106;
            p["b->s::c4"] = -0.08829686414;
            p["b->s::c5"] = 0.0003601965805;
            p["b->s::c6"] = 0.001020749573;
            p["b->s::Re{c7}"] = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"] = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"] = -0.1827530948;
            p["b->smumu::Re{c9}"] = 4.294489364 + 1;
            p["b->smumu::Im{c9}"] = 0.5;
            p["b->smumu::Re{c9'}"] = 2;
            p["b->smumu::Im{c9'}"] = 1.5;
            p["b->smumu::Re{c10}"] = -4.196294696 + 3;
            p["b->smumu::Im{c10}"] = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["b->smumu::Re{cS}"] = 0.5;
            p["b->smumu::Im{cS}"] = 1;
            p["b->smumu::Re{cS'}"] = 0.6;
            p["b->smumu::Im{cS'}"] = 1.1;
            p["b->smumu::Re{cP}"] = 0.7;
            p["b->smumu::Im{cP}"] = 1.2;
            p["b->smumu::Re{cP'}"] = 0.8;
            p["b->smumu::Im{cP'}"] = 1.3;
            p["b->smumu::Re{cT}"] = 0.9;
            p["b->smumu::Im{cT}"] = 1.4;
            p["b->smumu::Re{cT5}"] = 1.0;
            p["b->smumu::Im{cT5}"] = 1.5;
            p["K::a_1@1GeV"] = 0.1;
            p["K::a_2@1GeV"] = 0.1;
            p["B::1/lambda_B_p"] = 1.0 / 0.485;

            Options oo;
            oo.set("model", "WilsonScan");
            oo.set("scan-mode", "cartesian");
            oo.set("tag", "BFS2004");
            oo.set("qcdf-integrals", "mixed");
            oo.set("form-factors", "KMPW2010");
            oo.set("l", "mu");
            oo.set("q", "u");

            double eps = 1e-3;
            static const double s = 6.0;

            BToKDilepton d(p, oo);
            auto amps = d.amplitudes(s);
            std::array<double, 3> a = d.angular_coefficients(s);

            TEST_CHECK_RELATIVE_ERROR_C(amps.F_A, complex<double>(2.803705304, 6), 1e-14);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_S, complex<double>(3.277235546, 6.256540588), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_T, complex<double>(7.695315895, 11.97049139), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_T5, complex<double>(8.550350995, 12.82552649), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.F_P, complex<double>(4.010492477, 6.467135768), eps);

            /* difference comes from cal_T, F_V affects everything below */
            TEST_CHECK_RELATIVE_ERROR(std::real(amps.F_V), 7.756362368, eps);
            TEST_CHECK_RELATIVE_ERROR(std::imag(amps.F_V), 3.191642172, 6 * eps);

            eps *= 2.5;
            TEST_CHECK_RELATIVE_ERROR(a[0],  3.92053702e-20, eps);
            TEST_CHECK_RELATIVE_ERROR(a[1],  9.694697008e-21, eps);
            TEST_CHECK_RELATIVE_ERROR(a[2], -2.756810607e-20, eps);

            const double tau_over_hbar = p["life_time::B_u"] / p["QM::hbar"];
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(1, 6),
                                      2.898727023e-19 * tau_over_hbar, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio_cp_averaged(1, 6),
                                      2.8855929e-19 * tau_over_hbar, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_forward_backward_asymmetry(1, 6), 0.1097985735, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_flat_term(1, 6), 0.2788261376, eps);
            TEST_CHECK_RELATIVE_ERROR(d.integrated_cp_asymmetry(1, 6), 0.00455162022, 8 * eps);
        }
} b_to_k_dilepton_BFS2004_bobeth_compatibility_test;
