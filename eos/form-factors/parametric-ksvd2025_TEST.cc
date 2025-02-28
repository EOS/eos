/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Matthew Kirk
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
#include <eos/form-factors/parametric-ksvd2025.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricKSvD2025Test :
    public TestCase
{
    public:
        ParametricKSvD2025Test() :
            TestCase("parametric_KSvD2025_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            // t0 = -1
            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                    =  0.13957;
                p["mass::K_d"]                     =  0.497611;
                p["0->Kpi::t_0@KSvD2025"]         = -1.0;
                p["0->Kpi::b_+^1@KSvD2025"]   =  0.0;
                p["0->Kpi::b_+^2@KSvD2025"]   =  0.0;
                p["0->Kpi::b_+^3@KSvD2025"]   =  0.0;
                p["0->Kpi::b_+^4@KSvD2025"]   =  0.0;
                p["0->Kpi::b_+^5@KSvD2025"]   =  0.0;
                p["0->Kpi::b_+^6@KSvD2025"]   =  0;
                p["0->Kpi::b_+^7@KSvD2025"]   =  0;
                p["0->Kpi::b_+^8@KSvD2025"]   =  0;
                p["0->Kpi::b_+^9@KSvD2025"]   =  0;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->Kpi::KSvD2025", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }

                /* z mapping */
                {
                    KSvD2025FormFactors<VacuumToKPi> ff(p, Options{ });

                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(-2.0)),                      0.133502508, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(-2.0)),                      0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(-1.0)),                      0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(-1.0)),                      0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),                     -0.300926358,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),                      0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)),                     -0.363775160,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)),                     -0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)),                     -0.874666169,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)),                     -0.484725791,   eps);
                }

                /* f_+ at timelike q2 > 0.0 */
                {
                    KSvD2025FormFactors<VacuumToKPi> ff(p, Options{ });

                    const auto chi_1m = 0.5;

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z( 0.0), chi_1m)),  0.00221304637,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z( 0.0), chi_1m)),  0.0,            eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.1), chi_1m)),  0.000723952644, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.1), chi_1m)),  0.0,            eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.5), chi_1m)), -0.00512638515,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.5), chi_1m)),  0.00294207286,  eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z( 0.0), chi_1m)),  0.0225064602,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z( 0.0), chi_1m)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z(+0.1), chi_1m)),  0.0269438013,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z(+0.1), chi_1m)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z(+0.5), chi_1m)), -0.00056134248, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z(+0.5), chi_1m)), -0.00284024293, eps);
                }

                /* f_0 at timelike q2 > 0.0 */
                {
                    KSvD2025FormFactors<VacuumToKPi> ff(p, Options{ });

                    const auto chi_0p = 0.3;

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z( 0.0), chi_0p)),  0.00484747542, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z( 0.0), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z(+0.1), chi_0p)),  0.00322555580, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z(+0.1), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z(+0.5), chi_0p)),  0.00357029497, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z(+0.5), chi_0p)),  0.00365234921, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z( 0.0), chi_0p)),  0.0192623633,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z( 0.0), chi_0p)),  0.0,            eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z(+0.1), chi_0p)),  0.0417161044,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z(+0.1), chi_0p)),  0.0,            eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z(+0.5), chi_0p)),  0.000766164476, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z(+0.5), chi_0p)), -0.000103553419, eps);
                }
            }
        }
} parametric_KSvD2025_test;
