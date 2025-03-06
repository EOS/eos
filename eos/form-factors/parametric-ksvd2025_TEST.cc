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
                // K*(892)
                p["0->Kpi::M_(+,0)@KSvD2025"]     =  0.890;
                p["0->Kpi::Gamma_(+,0)@KSvD2025"] =  0.026;
                // K*(1410)
                p["0->Kpi::M_(+,1)@KSvD2025"]     =  1.368;
                p["0->Kpi::Gamma_(+,1)@KSvD2025"] =  0.106;
                // K0*(700)
                p["0->Kpi::M_(0,0)@KSvD2025"]     =  0.680;
                p["0->Kpi::Gamma_(0,0)@KSvD2025"] =  0.300;
                 // K0*(1430)
                // p["0->Kpi::M_(0,1)@KSvD2025"]     =  1.431;
                // p["0->Kpi::Gamma_(0,1)@KSvD2025"] =  0.110;
                p["0->Kpi::b_+^1@KSvD2025"]   =  0.5;
                p["0->Kpi::b_+^2@KSvD2025"]   =  0.1;
                p["0->Kpi::b_+^3@KSvD2025"]   =  0.05;
                p["0->Kpi::b_0^1@KSvD2025"]   =  0.3;
                p["0->Kpi::b_0^2@KSvD2025"]   =  0.08;
                p["0->Kpi::b_0^3@KSvD2025"]   =  0.01;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->Kpi::KSvD2025", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }

                /* z mapping */
                {
                    KSvD2025FormFactors<VacuumToKPi> ff(p, Options{ });

                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(-2.0)),                      0.133502508,   eps);
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

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z( 0.0), chi_1m)),  0.0129766667,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z( 0.0), chi_1m)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.1), chi_1m)),  0.0118136429,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.1), chi_1m)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.5), chi_1m)),  0.00511197788,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.5), chi_1m)), -0.00296703547,  eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z( 0.0), chi_1m)),  0.0198271348,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z( 0.0), chi_1m)),  0.0,          eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z(+0.1), chi_1m)),  0.0172636701,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z(+0.1), chi_1m)),  0.0,          eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_p(ff.z(+0.5), chi_1m)),  0.00465459238, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_p(ff.z(+0.5), chi_1m)), -0.00398567803, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_p(ff.z( 0.0))),  1.04722073,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_p(ff.z( 0.0))),  0.0,         eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_p(ff.z( 0.1))),  1.05179947,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_p(ff.z( 0.1))),  0.0,         eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_p(ff.z( 0.5))),  0.493194311, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_p(ff.z( 0.5))), -0.684494289, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_productprime_p(ff.z( 0.0))), -0.185712601,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_productprime_p(ff.z( 0.0))),  0.0,          eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_productprime_p(ff.z( 0.1))),  0.0419145052, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_productprime_p(ff.z( 0.1))),  0.0,          eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_productprime_p(ff.z( 0.5))),  1.95390159,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_productprime_p(ff.z( 0.5))), -1.30168386,   eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.b0_fp(), -0.0999435, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.0)), -19.5870533, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.0)),   0.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.1)), -24.1282649, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.1)),   0.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.5)), -75.5346359, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.5)),   1.3358251, eps);
                }

                /* f_0 at timelike q2 > 0.0 */
                {
                    KSvD2025FormFactors<VacuumToKPi> ff(p, Options{ {"n-resonances-0p", "1"} });

                    const auto chi_0p = 0.3;

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z( 0.0), chi_0p)),  0.00874124234, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z( 0.0), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z(+0.1), chi_0p)),  0.00818142811, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z(+0.1), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_z(ff.z(+0.5), chi_0p)),  0.00494611070, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_z(ff.z(+0.5), chi_0p)), -0.00127383276, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z( 0.0), chi_0p)),  0.00955439874, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z( 0.0), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z(+0.1), chi_0p)),  0.00829786553, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z(+0.1), chi_0p)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitildeprime_z(ff.z(+0.5), chi_0p)),  0.00195495683, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitildeprime_z(ff.z(+0.5), chi_0p)), -0.00229014562, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_z(ff.z( 0.0))),  0.547815118, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_z(ff.z( 0.0))),  0.0,         eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_z(ff.z( 0.1))),  0.587349411, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_z(ff.z( 0.1))),  0.0,         eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.resonance_product_z(ff.z( 0.5))),  1.004483094, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.resonance_product_z(ff.z( 0.5))),  0.556192495, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.b0_f0(), -0.229236, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_0(0.0)), -19.5870533, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_0(0.0)),   0.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_0(0.1)), -23.5662156, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_0(0.1)),   0.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_0(0.5)), -59.9125015, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_0(0.5)), -83.7890488, eps);
                }
            }
        }
} parametric_KSvD2025_test;
