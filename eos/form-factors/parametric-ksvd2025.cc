/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/form-factors/parametric-ksvd2025.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/maths/integrate.hh>

#include <functional>
#include <numeric>

namespace eos
{
    /* Vacuum -> pi pi */
    KSvD2025FormFactors<VacuumToKPi>::KSvD2025FormFactors(const Parameters & p, const Options & /*o*/) :
        _b_fp{{
            UsedParameter(p[_coeff_name("+", "1")], *this),
            UsedParameter(p[_coeff_name("+", "2")], *this),
            UsedParameter(p[_coeff_name("+", "3")], *this),
            UsedParameter(p[_coeff_name("+", "4")], *this),
            UsedParameter(p[_coeff_name("+", "5")], *this),
            UsedParameter(p[_coeff_name("+", "6")], *this),
            UsedParameter(p[_coeff_name("+", "7")], *this),
            UsedParameter(p[_coeff_name("+", "8")], *this),
            UsedParameter(p[_coeff_name("+", "9")], *this)
        }},
        _M_fp{{
            UsedParameter(p["0->Kpi::M_(+,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::M_(+,1)@KSvD2025"], *this)
        }},
        _G_fp{{
            UsedParameter(p["0->Kpi::Gamma_(+,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::Gamma_(+,1)@KSvD2025"], *this)
        }},
        _b_fz{{
            UsedParameter(p[_coeff_name("0", "1")], *this),
            UsedParameter(p[_coeff_name("0", "2")], *this),
            UsedParameter(p[_coeff_name("0", "3")], *this),
            UsedParameter(p[_coeff_name("0", "4")], *this),
            UsedParameter(p[_coeff_name("0", "5")], *this),
            UsedParameter(p[_coeff_name("0", "6")], *this),
            UsedParameter(p[_coeff_name("0", "7")], *this),
            UsedParameter(p[_coeff_name("0", "8")], *this),
            UsedParameter(p[_coeff_name("0", "9")], *this)
        }},
        _M_fz{{
            UsedParameter(p["0->Kpi::M_(0,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::M_(0,1)@KSvD2025"], *this)
        }},
        _G_fz{{
            UsedParameter(p["0->Kpi::Gamma_(0,0)@KSvD2025"], *this),
            UsedParameter(p["0->Kpi::Gamma_(0,1)@KSvD2025"], *this)
        }},
        _m_K(p["mass::K_d"], *this),
        _m_pi(p["mass::pi^-"], *this),
        _t_0(p["0->Kpi::t_0@KSvD2025"], *this),
        _hbar(p["QM::hbar"], *this)
    {
    }

    KSvD2025FormFactors<VacuumToKPi>::~KSvD2025FormFactors() = default;

    FormFactors<VacuumToPP> *
    KSvD2025FormFactors<VacuumToKPi>::make(const Parameters & p, const Options & o)
    {
        return new KSvD2025FormFactors<VacuumToKPi>(p, o);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::z(const complex<double> & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::dzdq2(const complex<double> & q2) const
    {
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return -sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) + sqrt(t_p - t_0)));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::dzdq2_II(const complex<double> & q2) const
    {
        // dz/dq2 on the second Riemann sheet
        const double t_p = this->_t_p();
        const double t_0 = this->_t_0();
        return sqrt(t_p - t_0) / (sqrt(t_p - q2) * power_of<2>(sqrt(t_p - q2) - sqrt(t_p - t_0)));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::series_m(const complex<double> & z, const std::array<double, 10u> & c) const
    {
        std::array<complex<double>, 10> zvalues;
        complex<double> current_zv = 1.0;
        for (complex<double> & zv: zvalues)
        {
            zv = current_zv;
            current_zv *= z;
        }

        return std::inner_product(c.cbegin(), c.cend(), zvalues.cbegin(), complex<double>(0.0, 0.0));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::w_p(const complex<double> & z) const
    {
        return power_of<2>(1.0 + z) * pow((1.0 - z), 5.0 / 2.0);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitilde_p(const complex<double> & z, const double & chi_1m) const
    {
        // The weight function ``(1.0 + z)^2 * (1.0 - z)^(+5/2)`` has been cancelled against the outer function
        // to remove unphysical singularities and correct the asymptotic behaviour
        const double t_p      = this->_t_p();
        const double t_0      = this->_t_0();
        const double t_m      = this->_t_m();
        const double tfactor  = 1.0 - t_0 / t_p;
        const double Q2       = 1.0;
        const double Q2factor = 1.0 + Q2 / t_p;
        const complex<double> zfactor = (1.0 + z) / (1.0 - z);

        return pow(t_p, -5.0/4.0)
            / sqrt(32.0 * M_PI * chi_1m)
            * pow(1.0 - z, -3.0) * pow(1.0 + z, -3.0/2.0) * sqrt(tfactor) * (1.0 + zfactor)
            * pow( power_of<2>(1.0 + z) * tfactor / power_of<4>(1.0 - z)
                   * (4.0 * t_p * z + t_m * power_of<2>(1.0 - z) - t_0 * power_of<2>(1.0 + z)), 3.0/4.0)
            / power_of<2>(zfactor * sqrt(tfactor) + 1.0)
            / power_of<3>(zfactor * sqrt(tfactor) + sqrt(Q2factor));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::resonance_product_p(const complex<double> & z) const
    {
        // Product of the resonance factors for f+
        // TODO
        return 1.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitildeprime_p(const complex<double> & z, const double & chi_1m) const
    {
        // Derivative of the modified outer function phitilde_+ with respect to z
        const double t_p      = this->_t_p();
        const double t_0      = this->_t_0();
        const double t_m      = this->_t_m();
        const double tfactor  = 1.0 - t_0 / t_p;
        const double Q2       = 1.0;
        const double Q2factor = 1.0 + Q2 / t_p;
        const complex<double> zfactor = (1.0 + z) / (1.0 - z);

        return -1.0 *(
            ( pow(tfactor, 1.5)*sqrt(zfactor) * (
                sqrt(tfactor) * ( t_0*power_of<2>(1.0 + z)*(5.0 + 6.0*z - 11.0*power_of<2>(z)
                    + sqrt(tfactor)*(-3.0 + 8.0*z + 11.0*power_of<2>(z)))
                    - 2*t_p*( 3.0 + z + 21.0*power_of<2>(z) - 25.0*power_of<3>(z)
                                + sqrt(tfactor)*(3.0 - 9.0*z + 13.0*power_of<2>(z) + 25.0*power_of<3>(z))
                            )
                    + t_m * (1.0 - 12.0*z + 11.0*power_of<2>(z) + sqrt(tfactor)*(9.0 - 2.0*z - 11.0*power_of<2>(z)))
                        * power_of<2>(1.0 - z)
                    )
                + sqrt(Q2factor)*(1.0 - z)
                    * ( t_0*(1.0 + z)*(17.0 - 6.0*z - 11.0*power_of<2>(z)
                                    + sqrt(tfactor)*(9.0 + 20.0*z + 11.0*power_of<2>(z)))
                        - 2*t_p*(3.0 + 22.0*z - 25.0*power_of<2>(z) + sqrt(tfactor)*(3.0 + 12.0*z + 25.0*power_of<2>(z)))
                        - t_m*power_of<2>(1.0 - z)*(sqrt(tfactor)*(3.0 + 11.0*z) + 11.0*(1.0 - z))
                    )
               )
            )
            / ( t_p*pow(1.0 - z, 3.5) * power_of<3>(1.0 - z + sqrt(tfactor)*(1.0 + z))
              * power_of<4>(sqrt(tfactor)*(1.0 + z) + sqrt(Q2factor)*(1.0 - z))
              * pow((tfactor*t_p*power_of<2>(1.0 + z)*(4*t_p*z - t_0*power_of<2>(1.0 + z) + t_m*power_of<2>(1.0 - z)))/power_of<4>(1.0 - z), 0.25)
              * sqrt(32*M_PI*chi_1m)
            )
        );
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::w_z(const complex<double> & z) const
    {
        return (1.0 + z) * pow((1.0 - z), 7.0 / 2.0);
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitilde_z(const complex<double> & z, const double & chi_0p) const
    {
        // The weight function ``(1.0 + z) * (1.0 - z)^(+7/2)`` has been cancelled against the outer function
        // to remove unphysical singularities and correct the asymptotic behaviour
        const double t_p      = this->_t_p();
        const double t_0      = this->_t_0();
        const double t_m      = this->_t_m();
        const double tfactor  = 1.0 - t_0 / t_p;
        const double Q2       = 1.0;
        const double Q2factor = 1.0 + Q2 / t_p;
        const complex<double> zfactor = (1.0 + z) / (1.0 - z);

        return pow(t_p, -3.0/4.0) * sqrt(t_m)
        / sqrt(32.0 * M_PI * chi_0p / 3.0)
        * pow(1.0-z, -4.0) * pow(1.0+z, -1.0/2.0) * sqrt(tfactor) * (1.0 + zfactor)
        * pow( power_of<2>(1.0+z) * tfactor / power_of<4>(1.0-z)
               * (4.0 * t_p * z + t_m * power_of<2>(1.0-z) - t_0 * power_of<2>(1.0+z)), 1.0/4.0)
        / power_of<2>(zfactor * sqrt(tfactor) + 1.0)
        / power_of<2>(zfactor * sqrt(tfactor) + sqrt(Q2factor));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::phitildeprime_z(const complex<double> & z, const double & chi_0p) const
    {
        // Derivative of the modified outer function phitilde_0 with respect to z
        const double t_p      = this->_t_p();
        const double t_0      = this->_t_0();
        const double t_m      = this->_t_m();
        const double tfactor  = 1.0 - t_0 / t_p;
        const double Q2       = 1.0;
        const double Q2factor = 1.0 + Q2 / t_p;
        const complex<double> zfactor = (1.0 + z) / (1.0 - z);


        return (pow(tfactor,1.5)*sqrt(t_m*t_p)*pow(zfactor,1.5)*
        (sqrt(tfactor)*(-(t_0*power_of<2>(1.0 + z)*(5.0 + 6.0*z - 11.0*power_of<2>(z) +
                  sqrt(tfactor)*(-3.0 + 8.0*z + 11.0*power_of<2>(z)))) +
             2*t_p*(1.0 + 7.0*z + 15.0*power_of<2>(z) - 23.0*power_of<3>(z) +
                sqrt(tfactor)*(1.0 - 7.0*z + 15.0*power_of<2>(z) + 23.0*power_of<3>(z))) +
             t_m*(3.0 + 8.0*z - 11.0*power_of<2>(z) + sqrt(tfactor)*(-5.0 + 6.0*z + 11.0*power_of<2>(z)))*power_of<2>(1.0 - z)) -
          sqrt(Q2factor)*(1.0 - z)*(t_0*(1.0 + z)*(13.0 - 2.0*z - 11.0*power_of<2>(z) +
                sqrt(tfactor)*(5.0 + 16.0*z + 11.0*power_of<2>(z))) -
             2*t_p*(1.0 + 22.0*z - 23.0*power_of<2>(z) + sqrt(tfactor)*(1.0 + 8.0*z + 23.0*power_of<2>(z))) -
             t_m*power_of<2>(1.0 - z)*(sqrt(tfactor)*(3.0 + 11.0*z) + 11.0*(1.0 - z)))))/
      (sqrt(t_p)*pow(1.0 - z,4.5)*pow(1.0 - z + sqrt(tfactor)*(1.0 + z),3)*
        pow(sqrt(tfactor)*(1.0 + z) + sqrt(Q2factor)*(1.0 - z),3)*
        pow((tfactor*t_p*power_of<2>(1.0 + z)*(4*t_p*z - t_0*power_of<2>(1.0 + z) + t_m*power_of<2>(1.0 - z)))/power_of<4>(1.0 - z),0.75)*
        sqrt(32*M_PI*chi_0p/3.));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::resonance_product_z(const complex<double> & z) const
    {
        // Product of the resonance factors for f0
        // TODO
        return 1.0;
    }

    double
    KSvD2025FormFactors<VacuumToKPi>::_b0_fp(const double & chi_1m) const
    {
        // determine the coefficient b^+_0 of f_+(q^2) by imposing that Im f_+(q^2) ~ sqrt(q^2 - t_+)^3
        const complex<double> zm1 = this->z(-1.0);

        const complex<double> phitilde_m1 = this->phitilde_p(zm1, chi_1m);
        const complex<double> phitildeprime_m1 = this->phitildeprime_p(zm1, chi_1m);

        const complex<double> resonance_productprime = this->resonance_productprime_p(zm1);

        const complex<double> X = phitilde_m1;
        const complex<double> Xprime = -1.0*phitildeprime_m1 / power_of<2>(phitilde_m1) + resonance_productprime / phitilde_m1;

        const double b0 = -1.0 / Xprime * ();

        return b0;
    }

    double
    KSvD2025FormFactors<VacuumToKPi>::_b0_fz(const double & chi_1m, const double & chi_0p) const
    {
        // determine the coefficient b^0_0 of f_0(q^2) by imposing that f_+(0) = f_0(0)
        const complex<double> z0 = this->z(0.0);

        std::array<double, 10> bp;
        bp[0] = this->_b0_fp(chi_1m);
        std::copy(_b_fp.cbegin(), _b_fp.cend(), bp.begin()+1);

        std::array<double, 10> bz;
        bz[0] = 0.0;
        std::copy(_b_fz.cbegin(), _b_fz.cend(), bz.begin()+1);

        const complex<double> bp_sum = this->series_m(z0, bp);
        const complex<double> b0_sum = this->series_m(z0, bz);

        const auto Pi_p = this->resonance_product_p(z0);
        const auto Pi_z = this->resonance_product_z(z0);

        const auto phitilde_p_z0 = this->phitilde_p(z0, chi_1m);
        const auto phitilde_z_z0 = this->phitilde_z(z0, chi_0p);


        const double b0 = std::real( (phitilde_z_z0 / phitilde_p_z0) * (Pi_p / Pi_z) * bp_sum - b0_sum );
        return b0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_p(const double & q2) const
    {
        static const double eps = 1.0e-12;
        return f_p(complex<double>(q2, eps));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_p(const complex<double> & q2) const
    {
        const auto z        = this->z(q2);
        const auto chi_1m   = 0.5; // TODO
        const auto phitilde = this->phitilde_p(z, chi_1m);

        const auto Pi_p = this->resonance_product_p(z);

        // prepare expansion coefficients
        std::array<double, 10> bp;
        std::copy(_b_fp.cbegin(), _b_fp.cend(), bp.begin()+1);
        // Fix b[0] to enforce f_+(q2=0) = 1
        bp[0] = _b0_fp(chi_1m);

        const auto series = this->series_m(z, bp);

        return series * Pi_p /  phitilde;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_0(const double & q2) const
    {
        static const double eps = 1.0e-12;
        return f_0(complex<double>(q2, eps));
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_0(const complex<double> & q2) const
    {
        const auto z        = this->z(q2);
        const auto chi_1m   = 0.5; // TODO
        const auto chi_0p   = 0.3; // TODO
        const auto phitilde = this->phitilde_z(z, chi_0p);

        const auto Pi_z = this->resonance_product_z(z);

        // prepare expansion coefficients
        std::array<double, 10> bz;
        std::copy(_b_fp.cbegin(), _b_fp.cend(), bz.begin()+1);
        // Fix b[0] to enforce f_0(0) = f_+(0)
        bz[0] = _b0_fz(chi_1m, chi_0p);

        const auto series = this->series_m(z, bz);

        return series * Pi_z /  phitilde;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    KSvD2025FormFactors<VacuumToKPi>::f_t(const complex<double> & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    double KSvD2025FormFactors<VacuumToKPi>::dispersive_integrand_p(const double & alpha) const
    {
        const complex<double> z = std::polar(1.0, alpha);
        const complex<double> w = this->w_p(z);

        const auto chi_1m = 1.0; // TODO

        const auto resonance_product = this->resonance_product_p(z);

        // prepare expansion coefficients
        std::array<double, 10> bp;
        std::copy(_b_fp.cbegin(), _b_fp.cend(), bp.begin()+1);
        // Fix b[0] from f_+'(q2=t+) = 0
        bp[0] = _b0_fp(chi_1m);

        const auto series = this->series_m(z, bp);

        return std::norm(w * resonance_product * series);
    }

    double KSvD2025FormFactors<VacuumToKPi>::saturation_p() const
    {
        std::function<double (const double &)> f = [this](const double & alpha) -> double { return this->dispersive_integrand_p(alpha); };
        return integrate<GSL::QAGS>(f, -M_PI, M_PI) / (2.0 * M_PI);
    }

    double KSvD2025FormFactors<VacuumToKPi>::dispersive_integrand_z(const double & alpha) const
    {
        const complex<double> z = std::polar(1.0, alpha);
        const complex<double> w = this->w_z(z);

        const auto chi_1m = 1.0; // TODO
        const auto chi_0p = 1.0; // TODO

        const auto resonance_product = this->resonance_product_z(z);

        // prepare expansion coefficients
        std::array<double, 10> bz;
        std::copy(_b_fz.cbegin(), _b_fz.cend(), bz.begin()+1);
        // Fix b[0] from f_0(0) = f_+(0)
        bz[0] = _b0_fz(chi_1m, chi_0p);

        const auto series = this->series_m(z, bz);

        return std::norm(w * resonance_product * series);
    }

    double KSvD2025FormFactors<VacuumToKPi>::saturation_z() const
    {
        std::function<double (const double &)> f = [this](const double & alpha) -> double { return this->dispersive_integrand_z(alpha); };
        return integrate<GSL::QAGS>(f, -M_PI, M_PI) / (2.0 * M_PI);
    }

    const std::set<ReferenceName>
    KSvD2025FormFactors<VacuumToKPi>::references
    {
    };

    const std::vector<OptionSpecification>
    KSvD2025FormFactors<VacuumToKPi>::option_specifications
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KSvD2025FormFactors<VacuumToKPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KSvD2025FormFactors<VacuumToKPi>::end_options()
    {
        return option_specifications.cend();
    }

    template class KSvD2025FormFactors<VacuumToKPi>;
}
