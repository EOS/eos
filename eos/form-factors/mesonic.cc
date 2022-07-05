/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014, 2015, 2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2019 Nico Gubernari
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

#include <eos/form-factors/analytic-b-to-pi.hh>
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/mesonic-hqet.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/qualified-name.hh>

#include <cmath>
#include <limits>
#include <map>

namespace eos
{
    /* P -> V Processes */

    /* B_{u,d} -> omega */
    // [BSZ2015]
    template class BSZ2015FormFactors<BToOmega, PToV>;

    /* B_{u,d} -> K^* */

    constexpr double BToDstar::mR2_0m;
    constexpr double BToDstar::mR2_1m;
    constexpr double BToDstar::mR2_1p;

    constexpr double BToKstar::mR2_0m;
    constexpr double BToKstar::mR2_1m;
    constexpr double BToKstar::mR2_1p;

    constexpr double BToRho::mR2_0m;
    constexpr double BToRho::mR2_1m;
    constexpr double BToRho::mR2_1p;

    constexpr double BsToPhi::mR2_0m;
    constexpr double BsToPhi::mR2_1m;
    constexpr double BsToPhi::mR2_1p;

    constexpr double BsToKstar::mR2_0m;
    constexpr double BsToKstar::mR2_1m;
    constexpr double BsToKstar::mR2_1p;

    constexpr double BsToDsstar::mR2_0m;
    constexpr double BsToDsstar::mR2_1m;
    constexpr double BsToDsstar::mR2_1p;


    // mass B_d, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_B       = 5.2795;
    // mass K^*0, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToV>::_m_Kstar   = 0.89594;
    // mass B_s (0-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_0m  = 5.366 * 5.366;
    // mass B_s (1-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1m  = 5.412 * 5.412;
    // mass B_s (1+), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToV>::_m_Bs2_1p  = 5.829 * 5.829;
    const double KMPW2010FormFactors<PToV>::_tau_p     = (_m_B + _m_Kstar) * (_m_B + _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_m     = (_m_B - _m_Kstar) * (_m_B - _m_Kstar);
    const double KMPW2010FormFactors<PToV>::_tau_0     = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToV>;


    FormFactors<PToV>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<PToV>::KeyType, FormFactorFactory<PToV>::ValueType>
    FormFactorFactory<PToV>::form_factors
    {
        { "B->omega::BSZ2015",    &BSZ2015FormFactors<BToOmega,   PToV>::make         },
        { "B->rho::BSZ2015",      &BSZ2015FormFactors<BToRho,     PToV>::make         },
        { "B->K^*::KMPW2010",     &KMPW2010FormFactors<PToV>::make                    },
        { "B->K^*::BSZ2015",      &BSZ2015FormFactors<BToKstar,   PToV>::make         },
        { "B->D^*::BSZ2015",      &BSZ2015FormFactors<BToDstar,   PToV>::make         },
        { "B->D^*::HQET",         &HQETFormFactors<BToDstar,      PToV>::make         },
        { "B->D^*::BGL1997",      &BGL1997FormFactors<BToDstar>::make                 },
        { "B_s->K^*::BSZ2015",    &BSZ2015FormFactors<BsToKstar,  PToV>::make         },
        { "B_s->D_s^*::BSZ2015",  &BSZ2015FormFactors<BsToDsstar, PToV>::make         },
        { "B_s->D_s^*::HQET",     &HQETFormFactors<BsToDsstar,    PToV>::make         },
        { "B_s->phi::BSZ2015",    &BSZ2015FormFactors<BsToPhi,    PToV>::make         },
        // analytic computations
        { "B->K^*::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::make   },
        { "B->D^*::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::make   },
        { "B->rho::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::make     },
        { "B_s->K^*::B-LCSR",     &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::make  },
        { "B_s->phi::B-LCSR",     &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::make    },
        { "B_s->D_s^*::B-LCSR",   &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::make }
    };

    complex<double>
    FormFactors<PToV>::v(const complex<double> &) const
    {
        throw InternalError("P->V form factor V for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_0(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_0 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_1(const complex<double> &) const
    {
        throw InternalError("P->V form factors A_1 complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_12(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_12 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_2(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_2 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_1(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_1 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_2(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_2 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_23(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_23 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    std::shared_ptr<FormFactors<PToV>>
    FormFactorFactory<PToV>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<PToV>> result;

        auto i = form_factors.find(name);
        if (form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToV>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToV>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToV>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : FormFactorFactory<PToV>::form_factors)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result { "form-factors", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

    /* P -> P Processes */

    /* B_{u,d} -> K */

    // [BCL2008]
    template class BCL2008FormFactors<BToK, 3u>;

    // mass B_u, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_B      = 5.27917;
    // mass K^+, cf. [PDG 2010]
    const double KMPW2010FormFactors<PToP>::_m_K      = 0.493677;
    // mass B_s^* (1-), cf. [KMPW2010]
    const double KMPW2010FormFactors<PToP>::_m_Bs2    = 5.412 * 5.412;
    const double KMPW2010FormFactors<PToP>::_tau_p    = (_m_B + _m_K) * (_m_B + _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_m    = (_m_B - _m_K) * (_m_B - _m_K);
    const double KMPW2010FormFactors<PToP>::_tau_0    = _tau_p - std::sqrt(_tau_p * _tau_p - _tau_m * _tau_p);
    template class KMPW2010FormFactors<PToP>;


    /* B_{u,d} -> pi */

    // [BCL2008]
    template class BCL2008FormFactors<BToPi, 3u>;
    template class BCL2008FormFactors<BToPi, 4u>;
    template class BCL2008FormFactors<BToPi, 5u>;

    /* B_{u,d} -> D */

    // [BCL2008]
    template class BCL2008FormFactors<BToD, 3u>;


    FormFactors<PToP>::~FormFactors() = default;

    double FormFactors<PToP>::f_m(const double & /*s*/) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double FormFactors<PToP>::f_p_d1(const double & s) const
    {
        using namespace std::placeholders;

        std::function<double (const double &)> f = [&, this](const double & q2) -> double { return this->f_p(q2); };

        return derivative<1u, deriv::TwoSided>(f, s);
    }

    double FormFactors<PToP>::f_p_d2(const double & s) const
    {
        using namespace std::placeholders;

        std::function<double (const double &)> f = [&, this](const double & q2) -> double { return this->f_p(q2); };

        return derivative<2u, deriv::TwoSided>(f, s);
    }

    const std::map<FormFactorFactory<PToP>::KeyType, FormFactorFactory<PToP>::ValueType>
    FormFactorFactory<PToP>::form_factors
    {
        // parametrizations
        // b -> s
        { "B->K::BCL2008",       &BCL2008FormFactors<BToK, 3u>::make              },
        { "B->K::KMPW2010",      &KMPW2010FormFactors<PToP>::make                 },
        { "B->K::BSZ2015",       &BSZ2015FormFactors<BToK,   PToP>::make          },
        // b -> u
        { "B->pi::BCL2008",      &BCL2008FormFactors<BToPi, 3u>::make             },
        { "B->pi::BCL2008-4",    &BCL2008FormFactors<BToPi, 4u>::make             },
        { "B->pi::BCL2008-5",    &BCL2008FormFactors<BToPi, 5u>::make             },
        { "B->pi::BSZ2015",      &BSZ2015FormFactors<BToPi,  PToP>::make          },
        { "B_s->K::BSZ2015",     &BSZ2015FormFactors<BsToK,  PToP>::make          },
        // b -> c
        { "B->D::BCL2008",       &BCL2008FormFactors<BToD, 3u>::make              },
        { "B->D::BSZ2015",       &BSZ2015FormFactors<BToD,   PToP>::make          },
        { "B->D::HQET",          &HQETFormFactors<BToD,      PToP>::make          },
        { "B->D::BGL1997",       &BGL1997FormFactors<BToD>::make                  },
        { "B_s->D_s::BSZ2015",   &BSZ2015FormFactors<BsToDs, PToP>::make          },
        { "B_s->D_s::HQET",      &HQETFormFactors<BsToDs,    PToP>::make          },
        // analytic computations
        { "B->pi::DKMMO2008",    &AnalyticFormFactorBToPiDKMMO2008::make          },
        { "B->pi::B-LCSR",       &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::make   },
        { "B->K::B-LCSR",        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::make    },
        { "B->D::B-LCSR",        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::make    },
        { "B_s->K::B-LCSR",      &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::make   },
        { "B_s->D_s::B-LCSR",    &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::make  }
    };

    complex<double>
    FormFactors<PToP>::f_p(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_+ for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToP>::f_0(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_0 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToP>::f_t(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_t for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    std::shared_ptr<FormFactors<PToP>>
    FormFactorFactory<PToP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<PToP>> result;

        auto i = FormFactorFactory<PToP>::form_factors.find(name);
        if (FormFactorFactory<PToP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToP>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : FormFactorFactory<PToP>::form_factors)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result { "form-factors", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

    /* P -> PP Processes */

    constexpr double BToPi::m2_Br1m;
    constexpr double BToPi::m2_Br0p;

    constexpr double BToK::m2_Br1m;
    constexpr double BToK::m2_Br0p;

    constexpr double BToD::m2_Br1m;
    constexpr double BToD::m2_Br0p;

    constexpr double BsToK::m2_Br1m;
    constexpr double BsToK::m2_Br0p;

    constexpr double BsToDs::m2_Br1m;
    constexpr double BsToDs::m2_Br0p;

    FormFactors<PToPP>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<PToPP>::KeyType, FormFactorFactory<PToPP>::ValueType>
    FormFactorFactory<PToPP>::form_factors
    {
        // analytic computations
        { "B->pipi::BFvD2016",            &AnalyticFormFactorBToPiPiBFvD2016::make   },
        { "B->pipi::FvDV2018-Dispersive", &AnalyticFormFactorBToPiPiFvDV2018::make   },
        { "B->pipi::FvDV2018",            &FvDV2018FormFactors<BToPiPi>::make        },
    };

    std::shared_ptr<FormFactors<PToPP>>
    FormFactorFactory<PToPP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<PToPP>> result;

        auto i = FormFactorFactory<PToPP>::form_factors.find(name);
        if (FormFactorFactory<PToPP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToPP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToPP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* V -> P Processes */

    FormFactors<VToP>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<VToP>::KeyType, FormFactorFactory<VToP>::ValueType>
    FormFactorFactory<VToP>::form_factors
    {
        // parametrizations
        // b -> c
        { "B^*->D::HQET",        &HQETFormFactors<BstarToD, VToP>::make          },
    };

    std::shared_ptr<FormFactors<VToP>>
    FormFactorFactory<VToP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<VToP>> result;

        auto i = FormFactorFactory<VToP>::form_factors.find(name);
        if (FormFactorFactory<VToP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<VToP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<VToP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return { result };
    }

    /* V -> V Processes */

    FormFactors<VToV>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<VToV>::KeyType, FormFactorFactory<VToV>::ValueType>
    FormFactorFactory<VToV>::form_factors
    {
        // parametrizations
        // b -> c
        // not yet supported
        { "B^*->D^*::HQET",      &HQETFormFactors<BstarToDstar, VToV>::make      },
    };

    std::shared_ptr<FormFactors<VToV>>
    FormFactorFactory<VToV>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        std::shared_ptr<FormFactors<VToV>> result;

        auto i = FormFactorFactory<VToV>::form_factors.find(name);
        if (FormFactorFactory<VToV>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<VToV>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<VToV>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }
}
