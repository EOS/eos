/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#include <eos/b-decays/b-to-vec-nu-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/quantum-numbers.hh>

#include <map>
#include <string>

namespace eos
{
    template <>
    struct Implementation<BToVectorDineutrino>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        QuarkFlavorOption opt_D;
        QuarkFlavorOption opt_q;
        SpecifiedOption opt_I;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_V;

        UsedParameter g_fermi;

        UsedParameter alpha_e;

        UsedParameter hbar;

        const double isospin_factor;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> ()> lambda_t;
        std::function<WilsonCoefficients<wc::SBNuNu> ()> wc;

        GSL::QAGS::Config int_config;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        // { D, q, I } -> { process, m_B, m_V, c_I }
        // D: s, d; the quark flavor in the weak transition
        // q: u, d, s: the spectar quark flavor
        // I: 1, 0, 1/2: the total isospin of the daughter meson
        // process: string that can be used to obtain the form factor
        // B: name of the B meson
        // V: name of the daughter meson
        // c_I: isospin factor by which the amplitudes are multiplied
        static const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, double>> process_map;

        inline std::string _process() const
        {
            const auto D = opt_D.value();
            const auto q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(D, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of D=" + opt_D.str() + ", q=" + opt_q.str() + ", I=" + I);

            return std::get<0>(p->second);
        }

        inline std::string _B() const
        {
            const auto D = opt_D.value();
            const auto q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(D, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of D=" + opt_D.str() + ", q=" + opt_q.str() + ", I=" + I);

            return std::get<1>(p->second);
        }

        inline std::string _V() const
        {
            const auto D = opt_D.value();
            const auto q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(D, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of D=" + opt_D.str() + ", q=" + opt_q.str() + ", I=" + I);

            return std::get<2>(p->second);
        }

        inline double _isospin_factor() const
        {
            const auto D = opt_D.value();
            const auto q = opt_q.value();
            const std::string I = opt_I.value();
            const auto p = process_map.find(std::make_tuple(D, q, I));

            if (p == process_map.end())
                throw InternalError("Unsupported combination of D=" + opt_D.str() + ", q=" + opt_q.str() + ", I=" + I);

            return std::get<3>(p->second);
        }


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_D(o, options, "D"),
            opt_q(o, options, "q"),
            opt_I(o, options, "I"),
            m_B(p["mass::" + _B()], u),
            tau_B(p["life_time::" + _B()], u),
            m_V(p["mass::" + _V()], u),
            g_fermi(p["WET::G_Fermi"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            hbar(p["QM::hbar"], u),
            isospin_factor(_isospin_factor()),
            mu(p[opt_D.str() + "b" + "nunu::mu"], u),
            int_config(GSL::QAGS::Config().epsrel(0.5e-3)),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            form_factors(FormFactorFactory<PToV>::create(_process() + "::" + o.get("form-factors", "BSZ2015"), p, o))
        {
            switch (opt_D.value())
            {
                case QuarkFlavor::strange:
                    lambda_t  = [*this] () { return model->ckm_tb() * std::conj(model->ckm_ts()); };
                    wc        = [*this] () { return model->wet_sbnunu(cp_conjugate); };
                    break;

                default:
                    throw InternalError("Unexpected quark flavor: '" + opt_D.str() + "'");
            }

            u.uses(*form_factors);
            u.uses(*model);
        }

        // differential decay width
        double differential_decay_width(const double & q2) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_V = this->m_V(), m_V2 = m_V * m_V;
            const double lambda = eos::lambda(m_B2, m_V2, q2), sqrt_lambda = std::sqrt(lambda);

            const auto wc = this->wc();

            const double V   = form_factors->v(q2);
            const double A1  = form_factors->a_1(q2);
            const double A12 = form_factors->a_12(q2);

            // using different normalization than [FLS:2021A], eq. (1)
            // note that eq. (1) is a Lagrangian, while we use the Hamiltonian definition
            const double norm = power_of<2>(4.0 * g_fermi() * alpha_e / (2.0 * M_PI)) / 2.0 * std::norm(lambda_t())
                // remainder as in [FLS:2021A], eq. (13), except for moving the q2 factor into the square brackets
                * sqrt_lambda / (power_of<3>(4.0 * M_PI * m_B));

            // first term in square brackets in [FLS:2021A], eq. (13)
            const double contr_vector = q2 * lambda / (12.0 * (m_B + m_V) * (m_B + m_V)) * V * V * std::norm(wc.cVL() + wc.cVR());
            const double contr_axial  = 8.0 * m_B2 * m_V2 / 3.0 * A12 * A12                      * std::norm(wc.cVL() - wc.cVR())
                                      + q2 * (m_B + m_V) * (m_B + m_V) / (12.0) * A1 * A1        * std::norm(wc.cVL() - wc.cVR());

            // assume the production of 3 diagonal neutrino flavors (nu_i nubar_i)
            return 3.0 * norm * (contr_vector + contr_axial);
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
        }
    };

    const std::map<std::tuple<QuarkFlavor, QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, double>>
    Implementation<BToVectorDineutrino>::Implementation::process_map
    {
        { { QuarkFlavor::strange, QuarkFlavor::up,   "1/2" }, { "B->K^*",     "B_u", "K_u^*",  1.0                  } },
        { { QuarkFlavor::strange, QuarkFlavor::down, "1/2" }, { "B->K^*",     "B_d", "K_d^*",  1.0                  } }
    };

    const std::vector<OptionSpecification>
    Implementation<BToVectorDineutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        { "D", { "s" },             "s" },
        { "q", { "u", "d" },        "u" },
        { "I", { "1", "0", "1/2" }, ""  },
    };

    BToVectorDineutrino::BToVectorDineutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToVectorDineutrino>(new Implementation<BToVectorDineutrino>(parameters, options, *this))
    {
    }

    BToVectorDineutrino::~BToVectorDineutrino()
    {
    }

    double
    BToVectorDineutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_branching_ratio(q2);
    }

    double
    BToVectorDineutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToVectorDineutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, q2_min, q2_max, _imp->int_config);
    }

    const std::string
    BToVectorDineutrino::description = "\
    The decay B->P nu nu, where both B=(B qbar) and P=(U qbar) are pseudoscalars, and l=e,mu,tau is a lepton.";

    const std::string
    BToVectorDineutrino::kinematics_description_q2 = "\
    The invariant mass of the nu-nubar pair in GeV^2.";

    const std::set<ReferenceName>
    BToVectorDineutrino::references
    {
        "BGNS:2014A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToVectorDineutrino::begin_options()
    {
        return Implementation<BToVectorDineutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToVectorDineutrino::end_options()
    {
        return Implementation<BToVectorDineutrino>::options.cend();
    }
}
