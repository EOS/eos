/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
 * Copyright (c) 2024 Méril Reboud
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

#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-nu-nu.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    template <>
    struct Implementation<BToDineutrino>
    {
        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        UsedParameter f_B;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter hbar;

        UsedParameter m_b;

        UsedParameter m_q;

        static const std::vector<OptionSpecification> options;

        std::function<complex<double> ()> lambda;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, options, "q"),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            mu(p[opt_q.str() + "bnunu" + "::mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            m_b(p["mass::b(MSbar)"], u),
            m_q(p["mass::" + opt_q.str() + "(2GeV)"], u)
        {
            Context ctx("When constructing B_q->nunu observables");

            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    lambda = [this]() { return this->model->ckm_tb() * conj(this->model->ckm_ts()); };
                    break;
                default:
                    // only neutral B mesons can decay in this channel
                    throw InternalError("ExclusiveBToDineutrino: q = '" + opt_q.str() + "' is not a valid option for a neutral decay channel");
            }

            u.uses(*model);
        }

        // cf. [BEKU2002], Eq. (3.6)
        double branching_ratio() const
        {
            double lambda_t = abs(lambda());

            WilsonCoefficients<wc::SBNuNu> wc = model->wet_sbnunu(false);

            return 0.0;
            //return power_of<2>(g_fermi() * alpha_e() * lambda_t * f_B()) / 64.0 / power_of<3>(M_PI) * tau_B / hbar
            //    * beta_l * power_of<3>(m_B()) * (
            //            power_of<2>(beta_l) * std::norm(m_B / (m_b + m_q) * (wc.cS() - wc.cSprime()))
            //            + std::norm(m_B / (m_b + m_q) * (wc.cP() - wc.cPprime()) + 2.0 * m_l / m_B * (wc.c10() - wc.c10prime())));
        }
    };

    BToDineutrino::BToDineutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDineutrino>(new Implementation<BToDineutrino>(parameters, options, *this))
    {
    }

    BToDineutrino::~BToDineutrino()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BToDineutrino>::options
    {
        Model::option_specification(),
        { "q", { "s" }, "s"}
    };

    double
    BToDineutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    const std::set<ReferenceName>
    BToDineutrino::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BToDineutrino::begin_options()
    {
        return Implementation<BToDineutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDineutrino::end_options()
    {
        return Implementation<BToDineutrino>::options.cend();
    }
}
