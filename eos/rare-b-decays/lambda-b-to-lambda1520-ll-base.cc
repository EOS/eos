/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Méril Reboud
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

#include <eos/rare-b-decays/lambda-b-to-lambda1520-ll-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    LambdaBToLambda1520Dilepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda(1520)::" + o.get("form-factors", "ABR2022"), p)),
        opt_l(o, "l", { "e", "mu", "tau" }, "mu"),
        mu(p["sb" + opt_l.value() + opt_l.value() + "::mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["WET::G_Fermi"], *this),
        hbar(p["QM::hbar"], *this),
        m_l(p["mass::" + opt_l.value()], *this),
        m_Lb(p["mass::Lambda_b"], *this),
        m_Lstar(p["mass::Lambda(1520)"], *this),
        cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
        lepton_flavor(opt_l.value())
    {
        if (! form_factors.get())
            throw InternalError("Form factors not found!");

        this->uses(*form_factors);
        this->uses(*model);
    }

    LambdaBToLambda1520Dilepton::AmplitudeGenerator::~AmplitudeGenerator()
    {
    }

    double
    LambdaBToLambda1520Dilepton::AmplitudeGenerator::lambda(const double & s) const
    {
        return eos::lambda(m_Lb() * m_Lb(), m_Lstar() * m_Lstar(), s);
    }

    double
    LambdaBToLambda1520Dilepton::AmplitudeGenerator::beta_l(const double & s) const
    {
        return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
    }
}
