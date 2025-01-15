/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2023-2024 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_USED_PARAMETER_READER_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_USED_PARAMETER_READER_HH 1

#include <eos/observable.hh>
#include <eos/utils/expression-fwd.hh>

namespace eos::exp
{
    // Visit the expression tree and return the set of used parameters.
    class ExpressionUsedParameterReader
    {
        public:
            std::set<Parameter::Id> parameter_ids;

            ExpressionUsedParameterReader() = default;
            ~ExpressionUsedParameterReader() = default;

            void visit(const BinaryExpression & e);

            void visit(const FunctionExpression & e);

            void visit(const ConstantExpression &);

            void visit(const ObservableNameExpression & e);

            void visit(const ObservableExpression & e);

            void visit(const ParameterNameExpression &);

            void visit(const ParameterExpression &);

            void visit(const KinematicVariableNameExpression & e);

            void visit(const KinematicVariableExpression & e);

            void visit(const CachedObservableExpression & e);
    };
}

#endif
