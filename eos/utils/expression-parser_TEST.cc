/*
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

#include <eos/utils/expression.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/expression-cloner.hh>
#include <eos/utils/expression-evaluator.hh>
#include <eos/utils/expression-kinematic-reader.hh>
#include <eos/utils/expression-maker.hh>
#include <eos/utils/expression-parser-impl.hh>
#include <eos/utils/expression-printer.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

#include <iostream>

using namespace test;
using namespace eos::exp;
using namespace eos;

class ExpressionTest
{
    private:
        const std::string _input;

    public:
        Expression e;
        bool completed;

        using It = std::string::const_iterator;
        ExpressionParser<It> parser;

        ExpressionTest(const std::string input) :
            _input(input)
        {
            It first(_input.begin()), last(_input.end());
            completed = qi::phrase_parse(first, last, parser, ascii::space, e) && (first == last);
        }
};


class ExpressionParserTest :
    public TestCase
{
    public:
        ExpressionParserTest() :
        TestCase("expression_parser_test")
        {
        }

        virtual void run() const
        {
            // testing parser failure
            {
                ExpressionTest test("1 /* 2");

                TEST_CHECK(! test.completed);
            }

            // testing basic parsing of constants and binary expressions
            {
                ExpressionTest test("1+2*3");

                std::stringstream out;
                ExpressionPrinter printer(out);
                test.e.accept(printer);

                ExpressionEvaluator evaluator;

                TEST_CHECK(test.completed);
                TEST_CHECK_EQUAL(test.e.accept_returning<double>(evaluator), 7.0);
                TEST_CHECK_EQUAL_STR("BinaryExpression(ConstantExpression(1) + BinaryExpression(ConstantExpression(2) * ConstantExpression(3)))", out.str());

                // Simple exponentiation
                ExpressionTest test2("1+2^2*3");
                TEST_CHECK(test2.completed);
                TEST_CHECK_EQUAL(test2.e.accept_returning<double>(evaluator), 13.0);

                // Non interger exponentiation
                ExpressionTest test3("2^(1+3.5)+3");
                TEST_CHECK(test3.completed);
                TEST_CHECK_RELATIVE_ERROR(test3.e.accept_returning<double>(evaluator), 25.627416998, 1e-5);
            }

            // testing parsing of an expression containing kinematic variables
            {
                ExpressionTest test("{q2_mu} - {q2_e}");

                std::stringstream out;
                ExpressionPrinter printer(out);
                test.e.accept(printer);

                ExpressionEvaluator evaluator;

                TEST_CHECK(test.completed);
                TEST_CHECK_EQUAL_STR(
                    "BinaryExpression(KinematicVariableNameExpression(q2_mu)"
                    " - "
                    "KinematicVariableNameExpression(q2_e))",
                    out.str()
                );

                // Cannot evaluate expression with KinematicVariableNameExpression objects
                TEST_CHECK_THROWS(InternalError, test.e.accept_returning<double>(evaluator));

                // Extract kinematic variables from an expression
                ExpressionKinematicReader kinematic_reader;
                std::set<std::string> kinematic_set = test.e.accept_returning<std::set<std::string>>(kinematic_reader);

                std::set<std::string> expected_kinematic{"q2_mu", "q2_e"};
                TEST_CHECK_EQUAL(expected_kinematic, kinematic_set);

                // Make and evaluate expression
                Kinematics k = Kinematics({{"q2_mu", 4.0},  {"q2_e", 3.0}});
                ExpressionMaker maker(Parameters::Defaults(), k, Options());
                Expression assessable_test = test.e.accept_returning<Expression>(maker);

                std::stringstream out2;
                ExpressionPrinter printer2(out2);
                assessable_test.accept(printer2);

                TEST_CHECK_EQUAL_STR(
                    "BinaryExpression(KinematicVariableExpression(q2_mu)"
                    " - "
                    "KinematicVariableExpression(q2_e))",
                    out2.str()
                );

                TEST_CHECK_EQUAL(assessable_test.accept_returning<double>(evaluator), 1.0);
            }

            // testing parsing and evaluation of observables
            {
                // test::obs1 is a test observable that requires two kinematic specifications, q2_min and q2_max
                // it returns p[mass::c] * multiplier * (q2_max - q2_min)
                ExpressionTest test("<<test::obs1;multiplier=2>>[q2_min=>q2_min_num] / <<test::obs1>>[q2_min=0.0]");

                std::stringstream out;
                ExpressionPrinter printer(out);
                test.e.accept(printer);

                ExpressionEvaluator evaluator;

                TEST_CHECK(test.completed);
                TEST_CHECK_EQUAL_STR(
                    "BinaryExpression(ObservableNameExpression(test::obs1;multiplier=2, aliases=[q2_min=>q2_min_num])"
                    " / "
                    "ObservableNameExpression(test::obs1, values=[q2_min=0]))",
                    out.str()
                );

                // Cannot evaluate expression with ObservableNameExpression objects
                TEST_CHECK_THROWS(InternalError, test.e.accept_returning<double>(evaluator));

                // Extract kinematic variables from an expression
                ExpressionKinematicReader kinematic_reader;
                std::set<std::string> kinematic_set = test.e.accept_returning<std::set<std::string>>(kinematic_reader);

                std::set<std::string> expected_kinematic{"q2_max", "q2_min_num"};
                TEST_CHECK_EQUAL(expected_kinematic, kinematic_set);
            }

            // Test numerical evaluation
            {
                ExpressionTest test("<<test::obs1;multiplier=2>>[q2_min=>q2_min_num] / <<test::obs1>>[q2_min=>q2_min_denom]");

                Parameters p = Parameters::Defaults();
                p["mass::c"] = 1.2;

                Kinematics k       = Kinematics({{"q2_min_num", 4.0}, {"q2_min_denom", 3.0}, {"q2_max", 10.0}});
                Kinematics k_num   = Kinematics({{"q2_min", 4.0},  {"q2_max", 10.0}});
                Kinematics k_denom = Kinematics({{"q2_min", 3.0},  {"q2_max", 10.0}});

                auto obs_num   = Observable::make("test::obs1;multiplier=2", p, k_num,   Options());
                auto obs_denom = Observable::make("test::obs1",              p, k_denom, Options());

                TEST_CHECK_RELATIVE_ERROR(obs_num->evaluate(),   14.4,  1e-10);
                TEST_CHECK_RELATIVE_ERROR(obs_denom->evaluate(),  8.4,  1e-10);

                // Make and evaluate expression
                ExpressionMaker maker(p, k, Options());
                Expression assessable_test = test.e.accept_returning<Expression>(maker);

                ExpressionEvaluator evaluator;

                TEST_CHECK_RELATIVE_ERROR(
                    assessable_test.accept_returning<double>(evaluator),
                    obs_num->evaluate() / obs_denom->evaluate(),
                    1e-3);


                // Observable with exponentiation
                ExpressionTest test2("<<mass::tau>>^2 - <<mass::mu>>^2");
                Expression assessable_test2 = test2.e.accept_returning<Expression>(maker);

                TEST_CHECK(test2.completed);
                TEST_CHECK_RELATIVE_ERROR(assessable_test2.accept_returning<double>(evaluator), 3.14592, 1e-3);
            }

            // testing cloning and usage of parameters
            {
                ExpressionTest test("{q2} - 4 * <<mass::mu>> * <<mass::tau>>");

                std::stringstream out;
                ExpressionPrinter printer(out);
                test.e.accept(printer);

                TEST_CHECK(test.completed);
                TEST_CHECK_EQUAL_STR(
                    "BinaryExpression("
                        "KinematicVariableNameExpression(q2) - "
                        "BinaryExpression("
                            "BinaryExpression(ConstantExpression(4) * ObservableNameExpression(mass::mu))"
                            " * ObservableNameExpression(mass::tau)"
                        ")"
                    ")",
                    out.str()
                );

                Options o;
                Kinematics k = Kinematics({{"q2", 10.0}});
                Parameters p = Parameters::Defaults();
                p["mass::mu"]     =  1.0;
                p["mass::tau"]    =  2.0;

                ExpressionMaker maker(p, k, o);
                Expression assessable_test = test.e.accept_returning<Expression>(maker);

                Kinematics k2 = Kinematics({{"q2", 20.0}});
                Parameters p2 = Parameters::Defaults();
                p2["mass::mu"]     =  2.0;
                p2["mass::tau"]    =  2.0;

                ExpressionCloner cloner(p2, k2, o);
                Expression cloned_test = assessable_test.accept_returning<Expression>(cloner);

                ExpressionEvaluator evaluator;

                TEST_CHECK_EQUAL(assessable_test.accept_returning<double>(evaluator),    2.0);
                TEST_CHECK_EQUAL(cloned_test.accept_returning<double>(evaluator),        4.0);

                // Test that parameters are considered as used
                ParameterUser parameter_user;
                ExpressionMaker maker_user(p, k, o, &parameter_user);
                assessable_test = test.e.accept_returning<Expression>(maker_user);
                std::set<Parameter::Id> used_ids(parameter_user.begin(), parameter_user.end());
                std::set<Parameter::Id> expected_ids{p["mass::mu"].id(), p["mass::tau"].id()};

                TEST_CHECK_EQUAL(used_ids, expected_ids);
            }

            // testing mixing of kinematic variables and observables
            {
                // { } are not allowed in the prefix of QualifiedNames
                TEST_CHECK_THROWS(QualifiedNameSyntaxError, ExpressionTest test("<<{test::obs1}>>"));

                // { } are allowed in the prefix of QualifiedNames (but test::obs1{} is not an existing obervable)
                ExpressionTest test("<<test::obs1{}>>");
                std::stringstream out;
                ExpressionPrinter printer(out);
                test.e.accept(printer);

                TEST_CHECK_EQUAL_STR("ObservableNameExpression(test::obs1{})", out.str());

                // Names of kinematic variables are not restricted
                ExpressionTest test2("{<<test::obs1>>}");
                TEST_CHECK(test2.completed);

                std::stringstream out2;
                ExpressionPrinter printer2(out2);
                test2.e.accept(printer2);
                TEST_CHECK_EQUAL_STR("KinematicVariableNameExpression(<<test::obs1>>)", out2.str());
            }

            // testing the compatibility of kinematic variables and aliases
            {
                // Simple case, no conflict
                ExpressionTest test("<<test::obs1>>[q2_min=>q2_min_num] * {q2_min_num}");
                ExpressionKinematicReader kinematic_reader;
                std::set<std::string> kinematic_set = test.e.accept_returning<std::set<std::string>>(kinematic_reader);
                std::set<std::string> expected_kinematic{"q2_max", "q2_min_num"};

                TEST_CHECK_EQUAL(expected_kinematic, kinematic_set);

                // Problematic case, conflict between the alias and the kinematic variable
                ExpressionTest test2("<<test::obs1>>[q2_min=>q2_min_num] * {q2_min}");
                std::set<std::string> kinematic_set2 = test2.e.accept_returning<std::set<std::string>>(kinematic_reader);
                std::set<std::string> expected_kinematic2{"q2_max", "q2_min", "q2_min_num"};

                TEST_CHECK_EQUAL(expected_kinematic2, kinematic_set2);

                // Make and evaluate expression
                Parameters p = Parameters::Defaults();
                Kinematics k = Kinematics({{"q2_min", 3.0}, {"q2_max", 10.0}});
                double obs_value = Observable::make("test::obs1", p, k, Options())->evaluate();

                k = Kinematics({{"q2_min", 4.0}, {"q2_min_num", 3.0}, {"q2_max", 10.0}});
                ExpressionMaker maker(p, k, Options());
                Expression assessable_test = test2.e.accept_returning<Expression>(maker);
                ExpressionEvaluator evaluator;

                // TEST_CHECK_RELATIVE_ERROR(assessable_test.accept_returning<double>(evaluator), 4.0 * obs_value, 1e-3);
            }
        }
} expression_parser_test;
