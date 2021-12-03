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

#include <test/test.hh>
#include <eos/utils/szego-polynomial.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class SzegoPolynomialTest :
    public TestCase
{
    public:
        SzegoPolynomialTest() :
            TestCase("szego_polynomial_test")
        {
        }

        virtual void run() const
        {
			// test case
			{
				SzegoPolynomial<5u> p{
					2.47895, // norm of the measure
					{ 0.762914, -0.7988, 0.807686, -0.81062, 0.811894 }
				};

				{
					const auto [p0, p1, p2, p3, p4, p5] = p(-0.1);

					TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p1, -0.847745,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p2, +1.54489,            1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p3, -2.82388,            1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p4, +5.166081,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p5, -9.464102,           1.0e-5);
				}

				{
					const auto [p0, p1, p2, p3, p4, p5] = p(+0.0);

					TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p1, -0.749503,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p2, +1.304458,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p3, -2.237009,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p4, +3.834085,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p5, -6.577731,           1.0e-5);
				}

				{
					const auto [p0, p1, p2, p3, p4, p5] = p(+0.1);

					TEST_CHECK_RELATIVE_ERROR(p0, +0.6351351032391984, 1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p1, -0.651261,           1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p2, +1.09668,            1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p3, -1.76188,            1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p4, +2.82970,            1.0e-5);
					TEST_CHECK_RELATIVE_ERROR(p5, -4.54691,            1.0e-5);
				}
			}
        }
} szego_polynomial_test;