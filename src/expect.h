/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_EXPECT_H
#define SWIFT_EXPECT_H

#if defined(__GNUC__)

/**
 * @brief Tag a branch of an if-statement to indicate that
 * this is the more likely one.
 *
 * It should be read as "I expect 'x' to be true".
 *
 * Note that this turns into a no-op but gives information to the compiler.
 *
 * Note: This should only be used for branches where the
 * compiler suggested this annotation when running with
 * the runtime profiler (-fprofile-generate with GCC).
 *
 * @brief The expression in the test.
 */
#define swift_likely_branch(x) __builtin_expect(!!(x), 1)

/**
 * @brief Tag a branch of an if-statement to indicate that
 * this is the more unlikely one.
 *
 * It should be read as "I expect 'x' to be false".
 *
 * Note that this turns into a no-op but gives information to the compiler.
 *
 * Note: This should only be used for branches where the
 * compiler suggested this annotation when running with
 * the runtime profiler (-fprofile-generate with GCC).
 *
 * @brief The expression in the test.
 */
#define swift_unlikely_branch(x) __builtin_expect(!!(x), 0)
#else

/* Default version when not using a GNUC compiler */

#define swift_likely_branch(x) (x)
#define swift_unlikely_branch(x) (x)
#endif

#endif /* SWIFT_EXPECT_H */
