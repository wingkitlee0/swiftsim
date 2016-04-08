/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_PARSER_H
#define SWIFT_PARSER_H

/* Config parameters. */
#include "../config.h"

/* Some constants. */
#define PARSER_MAX_LINE_SIZE 256
#define PARSER_MAX_NO_OF_PARAMS 512

/* A parameter in the input file */
struct parameter {
  char name[PARSER_MAX_LINE_SIZE];
  char value[PARSER_MAX_LINE_SIZE];
};

/* The array of parameters read from a file */
struct swift_params {
  struct parameter data[PARSER_MAX_NO_OF_PARAMS];
  int count;
};

/* Public API. */
void parser_read_file(const char *file_name, struct swift_params *params);
void parser_print_params(const struct swift_params *params);
void parser_write_params_to_file(const struct swift_params *params,
                                 const char *file_name);
char parser_get_param_char(const struct swift_params *params, const char *name);
int parser_get_param_int(const struct swift_params *params, const char *name);
float parser_get_param_float(const struct swift_params *params,
                             const char *name);
double parser_get_param_double(const struct swift_params *params,
                               const char *name);
void parser_get_param_string(const struct swift_params *params,
                             const char *name, char *retParam);

#endif /* SWIFT_PARSER_H */