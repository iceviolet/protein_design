/* -*- c++ -*- -------------------------------------------------------------------------
 * Saven Lab Boilerplate
 * University of Pennsylvania, Department of Chemistry, Makineni Theoretical Laboratories
 * Jeffery G. Saven, saven@sas.upenn.edu, http://saven.chem.upenn.edu
 *
 * Copyright (c) 2015
 * The Board of Trustees of the University of Pennsylvania.
 * All Rights Reserved.
 ------------------------------------------------------------------------------------- */

/**
 * @internal @file
 * @brief
 * Implements class project_ns::Simple
 *
 * @date   Jun 1 2015
 * @author Author <email>
 */

#include "simple.h"

using namespace project_ns;

Simple::Simple()
    : example_member_(0.0) {
}

Simple::~Simple() {
}

/**
 * @internal
 * Implementation specific (optional) comment
 */
bool Simple::Sample(double x) const {
  return x >= example_member_;
}
