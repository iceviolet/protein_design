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
 * @file
 * @brief
 * Declares class project_ns::Simple
 *
 * @date   Jun 1 2015
 * @author Author <email>
 */

#ifndef PROJECT_SIMPLE_H
#define PROJECT_SIMPLE_H

namespace project_ns {

/**
 * @brief
 * Class description.
 */
class Simple {

 public:

  /**
   * @brief
   * Constructor
   */
  Simple();

  /**
   * @brief
   * Destructor description
   */
  virtual ~Simple();

  /**
   * @brief
   * Description of method
   *
   * @param x  desc of parameter x
   * @return true if condition, false otherwise
   */
  bool Sample(double x) const;

 private:

  /**
   * @brief
   * Description of member
   */
  double example_member_;

  /**
   * @brief
   * Disabled copy constructor
   */
  Simple(const Simple&);

  /**
   * @brief
   * Disabled assignment operator
   */
  void operator=(const Simple&);

};

}

#endif
