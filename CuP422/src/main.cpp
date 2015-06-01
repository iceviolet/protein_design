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
 * Description of the main class routines
 *
 * @date   Jun 1 2015
 * @author Author <email>
 */

#include "mpi.h"

#include "simple.h"

#include "vergil.h"
#include "vergil_files.h"
#include "vergil_config.h"

using namespace vergil_ns;
using namespace project_ns;

/**
 * @brief
 * Design routine description...
 *
 * @param vergil   instance of the vergil superclass
 */
void run_design(Vergil *vergil) {

  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();

  Simple s;
  std::printf("%d\n", s.Sample(6));

  //etc.
}

/**
 * @brief
 * Main routine to create a Vergil instance and run the design
 *
 * @param argc   the number of arguments
 * @param argv   the arguments
 * @return main return
 */
int main(int argc, char **argv) {

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Create a VERGIL instance, pass arguments;
  Vergil *vergil = new Vergil();

  // Run design routine
  run_design(vergil);

  // Delete the VERGIL instance and clean up
  delete vergil;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
