/* -*- c++ -*- -------------------------------------------------------------------------
 * This file is part of the VERGIL protein design package.
 *
 * VERGIL - Probabilistic Protein Design Suite
 * University of Pennsylvania, Department of Chemistry, Makineni Theoretical Laboratories
 * Jeffery G. Saven, saven@sas.upenn.edu, http://saven.chem.upenn.edu
 * and others as listed in the AUTHORS file in the top-level source directory
 *
 * Copyright (c) 2011-2013
 * The Board of Trustees of the University of Pennsylvania.
 * All Rights Reserved.
 *
 * See the README file in the top-level Vergil directory.
 * Please consult the development guide in the Vergil Manual, located in the
 * docs directory, when contributing to the Vergil project. All additions should
 * have accompanying tests to ensure the functionality of the suite.
 *
 * To any researchers, we ask that you cite the research papers on the Vergil package.
 ------------------------------------------------------------------------------------- */

/**
 * @file
 * @brief
 * Design the smoothened protein to be water soluble
 *
 * @date   Jul 28 2014
 * @author Jose Villegas <josevill@sas.upenn.edu>
 */

#include "mpi.h"

#include "compute_unfolded_freeenergy.h"
#include "domain.h"
#include "function_environmental_energy.h"
#include "function_protein_tools.h"
#include "function_freeenergy_with_reference.h"
#include "function_meanfieldenergy.h"
#include "input_pdb_scaffold.h"
#include "problem_domainprobability.h"
#include "timer.h"
#include "vergil.h"
#include "vergil_files.h"
#include "vergil_config.h"

#include <iostream>

using namespace vergil_ns;

/**
 * @brief
 * Desing the smoothened protein to be water soluble...
 *
 * @param vergil   the vergil instance
 */
void Design(Vergil* vergil);

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

  // Run the routine
  Design(vergil);

  // Delete the VERGIL instance and clean up
  delete vergil;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

void Design(Vergil* vergil) {

  //Variables for the design
  double design_beta = 0.5;
  unsigned int max_dunbrack_rots = 2;

  Timer timer;

  //Load force field and dunbrack library
  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();
  vergil->LoadDunbrackRotamerLibrary(2002, max_dunbrack_rots);
  vergil->LoadScaffoldFromPDB(VERGIL_PATH "/working/input/4dzk_clean.pdb", "name C CA N O");

  //Compute and time the reference energy
  timer.Stamp();
  Log->print("Computing reference energies....");
  ComputeUnfoldedFreeEnergy unfolded_compute(ComputeUnfoldedFreeEnergy::AMBER_KONO, vergil->forcefield_parameter_library(),
                                             vergil->topology_library(), vergil->conformer_library());
  unfolded_compute.set_beta(0.5);  //ROOM_TEMPERATURE_BETA
  unfolded_compute.set_grid_space(10);
  unfolded_compute.AddCompute(vergil->compute("dihedral"));
  unfolded_compute.AddCompute(vergil->compute("vanderwaals"));
  unfolded_compute.AddCompute(vergil->compute("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  int num_amino_acids = 20;
  std::string amino_acid_set[] = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
      "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End(); ++it)
    for (int j = 0; j < num_amino_acids; ++j) 
      it->Add(amino_acid_set[j]);

  // Build all conformers
  vergil->BuildDomain();

  //Trim rotamers that clash with the backbone
  vergil->TrimDomain(30.0);  //1000.0

  // Set up energy function
  timer.Stamp();
  Log->print("Setting up energy function...");
  FunctionMeanFieldEnergy energy(vergil->domain());
  energy.set_pairwise_energy_cap_values(30.0, 30.0);
  energy.AddCompute(vergil->compute("dihedral"));
  energy.AddCompute(vergil->compute("vanderwaals"));
  energy.AddCompute(vergil->compute("electrostatic"));
  energy.ComputeEnergies();
  Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

  // Initialize probabilities
  vergil->domain()->SetUniformConformerProbabilities();

  //Environmental Energy function
  //double Eenv_target = 0.0;  //-48.0489;
  //FunctionEnvironmentalEnergy environmental_energy(*vergil->domain(), ENVIRONMENTAL_ENENERGY_FILE);

  double percent_charged_minimum = 30.0;
  FunctionProteinTools percent_charged(*vergil->domain(), FunctionProteinTools::PERCENT_CHARGED);

  //double percent_hydrophobic_minimum = 40.0;
  FunctionProteinTools percent_hydrophobic(*vergil->domain(), FunctionProteinTools::PERCENT_HYDROPHOBIC);

  double hydropathy_minimum = -1.0;
  double hydropathy_maximum = 1.0;
  FunctionProteinTools hydropathy(*vergil->domain(), FunctionProteinTools::HYDROPATHY_SCORE);

  // Set up the problem
  FunctionFreeEnergyWithReference free_energy(&energy, design_beta, &unfolded_compute);
  vergil->InitProbabilityProblem(free_energy);
  //vergil->problem()->AddConstraint("Eenv", &environmental_energy, Eenv_target);
  vergil->problem()->AddConstraint("Chrg", &percent_charged, ">", percent_charged_minimum);
  //vergil->problem()->AddConstraint("Hphobic", &percent_hydrophobic, ">", percent_hydrophobic_minimum);
  vergil->problem()->AddConstraint("Hydropathy_min", &hydropathy, ">", hydropathy_minimum);
  vergil->problem()->AddConstraint("Hydropathy_max", &hydropathy, "<", hydropathy_maximum);
  vergil->SolveProbabilityProblem();

  VectorN w = vergil->ConformerProbabilityVector();
  Log->print(percent_charged.to_string(w));
  Log->print(percent_hydrophobic.to_string(w));
  Log->print(hydropathy.to_string(w));

  //Output
  vergil->StandardOutput(VERGIL_PATH "/working/output/charged_constraint_test");
}

