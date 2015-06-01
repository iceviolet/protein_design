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
 * @date   Jul 8 2014
 * @author Krishna Vijayendran <kvija@mail.med.upenn.edu>
 * @author Jose Villegas <josevill@sas.upenn.edu>
 * @author Chris Von Bargen <vonc@sas.upenn.edu>
 */

#include "mpi.h"

#include "compute_unfolded_freeenergy.h"
#include "domain.h"
#include "function_environmental_energy.h"
#include "function_freeenergy_with_reference.h"
#include "function_meanfieldenergy.h"
#include "input_pdb_scaffold.h"
#include "problem_domainprobability.h"
#include "timer.h"
#include "vergil.h"
#include "vergil_files.h"
#include "vergil_config.h"

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
  unsigned int max_dunbrack_rots = 10;
  std::string pdb_filename = VERGIL_PATH "/working/input/4JKV.pdb";
  std::string output_fileprefix = VERGIL_PATH "/working/output/example_output";

  //Designable sites
  int num_designable_sites = 22;
  int designable_sites[] = { 243, 246, 253, 316, 362, 363, 369, 370, 373, 410, 413, 420, 424, 454, 457, 458, 460, 464,
      471, 537, 542, 543 };

  //Construct scaffold selection text from the array of designable sites
  std::string scaffold_seltext = "chain A and resid 190 to 550";  //Only consider residues 190-550 on chain A
  scaffold_seltext += " and (backbone or not resid ";  //Fix the backbone atoms, or all the atoms in sites that are NOT designable
  for (int i = 0; i < num_designable_sites; ++i) {
    scaffold_seltext += Log->to_str(designable_sites[i]) + " ";
  }
  scaffold_seltext += ")";

  //Create a local Timer instance
  Timer timer;

  //Load force field, dunbrack library, and scaffold from PDB file
  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();
  vergil->LoadDunbrackRotamerLibrary(2002, max_dunbrack_rots);
  vergil->LoadScaffoldFromPDB(pdb_filename, scaffold_seltext);

  //The three usual computes are automatically loaded into Vergil: "dihedral", "vanderwaals", and "electrostatic"
  //  for the AMBER84 modified ff, where "vanderwaals" is the KonoModified compute with hbonding

  //Compute and time the reference energy
  timer.Stamp();
  Log->print("Computing reference energies....");
  ComputeUnfoldedFreeEnergy unfolded_compute(ComputeUnfoldedFreeEnergy::AMBER_KONO, vergil->forcefield_parameter_library(),
                                             vergil->topology_library(), vergil->conformer_library());
  unfolded_compute.set_beta(ROOM_TEMPERATURE_BETA);
  unfolded_compute.set_grid_space(90);
  unfolded_compute.AddCompute(vergil->compute("dihedral"));
  unfolded_compute.AddCompute(vergil->compute("vanderwaals"));
  unfolded_compute.AddCompute(vergil->compute("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  //Type each site in the domain
  Path* chain = &(*vergil->domain())["4JKV"]["A"][0];

  //Disulfide bridges
  (*chain)[193].reference_type()->set_name("CYX");
  (*chain)[213].reference_type()->set_name("CYX");
  (*chain)[193].ConnectNextSiteWithIterface((*chain)[213], "DISU");

  (*chain)[217].reference_type()->set_name("CYX");
  (*chain)[295].reference_type()->set_name("CYX");
  (*chain)[217].ConnectNextSiteWithIterface((*chain)[295], "DISU");

  (*chain)[314].reference_type()->set_name("CYX");
  (*chain)[390].reference_type()->set_name("CYX");
  (*chain)[314].ConnectNextSiteWithIterface((*chain)[390], "DISU");

  (*chain)[490].reference_type()->set_name("CYX");
  (*chain)[507].reference_type()->set_name("CYX");
  (*chain)[490].ConnectNextSiteWithIterface((*chain)[507], "DISU");

  int num_amino_acids = 20;  //All 20 natural amino acids
  std::string amino_acid_set[] = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
      "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
  for (int i = 0; i < num_designable_sites; ++i) {
    for (int j = 0; j < num_amino_acids; ++j) {
      (*chain)[designable_sites[i]].Add(amino_acid_set[j]);
    }
  }

  // Build all conformers
  vergil->BuildDomain();

  //Trim rotamers that clash with the backbone
  vergil->TrimDomain(30.0);

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

  // Set up the problem and solve the unconstrained problem
  FunctionFreeEnergyWithReference free_energy(&energy, design_beta, &unfolded_compute);
  vergil->InitProbabilityProblem(free_energy);
  vergil->SolveProbabilityProblem();

  //Environmental Energy function
  double Eenv_target = 0.0;
  FunctionEnvironmentalEnergy environmental_energy(*vergil->domain(), ENVIRONMENTAL_ENERGY_FILE);

  //Print out the Environmental Energy value
  VectorN w = vergil->ConformerProbabilityVector();
  Log->print(environmental_energy.to_string(w));

  //Output standard files (pdb, psf, seq, csv)
  vergil->StandardOutput(output_fileprefix + "_unconstrained");

  //Now add the Environmental Energy constraint and resolve the problem
  vergil->problem()->AddConstraint("Eenv", &environmental_energy, Eenv_target);
  vergil->SolveProbabilityProblem();

  //Print out the Environmental Energy value to check its been constrainted properly
  w = vergil->ConformerProbabilityVector();
  Log->print(environmental_energy.to_string(w));

  //Output standard files (pdb, psf, seq, csv)
  vergil->StandardOutput(output_fileprefix);
}

