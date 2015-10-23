/*
 * P422.cpp
 *
 *  Created on: Oct 21, 2015
 *      Author: Violet
 */

#include "mpi.h"
#include "constants.h"

#include "potential_unfolded_freeenergy.h"
#include "domain.h"
#include "symmetry_generator.h"
#include "function_meanfieldenergy_lattice.h"
#include "function_freeenergy_with_reference.h"
#include "function_meanfieldenergy.h"
#include "input_pdb_scaffold.h"
#include "output_pdb.h"
#include "output_psf.h"
#include "problem_domainprobability.h"
#include "conformer_states.h"
#include "timer.h"
#include "geometry.h"
#include "vergil.h"
#include "vergil_files.h"
#include "vergil_config.h"

using namespace vergil_ns;

/**
 * @brief design routine with P4212 symmetry
 * @param vergil
 */
void DesignP4212(Vergil* vergil);

/**
 * @brief design routine with P4 symmetry
 * @param vergil
 */
void DesignP4(Vergil* vergil);

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
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        if(provided < MPI_THREAD_FUNNELED)
        {
           printf("MPI cannot support funneled!\n");
           MPI_Abort(MPI_COMM_WORLD,0);
        }

  // Create a VERGIL instance, pass arguments;
  Vergil *vergil = new Vergil();

  // Run the routine
  //DesignP4212(vergil);
  DesignP4(vergil);

  // Delete the VERGIL instance and clean up
  delete vergil;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

void DesignP4212(Vergil* vergil) {
  //set design parameters
  double design_beta = 0.5;
  double ref_beta = ROOM_TEMPERATURE_BETA;
  double twobody_energy_cap = 30;

  std::string home_dir = "/scratch/03048/huixi";
  std::string input_directory = home_dir + "/peptide_design/fullseq_design_approved_by_Jeff";
  std::string output_directory = home_dir + "/peptide_design/P4212";
  std::string pdb_filename = input_directory + "/C4459a_chainA_fixed_adg_NH2cap_for_vergil.pdb";
  std::string output_fileprefix = output_directory + "/C4459_P4212";

  //Create a local Timer instance
  Timer timer;

  //Load force field and design topology
  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();
  unsigned int max_dunbrack_rots = 81;
  vergil->LoadDunbrackRotamerLibrary(2002, max_dunbrack_rots);

  //Compute and time the reference energy
  timer.Stamp();
  Log->print("Computing reference energies....");
  PotentialUnfoldedFreeEnergy unfolded_compute(PotentialUnfoldedFreeEnergy::AMBER_KONO,
                                             vergil->forcefield_parameter_library(), vergil->topology_library(),
                                             vergil->conformer_library());
  unfolded_compute.set_beta(ref_beta);
  unfolded_compute.AddPotential(vergil->potential("dihedral"));
  unfolded_compute.AddPotential(vergil->potential("vanderwaals"));
  unfolded_compute.AddPotential(vergil->potential("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  std::string outfile_energy = output_fileprefix + "_a3.73to4.34.csv";
  FILE* outfile = fopen(outfile_energy.c_str(), "w");
  fprintf(outfile, "UnitCell_length_a(nm), Internal_energy(Kcal/mol)\n");
  fclose(outfile);

  for (double a = 37.3; a < 43.5; a+=0.1) {
    std::string scaffold ="protein and backbone";

    InputPDBScaffold input(vergil->domain(), scaffold);
    input.KeepSitesWithMissingAtoms();
    input.Read(pdb_filename);
    Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

    //Type each site in the domain using the original P422_Cu_1 sequence
    //DQEIRQM AEWIKKM AQMIDKM AHRIDRE A-NH2
    std::string typelist[] = { "ASP", "GLN", "GLU", "ILE", "ARG", "GLN", "MET", "ALA", "GLU", "TRP", "ILE", "LYS", "LYS",
        "MET", "ALA", "GLN", "MET", "ILE", "GLU", "LYS", "MET", "ALA", "HIS", "ARG", "ILE", "ASP", "ARG", "GLU", "ALA",
        "NH2" };

    for (Domain::SectionIterator jt = vergil->domain()->SectionIterator_Begin(); jt != vergil->domain()->SectionIterator_End(); ++jt) {
      size_t i = 0;
      for (Section::SiteIterator it = jt->SiteIterator_Begin(); it != jt->SiteIterator_End(); ++it) {
        it->Add(typelist[i]);
        ++i;
      }
    }

    // Build all conformers
     vergil->BuildDomain();
     vergil->TrimDomain(30.0);

     //rotate the domain(ASU) 45 degree around the 2-fold symmetry axis in the unit cell
     vergil->domain()->TransformBy(Geometry::Transformation(Vector3(0, 0, 1), 45.0, 0.0));

     //Set up symmetry and build lattice
     vergil->domain()->set_unit_cell_parameters("P4212", a, a, 1000, 90.0, 90.0, 90.0);
     SymmetryGenerator<Domain> symexp("P4212", a, a, 1000, 90.0, 90.0, 90.0);
     symexp.BuildLattice(*vergil->domain(), -1, -1, 0, 1, 1, 0);
     std::vector<Domain>* symmetry_related_elements = symexp.lattice();

     //Trim rotamers that clash with the scaffold
     vergil->TrimDomain(30.0, symmetry_related_elements);

     // Set up energy function
     timer.Stamp();
     Log->print("Setting up energy function...");
     FunctionMeanFieldEnergyLattice energy(vergil->domain(), symmetry_related_elements);
     energy.set_pairwise_energy_cap_values(twobody_energy_cap, twobody_energy_cap);
     energy.AddPotential(vergil->potential("dihedral"));
     energy.AddPotential(vergil->potential("vanderwaals"));
     energy.AddPotential(vergil->potential("electrostatic"));
     energy.ComputeEnergies();
     Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

     // Initialize probabilities
      vergil->domain()->SetUniformConformerProbabilities();

      // Set up the problem and solve the unconstrained problem
      FunctionFreeEnergyWithReference free_energy(&energy, design_beta, &unfolded_compute);
      vergil->InitProbabilityProblem(free_energy);
      vergil->SolveProbabilityProblem();

      //Output Energy
      VectorN prob = vergil->ConformerProbabilityVector();
      double meanfield_energy = free_energy.internal_energy()->Value(prob);
      FILE* outfile = fopen(outfile_energy.c_str(), "a");
      fprintf(outfile, "%10.5f, %10.5f\n", a/10, meanfield_energy);
      fclose(outfile);

      //Output standard files (pdb, psf, seq, csv)
      std::string output_file = output_directory + "/C4459_P4212_nofixcore_a" + Log->to_str(a) + "_mp_type_wCap";
      //vergil->StandardOutput(output_file);

      for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End(); ++it) {
        it->SetAggregatedTypeProbability();
      }
      vergil->RenameSymmetrySegname(symmetry_related_elements);
      OutputPDB output1(output_file + "_lattice.pdb", vergil->domain());
      output1.WriteLattice(symmetry_related_elements);
      Log->print("Most Probable PDB -- " + output_file + "_lattice.pdb");

      vergil->domain()->Clear();
  }
}

void DesignP4(Vergil* vergil) {
  //set design parameters
  double design_beta = 0.5;
  double ref_beta = ROOM_TEMPERATURE_BETA;
  double twobody_energy_cap = 30;

  std::string home_dir = "/scratch/03048/huixi";
  //std::string home_dir = "/Users/Violet/data";
  std::string input_directory = home_dir + "/peptide_design/fullseq_design_approved_by_Jeff";
  std::string output_directory = home_dir + "/peptide_design/P4";
  std::string pdb_filename = input_directory + "/P222_9_X0_X6_X12_X17_NH2cap.pdb";
  std::string output_fileprefix = output_directory + "/C4459_P4212";

  //Create a local Timer instance
  Timer timer;

  //Load force field and design topology
  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();
  unsigned int max_dunbrack_rots = 81;
  vergil->LoadDunbrackRotamerLibrary(2002, max_dunbrack_rots);

  //Compute and time the reference energy
  timer.Stamp();
  Log->print("Computing reference energies....");
  PotentialUnfoldedFreeEnergy unfolded_compute(PotentialUnfoldedFreeEnergy::AMBER_KONO,
                                             vergil->forcefield_parameter_library(), vergil->topology_library(),
                                             vergil->conformer_library());
  unfolded_compute.set_beta(ref_beta);
  unfolded_compute.AddPotential(vergil->potential("dihedral"));
  unfolded_compute.AddPotential(vergil->potential("vanderwaals"));
  unfolded_compute.AddPotential(vergil->potential("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  std::string outfile_energy = output_fileprefix + "_a4.04_theta0to45_d9to22_alpha0to45.csv";
  FILE* outfile = fopen(outfile_energy.c_str(), "w");
  fprintf(outfile, "UnitCell_length_a(nm), theta(deg), d(A), alpha(deg), Internal_energy(Kcal/mol)\n");
  fclose(outfile);

  double a = 40.4;
  for (double theta = 0.0; theta < 50; theta+=5) {
    for (double d = 9; d < 23; ++d) {
      for (double alpha = 0; alpha < 50; alpha+=5) {

        std::string scaffold ="protein and backbone";

        InputPDBScaffold input(vergil->domain(), scaffold);
        input.KeepSitesWithMissingAtoms();
        input.Read(pdb_filename);
        Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

        //Type each site in the domain using the original P422_Cu_1 sequence
        //DQEIRQM AEWIKKM AQMIDKM AHRIDRE A-NH2
        std::string typelist[] = { "ASP", "GLN", "GLU", "ILE", "ARG", "GLN", "MET", "ALA", "GLU", "TRP", "ILE", "LYS", "LYS",
            "MET", "ALA", "GLN", "MET", "ILE", "GLU", "LYS", "MET", "ALA", "HIS", "ARG", "ILE", "ASP", "ARG", "GLU", "ALA",
            "NH2" };

        for (Domain::SectionIterator jt = vergil->domain()->SectionIterator_Begin(); jt != vergil->domain()->SectionIterator_End(); ++jt) {
          size_t i = 0;
          for (Section::SiteIterator it = jt->SiteIterator_Begin(); it != jt->SiteIterator_End(); ++it) {
            it->Add(typelist[i]);
            //it->Add("GLY");
            ++i;
          }
        }

        // Build all conformers
         vergil->BuildDomain();
         vergil->TrimDomain(30.0);

         //move the domain(ASU) to the correct position in the unit cell
         //rotate ASU around z(c) axis by theta degree
         const Matrix m_z = Geometry::Transformation(Vector3(0, 0, 1), theta * DEGREES_TO_RADIANS);
         vergil->domain()->TransformBy(m_z);
         //translate ASU by d along the vector that is alpha degree away from a in the unit cell
         const Vector3 offset = Vector3(d * cos(alpha * DEGREES_TO_RADIANS), d * sin(alpha * DEGREES_TO_RADIANS), 0);
         const Matrix m_move = Geometry::Transformation(Vector3(0, 0, 0), 0, offset);
         vergil->domain()->TransformBy(m_move);

         //Set up symmetry and build lattice
         vergil->domain()->set_unit_cell_parameters("P4", a, a, 1000, 90.0, 90.0, 90.0);
         SymmetryGenerator<Domain> symexp("P4", a, a, 1000, 90.0, 90.0, 90.0);
         symexp.BuildLattice(*vergil->domain(), -1, -1, 0, 1, 1, 0);
         std::vector<Domain>* symmetry_related_elements = symexp.lattice();

         //Trim rotamers that clash with the scaffold
         vergil->TrimDomain(30.0, symmetry_related_elements);

         // Set up energy function
         timer.Stamp();
         Log->print("Setting up energy function...");
         FunctionMeanFieldEnergyLattice energy(vergil->domain(), symmetry_related_elements);
         energy.set_pairwise_energy_cap_values(twobody_energy_cap, twobody_energy_cap);
         energy.AddPotential(vergil->potential("dihedral"));
         energy.AddPotential(vergil->potential("vanderwaals"));
         energy.AddPotential(vergil->potential("electrostatic"));
         energy.ComputeEnergies();
         Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

         // Initialize probabilities
          vergil->domain()->SetUniformConformerProbabilities();

          // Set up the problem and solve the unconstrained problem
          FunctionFreeEnergyWithReference free_energy(&energy, design_beta, &unfolded_compute);
          vergil->InitProbabilityProblem(free_energy);
          vergil->SolveProbabilityProblem();

          //Output Energy
          VectorN prob = vergil->ConformerProbabilityVector();
          double meanfield_energy = free_energy.internal_energy()->Value(prob);
          FILE* outfile = fopen(outfile_energy.c_str(), "a");
          fprintf(outfile, "%10.5f, %10.5f, %10.5f, %10.5f, %10.5f\n", a/10, theta, d, alpha, meanfield_energy);
          fclose(outfile);

          //Output standard files (pdb, psf, seq, csv)
          std::string output_file = output_directory + "/C4459_P4_bundle_nofixcore_a" + Log->to_str(a) + "_theta" + Log->to_str(theta) + "_d" + Log->to_str(d) + "_alpha" + Log->to_str(alpha) + "_mp_type_wCap";
          //vergil->StandardOutput(output_file);

          for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End(); ++it) {
            it->SetAggregatedTypeProbability();
          }
          vergil->RenameSymmetrySegname(symmetry_related_elements);
          if(meanfield_energy < 0) {
            OutputPDB output1(output_file + "_lattice.pdb", vergil->domain());
            output1.WriteLattice(symmetry_related_elements);
            Log->print("Most Probable PDB -- " + output_file + "_lattice.pdb");
          }

          vergil->domain()->Clear();
      }
    }
  }
}
