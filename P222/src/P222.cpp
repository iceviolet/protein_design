/*
 * P222.cpp
 *
 *  Created on: Aug 28, 2015
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
 * @brief design routine
 * @param vergil
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
	Design(vergil);

	// Delete the VERGIL instance and clean up
	delete vergil;

	// Finalize MPI
	MPI_Finalize();

	return 0;
}

void Design(Vergil* vergil) {

	//set design parameters
	double design_beta = 0.5;
	double ref_beta = ROOM_TEMPERATURE_BETA;
	double twobody_energy_cap = 30;

	//file prefix
	std::string home_dir = "/scratch/03048/huixi";
	//std::string home_dir = "/Users/Violet/data";
	std::string input_directory =
			home_dir + "/peptide_design/fullseq_design_approved_by_Jeff";
	std::string output_directory =
			home_dir + "/peptide_design/P2_9";
	std::string pdb_filename = input_directory
			+ "/P222_9_X0_X6_X12_X17_NH2cap.pdb";
//	std::string pdb_filename = input_directory
//			+ "/P222_9_X0_X12_NH2cap.pdb";
  std::string output_fileprefix = output_directory + "/C4459_P2_9";


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
	PotentialUnfoldedFreeEnergy unfolded_compute(
			PotentialUnfoldedFreeEnergy::AMBER_KONO,
			vergil->forcefield_parameter_library(), vergil->topology_library(),
			vergil->conformer_library());
	unfolded_compute.set_beta(ref_beta);
	unfolded_compute.AddPotential(vergil->potential("dihedral"));
	unfolded_compute.AddPotential(vergil->potential("vanderwaals"));
	unfolded_compute.AddPotential(vergil->potential("electrostatic"));
	unfolded_compute.GenerateUnfoldedFreeEnergies();
	Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  std::string outfile_energy = output_fileprefix + "_bundle_a3.17_c3.25_beta78.1_theta35to50_d9.5to10.5.csv";
  FILE* outfile = fopen(outfile_energy.c_str(), "w");
  fprintf(outfile, "UnitCell_length_a(nm), UnitCell_length_c(nm), beta(deg), theta(deg), d(A), Internal_energy(Kcal/mol)\n");
  fclose(outfile);

	//for (double a = 28.6; a < 29.1; a+=0.1) {
	//	for (double c = 32.8; c < 33.9; c+=0.1) {
  double a = 31.7;
  double c = 32.5;
		 // for (double beta = 95.5; beta < 97.6; beta += 0.1) {
		    for (double theta = 35; theta < 51; ++theta ) {
		      for (double d = 9.5; d < 10.6; d += 0.1) {
          	double beta = 78.1;
		//double theta = 44;
		//double d = 10.1;
		// Load scaffold
          std::string scaffold ="protein and backbone";

          InputPDBScaffold input(vergil->domain(), scaffold);
          input.KeepSitesWithMissingAtoms();
          input.Read(pdb_filename);
          Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

          //Type each site in the domain using the original P222_9 sequence
          std::string typelist[] = { "ASP", "GLY", "ARG", "ILE", "GLU", "GLY", "MET", "ALA", "GLU", "ALA", "ILE", "LYS", "LYS",
              "MET", "ALA", "TYR", "ASN", "ILE", "ALA", "ASP", "MET", "ALA", "GLY", "ARG", "ILE", "TRP", "GLY", "GLU", "ALA",
              "NH2" };

          for (Domain::SectionIterator jt = vergil->domain()->SectionIterator_Begin(); jt != vergil->domain()->SectionIterator_End(); ++jt) {
            size_t i = 0;
            for (Section::SiteIterator it = jt->SiteIterator_Begin(); it != jt->SiteIterator_End(); ++it) {
              it->Add(typelist[i]);
              ++i;
              //it->Add("GLY");
            }
          }

          // Build all conformers
          vergil->BuildDomain();
          vergil->TrimDomain(30.0);

          //move the domain(ASU) to the correct position in the unit cell
          //rotate ASU around x axis by 90 degree
          const Matrix m = Geometry::Transformation(Vector3(1, 0, 0), 90 * DEGREES_TO_RADIANS);
          vergil->domain()->TransformBy(m);
          //rotate ASU around y(b) axis by theta degree
          const Matrix m_y = Geometry::Transformation(Vector3(0, 1, 0), theta * DEGREES_TO_RADIANS);
          vergil->domain()->TransformBy(m_y);
          //translate ASU by d A along the bisect of a and c in the unit cell
          const Vector3 offset = Vector3(d * cos(beta/2 * DEGREES_TO_RADIANS), 0, d * sin(beta/2 * DEGREES_TO_RADIANS));
          const Matrix m_move = Geometry::Transformation(Vector3(0, 0, 0), 0, offset);
          vergil->domain()->TransformBy(m_move);

          //Set up symmetry and build lattice
          vergil->domain()->set_unit_cell_parameters("P2", a, 1000, c, 90.0, beta, 90.0);
          SymmetryGenerator<Domain> symexp("P2", a, 1000, c, 90.0, beta, 90.0);
          symexp.BuildLattice(*vergil->domain(), -1, 0, -1, 1, 0, 1);
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
          fprintf(outfile, "%10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f\n", a/10, c/10, beta, theta, d, meanfield_energy);
          fclose(outfile);

          //Output standard files (pdb, psf, seq, csv)
          std::string output_file = output_directory + "/C4459_P2_9_bundle_nofixcore_a" + Log->to_str(a) + "_c" + Log->to_str(c) + "_beta" + Log->to_str(beta) + "_theta" + Log->to_str(theta) + "_d" + Log->to_str(d) + "_mp_type_wCap";
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
	//	  }
	//	}
//	}
}
