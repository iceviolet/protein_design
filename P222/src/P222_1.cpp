/*
 * P222_1.cpp
 *
 *  Created on: Oct 7, 2015
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
void DesignP2(Vergil* vergil, std::vector<double> parameters);

/**
 * @brief
 * Main routine to create a Vergil instance and run the design
 *
 * @param argc   the number of arguments, including upper and lower boundaries for 6 parameters
 * @param argv   the arguments
 * @return main return
 */
int main(int argc, char **argv) {

	if (argc != 13) {
		Log->error(FLERR, "The number of arguments passed to the program must be 12");
	}

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

	// Set up the unit cell and ASU placement parameters
	std::vector<double> parameters;
	for (int i = 1; i < argc; ++i) {
		parameters.push_back(atof(argv[i]));
	}

	//Run the routine
	DesignP2(vergil, parameters);

	// Delete the VERGIL instance and clean up
	delete vergil;

	// Finalize MPI
	MPI_Finalize();

	return 0;
}

void DesignP2(Vergil* vergil, std::vector<double> parameters) {

	//set design parameters
	double design_beta = 0.5;
	double ref_beta = ROOM_TEMPERATURE_BETA;
	double twobody_energy_cap = 30;

	//file prefix
	std::string home_dir = "/scratch/03048/huixi";
	//std::string home_dir = "/Users/Violet/data";
	std::string input_directory = home_dir
			+ "/peptide_design/fullseq_design_approved_by_Jeff";
	std::string output_directory = home_dir + "/peptide_design/P2_1";
	std::string pdb_filename = input_directory
			+ "/P222_9_X0_X6_X12_X17_NH2cap.pdb";
	std::string output_fileprefix = output_directory + "/C4459_P2_1";

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

	double a_min = parameters[0];
	double a_max = parameters[1];
	double c_min = parameters[2];
	double c_max = parameters[3];
	double beta_min = parameters[4];
	double beta_max = parameters[5];
	double theta_min = parameters[6];
	double theta_max = parameters[7];
	double d_min = parameters[8];
	double d_max = parameters[9];
	double alpha_min = parameters[10];
	double alpha_max = parameters[11];

	std::string outfile_energy = output_fileprefix + "_bundle_a" + Log->to_str(a_min) + "to" + Log->to_str(a_max) +
	"_c" + Log->to_str(c_min) + "to" + Log->to_str(c_max) + "_beta" + Log->to_str(beta_min) + "to" + Log->to_str(beta_max) +
	"_theta" + Log->to_str(theta_min) + "to" + Log->to_str(theta_max) + "_d" + Log->to_str(d_min) + "to" + Log->to_str(d_max) +
	"_alpha" + Log->to_str(alpha_min) + "to" + Log->to_str(alpha_max) + ".csv";
	FILE* outfile = fopen(outfile_energy.c_str(), "w");
	fprintf(outfile,
			"UnitCell_length_a(nm), UnitCell_length_c(nm), beta(deg), theta(deg), d(A), alpha(deg), Internal_energy(Kcal/mol)\n");
	fclose(outfile);

	for (double a = a_min; a < a_max + 0.1; a += 0.1) {
		for (double c = c_min; c < c_max + 0.1; c += 0.1) {
			for (double beta = beta_min; beta < beta_max + 0.1; beta += 0.1) {
				for (double theta = theta_min; theta < theta_max + 1; ++theta) {
					for (double d = d_min; d < d_max + 0.1; d += 0.1) {
						for (double alpha = alpha_min; alpha < alpha_max + 0.1;
								++alpha) {
							// Load scaffold
							std::string scaffold =
									"protein and backbone or site 30";

							InputPDBScaffold input(vergil->domain(), scaffold);
							input.KeepSitesWithMissingAtoms();
							input.Read(pdb_filename);
							Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

							//Type each site in the domain using the original P222_9 sequence
							std::string typelist[] = { "ASP", "GLY", "LYS",
									"ILE", "GLU", "GLY", "MET", "ALA", "GLU",
									"ALA", "ILE", "LYS", "LYS", "MET", "ALA",
									"ASN", "ASN", "ILE", "GLU", "GLN", "MET",
									"ALA", "GLY", "TRP", "ILE", "TRP", "GLY",
									"GLU", "ALA", "NH2" };

							for (Domain::SectionIterator jt =
									vergil->domain()->SectionIterator_Begin();
									jt
											!= vergil->domain()->SectionIterator_End();
									++jt) {
								size_t i = 0;
								for (Section::SiteIterator it =
										jt->SiteIterator_Begin();
										it != jt->SiteIterator_End(); ++it) {
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
							const Matrix m = Geometry::Transformation(
									Vector3(1, 0, 0), 90 * DEGREES_TO_RADIANS);
							vergil->domain()->TransformBy(m);
							//rotate ASU around y(b) axis by theta degree
							const Matrix m_y = Geometry::Transformation(
									Vector3(0, 1, 0),
									theta * DEGREES_TO_RADIANS);
							vergil->domain()->TransformBy(m_y);
							//translate ASU by d along the vector that is alpha degree away from a in the unit cell
							const Vector3 offset = Vector3(
									d * cos(alpha * DEGREES_TO_RADIANS), 0,
									d * sin(alpha * DEGREES_TO_RADIANS));
							const Matrix m_move = Geometry::Transformation(
									Vector3(0, 0, 0), 0, offset);
							vergil->domain()->TransformBy(m_move);

							//Set up symmetry and build lattice
							vergil->domain()->set_unit_cell_parameters("P2", a,
									1000, c, 90.0, beta, 90.0);
							SymmetryGenerator<Domain> symexp("P2", a, 1000, c,
									90.0, beta, 90.0);
							symexp.BuildLattice(*vergil->domain(), -1, 0, -1, 1,
									0, 1);
							std::vector<Domain>* symmetry_related_elements =
									symexp.lattice();

							//Trim rotamers that clash with the scaffold
							vergil->TrimDomain(30.0, symmetry_related_elements);

							//skip the case when all conformers are removed on any of the sites
							bool skip = 0;
							for (Domain::SectionIterator jt =
									vergil->domain()->SectionIterator_Begin();
									jt
											!= vergil->domain()->SectionIterator_End();
									++jt) {
								for (Section::SiteIterator it =
										jt->SiteIterator_Begin();
										it != jt->SiteIterator_End(); ++it) {
									if (it->NumberConformers() == 0) {
										skip = 1;
										break;
									}
								}
								if (skip)
									break;
							}

							if (skip) {
								vergil->domain()->Clear();
								continue;
							}

							// Set up energy function
							timer.Stamp();
							Log->print("Setting up energy function...");
							FunctionMeanFieldEnergyLattice energy(
									vergil->domain(),
									symmetry_related_elements);
							energy.set_pairwise_energy_cap_values(
									twobody_energy_cap, twobody_energy_cap);
							energy.AddPotential(vergil->potential("dihedral"));
							energy.AddPotential(
									vergil->potential("vanderwaals"));
							energy.AddPotential(
									vergil->potential("electrostatic"));
							energy.ComputeEnergies();
							Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

							// Initialize probabilities
							vergil->domain()->SetUniformConformerProbabilities();

							// Set up the problem and solve the unconstrained problem
							FunctionFreeEnergyWithReference free_energy(&energy,
									design_beta, &unfolded_compute);
							vergil->InitProbabilityProblem(free_energy);
							vergil->SolveProbabilityProblem();

							//Output Energy
							VectorN prob = vergil->ConformerProbabilityVector();
							double meanfield_energy =
									free_energy.internal_energy()->Value(prob);
							FILE* outfile = fopen(outfile_energy.c_str(), "a");
							fprintf(outfile,
									"%10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f\n",
									a / 10, c / 10, beta, theta, d, alpha,
									meanfield_energy);
							fclose(outfile);

							//Output standard files (pdb, psf, seq, csv)
							std::string output_file = output_directory
									+ "/C4459_P2_1_bundle_nofixcore_a" + Log->to_str(a) + "_c" + Log->to_str(c) + "_beta" + Log->to_str(beta) +
							"_theta" + Log->to_str(theta) + "_d" + Log->to_str(d) + "_alpha" + Log->to_str(alpha) + "_mp_type_wCap";
							//vergil->StandardOutput(output_file);

							if (meanfield_energy < -780) {
								for (Domain::SiteIterator it =
										vergil->domain()->SiteIterator_Begin();
										it
												!= vergil->domain()->SiteIterator_End();
										++it) {
									it->SetAggregatedTypeProbability();
								}
								vergil->RenameSymmetrySegname(
										symmetry_related_elements);
								OutputPDB output1(output_file + "_lattice.pdb",
										vergil->domain());
								output1.WriteLattice(symmetry_related_elements);
								Log->print("Most Probable PDB -- " + output_file + "_lattice.pdb");
							}

							vergil->domain()->Clear();
						}
					}
				}
			}
		}
	}
}

