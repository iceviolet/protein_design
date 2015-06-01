/*
 * CuP422.cpp
 *
 *  Created on: Nov 20, 2014
 *      Author: Violet
 */

#include "constants.h"
#include "mpi.h"

#include "compute_unfolded_freeenergy.h"
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
 * @brief
 *
 * @param vergil   the vergil instance
 */
void Design(Vergil* vergil);

/**
 * @brief non-design routine
 * @param vergil
 */
void Run(Vergil* vergil);

/**
 * @brief loop through the chi angles of the ligand and find the optimal orientation of the ligand
 *
 * @param ligand
 * @param chi_increment
 * @param chi_stdev
 * @param lp_atom
 * @param output_dir
 * @param angle_tolerance
 */
void OptimizeLigandGeomety(Conformer& ligand, double chi_increment, std::vector<double> chi_stdev, std::string lp_atom,
                           double unit_cell_length, double angle_tolerance, std::string output_dir);

/**
 * @brief get degrees of freedom and the chi angles of the ligand
 * @param ligand
 * @return a vector of the chi values
 */
std::vector<double> GetChiValues(Conformer& ligand);

/**
 * @brief
 * Construct a vector v from the point bisecting the lone pair angle to the lone pair N
 * Calculate the dot product of v and z axis unit vector,
 * get the angle between the lone-pair direction and the z axis, 90 is the optimal value
 * @param ligand
 * @param LP_atom name of the atom that bears the lone pair
 * @return the vector representing the lone pair direction (either on NE2 or ND1)
 */
Vector3 CalculateLonePairDirection(Conformer& ligand, std::string LP_atom);

/**
 * @brief output the results as a table containing chi angles and the lone pair direction angle
 *
 * @param ligand_geometry a map with key as sidechain chi angles, values of angles representing lone pair direction
 */
void OutputStat(std::map<std::vector<double>, std::vector<double> > ligand_geometry, std::string lp_atom,
                std::string output_directroy);

void OutputStructure(Domain* domain, std::vector<double> chis, std::string lp_atom, std::string output_directroy);

void OutputLatticeStructure(Domain* domain, std::vector<double> chis, double unit_cell_length, std::string lp_atom,
                            std::string output_directroy);

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
  //Design(vergil);
  Run(vergil);

  // Delete the VERGIL instance and clean up
  delete vergil;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

void Design(Vergil* vergil) {

  //set design parameters
  double beta = 0.5;
  double ref_beta = ROOM_TEMPERATURE_BETA;
  double twobody_energy_cap = 30;
  double a_list[] = { 31.7 };

  std::string input_directory = "/Users/Violet/data/peptide_design/tetramer_29/fullseq_design_approved_by_Jeff";
  //std::string input_directory = "/Users/Violet/data/peptide_design/2D_lattice/C4459/P422/structure";
  std::string output_directory = "/Users/Violet/data/peptide_design/2D_lattice/C4459/P422/vergil_design";
  std::string pdb_filename = input_directory + "/C4459a_chainA_fixed_adg_NH2cap_for_vergil.pdb";
  //std::string pdb_filename = input_directory + "/C4459_P422_a31.2_mP_type_wCap_W_H_segX0.pdb";

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
  ComputeUnfoldedFreeEnergy unfolded_compute(ComputeUnfoldedFreeEnergy::AMBER_KONO,
                                             vergil->forcefield_parameter_library(), vergil->topology_library(),
                                             vergil->conformer_library());
  unfolded_compute.set_beta(ref_beta);
  unfolded_compute.AddCompute(vergil->compute("dihedral"));
  unfolded_compute.AddCompute(vergil->compute("vanderwaals"));
  unfolded_compute.AddCompute(vergil->compute("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  for (int idx = 0; idx < 1; ++idx) {
    std::string output_fileprefix = output_directory + "/C4459_P422_a" + Log->to_str(a_list[idx]) + "fix_19W_27H";
    // Load scaffold
    std::string scaffold = "protein and backbone or site 1 4 7 8 11 14 15 18 21 22 25 28 29 30";
//    std::string scaffold = "protein and backbone or site 1 to 15 or site 17 to 30";
    InputPDBScaffold input(vergil->domain(), scaffold);
    input.KeepSitesWithMissingAtoms();
    input.Read(pdb_filename);
    Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

// Type each site in the domain
    int num_amino_acids = 18;
    std::string amino_acid_set[] = {"ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS","MET",
      "PHE", "SER", "THR", "TRP", "TYR", "VAL"};

    std::string typelist[] = {"ASP", "ILE", "MET", "ALA", "ILE", "MET", "ALA", "ILE", "MET", "ALA", "ILE", "GLU", "ALA"};
    //std::string typelist_ex[] = { "GLN", "GLU", "ARG", "GLN", "GLU", "TRP", "LYS", "LYS", "LEU", "MET", "ASP", "LYS", "HIS", "ARG", "ASP", "ARG", "NH2" };

    size_t i = 0;
//    size_t j = 0;
    for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
        ++it) {
      size_t heptad_id = it->name() % 7;
      if (heptad_id == 1 || heptad_id == 4 || heptad_id == 0) {
        it->Add(typelist[i]);
        ++i;
//      } else {
//        //must only add the reference type to the fixed sites
//        it->Add(typelist_ex[j]);
//        ++j;
//      }
      } else if (it->name() == 27) {
        it->Add("HIS");
      } else if (it->name() == 19) {
        it->Add("TRP");
      } else if (it->name() == 30) {
        it->Add("NH2");
      } else {
        for (int j = 0; j < num_amino_acids; ++j)
        it->Add(amino_acid_set[j]);
      }
    }

    // Build all conformers
    vergil->BuildDomain();
    vergil->TrimDomain(30.0);

    //Offset the domain(ASU) to the correct position in the unit cell
    const Vector3 offset(a_list[idx] / 2, 0, 0);
    vergil->domain()->TransformBy(Geometry::Transformation(Vector3(0, 0, 0), 0.0, offset));

    //Set up symmetry and build lattice
    vergil->domain()->set_unit_cell_parameters("P422", a_list[idx], a_list[idx], 1000, 90.0, 90.0, 90.0);
    SymmetryGenerator<Domain> symexp("P422", a_list[idx], a_list[idx], 1000, 90.0, 90.0, 90.0);
    symexp.BuildLattice(*vergil->domain(), -1, -1, 0, 1, 1, 0);
    std::vector<Domain>* symmetry_related_elements = symexp.lattice();

    //Trim rotamers that clash with the scaffold
    vergil->TrimDomain(30.0, symmetry_related_elements);

    // Set up energy function
    timer.Stamp();
    Log->print("Setting up energy function...");
    FunctionMeanFieldEnergyLattice energy(vergil->domain(), symmetry_related_elements);
    energy.set_pairwise_energy_cap_values(twobody_energy_cap, twobody_energy_cap);
    energy.AddCompute(vergil->compute("dihedral"));
    energy.AddCompute(vergil->compute("vanderwaals"));
    energy.AddCompute(vergil->compute("electrostatic"));
    energy.ComputeEnergies();
    Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

    // Initialize probabilities
    vergil->domain()->SetUniformConformerProbabilities();

    // Set up the problem and solve the unconstrained problem
    FunctionFreeEnergyWithReference free_energy(&energy, beta, &unfolded_compute);
    vergil->InitProbabilityProblem(free_energy);
    vergil->SolveProbabilityProblem();

    //Output standard files (pdb, psf, seq, csv)
    vergil->StandardOutput(output_fileprefix);

    for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End(); ++it) {
      it->SetAggregatedTypeProbability();
    }
    vergil->RenameSymmetrySegname(symmetry_related_elements);
    OutputPDB output1(output_fileprefix + "_lattice.pdb", vergil->domain());
    output1.WriteLattice(symmetry_related_elements);
    Log->print("Most Probable PDB -- " + output_fileprefix + "_lattice.pdb");

    vergil->domain()->Clear();
  }
}

void Run(Vergil* vergil) {
  std::string input_directory = "/Users/Violet/data/peptide_design/2D_lattice/C4459/P422/structure";
  std::string output_directory = "/Users/Violet/data/peptide_design/2D_lattice/C4459/P422/structure/P422_a31.2_Cu/";
  std::string pdb_filename = input_directory + "/C4459_P422_a31.2_mP_type_wCap_W_H_segX0.pdb";
  //std::string pdb_filename = input_directory + "/C4459_P422_a31.7fix_19W_27H.pdb";

  //Geometry variables
  double chi_incr = 1.0;
  std::vector<double> chi_stdev;
  chi_stdev.push_back(180.0);
  chi_stdev.push_back(180.0);
  std::string lp_atom = "NE2";
  double unit_cell_length = 31.2;
  double angle_tolerance = 5.0;

  //Create a local Timer instance
  Timer timer;

  //Load force field and design topology
  vergil->LoadAMBERUnitedAtom();
  vergil->LoadSavenDesignStandards();
  unsigned int max_dunbrack_rots = 9;
  vergil->LoadDunbrackRotamerLibrary(2002, max_dunbrack_rots);

  std::string scaffold = "protein and backbone or not site 23";
  InputPDBScaffold input(vergil->domain(), scaffold);
  input.KeepSitesWithMissingAtoms();
  input.Read(pdb_filename);
  Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

  std::string typelist[] = { "ASP", "GLN", "GLU", "ILE", "ARG", "GLN", "MET", "ALA", "GLU", "TRP", "ILE", "LYS", "LYS",
      "MET", "ALA", "GLN", "MET", "ILE", "ASP", "LYS", "MET", "ALA", "HIS", "ARG", "ILE", "ASP", "ARG", "GLU", "ALA",
      "NH2" };
//  std::string typelist[] = { "ASP", "GLN", "GLU", "ILE", "ARG", "GLN", "MET", "ALA", "GLU", "GLU", "ILE", "ARG", "LYS",
//      "MET", "ALA", "GLN", "GLN", "ILE", "TRP", "GLU", "MET", "ALA", "LYS", "ARG", "ILE", "ASP", "HIS", "GLU", "ALA", "NH2"};
  size_t i = 0;
  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
      ++it) {
    it->Add(typelist[i]);
    ++i;
  }

  vergil->BuildDomain();

  timer.Stamp();
  Conformer* ligand = &vergil->domain()[0]["X0"]["A"][0][23]["HIS"][0];
  OptimizeLigandGeomety(*ligand, chi_incr, chi_stdev, lp_atom, unit_cell_length, angle_tolerance, output_directory);
  Log->print_tag("RUNTIME", "Optimize ligand geometry time: " + timer.ElapsedToString());

  vergil->domain()->Clear();
}

void OptimizeLigandGeomety(Conformer& ligand, double chi_increment, std::vector<double> chi_stdev, std::string lp_atom,
                           double unit_cell_length, double angle_tolerance, std::string output_dir) {
  if (ligand.parent()->conformer_states() == NULL)
    Log->error(FLERR, "There are no conformer states for the ligand!");
  const std::vector<Linkage<std::string> >& chi_names = ligand.parent()->conformer_states()->degrees_of_freedom();

  std::vector<double> chi_values = GetChiValues(ligand);
  //adjust chi to values between [0, 360]
  for (size_t i = 0; i < chi_values.size(); ++i) {
    chi_values[i] = (chi_values[i] < 0 ? chi_values[i] += 360 : chi_values[i]);
  }

  //loop through the possible chi angles, filling the ligand geometry map
  std::map<std::vector<double>, std::vector<double> > ligand_geometry;
  for(double chi1 = floor(chi_values[0] - chi_stdev[0]); chi1 < ceil(chi_values[0] + chi_stdev[0]); chi1 += chi_increment) {
    for(double chi2 = floor(chi_values[1] - chi_stdev[1]); chi2 < ceil(chi_values[1] + chi_stdev[1]); chi2 += chi_increment) {
      ligand.AdjustLinkage(chi_names[0], chi1);
      ligand.AdjustLinkage(chi_names[1], chi2);

      std::vector<double> chis;
      chis.push_back(chi1);
      chis.push_back(chi2);
      Vector3 v = CalculateLonePairDirection(ligand, lp_atom);

      //Calculate the dot product of v and z axis unit vector
      //cos(q) = v*w / (|v|*|w|)
      double angle_z = RADIANS_TO_DEGREES * acos(v.DotProduct(Vector3(0, 0, 1)) / v.Length());
      //Calculate the dot product of v and the vector from the N to the metal (where the C4 symmetry axis is)
      Atom* N = ligand.Find(lp_atom);
      Vector3 N_metal(unit_cell_length - N->x(), - N->y(), 0);
      double angle_metal = RADIANS_TO_DEGREES * acos(v.DotProduct(N_metal) / (v.Length() * N_metal.Length()));

      if ((angle_z < (90.0 + angle_tolerance) && angle_z > (90.0 - angle_tolerance)) && (angle_metal < angle_tolerance)) {
        ligand_geometry[chis].push_back(angle_z);
        ligand_geometry[chis].push_back(angle_metal);

        //OutputStructure(ligand.parent()->parent()->parent()->parent()->parent()->parent(), chis, lp_atom, output_dir);
        Log->print("Optimized Cu binding PDB with chi1 " + Log->to_str(chi1) + ", chi2 " + Log->to_str(chi2));

        OutputLatticeStructure(ligand.parent()->parent()->parent()->parent()->parent()->parent(), chis, unit_cell_length, lp_atom, output_dir);
      }
    }
  }

  OutputStat(ligand_geometry, lp_atom, output_dir);
}

std::vector<double> GetChiValues(Conformer& ligand) {
  //get degrees of freedom and the chi angles of the ligand
  size_t number_of_chis = ligand.parent()->conformer_states()->NumberDegreesOfFreedom();
  if (ligand.parent()->conformer_states() == NULL)
    Log->error(FLERR, "There are no conformer states for the ligand!");
  const std::vector<Linkage<std::string> >& chi_names = ligand.parent()->conformer_states()->degrees_of_freedom();
  std::vector<double> chi_values;

  for (size_t i = 0; i < number_of_chis; ++i) {
    Atom* a = ligand.Find(chi_names[i][0]);
    Atom* b = ligand.Find(chi_names[i][1]);
    Atom* c = ligand.Find(chi_names[i][2]);
    Atom* d = ligand.Find(chi_names[i][3]);

    if (a == NULL || b == NULL || c == NULL || d == NULL) {
      Log->error(FLERR, "Check if there are missing atoms in the ligand");
    }

    chi_values.push_back(RADIANS_TO_DEGREES * Geometry::Dihedral(*a, *b, *c, *d));
  }

  return chi_values;
}

Vector3 CalculateLonePairDirection(Conformer& ligand, std::string LP_atom) {
  Atom* NE = ligand.Find("NE2");
  Atom* ND = ligand.Find("ND1");
  Atom* CG = ligand.Find("CG");
  Atom* CD = ligand.Find("CD2");
  Atom* CE = ligand.Find("CE1");
  if (NE == NULL || ND == NULL || CG == NULL || CD == NULL || CE == NULL)
    Log->error(FLERR, "Check the atom names in the ligand");

  Vector3 v;
  if (LP_atom == "NE2") {
    //define theta
    Vector3 r1(*CE, *NE);
    Vector3 r2(*CD, *NE);
    double theta = acos(r1.DotProduct(r2) / (r1.Length() * r2.Length()));

    //define gamma, gamma = 0.5 * theta
    Point middle;
    double a = 0.49969;
    middle.set_x(a * ND->x() + (1 - a) * CG->x());
    middle.set_y(a * ND->y() + (1 - a) * CG->y());
    middle.set_z(a * ND->z() + (1 - a) * CG->z());
    v = Vector3(middle, *NE);
    double gamma = acos(v.DotProduct(r2) / (v.Length() * r2.Length()));
  } else if (LP_atom == "ND1") {
    //define theta
    Vector3 r1(*CE, *ND);
    Vector3 r2(*CG, *ND);
    double theta = acos(r1.DotProduct(r2) / (r1.Length() * r2.Length()));

    //define gamma, gamma = 0.5 * theta
    Point middle;
    double a = 0.495;
    middle.set_x(a * NE->x() + (1 - a) * CD->x());
    middle.set_y(a * NE->y() + (1 - a) * CD->y());
    middle.set_z(a * NE->z() + (1 - a) * CD->z());
    v = Vector3(middle, *ND);
    double gamma = acos(v.DotProduct(r2) / (v.Length() * r2.Length()));
  }

  return v;
}

void OutputStat(std::map<std::vector<double>, std::vector<double> > ligand_geometry, std::string lp_atom, std::string output_directory) {
  std::string filename = output_directory + "/CuP422_his" + lp_atom + "_geometry_optimization.csv";
  FILE* outfile = fopen(filename.c_str(), "w");
  fprintf(outfile, "chi1, chi2, His_orienation_z, His_orienation_x\n");
  fclose(outfile);

  outfile = fopen(filename.c_str(), "a");
  for (std::map<std::vector<double>, std::vector<double> >::iterator it = ligand_geometry.begin(); it != ligand_geometry.end();
      ++it) {
    double chi1 = (it->first)[0] > 180 ? (it->first)[0] - 360 : (it->first)[0];
    double chi2 = (it->first)[1] > 180 ? (it->first)[1] - 360 : (it->first)[1];
    fprintf(outfile, "%10.5f, %10.5f, %10.5f, %10.5f\n", chi1, chi2, it->second[0], it->second[1]);
  }
  fclose(outfile);
}

void OutputStructure(Domain* domain, std::vector<double> chis, std::string lp_atom, std::string output_directory) {
  double chi1 = chis[0] > 180 ? chis[0] - 360 : chis[0];
  double chi2 = chis[1] > 180 ? chis[1] - 360 : chis[1];
  std::string filename = output_directory + "/CuP422_" + lp_atom + "_" + Log->to_str(chi1) + "_" + Log->to_str(chi2);
  OutputPDB output(filename + ".pdb", domain);
  output.Write();
}

void OutputLatticeStructure(Domain* domain, std::vector<double> chis, double unit_cell_length, std::string lp_atom, std::string output_directory) {
  double chi1 = chis[0] > 180 ? chis[0] - 360 : chis[0];
  double chi2 = chis[1] > 180 ? chis[1] - 360 : chis[1];

  //Set up symmetry and build lattice
  domain->set_unit_cell_parameters("P422", unit_cell_length, unit_cell_length, 1000, 90.0, 90.0, 90.0);
  SymmetryGenerator<Domain> symexp("P422", unit_cell_length, unit_cell_length, 1000, 90.0, 90.0, 90.0);
  symexp.BuildLattice(*domain, -1, -1, 0, 1, 1, 0);
  std::vector<Domain>* symmetry_related_elements = symexp.lattice();

  //Rename all motifs to be X- according to the asymmetric unit name
  int unit_count = 1;
  for (std::vector<Domain>::iterator it = symmetry_related_elements->begin(); it != symmetry_related_elements->end();
      ++it) {
    for (Domain::MotifIterator motif_it = it->MotifIterator_Begin(); motif_it != it->MotifIterator_End(); ++motif_it) {
      motif_it->set_name("X" + Log->to_str(unit_count));
    }
    unit_count++;
  }

  std::string filename = output_directory + "/CuP422_" + lp_atom + "_"+ Log->to_str(chi1) + "_" + Log->to_str(chi2) + "_lattice";
  OutputPDB output(filename + ".pdb", domain);
  output.WriteLattice(symmetry_related_elements);
}
