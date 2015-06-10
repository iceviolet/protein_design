/*
 * crosslinker.cpp
 *
 *  Created on: May 29, 2015
 *      Author: Violet
 */

#include "mpi.h"

#include "library_topology.h"
#include "input_pdb_topology.h"
#include "input_parameter_charmm.h"
#include "input_topology_charmm.h"
#include "input_pdb_scaffold.h"
#include "input_dunbrack.h"
#include "domain.h"
#include "potential_unfolded_freeenergy.h"
#include "function_meanfieldenergy_lattice.h"
#include "function_freeenergy_with_reference.h"
#include "output_pdb.h"
#include "vergil.h"
#include "timer.h"
#include "geometry.h"

using namespace vergil_ns;

/**
 * @brief non-design test routine
 * @param vergil
 */
void Test(Vergil* vergil);

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
  MPI_Init(&argc, &argv);

  // Create a VERGIL instance, pass arguments;
  Vergil *vergil = new Vergil();

  // Run the routine
  Design(vergil);
  //Test(vergil);

  // Delete the VERGIL instance and clean up
  delete vergil;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

void Design(Vergil *vergil) {
  std::string input_directory = "/Users/Violet/data/peptide_design/tetramer_29/fullseq_design_approved_by_Jeff";
  std::string pdb_filename = input_directory + "/C4459a_chainA_fixed_adg_NH2cap_for_vergil.pdb";
  std::string uaa_filename = "/Users/Violet/data/cross_linkers/geometry/PDBeChem/uaa_template.pdb";
  std::string uaa_param = "/Users/Violet/data/cross_linkers/toppar/par_uaa.inp";
  std::string uaa_top = "/Users/Violet/data/cross_linkers/toppar/top_uaa.inp";
  std::string uaa_rot_dir = "/Users/Violet/data/cross_linkers/parameterization_tools/SwissSideChain/rotamer_lib";
  std::string output_fileprefix = "/Users/Violet/data/cross_linkers/test";

  Timer timer;

  //Load force field and design topology
  vergil->LoadCHARMM(22);
  InputParameterCharmm input_uaa_param(vergil->forcefield_parameter_library());
  input_uaa_param.Read(uaa_param);
  Log->print_tag("INPUT", "uaa parameter -- " + uaa_param);
  vergil->LoadSavenDesignStandards();
  InputTopologyCharmm input_uaa_top(vergil->forcefield_parameter_library(), vergil->topology_library());
  input_uaa_top.Read(uaa_top);
  Log->print_tag("INPUT", "uaa topology -- " + uaa_top);
  InputPDBTopology inputuaa(vergil->topology_library(), "N CA C O");
  inputuaa.Read(uaa_filename);
  Log->print_tag("INPUT", "uaa PBD template -- " + uaa_filename);

  //Load rotamer libraries
  unsigned int max_dunbrack_rots = 2;
  vergil->LoadDunbrackRotamerLibrary(2010, max_dunbrack_rots);
  InputDunbrack input_lvg_rot(vergil->conformer_library(), max_dunbrack_rots);
  input_lvg_rot.Read(uaa_rot_dir + "/LVG_bbdep_Gfeller.lib");
  Log->print_tag("INPUT", "uaa LVG rotatmer library");

  // Calculate Reference Energy
  timer.Stamp();
  Log->print("Computing reference energies....");
  PotentialUnfoldedFreeEnergy unfolded_compute(PotentialUnfoldedFreeEnergy::CHARMM22,
                                             vergil->forcefield_parameter_library(), vergil->topology_library(),
                                             vergil->conformer_library());
  unfolded_compute.set_beta(ROOM_TEMPERATURE_BETA);
  unfolded_compute.AddPotential(vergil->potential("dihedral"));
  unfolded_compute.AddPotential(vergil->potential("vanderwaals"));
  unfolded_compute.AddPotential(vergil->potential("electrostatic"));
  unfolded_compute.GenerateUnfoldedFreeEnergies();
  Log->print_tag("RUNTIME", "Reference energy calculation: " + timer.ElapsedToString());

  // Load scaffold
  std::string scaffold = "backbone or site 1 4 7 8 11 14 15 18 21 22 25 28 29";
  InputPDBScaffold input(vergil->domain(), scaffold);
  input.KeepSitesWithMissingAtoms();
  input.Read(pdb_filename);
  Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

  // Type each site in the domain
  // todo: check fixed site 4, why cannot type it to Leu
  int num_amino_acids = 20;
  std::string amino_acid_set[] = { "ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
      "PHE", "SER", "THR", "TRP", "TYR", "VAL", "LVG", "CYS" };
  std::string typelist[] = { "ASP", "LEU", "MET", "ALA", "ILE", "MET", "ALA", "ILE", "MET", "ALA", "ILE", "GLU", "ALA" };
  size_t i = 0;
  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
      ++it) {
    size_t heptad_id = it->name() % 7;
    if (heptad_id == 1 || heptad_id == 4 || heptad_id == 0) {
      it->Add(typelist[i]);
      if (it->name() == 29) {
        // Cap the C term
        TopologyMonomer* monomer = vergil->topology_library()->FindMonomer("ALA");
        monomer->set_last_default_patch("CT2");
      }
      ++i;
    } else {
      for (int j = 0; j < num_amino_acids; ++j)
        it->Add(amino_acid_set[j]);
    }
  }

  // Build all conformers
  vergil->BuildDomain();
  vergil->TrimDomain(30.0);

  //Offset the domain(ASU) to the correct position in the unit cell
  const Vector3 offset(45/2, 0, 0);
  vergil->domain()->TransformBy(Geometry::Transformation(Vector3(0, 0, 0), 0.0, offset));

  //Set up symmetry and build lattice
  vergil->domain()->set_unit_cell_parameters("P622", 45, 45, 1000, 90.0, 90.0, 120.0);
  SymmetryGenerator<Domain> symexp("P622", 45, 45, 1000, 90.0, 90.0, 120.0);
  symexp.BuildLattice(*vergil->domain(), -1, -1, 0, 1, 1, 0);
  std::vector<Domain>* symmetry_related_elements = symexp.lattice();

  //Trim rotamers that clash with the scaffold
  vergil->TrimDomain(30.0, symmetry_related_elements);

  // Set up energy function
  timer.Stamp();
  Log->print("Setting up energy function...");
  FunctionMeanFieldEnergyLattice energy(vergil->domain(), symmetry_related_elements);
  energy.set_pairwise_energy_cap_values(30, 30);
  energy.AddPotential(vergil->potential("dihedral"));
  energy.AddPotential(vergil->potential("vanderwaals"));
  energy.AddPotential(vergil->potential("electrostatic"));
  energy.ComputeEnergies();
  Log->print_tag("RUNTIME", "Energy Matrix Fill Time: " + timer.ElapsedToString());

  // Initialize probabilities
  vergil->domain()->SetUniformConformerProbabilities();

  // Set up the problem
  FunctionFreeEnergyWithReference free_energy(&energy, 0.5, &unfolded_compute);
  vergil->InitProbabilityProblem(free_energy);

  //Todo: Set up composition constraints from the problem class

  //solve the unconstrained problem
  vergil->SolveProbabilityProblem();

  //Output standard files (pdb, psf, seq, csv)
  vergil->StandardOutput(output_fileprefix);

  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
      ++it) {
    it->SetAggregatedTypeProbability();
  }
  vergil->RenameSymmetrySegname(symmetry_related_elements);
  OutputPDB output1(output_fileprefix + "_lattice.pdb", vergil->domain());
  output1.WriteLattice(symmetry_related_elements);
  Log->print("Most Probable PDB -- " + output_fileprefix + "_lattice.pdb");

  vergil->domain()->Clear();
}

void Test(Vergil* vergil) {
  std::string input_directory = "/Users/Violet/data/peptide_design/tetramer_29/fullseq_design_approved_by_Jeff";
  std::string pdb_filename = input_directory + "/C4459a_chainA_fixed_adg_NH2cap_for_vergil.pdb";
  std::string uaa_filename = "/Users/Violet/data/cross_linkers/geometry/PDBeChem/uaa_template.pdb";
  std::string uaa_param = "/Users/Violet/data/cross_linkers/toppar/par_uaa.inp";
  std::string uaa_top = "/Users/Violet/data/cross_linkers/toppar/top_uaa.inp";
  std::string output_fileprefix = "/Users/Violet/data/cross_linkers/test";

  //Load force field and design topology
  vergil->LoadCHARMM(22);
  InputParameterCharmm input_uaa_param(vergil->forcefield_parameter_library());
  input_uaa_param.Read(uaa_param);
  Log->print_tag("INPUT", "uaa parameter -- " + uaa_param);
  vergil->LoadSavenDesignStandards();
  InputTopologyCharmm input_uaa_top(vergil->forcefield_parameter_library(), vergil->topology_library());
  input_uaa_top.Read(uaa_top);
  Log->print_tag("INPUT", "uaa topology -- " + uaa_top);
  InputPDBTopology inputuaa(vergil->topology_library(), "N CA C O");
  inputuaa.Read(uaa_filename);
  Log->print_tag("INPUT", "uaa PBD template -- " + uaa_filename);

  // Load scaffold
  std::string scaffold = "protein and backbone or site 1 4 7 8 11 14 15 18 21 22 25 28 29";
  //std::string scaffold = "all";
  InputPDBScaffold input(vergil->domain(), scaffold);
  input.KeepSitesWithMissingAtoms();
  input.Read(pdb_filename);
  Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

  //type the sites
  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
          ++it) {
    if(it->name() == 2)
      it->Add("CYS");
    else if(it->name() == 5)
      it->Add("HCS");
    else if(it->name() == 6)
      it->Add("LVG");
    else if(it->name() == 9)
      it->Add("2AG");
    else if(it->name() == 10)
      it->Add("AZD");
    else if(it->name() == 12)
      it->Add("AZH");
    else if(it->name() == 13)
      it->Add("NOI");
    else if(it->name() == 16)
      it->Add("NOT");
    else if(it->name() == 19)
      it->Add("LPG");
    else if(it->name() == 23)
      it->Add("GYA");
    else if(it->name() == 27)
      it->Add("UNK");
    else
      it->Add("GLY");
  }

  //still approach to connect on the site level, but need to figure our why uaas won't build
  //answer is need to input all, cause uaas are not included in protein
  //todo: need to solve how to add the additional hydrogen
  //Site& cys = (vergil->domain())[0]["X0"]["A"][0][13];
  //Site& lvg = (vergil->domain())[0]["X0"]["A"][0][16];
  //cys.ConnectNextSiteWithIterface(lvg, "UNK1");

  vergil->BuildDomain();

  //Output standard files (pdb, psf, seq, csv)
  vergil->StandardOutput(output_fileprefix);

  vergil->domain()->Clear();
}



