/*
 * crosslinker.cpp
 *
 *  Created on: May 29, 2015
 *      Author: Violet
 */

#include "mpi.h"

#include "library_topology.h"
#include "input_pdb_topology.h"
#include "input_topology_charmm.h"
#include "input_pdb_scaffold.h"
#include "domain.h"
#include "vergil.h"

using namespace vergil_ns;

/**
 * @brief non-design routine
 * @param vergil
 */
void Run(Vergil* vergil);

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

void Run(Vergil* vergil) {
  std::string input_directory = "/Users/Violet/data/peptide_design/tetramer_29/fullseq_design_approved_by_Jeff";
  std::string pdb_filename = input_directory + "/C4459a_chainA_fixed_adg_NH2cap_for_vergil.pdb";
  std::string uaa_filename = "/Users/Violet/data/cross_linkers/geometry/PDBeChem/uaa_template.pdb";
  std::string uaa_top = "/Users/Violet/data/cross_linkers/toppar/top_uaa.inp";
  std::string output_fileprefix = "/Users/Violet/data/cross_linkers/test";

  //Load force field and design topology
  vergil->LoadCHARMM(22);
  vergil->LoadSavenDesignStandards();
  InputTopologyCharmm input_uaa_top(vergil->forcefield_parameter_library(), vergil->topology_library());
  input_uaa_top.Read(uaa_top);
  Log->print_tag("INPUT", "uaa topology -- " + uaa_top);
  InputPDBTopology inputuaa(vergil->topology_library(), "N CA C O");
  inputuaa.Read(uaa_filename);
  Log->print_tag("INPUT", "uaa PBD template -- " + uaa_filename);

  // Load scaffold
  std::string scaffold = "protein and backbone or site 1 4 7 8 11 14 15 18 21 22 25 28 29";
  InputPDBScaffold input(vergil->domain(), scaffold);
  input.KeepSitesWithMissingAtoms();
  input.Read(pdb_filename);
  Log->print_tag("INPUT", "PDB scaffold -- " + pdb_filename);

  //type the sites
  for (Domain::SiteIterator it = vergil->domain()->SiteIterator_Begin(); it != vergil->domain()->SiteIterator_End();
          ++it) {
    if(it->name() == 2)
      it->Add("2AG");
    else if(it->name() == 9)
      it->Add("HCS");
    else if(it->name() == 16)
      it->Add("LVG");
    else
      it->Add("GLY");
  }

  vergil->BuildDomain();

  //Output standard files (pdb, psf, seq, csv)
  vergil->StandardOutput(output_fileprefix);

  vergil->domain()->Clear();
}



