/*
 * NSA_example.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Violet
 */

#include "mpi.h"
#include "vergil_config.h"
#include "vergil.h"
#include "input_topology_charmm.h"
#include "domain.h"
#include "section.h"
#include "site.h"
#include "surface_area_calculator.h"
#include "conformer_states_explicit.h"
#include "library_conformer_states.h"

#include "string_edit.h"
#include "string"

using namespace vergil_ns;

void NSA_calculation(Vergil* vergil);

int main(int argc, char **argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Create a VERGIL instance, pass arguments;
	Vergil* vergil = new Vergil();

	// Run the routine
	NSA_calculation(vergil);

	// Delete the VERGIL instance and clean up
	delete vergil;

	// Finalize MPI
	MPI_Finalize();

	return 0;
}

void NSA_calculation(Vergil* vergil) {
	//Force field
	vergil->LoadCHARMM(22);
	unsigned int max_rot = 2;
        vergil->LoadDunbrackRotamerLibrary(2002, max_rot);

	//load special topology files
	InputTopologyCharmm input_top(vergil->forcefield_parameter_library(),
			vergil->topology_library());
	input_top.Read(
			"/Users/Violet/data/Charmm/top_all27_prot_lipid_d-pro_nap.inp");
	Log->print_tag_color("INPUT", "green", "loaded special topology file containing NAP and D-PRO.");

	//Open the pdb list
	std::string output_directory =
			"/Users/Violet/data/beta_hairpin/max1/8mer/max1_anti_8mer/tip4p/NSC_SASA_analysis/";
	std::string list_file =
			output_directory + "100_random_pdbs_cross_angle_167.txt";
	FILE* pdblist_file = fopen(list_file.c_str(), "r");
	if (!pdblist_file) {
		Log->error(FLERR, "Cannot open the PDB list file: " + list_file);
	}

	//Grab the pdb file names
	std::vector<std::string> pdbs;
	const static int BLEN = 500;
	char buffer[BLEN];
	while (fgets(buffer, BLEN, pdblist_file)) {
		std::vector<std::string> words = StringEdit::Tokenize(std::string(buffer), "\n");
		for (std::vector<std::string>::iterator it = words.begin();
				it != words.end(); ++it) {
			*it = StringEdit::Trim(*it);
		}
		pdbs.push_back(words[0]);
	}

	//Output the header
//	std::string SASA_file = output_directory
//			+ "max1_pdbs_cross_angle_167_CDEF_NSC_SASA.csv";
//	FILE* SASA_outfile = fopen(SASA_file.c_str(), "w");
//	fprintf(SASA_outfile,
//			"pdb, total_SASA (A^2), hydrophobic_SASA (A^2), SASA_percent \n");
//	fclose(SASA_outfile);

	for (std::vector<std::string>::const_iterator it = pdbs.begin();
			it != pdbs.end(); ++it) {

		//Cap the Cterm as CT2
		vergil->topology_library()->SubstituteDefaultCap("CTER", "CT2");

		//Input pdb
		//std::string output_directory = "/Users/Violet/data/beta_hairpin/lnk1/8mer/lnk1_8mer_anti/tip4p/rep3/";
		vergil->LoadScaffoldFromPDB(
				output_directory + "pdbs_cross_angle_167/" + *it,
				"name N CA C O");
		/*	vergil->LoadScaffoldFromPDB(
		 output_directory + "lnk1_anti_8mer_design_fine_tip4p_rep3_70ns.pdb",
		 "name N CA C O");*/

		std::string MAX1_sequence[20] = { "VAL", "LYS", "VAL", "LYS", "VAL",
				"LYS", "VAL", "LYS", "VAL", "DPR", "PRO", "THR", "LYS", "VAL",
				"LYS", "VAL", "LYS", "VAL", "LYS", "VAL" };
		//std::string LNK1_sequence[20] = { "NAP", "LYS", "NAP", "LYS", "ALA",
		//		"LYS", "ALA", "LYS", "VAL", "DPR", "PRO", "THR", "LYS", "ALA",
		//		"LYS", "ALA", "LYS", "NAP", "LYS", "NAP" };

		//Type the protein
		for (Domain::SectionIterator it =
				vergil->domain()->SectionIterator_Begin();
				it != vergil->domain()->SectionIterator_End(); ++it) {
			unsigned int resid = 0;
			for (Section::SiteIterator jt = it->SiteIterator_Begin();
					jt != it->SiteIterator_End(); ++jt) {
				jt->Add(MAX1_sequence[resid]);
				//jt->Add(LNK1_sequence[resid]);
				resid++;
			}
		}

		//Create conformer states for DPR
		ConformerStatesExplicit* DPR_States = new ConformerStatesExplicit(
				"DPR");
		Linkage<std::string> DPR_dihedral1("N", "CA", "CB", "CG");
		Linkage<std::string> DPR_dihedral2("CA", "CB", "CG", "CD");
		DPR_States->AddDegreeOfFreedom(DPR_dihedral1);
		DPR_States->AddDegreeOfFreedom(DPR_dihedral2);
		std::vector<double> DPR_values;
		DPR_values.push_back(31.6100);
		DPR_values.push_back(-34.5900);
		DPR_States->AddState(DPR_values);
		vergil->conformer_library()->AddConformerStates("DPR", DPR_States);

		// Set root atoms for each monomer type
		vergil->topology_library()->FindMonomer("DPR")->set_root("N");

		//Create conformer states for NAP
		ConformerStatesExplicit* NAP_States = new ConformerStatesExplicit(
				"NAP");
		Linkage<std::string> NAP_dihedral1("N", "CA", "CB", "CG");
		Linkage<std::string> NAP_dihedral2("CA", "CB", "CG", "CD1");
		NAP_States->AddDegreeOfFreedom(NAP_dihedral1);
		NAP_States->AddDegreeOfFreedom(NAP_dihedral2);
		std::vector<double> NAP_values;
		NAP_values.push_back(180.0000);
		NAP_values.push_back(90.0000);
		NAP_States->AddState(NAP_values);
		vergil->conformer_library()->AddConformerStates("NAP", NAP_States);

		// Set root atoms for each monomer type
		vergil->topology_library()->FindMonomer("NAP")->set_root("N");

		vergil->BuildDomain();

		/*	OutputPDB output1(output_directory + "test.pdb", vergil->domain());
		 output1.Write();

		 OutputPSF output2(output_directory + "test.psf", vergil->domain());
		 output2.Write();*/

		//Some variables set up;
		//Use 50% surface area threshold to distinguish exposed site
		//Use 20% surface area threshold to distinguish buried site
		//comparing to the Ala-X-Ala model
		std::vector<Site*> buried_sites;
		std::vector<Site*> exposed_sites;
		std::string buried_or_exposed;
		//double hydrophobic_SASA = 0.0;
		//double total_SASA = 0.0;
		double exposed_threshold = 0.5;
		double buried_threshold = 0.2;

		SurfaceAreaCalculator NSC;

		//Switch the cap back for the the Ala-X-Ala model
		vergil->topology_library()->SubstituteDefaultCap("CT2", "CTER");

//	std::string SASA_file = output_directory
//			+ "max1_anti_8mer_tip4p_100ps_CDEF_NSC_SASA.csv";
//	FILE* SASA_outfile = fopen(SASA_file.c_str(), "w");
//	fprintf(SASA_outfile, "Chain, Site, SASA (A^2), buried or exposed \n");
//	fclose(SASA_outfile);

		//For LNK1 and Linear MAX1 (cross-angle >= 170), only calculates the middle section (chain C D E F)
		for (Domain::SectionIterator it =
				vergil->domain()->SectionIterator_Begin();
				it != vergil->domain()->SectionIterator_End(); ++it) {
			if (it->name() == "C" || it->name() == "D" || it->name() == "E"
					|| it->name() == "F") {
				for (Section::SiteIterator jt = it->SiteIterator_Begin();
						jt != it->SiteIterator_End(); ++jt) {
					//Calculate SASA ratio
					double SASA_ratio = NSC.RatioSASA(vergil, *jt);

					//Sum up the SASA of residues on the hydrophobic face
//					if (jt->reference_type()->name() == "VAL"
//							|| jt->reference_type()->name() == "THR"
//							|| jt->reference_type()->name() == "PRO"
//							|| jt->reference_type()->name() == "DPR")
//						hydrophobic_SASA += area;

					/*				//Sum up the SASA of residues on the hydrophobic face
					 if (jt->reference_type()->name() == "VAL"
					 || jt->reference_type()->name() == "THR"
					 || jt->reference_type()->name() == "PRO"
					 || jt->reference_type()->name() == "DPR"
					 || jt->reference_type()->name() == "NAP"
					 || jt->reference_type()->name() == "ALA")
					 hydrophobic_SASA += area;*/

					//Sum up the SASA for the entire protein
					//total_SASA += area;

					if (SASA_ratio < buried_threshold) {
						buried_sites.push_back(&*jt);
						buried_or_exposed = "buried";
					} else if (SASA_ratio > exposed_threshold) {
						exposed_sites.push_back(&*jt);
						buried_or_exposed = "exposed";
					} else {
						buried_or_exposed = "";
					}

					//Using NSC, buried if SASA less than 20% of that of the same residue in extended Ala-X-Als
//				SASA_outfile = fopen(SASA_file.c_str(), "a");
//				fprintf(SASA_outfile, "%s, %5d, %12.5f, %s \n",
//						jt->parent()->parent()->parent()->name().c_str(),
//						jt->name(), area, buried_or_exposed.c_str());
//				fclose(SASA_outfile);
				}
			}
		}

//	SASA_outfile = fopen(SASA_file.c_str(), "a");
//	fprintf(SASA_outfile, "Hydrophobic SASA is: %f \n", hydrophobic_SASA);
//	fprintf(SASA_outfile, "Total SASA is: %f \n", total_SASA);
//	fprintf(SASA_outfile, "Hydrophobic SASA percentage is: %f \n",
//			hydrophobic_SASA / total_SASA);
//	fclose(SASA_outfile);

//		SASA_outfile = fopen(SASA_file.c_str(), "a");
//		fprintf(SASA_outfile, "%s, %12.5f, %12.5f, %12.5f \n", it->c_str(),
//				total_SASA, hydrophobic_SASA, hydrophobic_SASA / total_SASA);
//		fclose(SASA_outfile);

		//populate a string of buried site ids
		std::string buried_sites_str;
		if (!buried_sites.empty()) {
			for (size_t i = 0; i < buried_sites.size(); ++i) {
				buried_sites_str.append(Log->to_str(buried_sites[i]->name()));
				buried_sites_str.append(" ");
			}
			Log->print("Buried sites: " + buried_sites_str);
		} else {
			Log->print("No buried sites in this pdb!");
			buried_sites_str = "";
		}

		//populate a string of exposed site ids
		std::string exposed_sites_str;
		if (!exposed_sites.empty()) {
			for (size_t i = 0; i < exposed_sites.size(); ++i) {
				exposed_sites_str.append(Log->to_str(exposed_sites[i]->name()));
				exposed_sites_str.append(" ");
			}
			Log->print("Exposed sites: " + exposed_sites_str);
		} else {
			Log->print("No exposed sites in this pdb");
			exposed_sites_str = "";
		}

		//Clear the domain
		vergil->domain()->Clear();
	}
}

