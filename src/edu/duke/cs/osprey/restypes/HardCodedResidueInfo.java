/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 *
 * @author mhall44
 */
public class HardCodedResidueInfo {
	/*
	 * This version of OSPREY aims to make residue-type-dependent operations as
	 * generic as possible, appropriate for D- and L-amino acids (natural or
	 * not) as well as non-amino acid residues However, we sometimes require
	 * information that is not in the template files The idea of this class is
	 * to get that stuff all in one place to facilitate adding support for more
	 * residue types in the future Also we may want to transfer information from
	 * here into some kind of template files in the future
	 */

	public static String[] possibleAABBAtoms = new String[] { "N", "H", "CA", "C", "O", "OXT", "H1", "H2", "H3" };
	// BB atoms, which should stay still in mutations and should be moved in
	// perturbations.
	// We'll move HA with the sidechain, so it's not included here.

	public static Set<String> possibleAABBAtomsLookup; // DZ: for O(1) lookup
	public static Set<String> possibleNABBAtomsLookup; // DZ: for O(1) lookup

	public static String[] possibleNABBAtoms = new String[] { "P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'",
			"H4'", "C3'", "H3'", "O3'", "HO3'", "O2'", "HO2'", "C2'", "H2'", "C1'", "O4'", "HO5'" };
	// Nucleic Acid BB atoms. Note that O2' is unique to RNA
	// H1' is moved with the "sidechain" in the same way as HA in amino acids.

	public static LinkedHashMap<String, String> three2one = null;
	public static LinkedHashMap<String, String> one2three = null;// reverse
																	// lookup

	static {
		initThree2One();

		one2three = new LinkedHashMap<String, String>();
		for (String threeLet : three2one.keySet())
			one2three.put(three2one.get(threeLet), threeLet);

		possibleAABBAtomsLookup = new HashSet<String>();
		for (String name : possibleAABBAtoms) {
			possibleAABBAtomsLookup.add(name);
		}

		possibleNABBAtomsLookup = new HashSet<String>();
		for (String name : possibleNABBAtoms) {
			possibleNABBAtomsLookup.add(name);
		}
	}

	public static void initThree2One() {
		three2one = new LinkedHashMap<String, String>();
		three2one.put("ALA", "A");
		three2one.put("CYS", "C");
		three2one.put("ASP", "D");
		three2one.put("GLU", "E");
		three2one.put("PHE", "F");
		three2one.put("GLY", "G");
		three2one.put("HIS", "H");
		three2one.put("ILE", "I");
		three2one.put("LYS", "K");
		three2one.put("LEU", "L");
		three2one.put("MET", "M");
		three2one.put("ASN", "N");
		three2one.put("PRO", "P");
		three2one.put("GLN", "Q");
		three2one.put("ARG", "R");
		three2one.put("SER", "S");
		three2one.put("THR", "T");
		three2one.put("VAL", "V");
		three2one.put("TRP", "W");
		three2one.put("TYR", "Y");
	}

	public static String getOneLet(String aa3Name) {
		String res = three2one.get(aa3Name);
		if (res == null)
			res = "X";
		return res;
	}

	// Here's some stuff we need to mutate amino acid protein residues
	public static boolean canMutateTo(ResidueTemplate templ) {
		// do we currently support mutations to the given amino-acid type?
		if (templ.templateRes.coords == null) {
			return false;
		} else {
			return hasAminoAcidBB(templ.templateRes) || hasNucleicAcidBB(templ.templateRes);
			// can currently mutate to any amino acid (D or L, naturally
			// occurring sidechain or not) whose sidechain attaches only
			// to CA and for which we have template coords
			// DZ: now supports nucleic acids as well
		}
	}

	public static ArrayList<String> listBBAtomsForMut(ResidueTemplate newTemplate, ResidueTemplate oldTemplate) {
		// list the backbone atom names shared between the residue templates
		// (in the sense of atoms that don't move at all when we mutate--
		// basically backbone atoms that are present both before and after the
		// mutation)

		if (!canMutateTo(newTemplate)) {
			throw new UnsupportedOperationException("ERROR: Can't currently mutate to " + newTemplate.name);
		}

		ArrayList<String> ans = new ArrayList<>();

		if (hasAminoAcidBB(newTemplate.templateRes) && hasAminoAcidBB(oldTemplate.templateRes)) {
			for (String atName : possibleAABBAtoms) {
				if (newTemplate.templateRes.getAtomIndexByName(atName) != -1) {
					// atom name appears in first template
					if (oldTemplate.templateRes.getAtomIndexByName(atName) != -1) {
						ans.add(atName);
					}
				}
			}
		} else if (hasNucleicAcidBB(newTemplate.templateRes) && hasNucleicAcidBB(oldTemplate.templateRes)) {
			// DZ: NA case
			for (String atName : possibleNABBAtoms) {
				if (newTemplate.templateRes.getAtomIndexByName(atName) != -1) {
					// atom name appears in first template
					if (oldTemplate.templateRes.getAtomIndexByName(atName) != -1) {
						ans.add(atName);
					}
				}
			}
		} else {
			throw new UnsupportedOperationException(
					"ERROR: Can't currently mutate from " + oldTemplate.name + " to " + newTemplate.name);
			// DZ: do not allow mutation between amino acids and nucleic acids
		}

		return ans;
	}
	
	// DZ: Get the key coordinates for a residue with a backbone.
	// These coordinates are used to pivot the sidechain during sidechain idealization  
	// They are analogous to the nitrogen, alpha carbon, and beta carbon in amino acids.
	public static double[][] getKeyCoords(Residue res) {
		double[] NCoords, CACoords, CBCoords, CCoords;
		if (hasNucleicAcidBB(res)) {
			// DZ: analogous positions in nucleic acids
			NCoords = res.getCoordsByAtomName("O4'");
			CACoords = res.getCoordsByAtomName("C1'");
			CBCoords = findPivotCoord(res);
			CCoords = res.getCoordsByAtomName("C2'");
		} else {
			NCoords = res.getCoordsByAtomName("N");
			CACoords = res.getCoordsByAtomName("CA");
			CBCoords = res.getCoordsByAtomName("CB");
			CCoords = res.getCoordsByAtomName("C");
		}
		return new double[][] { NCoords, CACoords, CBCoords, CCoords };
	}

	public static double[] findPivotCoord(Residue res) {
	    return HardCodedResidueInfo.isPyrimidine(res.template.name)
                ? res.getCoordsByAtomName("N1" )
                : HardCodedResidueInfo.isPurine(res.template.name)
                    ? res.getCoordsByAtomName("N9")
                    : null;
	}

	public static boolean isPyrimidine(String resName) {
		return resName.equalsIgnoreCase("ru3")
				|| resName.equalsIgnoreCase("ru2")
				|| resName.equalsIgnoreCase("rt3")
				|| resName.equalsIgnoreCase("rt2")
				|| resName.equalsIgnoreCase("rc3")
				|| resName.equalsIgnoreCase("rc2");
	}

	public static boolean isPurine(String resName) {
		return resName.equalsIgnoreCase("ra3")
				|| resName.equalsIgnoreCase("ra2")
				|| resName.equalsIgnoreCase("rg3")
				|| resName.equalsIgnoreCase("rg2");
	}

	public static boolean isNAPivot(Atom at) {
		// Is the atom next to the C1' carbon? Note that this method assumes
		// that the atom is not part of the NA backbone. This should return true
		// for N1 for pyrimidines and N9 for purines
		if (at.name.equalsIgnoreCase("H1'")
                || at.name.equalsIgnoreCase("C2'")
                || at.name.equalsIgnoreCase("O4'")) {
			return false;
		}
		for (Atom neighbors : at.bonds) {
			if (neighbors.name.equalsIgnoreCase("C1'")) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Sets the pucker of an NA residue whose pucker is not specified
	 * @param res
	 */
	public static void setPucker(Residue res) {
		if (hasNucleicAcidBB(res)) { // only adjust pucker for nucleic acids.
			double[] C5 = res.getCoordsByAtomName("C5'");
			double[] C4 = res.getCoordsByAtomName("C4'");
			double[] C3 = res.getCoordsByAtomName("C3'");
			double[] O3 = res.getCoordsByAtomName("O3'");
			double delta = Protractor.measureDihedral(new double[][] { C5, C4, C3, O3 });
			// The following values came from Swati's thesis, which can be found here:
            // /usr/project/dlab/Users/swati/Thesis/Document/dissertation_0420_2.pdf
			if (delta >= 60 && delta <= 105) { // C3' endo
				res.fullName = res.fullName
						.replace("  U", "RU3")
						.replace("  G", "RG3")
						.replace("  C", "RC3")
						.replace("  A", "RA3");
			} else if (delta >= 125 && delta <= 165) { // C2' endo
				res.fullName = res.fullName
						.replace("  U", "RU2")
						.replace("  G", "RG2")
						.replace("  C", "RC2")
						.replace("  A", "RA2");
			}
		}
	}

	public static int[][] findMutAlignmentAtoms(ResidueTemplate template1, ResidueTemplate template2) {
		/*
		 * List indices of atoms in the two residue templates that can be
		 * aligned to perform a mutation indices in first template go in ans[0],
		 * for second template go in ans[1]
		 * 
		 * Starting off with amino acids (L- or D- both supported, natural
		 * sidechains or not), add more later (until then, trying to mutate a
		 * non-amino acid residue will cause an error though)
		 * 
		 * The alignment will use RigidBodyMotion.superimposingMotion, which
		 * prioritizes alignment of earlier atoms so order matters
		 * 
		 * DZ: Now recognizes RNA as well.
		 * 
		 */

		int mutAtoms1[] = getTemplateMutAtoms(template1);
		int mutAtoms2[] = getTemplateMutAtoms(template2);

		if (template1.name.equalsIgnoreCase("PRO") || template2.name.equalsIgnoreCase("PRO")) {
			mutAtoms1 = getTemplateProMutAtoms(template1);
			mutAtoms2 = getTemplateProMutAtoms(template2);
		}

		for (int ma[] : new int[][] { mutAtoms1, mutAtoms2 }) {
			for (int atNum : ma) {
				if (atNum == -1) {// some atom(s) not found
					throw new UnsupportedOperationException(
							"ERROR: Mutation from " + template1.name + " to " + template2.name + " not supported yet");
				}
			}
		}

		return new int[][] { mutAtoms1, mutAtoms2 };
	}

	private static int[] getTemplateMutAtoms(ResidueTemplate template) {
		// Get atoms to use for alignment
		if (hasAminoAcidBB(template.templateRes)) {
			int N = template.templateRes.getAtomIndexByName("N");
			int CA = template.templateRes.getAtomIndexByName("CA");
			int C = template.templateRes.getAtomIndexByName("C");
			return new int[] { CA, N, C };
		} else if (hasNucleicAcidBB(template.templateRes)) {
			int O4 = template.templateRes.getAtomIndexByName("O4'");
			int C1 = template.templateRes.getAtomIndexByName("C1'");
			int C2 = template.templateRes.getAtomIndexByName("C2'");
			return new int[] { C1, O4, C2 };
		}
		return new int[] { -1 };
	}

	// If mutating to or from PRO need to make sure CD is aligned to H (its
	// non-Pro equivalent)
	// C and O will be OK...they are copied exactly
	private static int[] getTemplateProMutAtoms(ResidueTemplate template) {
		// Get atoms to use for alignment
		int N = template.templateRes.getAtomIndexByName("N");
		int CA = template.templateRes.getAtomIndexByName("CA");

		if (template.name.equalsIgnoreCase("PRO")) {// the PRO itself...use CD
			int CD = template.templateRes.getAtomIndexByName("CD");
			return new int[] { CA, N, CD };
		} else {// get H to align to CD
			int H = template.templateRes.getAtomIndexByName("H");
			return new int[] { CA, N, H };
		}
	}

	// The following functions are intended to mark inter-residue bonds between
	// amino acids
	// We're assuming for now that non-amino acid residues aren't bonded to
	// anything,
	// but could change that by changing just these functions
	// we might want a template file to tell this stuff eventually

	static double maxBondDist = 3;// liberal upper-bound on inter-residue bond
									// length. Given in angstroms

	// this version marks all the inter-residue bonds in a molecule
	public static void markInterResBonds(Molecule molec) {

		// first, peptide and phosphodiester bonds. Must be between consecutive
		// residues
		for (int resNum = 0; resNum < molec.residues.size() - 1; resNum++) {
			Residue res1 = molec.residues.get(resNum);
			Residue res2 = molec.residues.get(resNum + 1);
			if (!tryToMakeBond(res1, res2, InterResBondType.PEPTIDE)) {
				// DZ: only try to make a phosphodiester bond if the peptide
				// bond failed.
				tryToMakeBond(res1, res2, InterResBondType.PHOSPHODIESTER);
			}
		}

		// there could also be disulfide bonds. Look for any cysteines that are
		// close enough
		for (Residue res1 : molec.residues) {
			if (res1.fullName.startsWith("CYX")) {
				// unprotonated cysteine...could be in a disulfide bond
				for (Residue res2 : molec.residues) {
					if (res2.indexInMolecule < res1.indexInMolecule) {
						// don't double-count
						if (res2.fullName.startsWith("CYX")) {
							tryToMakeBond(res1, res2, InterResBondType.DISULFIDE);
						}
					}
				}
			}
		}

		// now all inter-res bonds are assumed to be made
		for (Residue res : molec.residues)
			res.interResBondsMarked = true;
	}

	// this version just tries to make inter-res bonds involving the specified
	// residue
	// it's intended for use right after res is mutated, to reconnect it to the
	// rest of the molecule
	public static void reconnectInterResBonds(Residue res) {

		// first try to make peptide bonds
		if (res.indexInMolecule > 0) {// res is not the first residue
			Residue prevRes = res.molec.residues.get(res.indexInMolecule - 1);
			if (!tryToMakeBond(prevRes, res, InterResBondType.PEPTIDE)) {
				// DZ: only try to make a phosphodiester bond if the peptide
				// bond
				// failed.
				tryToMakeBond(prevRes, res, InterResBondType.PHOSPHODIESTER);
			}
		}
		if (res.indexInMolecule < res.molec.residues.size() - 1) {// res is not
																	// the last
																	// residue
			Residue nextRes = res.molec.residues.get(res.indexInMolecule + 1);
			if (!tryToMakeBond(res, nextRes, InterResBondType.PEPTIDE)) {
				// DZ: only try to make a phosphodiester bond if the peptide
				// bond
				// failed.
				tryToMakeBond(res, nextRes, InterResBondType.PHOSPHODIESTER);
			}
		}

		// try to make any disulfide bonds, if res is an unprotonated cysteine
		if (res.fullName.startsWith("CYX")) {

			for (Residue res2 : res.molec.residues) {
				if (res2 != res) {
					if (res2.fullName.startsWith("CYX")) {
						tryToMakeBond(res, res2, InterResBondType.DISULFIDE);
					}
				}
			}
		}

		// OK that should be all for now
		res.interResBondsMarked = true;
	}

	/**
	 * Replaces tryToMakePeptideBond and tryToMakeDisulfideBond Given
	 * consecutive residues res1 and res2, make a peptide bond between them if
	 * appropriate
	 * 
	 * @param res1
	 *            first residue
	 * @param res2
	 *            second residue
	 * @param type
	 *            bond type
	 * @return true if bond was formed, false otherwise
	 */
	public static boolean tryToMakeBond(Residue res1, Residue res2, InterResBondType type) {

		if (res1.template == null || res2.template == null) {
			throw new RuntimeException("ERROR: Trying to bond residues without template");
		}
		int index1;
		int index2;
		switch (type) {
		case DISULFIDE:
			index1 = res1.getAtomIndexByName("SG");
			index2 = res2.getAtomIndexByName("SG");
			if (!(res1.fullName.startsWith("CYX") && res2.fullName.startsWith("CYX"))) {
				return false;
			}
			break;
		case PHOSPHODIESTER:
			index1 = res1.getAtomIndexByName("O3'");
			index2 = res2.getAtomIndexByName("P");
			if (!(hasNucleicAcidBB(res1) && hasNucleicAcidBB(res2))) {
				return false;
			}
			break;
		case PEPTIDE:
			index1 = res1.getAtomIndexByName("C");
			index2 = res2.getAtomIndexByName("N");
			if (!(hasAminoAcidBB(res1) && hasAminoAcidBB(res2))) {
				return false;
			}
			break;
		default:
			throw new RuntimeException("Bond type not recognized.");
		}

		if (index1 == -1 || index2 == -1) {
			return false;
			// atoms not found.
		}

		// Get distance between these atoms
		double dist = VectorAlgebra.distance(res1.coords, index1, res2.coords, index2);
		if (dist < maxBondDist) {
			Atom S1 = res1.atoms.get(index1);
			Atom S2 = res2.atoms.get(index2);
			S1.addBond(S2);
			return true;
		}
		return false;
	}

	// // DZ: now returns true if a bond was formed and false otherwise
	// public static boolean tryToMakeDisulfideBond(Residue res1, Residue res2)
	// {
	// // Given CYX residues res1 and res2, make a disulfide bond between them
	// // if appropriate
	// if (res1.template == null || res2.template == null)
	// throw new RuntimeException("ERROR: Trying to disulfide-bond residue
	// without template");
	//
	// int SIndex1 = res1.getAtomIndexByName("SG");
	// int SIndex2 = res2.getAtomIndexByName("SG");
	//
	// if (SIndex1 == -1 || SIndex2 == -1)
	// throw new RuntimeException("ERROR: Trying to disulfide-bond residue
	// without SG");
	//
	// // Get distance between these atoms
	// double SSDist = VectorAlgebra.distance(res1.coords, SIndex1, res2.coords,
	// SIndex2);
	// if (SSDist < maxBondDist) {
	// Atom S1 = res1.atoms.get(SIndex1);
	// Atom S2 = res2.atoms.get(SIndex2);
	// S1.addBond(S2);
	// return true;
	// }
	// return false;
	// }

	// public static boolean tryToMakePeptideBond(Residue res1, Residue res2) {
	// // Given consecutive residues res1 and res2, make a peptide bond between
	// // them if appropriate
	// if (res1.template == null || res2.template == null) {
	// throw new RuntimeException("ERROR: Trying to peptide-bond residue without
	// template");
	// }
	// if ((hasAminoAcidBB(res1) && hasAminoAcidBB(res2))) {// can only
	// // peptide-bond
	// // amino acids
	//
	// int CIndex = res1.getAtomIndexByName("C");
	// int NIndex = res2.getAtomIndexByName("N");
	//
	// // Get distance between these atoms
	// double CNDist = VectorAlgebra.distance(res1.coords, CIndex, res2.coords,
	// NIndex);
	// if (CNDist < maxBondDist) {
	// Atom C = res1.atoms.get(CIndex);
	// Atom N = res2.atoms.get(NIndex);
	// C.addBond(N);
	// return true;
	// }
	// }
	// return false;
	// }

	public static boolean hasAminoAcidBB(Residue res) {
		// does res have the three main amino-acid backbone atoms (N,CA, and C)?
		// This would allow it to peptide-bond
		return (res.getAtomIndexByName("N") >= 0 && res.getAtomIndexByName("CA") >= 0
				&& res.getAtomIndexByName("C") >= 0); // return true only if we
														// found them all.
	}

	public static boolean hasNucleicAcidBB(Residue res) {
		// does res have an NA backbone
		// This would allow it to form phosphodiester bonds
		// we check the carbons and oxygens found in both ribose and deoxyribose
		return (res.getAtomIndexByName("O5'") >= 0 && res.getAtomIndexByName("O4'") >= 0
				&& res.getAtomIndexByName("O3'") >= 0 && res.getAtomIndexByName("C1'") >= 0
				&& res.getAtomIndexByName("C2'") >= 0 && res.getAtomIndexByName("C3'") >= 0
				&& res.getAtomIndexByName("C4'") >= 0 && res.getAtomIndexByName("C5'") >= 0);
	}

	public static String getTemplateName(Residue res) {
		// some residues may have template names different than what's in the
		// full name for the residue
		// we identify those here
		// we cover HIS and CYS; we assume templates follow the HID/HIE/HIP and
		// CYS/CYX conventions as in
		// all_amino94.in

		String resName = res.fullName.substring(0, 3);

		if (resName.equalsIgnoreCase("HIS")) {

			int HDCount = 0;// Count delta and epsilon hydrogens
			int HECount = 0;

			for (Atom at : res.atoms) {
				if (at.name.contains("HD"))
					HDCount++;
				else if (at.name.contains("HE"))
					HECount++;
			}

			// Determine protonation state, assign template name accordingly
			if ((HDCount == 1) && (HECount == 2))
				return "HIE";
			else if ((HDCount == 2) && (HECount == 1))
				return "HID";
			else if ((HDCount == 2) && (HECount == 2))
				return "HIP";
			else {
				throw new RuntimeException("ERROR: Invalid protonation state for " + res.fullName);
			}
		} else if (resName.equalsIgnoreCase("CYS")) {

			if (res.getAtomIndexByName("HG") >= 0) // there is a gamma hydrogen
				return "CYS";
			else
				return "CYX";
		} else// by default, the template name is just the first three letters
				// of the full name
			return resName;
	}

	// Stuff about protonation-dependent rotamers for K* (Hie,Hid,etc.) could go
	// in this class too...

}
