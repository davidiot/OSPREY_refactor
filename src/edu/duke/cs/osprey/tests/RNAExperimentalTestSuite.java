package edu.duke.cs.osprey.tests;

import java.util.ArrayList;
import edu.duke.cs.osprey.dof.deeper.perts.RNABackrub;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/*
 * A bunch of test for determine the optimal scaling factor
 */
public class RNAExperimentalTestSuite {

	private static final String[] atomNamesForAngle1 = new String[] { "O3'", "C3'", "C2'" };
	private static final String[] atomNamesForAngle2 = new String[] { "C2'", "C1'", "O4'" };
	private static final String[] atomNamesForAngle3 = new String[] { "O4'", "C4'", "C5'" };
	private static final String[][] atomNamesForAngles = new String[][] { atomNamesForAngle1, atomNamesForAngle2,
			atomNamesForAngle3 }; // tau atom analogues

	private static final String[] files = new String[] { "1CC8.ss.pdb", "1CC8hel.pdb", "1CC8sheet.pdb", "1cslH.pdb",
			"1cslH.renum.pdb", "1fxlFH.pdb", "1igdFH.pdb", "2oeuH.pdb", "354dH.pdb", "3IWNFH.pdb", "3MXHFH.pdb",
			"3P49FH.pdb" };

	public static void runAllTests() {
		investigateRNABackrub();
	}

	public static void investigateRNABackrub() {

		int min = -101;
		double minStrain = Double.MAX_VALUE;
		for (int x = -100; x < 100; x++) {
			int resCount = 0;
			double totalStrain = 0;
			double endStrain = 0;
			double scaling = x / 100.0;
			for (String file : files) {
				Molecule m1 = PDBFileReader.readPDBFile(file);
				for (Residue r : m1.residues) {
					if (HardCodedResidueInfo.hasNucleicAcidBB(r)) {
						resCount++;
						double[] originalAngles = measureAngles(r);
						ArrayList<Residue> affected = new ArrayList<Residue>();
						affected.add(r);
						RNABackrub rub = new RNABackrub(affected, scaling);
						rub.doPerturbationMotion(2.5);
						double[] newAngles = measureAngles(r);
						double[] differences = new double[3];
						for (int i = 0; i < differences.length; i++) {
							differences[i] = Math.abs(newAngles[i] - originalAngles[i]);
							totalStrain += differences[i];
						}
						endStrain += differences[0] + differences[2];
					}
				}
			}
			double averageTotalStrain = totalStrain / resCount / 3;
			System.out.println(averageTotalStrain);
			double averageEndStrain = endStrain / resCount / 2;
			System.out.println(averageEndStrain);
			if (averageTotalStrain < minStrain) {
				minStrain = averageTotalStrain;
				min = x;
			}
		}
		System.out.println(min);
	}

	private static double[] measureAngles(Residue res) {
		double[] results = new double[3];
		for (int i = 0; i < atomNamesForAngles.length; i++) {
			String[] atomNames = atomNamesForAngles[i];
			double[] coordsA = res.getCoordsByAtomName(atomNames[0]);
			double[] coordsB = res.getCoordsByAtomName(atomNames[1]);
			double[] coordsC = res.getCoordsByAtomName(atomNames[2]);
			results[i] = Protractor.getAngleDegrees(coordsA, coordsB, coordsC);
		}
		return results;
	}

}
