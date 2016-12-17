package edu.duke.cs.osprey.tests;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.dof.deeper.perts.Backrub;
import edu.duke.cs.osprey.dof.deeper.perts.RNABackrub;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/*
 * A bunch of test for determine the optimal scaling factor
 */
public class RNAExperimentalTestSuite {

	private static final String[] tauAngles = new String[] { "N", "CA", "C" };
	private static final String[] atomNamesForAngle1 = new String[] { "O3'", "C3'", "C2'" };
	private static final String[] atomNamesForAngle2 = new String[] { "C2'", "C1'", "O4'" };
	private static final String[] atomNamesForAngle3 = new String[] { "O4'", "C4'", "C5'" };
	private static final String[][] atomNamesForAngles = new String[][] { atomNamesForAngle1, atomNamesForAngle2,
			atomNamesForAngle3 }; // tau atom analogues

	private static final String[] files = new String[] { "1cslH.pdb" };

	public static void runAllTests() {
		// printBackrub();
		// investigateBackrub();
		investigateRNABackrub();
	}

	public static void printBackrub() {

		Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");

		Molecule m1 = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		ArrayList<Residue> affected1 = new ArrayList<Residue>();
		affected1.add(m1.residues.get(55));
		affected1.add(m1.residues.get(56));
		affected1.add(m1.residues.get(57));
		Backrub rub = new Backrub(affected1, 0);
		rub.doPerturbationMotion(10);
		affected1.get(0).indexInMolecule = 90;
		affected1.get(0).fullName = "LEU A  90 ";
		affected1.get(1).indexInMolecule = 91;
		affected1.get(1).fullName = "GLU A  91 ";
		affected1.get(2).indexInMolecule = 92;
		affected1.get(2).fullName = "LYS A  92 ";
		m.appendResidue(affected1.get(0));
		m.appendResidue(affected1.get(1));
		m.appendResidue(affected1.get(2));

		Molecule m2 = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		ArrayList<Residue> affected2 = new ArrayList<Residue>();
		affected2.add(m2.residues.get(55));
		affected2.add(m2.residues.get(56));
		affected2.add(m2.residues.get(57));
		Backrub rub2 = new Backrub(affected2, 0.7);
		rub2.doPerturbationMotion(10);
		affected2.get(0).indexInMolecule = 93;
		affected2.get(0).fullName = "LEU A  93 ";
		affected2.get(1).indexInMolecule = 94;
		affected2.get(1).fullName = "GLU A  94 ";
		affected2.get(2).indexInMolecule = 95;
		affected2.get(2).fullName = "LYS A  95 ";
		m.appendResidue(affected2.get(0));
		m.appendResidue(affected2.get(1));
		m.appendResidue(affected2.get(2));

		Molecule m3 = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		ArrayList<Residue> affected3 = new ArrayList<Residue>();
		affected3.add(m3.residues.get(55));
		affected3.add(m3.residues.get(56));
		affected3.add(m3.residues.get(57));
		Backrub rub3 = new Backrub(affected3, 1);
		rub3.doPerturbationMotion(10);
		affected3.get(0).indexInMolecule = 96;
		affected3.get(0).fullName = "LEU A  96 ";
		affected3.get(1).indexInMolecule = 97;
		affected3.get(1).fullName = "GLU A  97 ";
		affected3.get(2).indexInMolecule = 98;
		affected3.get(2).fullName = "LYS A  98 ";
		m.appendResidue(affected3.get(0));
		m.appendResidue(affected3.get(1));
		m.appendResidue(affected3.get(2));

		PDBFileWriter.writePDBFile(m, "testResults/1CC8.backrubcomparison.pdb");
	}

	public static void investigateBackrub() {
		for (int x = 0; x < 60; x++) {
			Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");
			ArrayList<Residue> affected = new ArrayList<Residue>();
			affected.add(m.residues.get(x));
			affected.add(m.residues.get(x + 1));
			affected.add(m.residues.get(x + 2));
			double[] originalAngles = measureTauAngles(affected);
			Backrub rub = new Backrub(affected, 2);
			rub.doPerturbationMotion(2.5);
			double[] newAngles = measureTauAngles(affected);
			double[] differences = new double[3];
			for (int i = 0; i < differences.length; i++) {
				differences[i] = Math.abs(newAngles[i] - originalAngles[i]);
			}
			System.out.println(Arrays.toString(differences));
		}
	}

	public static void investigateRNABackrub() {

		int min = -101;
		double minStrain = Double.MAX_VALUE;
		for (int x = -10; x < 10; x++) {
			int resCount = 0;
			double scaling = x / 10.0;
			for (String file : files) {
				Molecule m1 = PDBFileReader.readPDBFile(file);
				int count = 0;
				for (Residue r : m1.residues) {
					if (HardCodedResidueInfo.hasNucleicAcidBB(r)) {
						resCount++;
						double[] originalAngles = measureRNAAngles(r);
						ArrayList<Residue> affected = new ArrayList<Residue>();
						affected.add(r);
						RNABackrub rub = new RNABackrub(affected, scaling);
						rub.doPerturbationMotion(2.5);
						double[] newAngles = measureRNAAngles(r);
						double[] differences = new double[3];
						for (int i = 0; i < differences.length; i++) {
							differences[i] = Math.abs(newAngles[i] - originalAngles[i]);
							if (newAngles[i] > 116.5 || newAngles[i] < 105.5) {
								count++;
								System.out.println(r.fullName);
							}
						}
						System.out.println(scaling);
						System.out.println(Arrays.toString(newAngles));
					}
				}
				System.out.println("Invalid perturbations: " + count);
			}
		}
	}

	private static double[] measureTauAngles(ArrayList<Residue> res) {
		double[] results = new double[3];
		for (int i = 0; i < atomNamesForAngles.length; i++) {
			double[] coordsA = res.get(i).getCoordsByAtomName(tauAngles[0]);
			double[] coordsB = res.get(i).getCoordsByAtomName(tauAngles[1]);
			double[] coordsC = res.get(i).getCoordsByAtomName(tauAngles[2]);
			results[i] = Protractor.getAngleDegrees(coordsA, coordsB, coordsC);
		}
		return results;
	}

	private static double[] measureRNAAngles(Residue res) {
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
