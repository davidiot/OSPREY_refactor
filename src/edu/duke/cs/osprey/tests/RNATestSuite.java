package edu.duke.cs.osprey.tests;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.ResBBState;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.dof.deeper.perts.Backrub;
import edu.duke.cs.osprey.dof.deeper.perts.PartialStructureSwitch;
import edu.duke.cs.osprey.dof.deeper.perts.RingPucker;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

public class RNATestSuite {
	public static void runAllTests() {
		// test1cslHEnergy();
		// testMutation();
		// testDihedral();
		// measureAngles();
		// testProteinSwitch();
		// testRNASwitch();
		// testBackrub();
		// testRingPucker();
		try {
			analyze3p49tau();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static final String[] tauAngles = new String[] { "N", "CA", "C" };
	private static final String[] rnaAngles = new String[] { "C1'", "C2'", "O2'" };

	private static double measureTauAngles(Residue res) {
		double[] coordsA = res.getCoordsByAtomName(tauAngles[0]);
		double[] coordsB = res.getCoordsByAtomName(tauAngles[1]);
		double[] coordsC = res.getCoordsByAtomName(tauAngles[2]);
		return Protractor.getAngleDegrees(coordsA, coordsB, coordsC);
	}

	private static double measureRNAAngles(Residue res) {
		double[] coordsA = res.getCoordsByAtomName(rnaAngles[0]);
		double[] coordsB = res.getCoordsByAtomName(rnaAngles[1]);
		double[] coordsC = res.getCoordsByAtomName(rnaAngles[2]);
		return Protractor.getAngleDegrees(coordsA, coordsB, coordsC);
	}

	public static void analyze3p49() throws FileNotFoundException {
		PrintWriter output = new PrintWriter("testResults/3p49.txt");
		for (int i = 723; i < 730; i++) {
			ArrayList<Double> in = new ArrayList<Double>();
			ArrayList<Double> out = new ArrayList<Double>();
			for (int j = -1000; j <= 1000; j++) {
				Molecule m = PDBFileReader.readPDBFile("3P49FH.pdb");
				Residue r = m.getResByPDBResNumber(i + "");
				ArrayList<Residue> affected = new ArrayList<Residue>();
				affected.add(r);
				double originalAngle = measureRNAAngles(r);
				RingPucker rp = new RingPucker(affected);
				rp.doPerturbationMotion(j / 100.0);
				double newAngle = measureRNAAngles(r);
				in.add(j / 100.0);
				out.add(newAngle - originalAngle);
			}
			output.println(i);
			output.println(out.toString());
		}
		output.close();
	}

	public static void analyze3p49tau() throws FileNotFoundException {
		PrintWriter output = new PrintWriter("testResults/3p49tau.txt");
		for (int i : new int[] { 79, 83, 85, 88, 89, 91 }) {
			ArrayList<Double> in = new ArrayList<Double>();
			ArrayList<Double> out = new ArrayList<Double>();
			for (int j = -250; j <= 250; j++) {
				Molecule m = PDBFileReader.readPDBFile("3P49FH.pdb");
				Residue r = m.getResByPDBResNumber(i + "");
				ArrayList<Residue> affected = new ArrayList<Residue>();
				affected.add(m.getResByPDBResNumber((i - 1) + ""));
				affected.add(r);
				affected.add(m.getResByPDBResNumber((i + 1) + ""));
				double originalAngle = measureTauAngles(m.getResByPDBResNumber((i + 1) + ""));
				Backrub br = new Backrub(affected);
				br.doPerturbationMotion(j / 100.0);
				double newAngle = measureTauAngles(m.getResByPDBResNumber((i + 1) + ""));
				in.add(j / 100.0);
				out.add(newAngle - originalAngle);
			}
			output.println(i);
			output.println(out.toString());
		}
		output.close();
	}

	public static void testBackrub() {
		Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		ArrayList<Residue> affected = new ArrayList<Residue>();
		affected.add(m.residues.get(55));
		affected.add(m.residues.get(56));
		affected.add(m.residues.get(57));
		Backrub rub = new Backrub(affected);
		rub.doPerturbationMotion(10);
		// superimpose old and new PDBs
		Molecule m2 = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		affected.get(0).indexInMolecule = 90;
		affected.get(0).fullName = "LEU A  90 ";
		affected.get(1).indexInMolecule = 91;
		affected.get(1).fullName = "GLU A  91 ";
		affected.get(2).indexInMolecule = 92;
		affected.get(2).fullName = "LYS A  92 ";
		m2.appendResidue(affected.get(0));
		m2.appendResidue(affected.get(1));
		m2.appendResidue(affected.get(2));
		PDBFileWriter.writePDBFile(m2, "testResults/1CC8rub.pdb");
	}

	public static void testRingPucker() {
		Molecule molecule = PDBFileReader.readPDBFile("354dH.pdb");
		for (int i : new int[] { -20, -10, 10, 20 }) {
			Molecule m = PDBFileReader.readPDBFile("354dH.pdb");
			ArrayList<Residue> affected = new ArrayList<Residue>();
			affected.add(m.residues.get(0));
			RingPucker rp = new RingPucker(affected);
			rp.doPerturbationMotion(i);
			affected.get(0).indexInMolecule = 200 + i;
			affected.get(0).fullName = "RC3 C  " + affected.get(0).indexInMolecule + " ";
			molecule.appendResidue(affected.get(0));
		}
		PDBFileWriter.writePDBFile(molecule, "testResults/354dHpuck.pdb");
	}

	public static void testProteinSwitch() {
		Molecule m = PDBFileReader.readPDBFile("1CC8hel.pdb");
		ArrayList<String> altConfPDBFiles = new ArrayList<String>();
		altConfPDBFiles.add("1CC8sheet.pdb");
		PartialStructureSwitch switcheroo = new PartialStructureSwitch(m.residues, altConfPDBFiles);
		switcheroo.doPerturbationMotion(1);
		PDBFileWriter.writePDBFile(m, "testResults/1CC8switch.pdb");
	}

	public static void testRNASwitch() {
		Molecule m = PDBFileReader.readPDBFile("354dH.pdb");
		ArrayList<String> altConfPDBFiles = new ArrayList<String>();
		altConfPDBFiles.add("1cslH.renum.pdb");
		ArrayList<Double> dependentGenChi1 = new ArrayList<>();
		for (Residue res : m.residues) {
			dependentGenChi1.add(GenChi1Calc.getGenChi1(res));
			// record gen chi1 so we can restore it later
		}
		PartialStructureSwitch switcheroo = new PartialStructureSwitch(m.residues, altConfPDBFiles);
		switcheroo.doPerturbationMotion(1);
		for (int resNum = 0; resNum < m.residues.size(); resNum++) {
			Residue res = m.residues.get(resNum);
			SidechainIdealizer.idealizeSidechain(res);
			GenChi1Calc.setGenChi1(res, dependentGenChi1.get(resNum));
		}
		PDBFileWriter.writePDBFile(m, "testResults/354dHswitch.pdb");
	}

	public static void measureAngles() {
		Molecule m = PDBFileReader.readPDBFile("1cslH.pdb");
		for (Residue res : m.residues) {
			double C5[] = res.getCoordsByAtomName("C5'");
			double C4[] = res.getCoordsByAtomName("C4'");
			double C3[] = res.getCoordsByAtomName("C3'");
			double O3[] = res.getCoordsByAtomName("O3'");

			double chi1Measured = Protractor.measureDihedral(new double[][] { C5, C4, C3, O3 });

			System.out.println(res.fullName + ": " + chi1Measured);
			System.out.println(res.template.name);
		}
	}

	public static void test1cslHEnergy() {
		Molecule m = PDBFileReader.readPDBFile("1cslH.pdb");
		EnergyFunction fullEFunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(m);
		double fullE = fullEFunc.getEnergy();
		System.out.println("1cslH.pdb full energy: " + fullE);
	}

	public static void testMutation() {
		Molecule m = PDBFileReader.readPDBFile("1cslH.pdb");
		PDBFileWriter.writePDBFile(m, "testResults/1cslH.pdb");
		String oldS = m.residues.get(5).fullName;
		ResidueTypeDOF mutDOF = new ResidueTypeDOF(m.residues.get(5)); // Originally
																		// RG3 A
																		// 48
		mutDOF.mutateTo("RU3");
		String newS = m.residues.get(5).fullName;
		PDBFileWriter.writePDBFile(m, "testResults/1cslH.G48U.pdb");
		System.out.println("Mutated " + oldS + " to " + newS);
		System.out.println("Wrote original file: testResults/1cslH.pdb");
		System.out.println("Wrote mutation test output: testResults/1cslH.G48U.pdb");
	}

	public static void testDihedral() {

		Molecule m = PDBFileReader.readPDBFile("1cslH.pdb");
		Residue res = m.residues.get(5); // RG3 A 48

		FreeDihedral chi1 = new FreeDihedral(res, 0);

		double O4[] = res.getCoordsByAtomName("O4'");
		double C1[] = res.getCoordsByAtomName("C1'");
		double N9[] = res.getCoordsByAtomName("N9");
		double C4[] = res.getCoordsByAtomName("C4");

		double chi1Measured = Protractor.measureDihedral(new double[][] { O4, C1, N9, C4 });

		System.out.println("chi1 originally measured as " + chi1Measured);

		chi1.apply(55.);

		// measure dihedrals. Start by collecting coordinates
		O4 = res.getCoordsByAtomName("O4'");
		C1 = res.getCoordsByAtomName("C1'");
		N9 = res.getCoordsByAtomName("N9");
		C4 = res.getCoordsByAtomName("C4");

		chi1Measured = Protractor.measureDihedral(new double[][] { O4, C1, N9, C4 });

		System.out.println("chi1 applied as 55, measured as " + chi1Measured);
		System.out.println("Outputting dihedral-adjusted structure as testResults/1cslH.dih.pdb");
		PDBFileWriter.writePDBFile(m, "testResults/1cslH.dih.pdb");
	}
}
