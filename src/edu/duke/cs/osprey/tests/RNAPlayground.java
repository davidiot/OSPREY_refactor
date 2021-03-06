package edu.duke.cs.osprey.tests;

import java.util.ArrayList;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.deeper.perts.Backrub;
import edu.duke.cs.osprey.dof.deeper.perts.PartialStructureSwitch;
import edu.duke.cs.osprey.dof.deeper.perts.RingPucker;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

public class RNAPlayground {
	public static void runAllTests() {
		testDihedral();
		measureAngles();
		testProteinSwitch();
		testBackrub();
		showRingPucker();
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

	public static void showRingPucker() {
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
