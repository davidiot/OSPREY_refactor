package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

public class RNATestSuite {
	public static void runAllTests() {
		// UnitTestSuite.runAllTests();
		test1cslHEnergy();
		testMutation();
		testDihedral();
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
