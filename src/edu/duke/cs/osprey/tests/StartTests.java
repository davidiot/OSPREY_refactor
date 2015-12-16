package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;

public class StartTests {
	public static void test() {
		System.out.println("START");
		Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");
		EnergyFunction fullEFunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(m);
		double fullE = fullEFunc.getEnergy();

		System.out.println("1CC8.ss.pdb full energy: " + fullE);

		for (Residue res : m.residues) {
			System.out.println(res.fullName);
			for (int i = 0; i < res.template.dihedral4Atoms.length; i++) {
				for (int j = 0; j < res.template.dihedral4Atoms[0].length; j++) {
					System.out.print(res.template.dihedral4Atoms[i][j]);
				}
			}
			System.out.println();
		}

		System.out.println("END");
	}
}
