package edu.duke.cs.osprey.rna;

import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.dof.deeper.perts.PartialStructureSwitch;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.*;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;

/**
 * Unit tests pertaining to RNA structure functionality
 *
 * @author dzhou
 */
public class TestRNAStructure extends TestBase {
	
    @BeforeClass
    public static void before() {
    	initRNAEnvironment();
    }

	/**
	 * Tests pucker detection functionality.
	 * Loading a file with unspecified puckers will automatically have them encoded.
	 */
	@Test
	public void test354dHPuckerDetection() {
		String folder = "test/354dH.junit/";
		Molecule m1 = PDBFileReader.readPDBFile(folder + "354dH.pdb");
		Molecule m2 = PDBFileReader.readPDBFile(folder + "354dH_unspecified_pucker.pdb");
		String m1FileName = folder + "testResults/354dH1.pdb";
		String m2FileName = folder + "testResults/354dH2.pdb";
		PDBFileWriter.writePDBFile(m1, m1FileName);
		PDBFileWriter.writePDBFile(m2, m2FileName);
		comparePDB(m1FileName, m2FileName);
	}

	/**
	 * If the energies are ridiculous, chances are that bonds are not being properly formed.
	 */
	@Test
	public void test1cslHEnergy() {
		String folder = "test/1CSLH.junit/";
		Molecule m = PDBFileReader.readPDBFile(folder + "1cslH.pdb");
		EnergyFunction fullEFunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(m);
		assertThat(fullEFunc.getEnergy(), isRelatively(-1011.018, 1e-3));
	}

	@Test
	public void testRNASwitch() {
		String folder = "test/354dH.junit/";
		Molecule m = PDBFileReader.readPDBFile(folder + "354dH.pdb");
		ArrayList<String> altConfPDBFiles = new ArrayList<String>();
		altConfPDBFiles.add(folder + "1cslH.renum.pdb");
		ArrayList<Double> dependentGenChi1 = new ArrayList<>();
		for (Residue res : m.residues) {
			dependentGenChi1.add(GenChi1Calc.getGenChi1(res));
			// record gen chi1 so we can restore it later
		}
		PartialStructureSwitch switcheroo = new PartialStructureSwitch(m.residues, altConfPDBFiles);
		switcheroo.doPerturbationMotion(1);

        Molecule m2 = PDBFileReader.readPDBFile(folder + "1cslH.renum.pdb");
		for (int resNum = 0; resNum < m.residues.size(); resNum++) {
			Residue res1 = m.residues.get(resNum);
            Residue res2 = m2.getResByPDBResNumber(res1.getPDBResNumber());
            for (String bbAtom : HardCodedResidueInfo.possibleNABBAtoms) {
                double[] c1 = res1.getCoordsByAtomName(bbAtom);
                double[] c2 = res2.getCoordsByAtomName(bbAtom);
                if (c1 != null && c2 != null && !res1.getPDBResNumber().equals("96")) {
                    for (int i = 0; i < c1.length; i++) {
                        assertThat(c1[i], isRelatively(c2[i], 1e-9));
                    }
                }
            }

            SidechainIdealizer.idealizeSidechain(res1);
            GenChi1Calc.setGenChi1(res1, dependentGenChi1.get(resNum));
            assertThat(GenChi1Calc.getGenChi1(res1), isRelatively(dependentGenChi1.get(resNum), 1e-9));

		}

		PDBFileWriter.writePDBFile(m, folder + "testResults/354dHswitch.pdb");
	}

}
