package edu.duke.cs.osprey.rna;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;

import java.util.ArrayList;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
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
     * Tests proper bonding and energy
	 * Note: If the energies are ridiculous, chances are that bonds are not being properly formed.
	 */
	@Test
	public void test1cslHEnergy() {
		String folder = "test/1CSLH.junit/";
		Molecule m = PDBFileReader.readPDBFile(folder + "1cslH.pdb");
		EnergyFunction fullEFunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(m);
		assertThat(fullEFunc.getEnergy(), isRelatively(-1011.018, 1e-3));
	}

	@Test
    /**
     * Tests full and partial structure switches
     * Note: The HO3' hydrogens do not move because they do not exist in the switched 1cslH structure
     */
	public void testRNASwitch() {
		String folder = "test/354dH.junit/";
		Molecule m = PDBFileReader.readPDBFile(folder + "354dH.pdb");
		ArrayList<String> altConfPDBFiles = new ArrayList<String>();
		altConfPDBFiles.add(folder + "1cslH.renum.pdb");
		PartialStructureSwitch switcheroo = new PartialStructureSwitch(m.residues, altConfPDBFiles);
		switcheroo.doPerturbationMotion(1);

        Molecule m2 = PDBFileReader.readPDBFile(folder + "1cslH.renum.pdb");
		for (int resNum = 0; resNum < m.residues.size(); resNum++) {
			Residue res1 = m.residues.get(resNum);
            Residue res2 = m2.getResByPDBResNumber(res1.getPDBResNumber());
            for (String bbAtom : HardCodedResidueInfo.possibleNABBAtoms) {
                double[] c1 = res1.getCoordsByAtomName(bbAtom);
                double[] c2 = res2.getCoordsByAtomName(bbAtom);
                if (c1 != null && c2 != null) {
                    for (int i = 0; i < c1.length; i++) {
                        assertThat(c1[i], isRelatively(c2[i], 1e-9));
                    }
                }
            }
            SidechainIdealizer.idealizeSidechain(res1);
		}

		PDBFileWriter.writePDBFile(m, folder + "testResults/354dHswitch.pdb");
	}

    @Test
    /**
     * Tests Chi1 angle modification for DEEPer
     */
    public void testRNAChi1() {
        String folder = "test/2OEUH.junit/";
        Molecule m = PDBFileReader.readPDBFile(folder + "2oeuH.pdb");
        for (Residue res : m.residues) {
            double chi1 = GenChi1Calc.getGenChi1(res);
            GenChi1Calc.setGenChi1(res, chi1 / 2);
            assertThat(GenChi1Calc.getGenChi1(res), isRelatively(chi1 / 2, 1e-9));
        }
        PDBFileWriter.writePDBFile(m, folder + "testResults/2oeuH.chi.pdb");
    }

    @Test
    /**
     * Tests mutations
     */
    public void testMutation() {
        String folder = "test/1CSLH.junit/";
        Molecule m = PDBFileReader.readPDBFile(folder + "1cslH.pdb");

        Residue res1 = m.getResByPDBResNumber("68");
        assertEquals(res1.fullName, "RA3 B  68 ");
        ResidueTypeDOF mutDOF1 = new ResidueTypeDOF(res1);
        mutDOF1.mutateTo("RG3");
        assertEquals(res1.fullName, "RG3 B  68 ");

        Residue res2 = m.getResByPDBResNumber("72");
        assertEquals(res2.fullName, "RU2 B  72 ");
        ResidueTypeDOF mutDOF2 = new ResidueTypeDOF(res2);
        mutDOF2.mutateTo("RC2");
        assertEquals(res2.fullName, "RC2 B  72 ");

        PDBFileWriter.writePDBFile(m, folder + "testResults/1cslH.A68G.U72C.pdb");
    }

}
