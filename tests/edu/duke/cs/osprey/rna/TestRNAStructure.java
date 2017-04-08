package edu.duke.cs.osprey.rna;

import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;

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
    
	@Test
	public void testPuckerDetection() {
		
		String folder = "test/354dH.junit/";
		Molecule m1 = PDBFileReader.readPDBFile(folder + "354dH.pdb");
		Molecule m2 = PDBFileReader.readPDBFile(folder + "354dH_unspecified_pucker.pdb");
		String m1FileName = folder + "testResults/354dH1.pdb";
		String m2FileName = folder + "testResults/354dH2.pdb";
		PDBFileWriter.writePDBFile(m1, m1FileName);
		PDBFileWriter.writePDBFile(m2, m2FileName);
		
		try {
			byte[] m1Bytes = Files.readAllBytes(Paths.get(m1FileName));
			byte[] m2Bytes = Files.readAllBytes(Paths.get(m2FileName));

			String m1String = new String(m1Bytes, StandardCharsets.UTF_8);
			String m2String = new String(m2Bytes, StandardCharsets.UTF_8);

			assertThat(m1String, is(m2String));
		} catch (IOException e) {
			throw new RuntimeException("FILES DIFFERED");
		}
	}
	
}
