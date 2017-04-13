package edu.duke.cs.osprey.rna;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;
import org.junit.Test;


import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

/**
 * Unit tests pertaining to RNA GMEC functionality
 *
 * @author dzhou
 */
public class TestRNAGMEC extends TestBase{

    /**
     * Find a GMEC for an RNA double strand with some base flips
     */
    @Test
    public void test1cslHGMEC() {

        String folder = "test/1CSLH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-48.751, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {0, 0, 0 , 38}));
    }

    /**
     * Find a GMEC when both RNA and protein are involved
     * Note that the stack size might need to be increased.
     */
    @Test
    public void test1fxlFHGMEC() {

        String folder = "test/1fxlFH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-61.737, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {33, 1, 14}));
    }

    /**
     * Tests GMEC with protein and RNA simultaneously.
     * Same as the DEEPer configuration only we do not allow perturbations
     * ran with -Xss10m
     */
    @Test
    public void test3P49FHGMEC(){
        String folder = "test/3P49FH.junit/";
        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively( -91.516, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {0, 20, 42, 0}));
    }

    /**
     * Tests DEEPer with RNA ring puckers
     * ran with -Xss10m
     */
    @Test
    public void test2OEUHDEEPer(){
        String folder = "test/2OEUH.deeper.junit/";
        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-64.117, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {4, 41, 30}));
    }

    /**
     * Tests DEEPer with protein and RNA simultaneously.  One backrub and two ring puckers.
     * ran with -Xss10m
     */
    @Test
    public void test3P49FHDEEPer(){
        String folder = "test/3P49FH.deeper.junit/";
        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-91.002, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {0, 20, 42, 0}));
    }

    private ConfSearch.EnergiedConf checkGMEC(String[] args) {
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        cfp.loadData();

        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        return gf.calcGMEC().get(0);
    }
}
